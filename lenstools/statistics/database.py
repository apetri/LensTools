"""

.. module:: database
	:platform: Unix
	:synopsis: This module is a convenient wrapper of utilities for I/O between SQL databases and Ensembles


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division,print_function,with_statement
from .ensemble import Ensemble
from ..simulations.logs import logdriver
from ..utils.decorators import Parallelize

import numpy as np
import pandas as pd

try:
	import sqlalchemy
	sqlalchemy=sqlalchemy
except ImportError:
	sqlalchemy=None

#########################
#####Database class######
#########################

class Database(object):

	_constructor_ensemble = Ensemble

	#Create a connection to a database
	def __init__(self,name,connection_string="sqlite:///{0}"):

		if sqlalchemy is None:
			raise ImportError("sqlalchemy is not installed!!")

		self.connection = sqlalchemy.create_engine(connection_string.format(name))


	#Set constructor for the query results
	def set_constructor(self,constructor):
		assert issubclass(constructor,Ensemble),"The constructor should be a sub-class of Ensemble"
		self._constructor_ensemble = constructor


	#For context manager
	def __enter__(self):
		return self

	def __exit__(self,type,value,tb):
		self.connection.dispose()

	#Insert records in the database
	def insert(self,df,table_name="data"):

		"""
		:param df: records to insert in the database, in Ensemble (or pandas DataFrame) format
		:type df: :py:class:`Ensemble`

		"""

		assert isinstance(df,pd.DataFrame)
		df.to_sql(table_name,self.connection,if_exists="append",index=False)

	#Query the database
	def query(self,sql):

		"""
		:param sql: sql query string
		:type sql: str.

		:returns: :py:class:`Ensemble`

		"""

		return self._constructor_ensemble.read_sql_query(sql,self.connection)

	#Query a list of databases and combine the results
	@classmethod
	def query_all(cls,db_names,sql):

		"""
		Perform the same SQL query on a list of databases and combine the results

		:param db_names: list of names of the databases to query
		:type db_names: list.

		:param sql: sql query string
		:type sql: str.

		:returns: :py:class:`Ensemble`

		"""

		all_results = list()

		#Query each database
		for db_name in db_names:
			with cls(db_name) as db:
				all_results.append(db.query(sql))

		#Combine and return
		return cls._constructor_ensemble.concat(all_results,axis=0,ignore_index=True)

	#Visualize information about a table in the database
	def info(self,table_name="data"):
		assert table_name in self.tables,"Table {0} does not exist!".format(table_name)
		return self.query("PRAGMA table_info({0})".format(table_name))

	#List tables in the database
	@property
	def tables(self):
		return self.connection.table_names()

	#Read table in a database
	def read_table(self,table_name):

		if hasattr(table_name,"format"):
			assert table_name in self.tables,"Table {0} does not exist!".format(table_name)
			return self._constructor_ensemble.read_sql_table(table_name,self.connection) 
		elif type(table_name)==int:
			return self._constructor_ensemble.read_sql_table(self.tables[table_name],self.connection)
		else:
			raise TypeError("table_name type not recognized")

	#Read table in a list of databases and combine the results
	@classmethod
	def read_table_all(cls,db_names,table_name):


		"""
		Read the same SQL table from a list of databases and combine the results

		:param db_names: list of names of the databases to query
		:type db_names: list.

		:param table: table to read
		:type table: str.

		:returns: :py:class:`Ensemble`

		"""

		all_results = list()

		#Query each database
		for db_name in db_names:
			with cls(db_name) as db:
				all_results.append(db.read_table(table_name))

		#Combine and return
		return cls._constructor_ensemble.concat(all_results,axis=0,ignore_index=True)


################################
#####ScoreDatabase class########
################################

class ScoreDatabase(Database):

	def __init__(self,*args,**kwargs):
		super(ScoreDatabase,self).__init__(*args,**kwargs)
		self._parameters = ["Om","w","sigma8"]

	def set_parameters(self,parameters):
		self._parameters = parameters

	@property 
	def parameters(self):
		return self._parameters

	def pull_features(self,feature_list,table_name="scores",score_type="likelihood"):

		"""
		Pull out the scores for a subset of features

		:param feature_list: feature list to pull out from the database
		:type feature_list: list.

		:param score_type: name of the column that contains the particular score you are considering
		:type score_type: str.

		"""
		if not len(feature_list):
			raise ValueError("The feature_list is empty!")

		quoted = ",".join(["'{0}'".format(f) for f in feature_list])

		#Query the score database
		query = "SELECT {0},feature_type,{1} FROM '{2}' WHERE feature_type IN ({3})".format(",".join(self.parameters),score_type,table_name,quoted)
		logdriver.info("Executing SQL query: {0}".format(query))
		scores = self.query(query)

		#Pivot the database so that each feature has its own column,rename the columns
		l,scores = scores.suppress_indices(by=self.parameters,suppress=["feature_type"],columns=[score_type])
		rename = lambda s:l.query('suppress_group_id=={0}'.format(s[-1]))['feature_type'].iloc[0] if type(s)==tuple else s
		scores.columns = map(rename,scores.columns)

		return scores


###################################################################
#######Compute scores of a grid of parameter combinations##########
###################################################################

def chi2score(emulator,parameters,data,data_covariance,nchunks,pool):

	#Score the data on each of the parameter combinations provided
	scores = emulator.score(parameters,data,features_covariance=data_covariance,split_chunks=nchunks,pool=pool)

	#Pop the parameter columns, compute the likelihoods out of the chi2
	for p in parameters.columns:
		scores.pop(p)

	return scores,scores.apply(lambda c:np.exp(-0.5*c),axis=0)

@Parallelize.masterworker
def chi2database(db_name,parameters,specs,table_name="scores",pool=None,nchunks=None):

	"""
	Populate an SQL database with the scores of different parameter sets with respect to the data; supports multiple features

	:param db_name: name of the database to populate
	:type db_name: str.

	:param parameters: parameter combinations to score
 	:type parameters: :py:class:`Ensemble`

 	:param specs: dictionary that should contain the emulator,data, and covariance matrix of each feature to consider; each value in this dictionary must be a dictionary with keys 'emulator', 'data' and 'data covariance'
 	:type specs: dict.

 	:param table_name: table name to populate in the database
 	:type table_name: str.

 	:param pool: MPIPool to spread the calculations over (pass None for automatic pool handling)
 	:type pool: MPIPool

 	:param nchunks: number of chunks to split the parameter score calculations in (one chunk per processor ideally) 
 	:type nchunks: int.

	"""

	#Each processor should have the same exact workload
	if nchunks is not None:
		assert not len(parameters)%nchunks

	#Database context manager
	logdriver.info("Populating table '{0}' of score database {1}...".format(table_name,db_name))
	with ScoreDatabase(db_name) as db:

		#Repeat the scoring for each key in the specs dictionary
		for feature_type in specs.keys():

			#Log
			logdriver.info("Processing feature_type: {0} ({1} feature dimensions, {2} parameter combinations)...".format(feature_type,len(specs[feature_type]["data"]),len(parameters)))
			
			#Score
			chi2,likelihood = chi2score(emulator=specs[feature_type]["emulator"],parameters=parameters,data=specs[feature_type]["data"],data_covariance=specs[feature_type]["data_covariance"],nchunks=nchunks,pool=pool)
			assert (chi2.columns==[feature_type]).all()

			#Add to the database
			db_chunk = parameters.copy()
			db_chunk["feature_type"] = feature_type
			db_chunk["chi2"] = chi2
			db_chunk["likelihood"] = likelihood

			db.insert(db_chunk,table_name)