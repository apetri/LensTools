"""

.. module:: database
	:platform: Unix
	:synopsis: This module is a convenient wrapper of utilities for I/O between SQL databases and Ensembles


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division,print_function,with_statement
from .ensemble import Ensemble

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

	#Create a connection to a database
	def __init__(self,name):

		if sqlalchemy is None:
			raise ImportError("sqlalchemy is not installed!!")

		self._constructor_ensemble = Ensemble
		self.connection = sqlalchemy.create_engine("sqlite:///"+name)


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
		:type df: Ensemble

		"""

		assert isinstance(df,pd.DataFrame)
		df.to_sql(table_name,self.connection,if_exists="append",index=False)

	#Query the database
	def query(self,sql):

		"""
		:param sql: sql query string
		:type sql: str.

		:returns: Ensemble

		"""

		return self._constructor_ensemble.read_sql_query(sql,self.connection)

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

		if type(table_name)==str:
			assert table_name in self.tables,"Table {0} does not exist!".format(table_name)
			return self._constructor_ensemble.read_sql_table(table_name,self.connection) 
		elif type(table_name)==int:
			return self._constructor_ensemble.read_sql_table(self.tables[table_name],self.connection)
		else:
			raise TypeError("table_name type not recognized")


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
		query = "SELECT {0},feature_type,{1} FROM {2} WHERE feature_type IN ({3})".format(",".join(self.parameters),score_type,table_name,quoted)
		print("[+] Executing SQL query: {0}".format(query))
		scores = self.query(query)

		#Pivot the database so that each feature has its own column,rename the columns
		l,scores = scores.suppress_indices(by=self.parameters,suppress=["feature_type"],columns=[score_type])
		rename = lambda s:l.query('suppress_group_id=={0}'.format(s[-1]))['feature_type'].iloc[0] if type(s)==tuple else s
		scores.columns = map(rename,scores.columns)

		return scores