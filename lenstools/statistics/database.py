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

		self.connection = sqlalchemy.create_engine("sqlite:///"+name)

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

		return Ensemble.read_sql_query(sql,self.connection)

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
		assert table_name in self.tables,"Table {0} does not exist!".format(table_name)
		return Ensemble.read_sql_table(table_name,self.connection) 