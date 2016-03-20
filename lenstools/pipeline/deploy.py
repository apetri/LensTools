from __future__ import division

import StringIO
from operator import add
from functools import reduce

from abc import ABCMeta,abstractproperty,abstractmethod

from .settings import JobSettings

############################################################
###########JobHandler abstract class########################
############################################################

class JobHandler(object):

	__metaclass__ = ABCMeta

	##################################
	######Abstract methods############
	##################################

	@abstractmethod
	def setDirectives(self,directives):
		pass

	@abstractmethod
	def setClusterSpecs(self,cluster_specs):
		pass


	########################################
	#####Default non--abstract methods######
	########################################

	def __init__(self):
		
		self.setDirectives()
		self.setClusterSpecs()

	def __repr__(self):
		directives = "\n".join([ "{0} : {1}".format(k,getattr(self.directives,k)) for k in self.directives._metadata])
		cluster = "\n".join([ "{0} : {1}".format(k,getattr(self.cluster_specs,k)) for k in self.cluster_specs._metadata])
		return "Directives:\n" + directives + "\n\nCluster specifications:\n" + cluster

	@property
	def directives(self):
		return self._directives

	@property
	def cluster_specs(self):
		return self._cluster_specs


	def writePreamble(self,settings,auto_num_nodes=True):

		"""
		Writes the preamble of the job script (resources request,job name, etc...)

		:param settings: job settings
		:type settings: JobSettings

		:param auto_num_nodes: if True, the number of requested nodes is computed automatically from the number of requested cores (knowing the cluster specifications)
		:type auto_num_nodes: bool.

		:returns: StringIO object

		"""

		#Type safety check
		assert isinstance(settings,JobSettings)

		#Write the preamble
		s = StringIO.StringIO()

		#Shell type
		s.write("{0}\n".format(self.cluster_specs.shell_prefix))
		
		#Write allocation-ID, if any
		if self.cluster_specs.charge_account_switch is not None:
			s.write("""
################################
######Allocation ID#############
################################

{0} {1}{2}
""".format(self.directives.directive_prefix,self.cluster_specs.charge_account_switch,settings.charge_account))

		#Write the rest of the preamble (except the resources allocations)
		s.write("""

##########################################
#############Directives###################
##########################################

{0} {1}{2}

{0} {3}{4}
{0} {5}{6}

""".format(self.directives.directive_prefix,self.directives.job_name_switch,settings.job_name,self.directives.stdout_switch,settings.redirect_stdout,self.directives.stderr_switch,settings.redirect_stderr))

		s.write("""
{0} {1}{2}
{0} {3}{4}

{0} {5}{6}
{0} {7}

""".format(self.directives.directive_prefix,self.directives.queue_type_switch,settings.queue,self.directives.wallclock_time_switch,settings.wallclock_time,self.directives.user_email_switch,settings.email,self.directives.user_email_type))

		#Write the resources requests
		s.write("""
##########################################
#############Resources####################
##########################################

{0} {1}{2}
""".format(self.directives.directive_prefix,self.directives.num_cores_switch,settings.num_cores))

		if auto_num_nodes:
			num_nodes = settings.num_cores//self.cluster_specs.cores_per_node
			if settings.num_cores%self.cluster_specs.cores_per_node:
				num_nodes+=1
		else:
			num_nodes = settings.num_nodes

		if self.directives.num_nodes_switch is not None:
			s.write("{0} {1}{2}\n".format(self.directives.directive_prefix,self.directives.num_nodes_switch,num_nodes))

		if (self.directives.tasks_per_node_switch is not None) and (hasattr(settings,"tasks_per_node")):
			s.write("{0} {1}{2}\n\n\n".format(self.directives.directive_prefix,self.directives.tasks_per_node_switch,settings.tasks_per_node))


		#Done, return to user
		s.seek(0)
		return s.read()


	def writeExecution(self,executables,cores,settings):

		"""
		Write the execution part of the script

		:param executables: list of executables to run on the compute nodes
		:type executables: list.

		:param cores: list of numbers of cores for each executable (must have the same length as executables)
		:type cores: list.

		:param settings: job settings
		:type settings: JobSettings

		:returns: StringIO object

		"""

		#Type safety check
		assert isinstance(settings,JobSettings)
		assert len(executables)==len(cores),"You must specify the number of cores for each executable!"

		#Check that the sum of the cores requested matches the job settings
		assert reduce(add,cores)==settings.num_cores,"The number of cores requested does not match the execution statement!"

		s = StringIO.StringIO()
		s.write("""
###################################################
#################Execution#########################
###################################################

""")

		if self.cluster_specs.execution_preamble is not None:
			s.write("{0}\n\n".format(self.cluster_specs.execution_preamble))

		offset = 0
		for n,executable in enumerate(executables):

			if self.cluster_specs.multiple_executables_on_node:
				s.write("{0} {1}{2} {3}{4} {5} &\n".format(self.cluster_specs.job_starter,self.cluster_specs.cores_at_execution_switch,cores[n],self.cluster_specs.offset_switch,offset,executable))
			else:
				
				if self.cluster_specs.offset_switch is not None:
					nodes = cores[n]//self.cluster_specs.cores_per_node + (cores[n]%self.cluster_specs.cores_per_node>0)
					s.write("{0} {1}{2} {3}{4} {5} &\n".format(self.cluster_specs.job_starter,self.cluster_specs.cores_at_execution_switch,cores[n],self.cluster_specs.offset_switch,nodes,executable))
				else:
					s.write("{0} {1}{2} {3} &\n".format(self.cluster_specs.job_starter,self.cluster_specs.cores_at_execution_switch,cores[n],executable))

			#Increase offset
			offset += cores[n]

		#wait statement
		s.write("{0}\n".format(self.cluster_specs.wait_switch))

		#Done, return to user
		s.seek(0)
		return s.read()


##########################################
########Directives class##################
##########################################

class Directives(object):

	def __init__(self,**kwargs):

		self._metadata = list()

		for key in kwargs:
			setattr(self,key,kwargs[key])
			self._metadata.append(key)


############################################
########ClusterSpecs class##################
############################################


class ClusterSpecs(object):

	def __init__(self,**kwargs):

		self._metadata = list()

		for key in kwargs:
			setattr(self,key,kwargs[key])
			self._metadata.append(key)



