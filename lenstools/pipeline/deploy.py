import os
import StringIO
from abc import ABCMeta,abstractproperty,abstractmethod

import astropy.units as u

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

	@property
	def directives(self):
		return self._directives

	@property
	def cluster_specs(self):
		return self._cluster_specs


	def writePreamble(self,settings):

		"""
		Writes the preamble of the job script (resources request,job name, etc...)

		:param settings: job settings
		:type settings: JobSettings

		:returns: StringIO object

		"""

		#Type safety check
		assert isinstance(settings,JobSettings)

		#Write the preamble
		s = StringIO.StringIO()

		#Shell type
		s.write("{0}\n".format(self.cluster_specs.shell_prefix))
		
		#Write allocation-ID, if any
		if settings.charge_account is not None:
			s.write("""
################################
######Allocation ID#############
################################

{0} {1} {2}
""".format(self.directives.directive_prefix,self.directives.charge_account_switch,settings.charge_account))

		#Write the rest of the preamble (except the resources allocations)
		s.write("""

##########################################
#############Directives###################
##########################################

{0} {1} {2}

{0} {3} {4}
{0} {5} {6}

""".format(self.directives.directive_prefix,self.directives.job_name_switch,settings.job_name,self.directives.stdout_switch,settings.redirect_stdout,self.directives.stderr_switch,settings.redirect_stderr))

		s.write("""
{0} {1} {2}
{0} {3} {4}

{0} {5}={6}
{0} {7}

""".format(self.directives.directive_prefix,self.directives.queue_type_switch,settings.queue,self.directives.wallclock_time_switch,settings.wallclock_time,self.directives.user_email_switch,settings.email,self.directives.user_email_type))

		#Write the resources requests
		s.write("""
##########################################
#############Resources####################
##########################################

{0} {1}{2}
""".format(self.directives.directive_prefix,self.directives.num_cores_switch,settings.num_cores))

		if self.directives.num_nodes_switch is not None:
			s.write("{0} {1}{2}\n".format(self.directives.directive_prefix,self.directives.num_nodes_switch,settings.num_nodes))

		if self.directives.tasks_per_node_switch is not None:
			s.write("{0} {1}{2}\n\n\n".format(self.directives.directive_prefix,self.directives.tasks_per_node_switch,settings.tasks_per_node))


		#Done, return to user
		s.seek(0)
		return s.read()


	def writeExecution(self,executables,cores,settings):

		"""
		Write the execution part of the script

		"""

		assert len(executables)==len(cores)

		s = StringIO.StringIO()
		s.write("""
###################################################
#################Execution#########################
###################################################

""")

		#TODO: handle the splitting between multiple executables
		for executable in executables:
			s.write("{0} {1}\n".format(self.cluster_specs.job_starter,executable))

		#Done, return to user
		s.seek(0)
		return s.read()



########################################################
###########StampedeHandler class########################
########################################################

SLURMspecs = {
"directive_prefix" : "#SBATCH",
"charge_account_switch" : "-A",
"job_name_switch" : "-J",
"stdout_switch" : "-o",
"stderr_switch" : "-e",
"num_cores_switch" : "-n ",
"num_nodes_switch" : "-N ",
"tasks_per_node_switch" : None,
"queue_type_switch" : "-p",
"wallclock_time_switch" : "-t",
"user_email_switch" : "--mail-user=",
"user_email_type" : "--mail-type=all",
}

StampedeClusterSpecs = {
"shell_prefix" : "#!/bin/bash",
"execution_preamble" : None,
"job_starter" : "ibrun",
"cores_per_node" : 16,
"memory_per_node" : 32.0*u.Gbyte,
"cores_at_execution_switch" : "-n",
"offset_switch" : "-o",
"wait_switch" : "wait",
"multiple_executables_on_node" : True
} 

class StampedeHandler(JobHandler):

	"""
	Handler class for SLURM jobs

	"""

	#############################################
	######Abstract method definitions############
	#############################################

	def setDirectives(self):
		self._directives = Directives(**SLURMspecs)

	def setClusterSpecs(self):
		self._cluster_specs = ClusterSpecs(**StampedeClusterSpecs)

########################################################################################################

##########################################
########Directives class##################
##########################################

class Directives(object):

	def __init__(self,**kwargs):

		for key in kwargs.keys():
			setattr(self,key,kwargs[key])


############################################
########ClusterSpecs class##################
############################################


class ClusterSpecs(object):

	def __init__(self,**kwargs):

		for key in kwargs.keys():
			setattr(self,key,kwargs[key])



