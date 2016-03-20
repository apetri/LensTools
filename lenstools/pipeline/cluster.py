from .deploy import JobHandler,Directives,ClusterSpecs

import astropy.units as u

########################################################
###########Job scheduler specs##########################
########################################################

_SLURMspecs = {
"directive_prefix" : "#SBATCH",
"job_name_switch" : "-J ",
"stdout_switch" : "-o ",
"stderr_switch" : "-e ",
"num_cores_switch" : "-n ",
"num_nodes_switch" : "-N ",
"tasks_per_node_switch" : None,
"queue_type_switch" : "-p ",
"wallclock_time_switch" : "-t ",
"user_email_switch" : "--mail-user=",
"user_email_type" : "--mail-type=all",
}

_PBSspecs = {
"directive_prefix" : "#PBS",
"job_name_switch" : "-N ",
"stdout_switch" : "-o ",
"stderr_switch" : "-e ",
"num_cores_switch" : "-l mppwidth=",
"num_nodes_switch" : None,
"tasks_per_node_switch" : "-lmppnppn=",
"queue_type_switch" : "-q ",
"wallclock_time_switch" : "-l walltime=",
"user_email_switch" : "-M ",
"user_email_type" : "-m abe",
}

########################################################
###########Cluster specs################################
########################################################

_StampedeClusterSpecs = {
"shell_prefix" : "#!/bin/bash",
"execution_preamble" : None,
"charge_account_switch" : "-A ",
"job_starter" : "ibrun",
"cores_per_node" : 16,
"memory_per_node" : 32.0*u.Gbyte,
"cores_at_execution_switch" : "-n ",
"offset_switch" : "-o ",
"wait_switch" : "wait",
"multiple_executables_on_node" : True
}

_EdisonClusterSpecs = {
"shell_prefix" : "#!/bin/bash -l",
"execution_preamble" : None,
"charge_account_switch" : None,
"job_starter" : "srun",
"cores_per_node" : 24,
"memory_per_node" : 64.0*u.Gbyte,
"cores_at_execution_switch" : "-n ",
"offset_switch" : "-N ",
"wait_switch" : "wait",
"multiple_executables_on_node" : False	
} 

_CoriClusterSpecs = {
"shell_prefix" : "#!/bin/bash -l",
"execution_preamble" : None,
"charge_account_switch" : None,
"job_starter" : "srun",
"cores_per_node" : 32,
"memory_per_node" : 128.0*u.Gbyte,
"cores_at_execution_switch" : "-n ",
"offset_switch" : "-N ",
"wait_switch" : "wait",
"multiple_executables_on_node" : False	
} 

########################################################
###########StampedeHandler class########################
########################################################

class StampedeHandler(JobHandler):

	"""
	Job handler for the TACC Stampede cluster

	"""

	#############################################
	######Abstract method definitions############
	#############################################

	def setDirectives(self):
		self._directives = Directives(**_SLURMspecs)

	def setClusterSpecs(self):
		self._cluster_specs = ClusterSpecs(**_StampedeClusterSpecs)


########################################################
###########EdisonHandler class##########################
########################################################

class EdisonHandler(JobHandler):

	"""
	Job handler for the NERSC Edison cluster

	"""

	#############################################
	######Abstract method definitions############
	#############################################

	def setDirectives(self):
		self._directives = Directives(**_SLURMspecs)

	def setClusterSpecs(self):
		self._cluster_specs = ClusterSpecs(**_EdisonClusterSpecs)


######################################################
###########CoriHandler class##########################
######################################################

class CoriHandler(JobHandler):

	"""
	Job handler for the NERSC Cori Phase 1 cluster

	"""

	#############################################
	######Abstract method definitions############
	#############################################

	def setDirectives(self):
		self._directives = Directives(**_SLURMspecs)

	def setClusterSpecs(self):
		self._cluster_specs = ClusterSpecs(**_CoriClusterSpecs)	