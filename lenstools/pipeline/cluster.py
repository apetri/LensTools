from .deploy import JobHandler,Directives,ClusterSpecs

import astropy.units as u

########################################################
###########Job scheduler specs##########################
########################################################

_SLURMspecs = {
"directive_prefix" : "#SBATCH",
"charge_account_switch" : "-A ",
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

########################################################
###########Cluster specs################################
########################################################

_StampedeClusterSpecs = {
"shell_prefix" : "#!/bin/bash",
"execution_preamble" : None,
"job_starter" : "ibrun",
"cores_per_node" : 16,
"memory_per_node" : 32.0*u.Gbyte,
"cores_at_execution_switch" : "-n ",
"offset_switch" : "-o ",
"wait_switch" : "wait",
"multiple_executables_on_node" : True
} 

########################################################
###########StampedeHandler class########################
########################################################

class StampedeHandler(JobHandler):

	"""
	Handler class for SLURM jobs

	"""

	#############################################
	######Abstract method definitions############
	#############################################

	def setDirectives(self):
		self._directives = Directives(**_SLURMspecs)

	def setClusterSpecs(self):
		self._cluster_specs = ClusterSpecs(**_StampedeClusterSpecs)

########################################################################################################