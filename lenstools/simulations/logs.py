import sys,platform
import resource
import logging

import numpy as np
import astropy.units as u

#########
#Logging#
#########

formatter = logging.Formatter("%(asctime)s.%(msecs)d:%(name)-12s:%(levelname)-4s: %(message)s",datefmt='%m-%d %H:%M:%S')

console = logging.StreamHandler(sys.stdout)
console_error = logging.StreamHandler(sys.stderr)

console.setFormatter(formatter)
console_error.setFormatter(formatter)

logpreamble = logging.getLogger("lenstools.preamble")
logdriver = logging.getLogger("lenstools.driver")
logplanes = logging.getLogger("lenstools.planes")
logray = logging.getLogger("lenstools.raytracing")
logcmb = logging.getLogger("lenstools.cmb")
logstderr = logging.getLogger("lenstools.stderr")

for logger in [logpreamble,logdriver,logplanes,logray,logcmb]:
	logger.addHandler(console)
	logger.propagate = False

logstderr.addHandler(console_error)
logstderr.propagate = False

#######################
#Peak memory usage log#
#######################

#Multiplicative factor to convert resource output into GB
ostype = platform.system()
if ostype in ["Darwin","darwin"]:
	to_gbyte = 1024.**3
elif ostype in ["Linux","linux"]:
	to_gbyte = 1024.**2
else:
	to_gbyte = np.nan

#Get the peak memory usage for a single task
def peakMemory():
	return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*u.Gbyte/to_gbyte

#Get the peak memory usage for all tasks
def peakMemoryAll(pool):
	memory_task = peakMemory()
	if pool is None:
		return memory_task,1
	memory_task_raw,unit = np.array([memory_task.value]),memory_task.unit
	memory_task_all = np.zeros(1)
	pool.comm.Reduce(memory_task_raw,memory_task_all)

	return memory_task_all[0]*unit,pool.size+1
