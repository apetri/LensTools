import sys
import logging

formatter = logging.Formatter("%(asctime)s.%(msecs)d:%(name)-12s:%(levelname)-4s: %(message)s",datefmt='%m-%d %H:%M:%S')

console = logging.StreamHandler(sys.stdout)
console_error = logging.StreamHandler(sys.stderr)

console.setFormatter(formatter)
console_error.setFormatter(formatter)

logpreamble = logging.getLogger("lenstools.preamble")
logdriver = logging.getLogger("lenstools.driver")
logplanes = logging.getLogger("lenstools.planes")
logray = logging.getLogger("lenstools.raytracing")
logstderr = logging.getLogger("lenstools.stderr")

for logger in [logpreamble,logdriver,logplanes,logray]:
	logger.addHandler(console)
	logger.propagate = False

logstderr.addHandler(console_error)
logstderr.propagate = False