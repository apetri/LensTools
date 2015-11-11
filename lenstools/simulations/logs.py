import sys
import logging

console = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("%(asctime)s.%(msecs)d:%(name)-12s:%(levelname)-4s: %(message)s",datefmt='%m-%d %H:%M:%S')
console.setFormatter(formatter)

logpreamble = logging.getLogger("lenstools.preamble")
logdriver = logging.getLogger("lenstools.driver")
logplanes = logging.getLogger("lenstools.planes")
logray = logging.getLogger("lenstools.raytracing")

for logger in [logpreamble,logdriver,logplanes,logray]:
	logger.addHandler(console)
	logger.propagate = False