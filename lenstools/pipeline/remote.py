import os,glob

###################################################
###########SystemHandler class#####################
###################################################

class SystemHandler(object):

	"""
	SystemHandler class
	
	"""


#############################################
#########Local filesystem ###################
#############################################

class LocalSystem(SystemHandler):

	def __init__(self):
		
		self.mkdir = os.mkdir
		self.exists = os.path.exists
		self.glob = glob.glob
		self.open = open