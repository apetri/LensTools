from abc import ABCMeta,abstractproperty,abstractmethod
import os,glob
import cPickle

try:
	from paramiko import SSHClient
	SSHClient = SSHClient
except ImportError:
	SSHClient = None

###################################################
###########SystemHandler class#####################
###################################################

class SystemHandler(object):

	__metaclass__ = ABCMeta

	##################################
	######Abstract methods############
	##################################

	@abstractmethod
	def __init__(self):
		pass

	@abstractmethod
	def init(self,d):
		pass

	@abstractmethod
	def mkdir(self,d):
		pass

	@abstractmethod
	def isbatch(self,d):
		pass

	@abstractmethod
	def exists(self,d):
		pass

	@abstractmethod
	def glob(self,n):
		pass

	@abstractmethod
	def open(self,f,mode):
		pass

	@abstractmethod
	def pickleload(self,fp):
		pass

	@abstractmethod
	def pickledump(self,obj,fp):
		pass


############################################
#########Local filesystem###################
############################################

class LocalSystem(SystemHandler):

	"""
	Local system handler

	"""

	#############################################
	######Abstract method definitions############
	#############################################

	def __init__(self,readonly=False):
		
		self.name = "localhost"
		self.readonly = readonly

	def init(self,d):
		self.mkdir(d)

	def mkdir(self,d):

		if self.readonly:
			raise IOError("Simulation batch is read only!")

		os.mkdir(d)

	def isbatch(self,d):
		return self.exists(d)

	def exists(self,d):
		return os.path.exists(d)

	def glob(self,n):
		return glob.glob(n)

	def open(self,f,mode):

		if (self.readonly) and ("w" in mode or "a" in mode):
			raise IOError("Simulation batch is read only!")

		return open(f,mode)

	def pickleload(self,fp):
		return cPickle.load(fp)

	def pickledump(self,obj,fp):
		cPickle.dump(obj,fp)



##########################################################
#########Remote Unix filesystem via SSH###################
##########################################################

class UnixSSH(SystemHandler):

	"""
	Unix SSH client

	"""

	#############################################
	######Abstract method definitions############
	#############################################

	def __init__(self,client,readonly=True):

		if SSHClient is None:
			raise ImportError("paramiko needs to be installed to use remote SSH functionality!")

		assert isinstance(client,SSHClient)
		self.client = client
		self.readonly = readonly

		#Assign IP address and port as name
		host_ip,port = self.client.get_transport().getpeername()
		print("[+] System handler will read/write from {0}:{1}".format(host_ip,port))
		self.name = "{0}:{1}".format(host_ip,port)

		#Open SFTP session
		self.sftp = self.client.open_sftp()

	def init(self,d):
		self.mkdir(d)

	def mkdir(self,d):
		stdin,stdout,stderr = self.client.exec_command("mkdir {0}".format(d))

	def isbatch(self,d):
		return self.exists(d)

	def exists(self,d):

		stdin,stdout,stderr = self.client.exec_command("ls -d {0}".format(d))
		result = stdout.read().strip("\n")

		if result=="":
			return False
		else:
			return True

	def glob(self,n):

		stdin,stdout,stderr = self.client.exec_command("ls -d {0}".format(n))
		return [ d.rstrip("\n").rstrip(":") for d in stdout.readlines("\n") ]

	def open(self,f,mode):

		if (self.readonly) and ("w" in mode or "a" in mode):
			raise IOError("Simulation batch is read only!") 

		return self.sftp.file(f,mode)

	def pickleload(self,fp):
		return cPickle.loads(fp.read())

	def pickledump(self,obj,fp):
		fp.write(cPickle.dumps(obj))
		

