from .remote import LocalSystem
from ..simulations import Nicaea
from ..simulations.gadget2 import Gadget2SnapshotDE

###################
#Default cosmology#
###################

class LensToolsCosmology(Nicaea):

	"""
	LensTools pipeline cosmology handler

	"""

CosmoDefault = LensToolsCosmology() 

##############################################
#Parametere name to attribute name dictionary#
##############################################

name2attr = dict()
name2attr["Om"] = "Om0"
name2attr["Ol"] = "Ode0"
name2attr["w"] = "w0"
name2attr["wa"] = "wa"
name2attr["h"] = "h"
name2attr["Ob"] = "Ob0"
name2attr["si"] = "sigma8"
name2attr["ns"] = "ns"

############################################
#Number of digits of precision for cosmo_id#
############################################

cosmo_id_digits = 3

########################
#Default system handler#
########################

syshandler = LocalSystem()

###################################
#Default simulation tree json file#
###################################

json_tree_file = ".tree.json"

##########################
#Default snapshot handler#
##########################

snapshot_handler = Gadget2SnapshotDE