try:
	
	from lenstools.pipeline import SimulationModel,SimulationCollection,SimulationIC
	from lenstools.pipeline.settings import EnvironmentSettings,NGenICSettings,PlaneSettings
	from lenstools.simulations import Nicaea,Gadget2Settings
	from lenstools import data

except ImportError:
	
	import sys
	sys.path.append("..")
	from lenstools.pipeline import SimulationModel,SimulationCollection,SimulationIC
	from lenstools.pipeline.settings import EnvironmentSettings,NGenICSettings,PlaneSettings
	from lenstools.simulations import Nicaea,Gadget2Settings
	from lenstools import data

home = "../SimTest/Home"
storage = "../SimTest/Storage"

cosmologies = [Nicaea(),Nicaea(Om0=0.4)]
box_sizes = [15.0,240.0]
part = [32,512]
seeds = [0,11,222]

#Check environment settings
env = EnvironmentSettings(home=home,storage=storage)

def test_directory_tree():

	#Create two simulation models
	for cosmo in cosmologies:
	
		simulation_model = SimulationModel(cosmology=cosmo,environment=env,parameters=["Om","Ol","w","si","ns"])

		#Create two different collections, three initial conditions per collection
		for i,box_size in enumerate(box_sizes):
			
			collection = simulation_model.newCollection(box_size=box_size*simulation_model.Mpc_over_h,nside=part[i])
			
			for seed in seeds:
				
				ic = collection.newRealization(seed)
				pln_settings = PlaneSettings()
				pln = ic.newPlaneSet(pln_settings)



def test_present():

	for model in SimulationModel.available(env):
		for collection in model.collections:
			for ic in collection.realizations:
				
				print(ic)
				pln = ic.getPlaneSet("Planes")
				print(pln)


def test_NGenICParam():

	#Create a parameter file for all the initial conditions present in the batch
	settings = NGenICSettings(GlassFile=data("dummy_glass_little_endian.dat"))

	for model in SimulationModel.available(env):
		for collection in model.collections:
			for ic in collection.realizations:
				ic.writeNGenIC(settings)


def test_Gadget2Param():

	#Create a parameter file for all the initial conditions present in the batch
	settings = Gadget2Settings()

	for model in SimulationModel.available(env):
		for collection in model.collections:
			for ic in collection.realizations:
				ic.writeGadget2(settings)
