try:
	
	from lenstools.pipeline import EnvironmentSettings,SimulationModel,SimulationCollection,SimulationIC
	from lenstools.simulations import Nicaea

except ImportError:
	
	import sys
	sys.path.append("..")
	from lenstools.pipeline import EnvironmentSettings,SimulationModel,SimulationCollection,SimulationIC
	from lenstools.simulations import Nicaea

home = "../SimTest/Home"
storage = "../SimTest/Storage"

cosmologies = [Nicaea(),Nicaea(wa=0.5)]
box_sizes = [15.0,240.0]
part = [32,512]
seeds = [0,11,222]

#Check environment settings
env = EnvironmentSettings(home=home,storage=storage)

def test_directoryTree():

	#Create two simulation models
	for cosmo in cosmologies:
	
		simulation_model = SimulationModel(cosmology=cosmo,environment=env,parameters=["Om","Ol","w","wa"])

		#Create two different collections, three initial conditions per collection
		for i,box_size in enumerate(box_sizes):
			
			collection = simulation_model.newCollection(box_size=box_size*simulation_model.Mpc_over_h,nside=part[i])
			
			for seed in seeds:
				ic = collection.newInitialCondition(seed)