from matplotlib import rc
import daft


def flow(ftype="png"):

	rc("font", family="serif", size=12)
	rc("text", usetex=False)

	r_color = {"ec" : "red"}
	g_color = {"ec" : "green"}

	#Instantiate PGM
	pgm = daft.PGM([16,7],origin=[0,0])

	#Nodes
	pgm.add_node(daft.Node("parameters","Parameters",1,2.5,aspect=3.,plot_params=r_color))
	pgm.add_node(daft.Node("geometry","geometry+seeds",4,2.5,aspect=4.))
	
	#ICS
	pgm.add_node(daft.Node("ic1","IC seed 1",6.5,3.25,aspect=2.5))
	pgm.add_node(daft.Node("ic2","IC seed 2",6.5,2.5,aspect=2.5))
	pgm.add_node(daft.Node("icN",r"IC seed $N$",6.5,1.75,aspect=2.5))

	#Evolution
	pgm.add_node(daft.Node("e1",r"$\delta^1(\mathbf{x},z)$",8.5,3.25,aspect=3.))
	pgm.add_node(daft.Node("e2",r"$\delta^2(\mathbf{x},z)$",8.5,2.5,aspect=3.))
	pgm.add_node(daft.Node("eN",r"$\delta^N(\mathbf{x},z)$",8.5,1.75,aspect=3.))

	#Lens Planes
	pgm.add_node(daft.Node("p1",r"Lens: $\sigma^1(\mathbf{x},z)$",10.5,3.25,aspect=3.5))
	pgm.add_node(daft.Node("p2",r"Lens: $\sigma^2(\mathbf{x},z)$",10.5,2.5,aspect=3.5))
	pgm.add_node(daft.Node("pN",r"Lens: $\sigma^N(\mathbf{x},z)$",10.5,1.75,aspect=3.5))

	#Mix planes
	pgm.add_plate(daft.Plate([9.4,1.0,2.1,3.0],label="Mix seeds"))

	#Lensing maps
	pgm.add_node(daft.Node("lens","Lensing maps " + r"$(\kappa,\gamma)$",13.0,2.5,aspect=4.5,plot_params=g_color))

	#Executables
	pgm.add_node(daft.Node("camb","CAMB+NGen-IC",4,0.5,aspect=4.5,observed=True))
	pgm.add_node(daft.Node("gadget","Gadget2",6.5,0.5,aspect=2.,observed=True))
	pgm.add_node(daft.Node("planes","lenstools.planes",8.5,0.5,aspect=4.,observed=True))
	pgm.add_node(daft.Node("ray","lenstools.raytracing",10.5,4.5,aspect=5.,observed=True))

	#Edges
	pgm.add_edge("parameters","geometry")
	pgm.add_edge("geometry","ic1")
	pgm.add_edge("geometry","ic2")
	pgm.add_edge("geometry","icN")
	pgm.add_edge("ic1","e1")
	pgm.add_edge("ic2","e2")
	pgm.add_edge("icN","eN")
	pgm.add_edge("e1","p1")
	pgm.add_edge("e2","p2")
	pgm.add_edge("eN","pN")
	pgm.add_edge("p2",'lens')
	pgm.add_edge("camb","geometry")
	pgm.add_edge("gadget","icN")
	pgm.add_edge("planes","eN")
	pgm.add_edge("ray","p1")

	#Render and save
	pgm.render()
	pgm.figure.savefig("figures/flow."+ftype)


if __name__=="__main__":
	flow()

