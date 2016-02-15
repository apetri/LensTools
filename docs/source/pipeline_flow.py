from matplotlib import rc
import daft

rc("font", family="serif", size=12)
rc("text", usetex=False)

r_color = {"ec" : "red"}
g_color = {"ec" : "green"}

#Instantiate PGM
pgm = daft.PGM([16,7],origin=[0,0])

#Nodes
pgm.add_node(daft.Node("parameters","Parameters",1,2.5,aspect=3.,plot_params=r_color))
pgm.add_node(daft.Node("geometry","geometry",3,2.5,aspect=3.))

#Seeds
pgm.add_node(daft.Node("seed1","seed 1",5,3.25,aspect=2.))
pgm.add_node(daft.Node("seed2","seed 2",5,2.5,aspect=2.))
pgm.add_node(daft.Node("seedN",r"seed $N$",5,1.75,aspect=2.))

#ICS
pgm.add_node(daft.Node("ic1","IC 1",6.5,3.25,aspect=2.))
pgm.add_node(daft.Node("ic2","IC 2",6.5,2.5,aspect=2.))
pgm.add_node(daft.Node("icN","IC N",6.5,1.75,aspect=2.))

#Evolution
pgm.add_node(daft.Node("e1",r"$\rho^1(\mathbf{x},z)$",8.5,3.25,aspect=3.))
pgm.add_node(daft.Node("e2",r"$\rho^2(\mathbf{x},z)$",8.5,2.5,aspect=3.))
pgm.add_node(daft.Node("eN",r"$\rho^N(\mathbf{x},z)$",8.5,1.75,aspect=3.))

#Lens Planes
pgm.add_node(daft.Node("p1","Lens Planes 1",10.5,3.25,aspect=3.5))
pgm.add_node(daft.Node("p2","Lens Planes 2",10.5,2.5,aspect=3.5))
pgm.add_node(daft.Node("pN","Lens Planes N",10.5,1.75,aspect=3.5))

#Mix planes
pgm.add_plate(daft.Plate([9.4,1.0,2.1,3.0],label="Mix realizations"))

#Lensing maps
pgm.add_node(daft.Node("lens","Lensing maps",13.0,2.5,aspect=3.5,plot_params=g_color))

#Executables
pgm.add_node(daft.Node("camb","CAMB",3,0.5,aspect=2.,observed=True))
pgm.add_node(daft.Node("ngenic","NGen-IC",5,0.5,aspect=2.,observed=True))
pgm.add_node(daft.Node("gadget","Gadget2",6.5,0.5,aspect=2.,observed=True))
pgm.add_node(daft.Node("planes","lenstools.planes",8.5,0.5,aspect=4.,observed=True))
pgm.add_node(daft.Node("ray","lenstools.raytracing",10.5,4.5,aspect=5.,observed=True))

#Edges
pgm.add_edge("parameters","geometry")
pgm.add_edge("geometry","seed1")
pgm.add_edge("geometry","seed2")
pgm.add_edge("geometry","seedN")
pgm.add_edge("seed1","ic1")
pgm.add_edge("seed2","ic2")
pgm.add_edge("seedN","icN")
pgm.add_edge("ic1","e1")
pgm.add_edge("ic2","e2")
pgm.add_edge("icN","eN")
pgm.add_edge("e1","p1")
pgm.add_edge("e2","p2")
pgm.add_edge("eN","pN")
pgm.add_edge("p2",'lens')
pgm.add_edge("camb","geometry")
pgm.add_edge("ngenic","seedN")
pgm.add_edge("gadget","icN")
pgm.add_edge("planes","eN")
pgm.add_edge("ray","p1")

#Render and save
pgm.render()
pgm.figure.savefig("figures/flow.png")

