from matplotlib import rc
import daft

rc("font", family="serif", size=12)
rc("text", usetex=False)

r_color = {"ec" : "red"}
g_color = {"ec" : "green"}

#Instantiate PGM
pgm = daft.PGM([13,7],origin=[0,0])

#Batch
pgm.add_node(daft.Node("batch","SimulationBatch",2,3.5,aspect=5.))

#Model
pgm.add_node(daft.Node("model","SimulationModel",2,2.5,aspect=5.,plot_params=r_color))

#Collection
pgm.add_node(daft.Node("collection","SimulationCollection",5,2.5,aspect=5.))

#IC
pgm.add_node(daft.Node("ic","SimulationIC",8,2.5,aspect=5.))

#Planes
pgm.add_node(daft.Node("planes","SimulationPlanes",11,2.5,aspect=5.))

#Maps and catalogs
pgm.add_node(daft.Node("maps","SimulationMaps",8,3.5,aspect=5.,plot_params=g_color))
pgm.add_node(daft.Node("catalog","SimulationCatalog",8,1.5,aspect=5.,plot_params=g_color))


#Edges
pgm.add_edge("model","collection")
pgm.add_edge("collection","ic")
pgm.add_edge("ic","planes")
pgm.add_edge("collection","maps")
pgm.add_edge("collection","catalog")

#Render and save
pgm.render()
pgm.figure.savefig("inheritance.png")

