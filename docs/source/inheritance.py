from matplotlib import rc
import daft

rc("font", family="serif", size=12)
rc("text", usetex=False)

r_color = {"ec" : "red"}
g_color = {"ec" : "green"}

#Instantiate PGM
pgm = daft.PGM([17,7],origin=[0,0])

#Batch
pgm.add_node(daft.Node("batch","SimulationBatch",2,2.5,aspect=5.))

#Model
pgm.add_node(daft.Node("model","SimulationModel",5,2.5,aspect=5.,plot_params=r_color))
pgm.add_node(daft.Node("telescopicmaps","SimulationTelescopicMaps",5,3.5,aspect=6.,plot_params=g_color))

#Collection
pgm.add_node(daft.Node("collection","SimulationCollection",8,2.5,aspect=5.))

#IC
pgm.add_node(daft.Node("ic","SimulationIC",11,2.5,aspect=5.))

#Planes
pgm.add_node(daft.Node("planes","SimulationPlanes",14,2.5,aspect=5.))

#Maps and catalogs
pgm.add_node(daft.Node("maps","SimulationMaps",11,3.5,aspect=5.,plot_params=g_color))
pgm.add_node(daft.Node("catalog","SimulationCatalog",11,1.5,aspect=5.,plot_params=g_color))
pgm.add_node(daft.Node("subcatalog","SimulationSubCatalog",14,1.5,aspect=5.,plot_params=g_color))

#Plate
pgm.add_plate(daft.Plate([3.4,1.0,12,3.0]))


#Edges
pgm.add_edge("model","collection")
pgm.add_edge("model","telescopicmaps")
pgm.add_edge("collection","ic")
pgm.add_edge("ic","planes")
pgm.add_edge("collection","maps")
pgm.add_edge("collection","catalog")
pgm.add_edge("catalog","subcatalog")

#Render and save
pgm.render()
pgm.figure.savefig("figures/inheritance.png")

