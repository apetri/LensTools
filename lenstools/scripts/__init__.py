"""
Command line scripts that come with lenstools

"""

##############################################
#######Line of sight integration##############
##############################################

integration_types = {

"born" : "Born approximation for the convergence (linear order in the lensing potential)",
"born-rt" : "Semi-Born approximation: integrate the density on the real ray trajectory",
"postBorn2" : "Convergence at second order in the lensing potential (lens-lens + geodesic perturbation)",
"postBorn2-ll" : "Convergence due to lens-lens coupling", 
"postBorn2-gp" : "Convergence due to geodesic perturbation",
"postBorn1+2" : "Convergence at Born + second post Born order in the lensing potential",
"postBorn1+2-ll" : "Convergence at Born + second post Born order (lens-lens only)",
"postBorn1+2-gp" : "Convergence at Born + second post Born order (geodesic perturbation only)",
"omega2" : "Rotation at second order in the lensing potential"

}