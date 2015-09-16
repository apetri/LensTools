import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lenstools.statistics.ensemble import Ensemble,Series
from lenstools.statistics.constraints import Emulator
from lenstools.statistics.contours import ContourPlot


p = np.outer(np.arange(10.0),np.ones(2))
f = np.outer(np.arange(10.0)+0.1*2,np.ones(5))
emulator = Emulator.from_features(f,p,parameter_index=[r"$\alpha$",r"$\beta$"],feature_index=[r"$f_{0}$".format(n) for n in range(5)])

test_feature = Series(np.ones(5)*4.0 + 0.05*np.random.randn(5),index=emulator[["features"]].columns)

features_covariance = Ensemble(np.eye(5)*0.5,index=test_feature.index,columns=test_feature.index)

emulator.train()

#Evaluate the score of each of these parameter combinations on the test feature
g = np.arange(0.0,10.0,0.5)
p = np.array(np.meshgrid(g,g,indexing="ij")).reshape(2,400).T
test_parameters = Ensemble(p,columns=emulator.parameter_names)

scores = emulator.score(parameters=test_parameters,observed_feature=test_feature,features_covariance=features_covariance)
scores['features'] = np.exp(-0.5*scores['features'])

contour = ContourPlot.from_scores(scores,parameters=emulator.parameter_names)
contour.show()
contour.getLikelihoodValues([0.684],precision=0.1)
contour.plotContours(colors=["red"])
contour.labels()

contour.savefig("contour_example.png")