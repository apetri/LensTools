#!/bin/bash

import re
import pandas as pd

plane_checkpoints = ("Snapshot Input","Gridding","MPI Communication","FFT operations","Plane Output")
raytracing_checkpoints = ("Lens input","Roll plane","Compute deflections","Compute shear matrix","Map output")

def date(line):
	datestr = re.match(r"([0-9]+-[0-9]+ [0-9]+:[0-9]+:[0-9]+.[0-9]+)",line).group()
	return pd.to_datetime(datestr,format="%m-%d %H:%M:%S")

def parse_plane_log(s):
	pass

def parse_raytracing_log(s):
	pass