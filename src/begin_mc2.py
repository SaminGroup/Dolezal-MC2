import numpy as np

step_limit = 2500 # max steps to run for
T = 1873.15 # simulation temperature
param = "begin" # should be 'begin' or 'continue'
supcomp_phrase = "vasp" # for Mustang runs
intra = False # True, perform intraswaps, False, do not
w = 0.60 # energetic penalty cost for acceptance criterion
C = np.loadtxt('concentration')

from src.mc2_global import mc2

mc2(step_limit, C, T, w, param, intra, supcomp_phrase)
