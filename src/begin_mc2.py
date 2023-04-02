
step_limit = 2500 # max steps to run for
Ncells = 3 # should match number of species
T = 1873.15 # simulation temperature
param = "begin" # should be 'begin' or 'continue'
supcomp_phrase = "vasp" # for Mustang runs
intra = False # True, perform intraswaps, False, do not

from src.mc2_global import mc2

mc2(step_limit, Ncells, T, param, intra, supcomp_phrase)
