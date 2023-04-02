
step_limit = 2500
Ncells = 2
T = 1750
param = 'begin'
supcomp_phrase = 'vasp'
intra = True

from src.mc2_global import mc2

mc2(step_limit, Ncells, T, param, intra, supcomp_phrase)
