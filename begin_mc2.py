from mc2 import mc2


"""
Ncells -->  Ncells should be equal to the number of elements being considered
               and listed in the POTCAR file. an example: NiCr requires 2 cells,
               NiCrCo requires 3 cells, etc.

T --> Temperature which appears in the kT acceptance term
"""

step_limit = 2500
Ncells = 5
T = 1273
param = "begin" # should be 'begin' or 'continue'
supcomp_phrase = "vasp" # for Mustang runs
#supcomp_phrase = "mpirun -np $SLURM_NTASKS /opt/packages/VASP/VASP5/INTEL/vasp_std"

mc2(step_limit, Ncells, T, param, supcomp_phrase)
