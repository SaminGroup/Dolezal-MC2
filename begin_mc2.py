"""
Ncells -->  Ncells should be equal to the number of elements being considered
               and listed in the POTCAR file. an example: NiCr requires 2 cells,
               NiCrCo requires 3 cells, etc.

T --> Temperature which appears in the kT acceptance term
"""

step_limit = 2500 # max steps to run for
Ncells = 3 # should match number of species
T = 1873.15 # simulation temperature
param = "begin" # should be 'begin' or 'continue'
supcomp_phrase = "vasp" # for Mustang runs
#supcomp_phrase = "mpirun -np $SLURM_NTASKS /opt/packages/VASP/VASP5/INTEL/vasp_std"
intra = False # True, perform intraswaps, False, do not


from mc2_global import mc2



mc2(step_limit, Ncells, T, param, intra, supcomp_phrase)
