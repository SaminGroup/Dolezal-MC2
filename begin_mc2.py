from mc2 import mc2


"""
all the data is saved to text files which are converted to a csv file once it
is off the supercomputer servers

mc2 takes 3 inputs
1. step_limit -->  temporary input while we iron out how well the algorithm
                   works. This input will be removed and a termination condition
                   will take its place (neglible dE between two steps)

3. Ncells -->  Ncells should be equal to the number of elements being considered
               and listed in the POTCAR file. an example: NiCr requires 2 cells,
               NiCrCo requires 3 cells, etc.

4. T --> Temperature which appears in the kT acceptance term

Before running mc2 make sure to have the correct number of POSCAR files created.
There needs to be one POSCAR{}_initial file per cell. Two cells will require a
POSCAR1_initial and POSCAR2_initial file. The user is required to provide the
necessary vectors and atoms per unit cell.

NOTE: please remember to list all elements alphabetically in the POTCAR file and
on line 6 in the POSCAR{}_initial file
"""

step_limit = 2500
Ncells = 5
T = 1273
param = "begin"
supcomp = "mustang"

mc2(step_limit, Ncells, T, param,supcomp)
