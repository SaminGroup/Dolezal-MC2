# Dolezal-MC2

This is my version of MC2 inspired by the released version (https://www.nature.com/articles/s41524-019-0259-z), but with the updated acceptance criterion derived
here (https://ui.adsabs.harvard.edu/abs/2020PhRvE.101f3306A/abstract). To use this method, please install CVXPY, "pip install cvxpy", "conda install -c conda-forge cvxpy", https://www.cvxpy.org/install/index.html. This package is required to solve Ax = b subject to the constraint that the molar fraction vector sums to 1 and the indices of f be in the domain [0,1]. 

MC2 is intitialized by calling the begin_mc2.py script, "python3 begin_mc2.py". Before calling this make sure to have the following in the parent directory:

        1. data/ directory
        2. potcar/ directory with each species POTCAR file inside labeled {}_POTCAR (i.e. Al_POTCAR, Nb_POTCAR, Ta_POTCAR, etc.)
        3. POSCAR1, POSCAR2, ..., POSCARm 
        4. Concatenated POTCAR1, POTCAR2,..., POTCARm files
        5. INCAR0 for the intitial relaxation VASP calculation and INCAR1 for the MC2-ISIF-3 calculations
        6. KPOINTS
        7. begin_mc2.py
        8. mc2.py
        9. mc2_global.py
        10. myfuncs.py


### Unsure how to make the MC2 simulation cells? Please check out MC2_simcells/ where you can easily make the simulation cells.  

For my work I used Bridges-2 at Pittsburgh Supercomputing Center and AFRL DSRC Mustang to execute vasp which is why I've included the supcomp_phrase param. Please update this parameter with a string containing the correct phrase to call VASP for your work. I usually did not implement the global flip and I would turn off intraswap (intra = False) until the system was in equilibrium. The best way to see equilibrium is to watch the molar fractions vs. steps.

mc2.py executes the algorithm with no global flips while mc2-global.py does include the global flip. I've generated two versions just in case. I've included INCAR0 and INCAR1 files as a sample of how I ran the mc2 algorithm for my research. I've also included data_to_csv.py as an example of how I analyzed the MC2 results. This file generates a CSV file which holds all of the MC data. It also generates plots of the data. There are a few user inputs and I've listed an example below,

Please list the species in order as found in the POSCAR: Al Nb Ta Ti Zr

Simulation Temperature: 1373

How many atoms per cell? 48 48 48 48 48

Equilibrium Begins: 500

Cell Shapes: BCC HCP FCC BCC BCC
