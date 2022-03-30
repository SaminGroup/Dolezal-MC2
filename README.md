# Dolezal-MC2

This is my version of MC2 inspired by the released version (https://www.nature.com/articles/s41524-019-0259-z), but with the updated acceptance criterion derived
here (https://ui.adsabs.harvard.edu/abs/2020PhRvE.101f3306A/abstract). If you use this version please cite our paper, https://pubs.acs.org/doi/10.1021/acs.langmuir.1c03191, thanks!

MC2 is intitialized by calling the begin_mc2.py script, "python3 begin_mc2.py". Before calling this make sure to have the following in the parent directory:

        1. data/ directory
        2. potcar/ directory with each species POTCAR file inside labeled {}_POTCAR (i.e. Al_POTCAR, Nb_POTCAR, Ta_POTCAR, etc.)
        3. POSCAR1, POSCAR2, ..., POSCARm 
        4. Concatenated POTCAR1, POTCAR2,..., POTCARm files
        5. INCAR0/1 for the initial volume/position relaxation VASP calculation
        6. INCAR2 for the MC2-ISIF-3 calculations
        7. INCAR-S for flipping away from singular states
        8. KPOINTS
        9. begin_mc2.py
        10. mc2.py
        11. mc2_global.py
        12. myfuncs.py


### Unsure how to make the MC2 simulation cells? Please check out SaminGroup/Dolezal-MC2_simcells where you can easily make the simulation cells.  

For my work I used Bridges-2 at Pittsburgh Supercomputing Center and AFRL DSRC Mustang to execute vasp which is why I've included the supcomp_phrase param. Please update this parameter with a string containing the correct phrase to call VASP for your work. I usually did not implement the global flip and I would turn off intraswap (intra = False) until the system was in equilibrium. The best way to see equilibrium is to watch the molar fractions vs. steps.

mc2.py executes the algorithm with no global flips while mc2-global.py does include the global flip. I've generated two versions just in case. I've included the INCAR files as a sample of how I ran the mc2 algorithm for my research. I've also included data_to_csv.py as an example of how I analyzed the MC2 results. This file generates a CSV file which holds all of the MC data. It also generates plots of the data. There are a few user inputs and I've listed an example below,

Please list the species in order as found in the POSCAR: Al Nb Ta Ti Zr

Simulation Temperature: 1373

How many atoms per cell? 48

Equilibrium Begins: 500

Cell Shapes: BCC HCP FCC BCC BCC
