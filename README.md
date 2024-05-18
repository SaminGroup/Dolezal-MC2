# Dolezal-MC2

This is our version of MC2 inspired by the released version (https://www.nature.com/articles/s41524-019-0259-z), but with the updated acceptance criterion derived
here (https://ui.adsabs.harvard.edu/abs/2020PhRvE.101f3306A/abstract). If you use this version please cite our paper, https://pubs.acs.org/doi/10.1021/acs.langmuir.1c03191, thanks!

I have updated all the preparation to be handled through a GUI, including the creation of simulation cells. Now, simply execute the prepare_mc2.py file and select your settings via the GUI. All necessary directories will be created if not found. EXCEPTION: Please remember to name the POTCAR files as follows, {}_POTCAR, where the atomic symbol should be provided as the prefix (i.e. Al_POTCAR, Nb_POTCAR, Ta_POTCAR, etc.) and place them in a directory named "potcar". Also, don't forget to provide a KPOINTS file and to edit the INCAR files appropriately for your system and HPC architecture (i.e., the NCORE tag should be updated depending on the nodes available for the MC2 run).

After executing the preparation script, execute the "begin.py" file to kick-off MC2.

A. Note: Because VASP is executed differently depending on the system, the GUI is equipped to receive the vasp execution command for your particular system.

B. Note: If you select the "continue" option from the GUI this will update the "begin.py" file accordingly. So if you want to launch MC2 from where it left off execute the prepare_mc2.py file, select continue, and then execute "begin.py" to kick-off MC2 in continue mode.

I've also included data_to_csv.py as an example of how I analyzed the MC2 results. This file generates a CSV file which holds all of the MC data. It also generates plots of the data. There are a few user inputs and I've listed an example below,

Please list the species in order as found in the POSCAR: Al Nb Ta Ti Zr

Simulation Temperature: 1373

How many atoms per cell? 48

Equilibrium Begins: 500

Cell Shapes: BCC HCP FCC BCC BCC
