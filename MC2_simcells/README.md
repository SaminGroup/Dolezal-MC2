## Dolezal MC2 Simulation Cell Routine

**Please make sure you have Atomic Simulation Environment (ASE) installed before trying to use this tool "pip install ase", "conda install -c conda-forge ase", https://github.com/rosswhitfield/ase

Before calling "write_poscars.py" make sure you have the right POTCAR files in the potcars/ directory and they are of the form {}_POTCAR (Au_POTCAR, Pt_POTCAR, etc)

The POSCAR1 file in potcars/ can be deleted. Once ready, call write_poscars.py, "python3 write_poscars.py" where a few prompts will appear, here are two examples,

------------- Using Default ----------------

Op. System (windows, linux, unix)? windows

Species Names (Au Pt Zr Ti etc)? Au Pt

Total system concentration? 0.25 0.75

Do you wish to create your own simcells or use default (y or d)? d

How many atoms total per sim cell (32 or 48)? 32

-------------------- End --------------------

--------------- Create Own ------------------

Op. System (windows, linux, unix)? windows

Species Names (Au Pt Zr Ti etc)? Au Pt

Total system concentration? 0.25 0.75

Do you wish to create your own simcells or use default (y or d)? y

How many atoms total per sim cell? 64

Lattice type for each cell (sc bcc fcc hcp) fcc bcc

Input cell dimensions for Cell 1, fcc, which has 4 atoms in the unit cell (x y z): 2 2 4

Input cell dimensions for Cell 2, bcc, which has 2 atoms in the unit cell (x y z): 2 4 4

-------------------- End --------------------

--------------- Larger System ---------------

Op. System (windows, linux, unix)? windows

Species Names (Au Pt Zr Ti etc)? Al Nb Ta Ti Zr

Total system concentration? 0.1 0.15 0.05 0.3 0.4

Do you wish to create your own simcells or use default (y or d)? d

How many atoms total per sim cell (32 or 48)? 48

-------------------- End --------------------

At the end you will have "m" poscars and "m" concatenated POTCAR files, where m is the number of species in the material you are investigating
