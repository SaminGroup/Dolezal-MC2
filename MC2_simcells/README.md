## Dolezal MC2 Simulation Cell Routine

**Please make sure you have Atomic Simulation Environment (ASE) installed before trying to use this tool "pip install ase", "conda install -c conda-forge ase", https://github.com/rosswhitfield/ase

Before calling "write_poscars.py" make sure you have the right POTCAR files in the potcars/ directory and they are of the form {}_POTCAR (Au_POTCAR, Pt_POTCAR, etc)

The POSCAR1 file in potcars/ can be deleted. Once ready, call write_poscars.py, "python3 write_poscars.py" where a few prompts will appear, here are two examples,

------------- Using Default ----------------

![image](https://user-images.githubusercontent.com/47109396/143620800-60f45a04-6130-4f81-b987-e67527c24984.png)

-------------------- End --------------------

--------------- Create Own ------------------

![image](https://user-images.githubusercontent.com/47109396/143621142-5c527c29-308d-4408-98c5-8f2dea23107f.png)

-------------------- End --------------------

Make sure to list the concenrations in the same order as the line right above it. In both examples I selected "n" for the POTCARs -- only do this if you already have generated all the concatentated POTCARs. Otherwise, select "y", and POTCAR1,...POTCARm will be generated.
