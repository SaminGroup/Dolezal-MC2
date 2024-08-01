import os
import re
import copy
import numpy as np
from math import log
from random import uniform
from shutil import copyfile


def acceptance(dG, cost, w):
    """
    plugs the enthalpy and volume change into the exponential and checks it
    against a random number between 0 and 1
    """
    
    
    # I take care of 1/kT in the dG calculation
    equation = exp(-(dG*(1-w) + (w/2)*cost))
    choice = min(equation,1)
    
    # Generate a random number between 0 and 1
    r = np.random.uniform(0, 1)
        
    if r < choice:
        return(1)
    else:
        return(0)


def atomic_percents(X, xdata):
    """
    xdata = [[cell 1], [cell 2], ..., [cell m]]
    cell 1 = [[type 1], [type 2], ...
    , [type m-1]]
    """
    m = len(X[0])
    for i in range(m):
        for j in range(m):
            xdata[j][i].append(X[i][j])

    return(xdata)


def concentrations(cells, N):
    """
    Returns the concentration of each species in the whole system

    total --> total number of atoms in the system

    counts --> total number of each species type in the system
    """

    m = len(cells)
    types = np.arange(1, m+1)
    counts = []
    X = []

    for atype in types:
        type_count = []
        for acell in cells:
            for apos in acell:
                if apos[3] == atype:
                    type_count.append(1)
        counts.append(sum(type_count))

    concen = np.array(counts)/N
    np.savetxt("data/concentration",concen)
    return(concen)


def create_cells(Ncells):
    cell_number = np.arange(1,Ncells+1)
    all_positions = []
    for acell in cell_number:
    # with the new POSCAR file written, we want to extract the atom positions
        with open('POSCAR{}'.format(acell)) as f:
            lines = f.readlines()
            f.close()
    # pull out the names and count of each species and total count
        names = lines[5].split()
        number_of_each_type = np.asarray(lines[6].split(),dtype=int)
        Natoms = sum(number_of_each_type)
    # now grab all positions in this POSCAR file
        newpos = []
        for i in range(8,Natoms+8):
            pos = np.asarray(lines[i].split(),dtype=float);pos=list(pos)
            newpos.append(pos)


        positions = [] # simcell_{}
        species = cell_number
        for i in range(len(number_of_each_type)):
            trials = 0
            while trials < number_of_each_type[i]:
                    newpos[0].append(species[i])
                    positions.append(newpos[0])
                    newpos.remove(newpos[0])
                    trials += 1

        all_positions.append(positions)

    return(np.array(all_positions), names, Natoms*Ncells)


def accepted_state_update(cells,r,m):
    with open('POSCAR{}'.format(r+1)) as f:
        lines = f.readlines()
        f.close()
    # pull out the names and count of each species and total count
    names = lines[5].split()
    number_of_each_type = np.asarray(lines[6].split(),dtype=int)
    Natoms = sum(number_of_each_type)
    # now grab all positions in this POSCAR file
    newpos = []
    for i in range(8,Natoms+8):
        pos = np.asarray(lines[i].split(),dtype=float);pos=list(pos)
        newpos.append(pos)


    positions = [] # simcell_{}
    species = np.arange(1,m+1)
    for i in range(len(number_of_each_type)):
        trials = 0
        while trials < number_of_each_type[i]:
                newpos[0].append(species[i])
                positions.append(newpos[0])
                newpos.remove(newpos[0])
                trials += 1


    cells[r] = np.array(positions)
    return(cells)


def global_accepted_state_update(cells,m):
    for k in range(m):
        with open('POSCAR{}'.format(k+1)) as f:
            lines = f.readlines()
            f.close()
        # pull out the names and count of each species and total count
        names = lines[5].split()
        number_of_each_type = np.asarray(lines[6].split(),dtype=int)
        Natoms = sum(number_of_each_type)
        # now grab all positions in this POSCAR file
        newpos = []
        for i in range(8,Natoms+8):
            pos = np.asarray(lines[i].split(),dtype=float);pos=list(pos)
            newpos.append(pos)


        positions = [] # simcell_{}
        species = np.arange(1,m+1)
        for i in range(len(number_of_each_type)):
            trials = 0
            while trials < number_of_each_type[i]:
                    newpos[0].append(species[i])
                    positions.append(newpos[0])
                    newpos.remove(newpos[0])
                    trials += 1


        cells[k] = np.array(positions)
    return(cells)


def data_text(dataset,stepcount,xdata):
    #-----------------------------------------------
    # Step 1: write the mofrac, energy, and volume data to separate
    # text files. Each row is one cell's data
    #-----------------------------------------------
    filenames = ['data/mofac', 'data/energy', 'data/volume']
    for i in range(len(filenames)):
        np.savetxt(filenames[i],dataset[i])
    #-----------------------------------------------
    # Step 2: write the MC steps to a text file
    #-----------------------------------------------
    np.savetxt("data/stepcount",[stepcount],fmt="%s")
    #-----------------------------------------------
    # Step 3: record the atomic percents
    #-----------------------------------------------
    xdata = np.array(xdata)
    np.savetxt("data/xdata",xdata.flatten())

def data_update(r, f, mofac, newvars, ev_data):
    m = len(mofac)
    #---------------------------------------
    # Step 1: molar fractions change regardless of which cell was chosen so
    # we need to handle that data differenlty
    #---------------------------------------
    for i in range(m):
        mofac[i].append(f[i])
    #---------------------------------------
    # Step 2: now we can update the energy and volume data, starting with the
    # selected cell's data
    #---------------------------------------
    for i in range(0, len(ev_data)):
        ev_data[i][r].append(newvars[i])
    #---------------------------------------
    # Step 2.1: update the cell's that were not chosen (repeat latest value)
    #---------------------------------------
    for i in range(len(ev_data)):
        for k in range(len(ev_data[i])):
            if k != r:
                ev_data[i][k].append(ev_data[i][k][-1])
    #---------------------------------------
    # Step 3: now we create dataset
    #---------------------------------------
    dataset = [mofac, ev_data[0], ev_data[1]]
    return(dataset)



def global_data_update(f, mofac, E_data, elist, V_data, vlist):
    m = len(mofac)
    #---------------------------------------
    # Step 1: molar fractions change regardless of which cell was chosen so
    # we need to handle that data differenlty
    #---------------------------------------
    for i in range(m):
        mofac[i].append(f[i])
        E_data[i].append(elist[i])
        V_data[i].append(vlist[i])
    return([mofac, E_data, V_data])


# High precision floating-point library
from mpmath import mp, log, exp
# Set the precision to 25 decimal places
mp.dps = 25


def enthalpy_and_volume(E, V, X, X_new, Flist, m, n, T):
    k = 8.617333262145e-5  # Boltzmann constant in eV/K
    kT = k * T
    B = 1 / kT
    N = n*m
    # E = [[e_i, e_f], [e_i, e_f]] for cell 1 and cell 2
    f0 = Flist[0]
    f = Flist[1]
    
    dH_before = 0.0
    dH_after = 0.0
    
    dH_cost_before = 0.0
    dH_cost_after = 0.0
    
    for i in range(m):
        dH_before += float(E[i][0]) * float(f0[i])
        dH_after += float(E[i][1]) * float(f[i])
        
        dH_cost_before += float(E[i][0])
        dH_cost_after += float(E[i][1])
    
    # Initialize the terms
    dV_before = mp.mpf('0.0')
    dV_after = mp.mpf('0.0')
    # Compute dV_before
    for i in range(m):  # Loop over phases (nu)
        if f0[i] > 0:
            dv_term = N * float(f0[i]) * log(float(V[i][0]) / n)
            dV_before += dv_term

    # Compute dV_after
    for i in range(m):  # Loop over phases (nu)
        if f[i] > 0:
            dv_term = N * float(f[i]) * log(float(V[i][1]) / n)
            dV_after += dv_term
    
    
    dX_before = 0.0
    dX_after = 0.0
    
    # Compute dX_before
    for j in range(m):  # Loop over phases (nu)
        N_f0 = float(N) * float(f0[j])
        for i in range(m):  # Loop over species (i)
            if X[i, j] > 0 and f0[j] > 0:
                dx_term_before = N_f0 * float(X[i, j]) * log(float(X[i, j])) - N_f0 * float(X[i, j])
                dX_before += dx_term_before
             

    # Compute dX_after
    for j in range(m):  # Loop over phases (nu)
        N_f = mp.mpf(float(N)) * mp.mpf(float(f[j]))
        for i in range(m):  # Loop over species (i)
            if X_new[i, j] > 0 and f[j] > 0:
                dx_term_after = N_f * float(X_new[i, j]) * log(float(X_new[i, j])) - N_f * float(X_new[i, j])
                dX_after += dx_term_after
    
    
    dX = dX_after - dX_before
    
    dV = dV_after - dV_before
        
    # reduced values comes to this expression (N = n.m and u = U/n)
    dH = (m * B) * (dH_after - dH_before)
    
    dH_cost = (m * B) * (dH_cost_after - dH_cost_before)
    # reduce the entropic contribution by a factor of dH_cost
    dG = dH - dV - dX
    
    return dG, dH_cost



def flip_and_lever(cells,C):

    m = len(cells)
    #----------------------------------
    # Step 1: perform a flip
    #----------------------------------
    singular = 1 ; cycle_count = 0
    while singular == 1:
        if cycle_count > 100:
            r = 0
            X = build_X(cells)
            F = calc_f(X,C)
            return(cells, r, X, F, 1)
        chosen = []
        copy_cells = copy.deepcopy(cells) # nested lists need deep copy
        r = int(uniform(0,len(copy_cells))) # the random cell
        r1 = int(uniform(0,len(copy_cells[r]))) # the random atom in cell r
        flip_choice = [x for x in range(1,m+1) if x != copy_cells[r][r1][3]]
        random_flip = flip_choice[int(uniform(0,len(flip_choice)))]
        copy_cells[r][r1][3] = random_flip
        chosen.append(r1)
        flipcount = len(chosen)
        #----------------------------------------
        # if we failed lever_check, flip 5 atoms
        #----------------------------------------

        X = build_X(copy_cells)
        F = calc_f(X,C)
        cycle_count += 1
        if mofac_test(F) == 1:
            singular = 0

    #---------------------------------
    # Step 2: update the state with the lever-approved flip
    #---------------------------------
    cells = copy_cells
    #---------------------------------
    # Step 3: edit the correct POSCAR file; edit the atomic positions
    #---------------------------------
    cell_number = r+1
    #---------------------------------
    with open('POSCAR{}'.format(cell_number)) as f:
        lines = f.readlines()

    # identify the correct spot in the file to add the new positions
    for i in range(len(lines)):
        if lines[i] == 'Direct\n':
            lines = lines[:i+1]
            break

    # update the atom positions and append to lines

    types = np.arange(1, m+1)
    counts = []
    for atype in types:
        count = []
        for apos in cells[r]:
            if apos[3] == atype:
                count.append(1)
                lines.append('  {}  {}  {}\n'.format(apos[0], apos[1], apos[2]))
        counts.append(sum(count))

    # now update the number of each atom type
    # since we save CONTCAR as latest POSCAR it is actually line 6 that we
    # need to update

    lines[6] = formatter(counts)

    # We will now update the POSCAR file, making sure to keep
    # POSCAR{} untouched until the acceptance has been met
    with open('POSCAR', 'w') as f:
        f.writelines(lines)
    #--------------------------------
    return(cells, r, X, F, 0)


def intraswap(cells):

    m = len(cells)
    #----------------------------------
    # Step 1: perform a flip
    #----------------------------------
    copy_cells = copy.deepcopy(cells) # nested lists need deep copy

    r = int(uniform(0,len(copy_cells))) # the random cell
    r1 = int(uniform(0,len(copy_cells[r]))) # the random atom in cell r
    r2 = r1
    # make sure we do not get stuck looking for a species that disappeared
    with open("POSCAR{}".format(r+1)) as f:
        lines = f.readlines()

    typecount = np.asarray(lines[6].split(),dtype=int)

    all_options = np.arange(1,m+1,dtype=int)
    only_options = np.array([all_options[i] for i in range(m) if typecount[i] != 0])
    if len(only_options) < 2:
        return(cells,r,0,0)

    # first pick what to flip atom 1 to
    flip_choice = [x for x in only_options if x != copy_cells[r][r1][3]]
    random_flip = flip_choice[int(uniform(0,len(flip_choice)))]
    # now find another atom that is the species that atom 1 will become
    while copy_cells[r][r2][3] != random_flip:
        r2 = int(uniform(0,len(copy_cells[r])))
    # once found, flip atom 2 to what atom 1 is
    copy_cells[r][r2][3] = copy_cells[r][r1][3]
    # now flip atom 1 to what atom 2 was
    copy_cells[r][r1][3] = random_flip



    X = build_X(copy_cells)


    #---------------------------------
    # Step 2: update the state with the lever-approved flip
    #---------------------------------
    cells = copy_cells
    #---------------------------------
    # Step 3: edit the correct POSCAR file; edit the atomic positions
    #---------------------------------
    cell_number = r+1
    #---------------------------------
    with open('POSCAR{}'.format(cell_number)) as f:
        lines = f.readlines()

    # identify the correct spot in the file to add the new positions
    for i in range(len(lines)):
        if lines[i] == 'Direct\n':
            lines = lines[:i+1]
            break

    # update the atom positions and append to lines

    types = np.arange(1, m+1)
    counts = []
    for atype in types:
        count = []
        for apos in cells[r]:
            if apos[3] == atype:
                count.append(1)
                lines.append('  {}  {}  {}\n'.format(apos[0], apos[1], apos[2]))
        counts.append(sum(count))

    # now update the number of each atom type
    # since we save CONTCAR as latest POSCAR it is actually line 6 that we
    # need to update
    lines[6] = formatter(counts)

    # We will now update the POSCAR file, making sure to keep
    # POSCAR{} untouched until the acceptance has been met
    with open('POSCAR', 'w') as f:
        f.writelines(lines)
    #--------------------------------
    return(cells, r, X, 1)


def formatter(alist):
    space1 = "   "
    newline = space1 + '  '.join(str(num) for num in alist)+"\n"
    return(newline)



def global_flip(cells,C):
    """
    randomly flip one or more atoms in all cells simultaneously
    -----------------------------------------------------------
    Natoms = how many atoms we have in each cell
    N_flip = max number of atoms that can be flipped at once
             i.e. flip is uniform(1,N_flip+1)
    """
    m = len(cells)
    Natoms = len(cells[0])

    copy_cells = copy.deepcopy(cells)
    number_to_flip = 2
    singular = 1
    while singular == 1:
        for i in range(m):
            flipped = 0
            options = [x for x in range(0,Natoms)]
            selected = []
            while flipped < number_to_flip:
                the_selection = int(uniform(0,len(options)))
                selected.append(options[the_selection])
                options.remove(options[the_selection])
                flipped += 1

            for atom in selected:
                flip_choice = [x for x in range(1,m+1) if x != copy_cells[i][atom][3]]
                random_flip = flip_choice[int(uniform(0,len(flip_choice)))]
                copy_cells[i][atom][3] = random_flip

        X = build_X(copy_cells)
        F = calc_f(X,C)
        if mofac_test(F) == 1:
            singular = 0


    cells = copy_cells

    for cell_number in range(0,m):
        r = cell_number + 1
        with open('POSCAR{}'.format(r)) as f:
            lines = f.readlines()

        ## identify the correct spot in the file to add the new positions
        for i in range(len(lines)):
            if lines[i] == 'Direct\n':
                lines = lines[:i+1]
                break

        ## update the atom positions and append to lines

        types = np.arange(1, m+1)
        counts = []
        for atype in types:
            count = []
            for apos in cells[cell_number]:
                if apos[3] == atype:
                    count.append(1)
                    lines.append('  {}  {}  {}\n'.format(apos[0], apos[1], apos[2]))
            counts.append(sum(count))

        ## now update the number of each atom type
        ## since we save CONTCAR as latest POSCAR it is actually line 6 that we
        ## need to update
        lines[6] = formatter(counts)
        ## We will now update the POSCAR file, making sure to keep
        ## POSCAR{} untouched until the acceptance has been met
        with open('global_POSCAR{}'.format(r), 'w') as f:
            f.writelines(lines)

    return(cells, X, F)

def import_data(m):

    # to avoid issues with appending to numpy array, convert to lists
    mofac = list(np.loadtxt("data/mofac",dtype=float))
    energy = list(np.loadtxt("data/energy",dtype=float))
    volume = list(np.loadtxt("data/volume",dtype=float))
    stepcount = np.loadtxt("data/stepcount",dtype=int)

    for i in range(m):
        mofac[i] = list(mofac[i])
        energy[i] = list(energy[i])
        volume[i] = list(volume[i])

    xdata = list(np.loadtxt("data/xdata").reshape(m,m,stepcount+1))
    for i in range(m):
        xdata[i] = list(xdata[i])
        for j in range(m):
            xdata[i][j] = list(xdata[i][j])
    #-------------------------------------
    return(mofac, energy, volume, stepcount, xdata)


def mofac_test(F):
    for f in F:
        if abs(f) < 1e-10:
            f = 0.0
        if f < 0.0:
            return(0)
        if f > 1.0:
            return(0)

    return(1)



def oszicar():
    """
    this function strips the latest energy value from the OSZICAR VASP
    output file and adds it to the recording of energy values. It also
    returns the new energy value so that it can be passed to the acceptance
    function.
    """
    # start by saving every line of the file in a list
    with open('OSZICAR') as f:
        lines = f.read().splitlines()
    # the very last line in the OSZICAR file is the one we are interested in
    new_energy = lines[-1]
    # now pull only the digits from the final string

    new_energy = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?"
                                                            , new_energy)

    # the format of the OSZICAR file will always be the same -- the value we
    # are interested in will always be new_energy[1]
    new_energy = float(new_energy[1])
    return(new_energy)

def outcar():
    """
    This function extracts the new cell volume, updates the volume record, and
    returns the new volume so that it may be passed to the acceptance function

    units are cubic angstroms
    """
    with open('OUTCAR') as f:
        lines = f.read().splitlines()
    # we need to run a search through the lines and find the last time the
    # volume was updated
    volumes = [aline for aline in lines if "volume" in aline]
    # the final volume update is the last entry in the volumes list
    new_volume = volumes[-1]
    match = "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?"
    new_volume = [float(x) for x in re.findall(match, new_volume)][0]

    return(new_volume)

def total_count(cells):
    total = []
    for acell in cells:
      total.append(len(acell))

    return(sum(total))



def update_poscar(r, cells, names):
    #-----------------------------------------------------
    # Step 1: since the step has been accepted, we can save this
    #         poscar file
    #-----------------------------------------------------
    m = len(cells)
    copyfile('CONTCAR{}'.format(r+1), 'POSCAR{}'.format(r+1))
    #-----------------------------------------------------
    # POSCAR{} will always contain the full names list and
    # all species counts, including 0s
    #-----------------------------------------------------
    cell_number = r+1
    with open('POSCAR{}'.format(cell_number)) as f:
        lines = f.readlines()


    types = np.arange(1, m+1)
    counts = []
    for atype in types:
        count = []
        for apos in cells[r]:
            if apos[3] == atype:
                count.append(1)
        counts.append(sum(count))

    lines[5] = '   ' + ' '.join(names)+'\n'
    lines[6] = formatter(counts)

    # We will now update the POSCAR file, making sure to keep
    # POSCAR{} untouched until the acceptance has been met
    with open('POSCAR{}'.format(cell_number), 'w') as f:
        f.writelines(lines)


def global_poscar_update(m,cells,names):
    """
    if the global flip was accepted then all poscar files need to be updated
    to their post-global-flip CONTCAR files
    """
    for i in range(m):
        copyfile('CONTCAR{}'.format(i+1), 'POSCAR{}'.format(i+1))


        with open('POSCAR{}'.format(i+1)) as f:
            lines = f.readlines()


        types = np.arange(1, m+1)
        counts = []
        for atype in types:
            count = []
            for apos in cells[i]:
                if apos[3] == atype:
                    count.append(1)
            counts.append(sum(count))


        lines[5] = '   ' + ' '.join(names)+'\n'
        lines[6] = formatter(counts)

        # We will now update the POSCAR file, making sure to keep
        # POSCAR{} untouched until the acceptance has been met
        with open('POSCAR{}'.format(i+1), 'w') as f:
            f.writelines(lines)


def local_potcar_update(m, names, r):
    '''
    if we do a local flip then we only need to update POSCAR{r}
    and POTCAR{r}
    '''
    updated_names = names.copy()
    with open('POSCAR') as f:
        lines = f.readlines()
    #----------------------
    # Step 1: read in counts from POSCAR
    #----------------------
    # Step 1a: if there is a 0 count (determined in flip_and_lever), then we
    #          need to eliminate that name from the POSCAR file
    #----------------------
    counts = [int(x) for x in lines[6].split()]

    updated_names = [updated_names[i] for i in range(m) if counts[i] != 0]

    lines[5] = '   ' + ' '.join(updated_names)+'\n'
    counts = [x for x in counts if x != 0]
    lines[6] = formatter(counts)
    #-----------------------
    # Step 1b: write new names and counts to the POSCAR file
    #-----------------------
    with open('POSCAR', 'w') as f:
        f.writelines(lines)
    #-----------------------
    # Step 2: update POTCAR{r} to reflect the new species in the cell
    #-----------------------
    new_potcar(updated_names, r+1)
    #-----------------------



def global_potcar_update(m, names):
    '''
    for a global flip we will need to read in and edit all global_POSCAR{} and
    POTCAR files
    '''
    filename = 'global_POSCAR{}'
    for i in range(m):
        #--------------------
        # Step 1: read in global_POSCAR{i} and extract the counts
        #--------------------
        # Step 1a: eliminate name of 0 count species and update global_POSCAR{i}
        updated_names = names.copy()
        with open(filename.format(i+1)) as f:
            lines = f.readlines()

        counts = [int(x) for x in lines[6].split()]
        updated_names = [updated_names[i] for i in range(m) if counts[i] != 0]

        lines[5] = '   ' + ' '.join(updated_names)+'\n'
        counts = [x for x in counts if x != 0]
        lines[6] = formatter(counts)

        with open(filename.format(i+1), 'w') as f:
            f.writelines(lines)
        #--------------------
        # Step 2: update POTCAR{i}
        new_potcar(updated_names, i+1)
        #--------------------

def new_potcar(names, cell_number):
    """
    names is a list of length m and each index is a string of the element name
    we will write to the POTCAR file in the order that is provided in the names
    list

    names = ['Type 1', 'Type 2',.., 'Type m']
    syntax of names should be captial letter followed by lowercase:

    He, Li, Au, etc
    """
    if len(names) > 1:
        file_name = ' potcar/{}_POTCAR '
        beginning = 'cat'
        ending = ' > POTCAR{}'.format(cell_number)
        file_names = []
        for k in range(len(names)):
            file_names.append(file_name.format(names[k]))
            middle = "".join(file_names)
            phrase = beginning+middle+ending

        os.system(phrase)

    else:
        with open('potcar/{}_POTCAR'.format(names[0])) as f:
            lines = f.readlines()
        with open('POTCAR{}'.format(cell_number), 'w') as f:
            f.writelines(lines)

def initial_vasp_run(m, names, cells, supcomp_phrase):
    """
    will run vasp on all initial cells and return E,V
    """
    E = []
    V = []
    for i in range(m):
        #-----------------------------------------------------
        # Step 1: the initial file to POSCAR and POTCAR for the VASP run
        #-----------------------------------------------------
        copyfile("INCAR0", "INCAR")
        copyfile('POSCAR{}'.format(i+1), 'POSCAR')
        local_potcar_update(m,names,i)
        copyfile('POTCAR{}'.format(i+1), 'POTCAR')
        #-----------------------------------------------------
        # Step 2: run VASP to relax volume and positions
        #-----------------------------------------------------
        os.system(supcomp_phrase) # VASP on Mustang
        copyfile("CONTCAR","POSCAR")
        copyfile("INCAR1", "INCAR")
        os.system(supcomp_phrase)
        #-----------------------------------------------------
        # Step 3: record the energy and volume from oszicar and outcar
        #-----------------------------------------------------
        E.append(oszicar())
        V.append(outcar())
        #-----------------------------------------------------
        # Step 4: CONTCAR becomes the new POSCAR{i}
        #-----------------------------------------------------
        copyfile("CONTCAR","CONTCAR{}".format(i+1))
        update_poscar(i, cells, names)
        #-----------------------------------------------------
    copyfile("INCAR2", "INCAR")
    return(E,V)


def vasp_run(r,supcomp_phrase):
    """
    runs vasp on the cell that was chosen at random and flipped

    returns the new energy and volume which will be used to decide the new state
    will remain

    """
    cellnum = r+1
    #----------------------------------------------
    # Step 1: copy POTCAR{r} to POTCAR
    #----------------------------------------------
    copyfile("INCAR2", "INCAR")
    copyfile('POTCAR{}'.format(cellnum), 'POTCAR')
    #----------------------------------------------
    # Step 2: run VASP
    #----------------------------------------------
    os.system(supcomp_phrase) # VASP on Mustang
    #os.system("rm CHG CHGCAR DOSCAR EIGENVAL IBZKPT PCDAT vasprun.xml REPORT XDATCAR")
    #-----------------------------------------------------
    # Step 3: record the energy and volume from oszicar and outcar
    #-----------------------------------------------------
    E = oszicar()
    V = outcar()
    #----------------------------------------------
    copyfile("CONTCAR", "CONTCAR{}".format(cellnum))

    return(E,V)


def global_vasp_run(m,supcomp_phrase,step):
    """
    for global flip we need to run and store all outputs until acceptance has
    been determined
    """
    E = []
    V = []
    copyfile("INCAR2", "INCAR")
    for i in range(1,m+1):
        #-----------------------------------------------------
        # Step 1: the initial file to POSCAR and POTCAR for the VASP run
        #-----------------------------------------------------
        copyfile('global_POSCAR{}'.format(i), 'POSCAR')
        copyfile('POTCAR{}'.format(i), 'POTCAR')
        #-----------------------------------------------------
        # Step 2: run VASP
        #-----------------------------------------------------
        os.system(supcomp_phrase) # VASP on Mustang
        #os.system("rm CHG CHGCAR DOSCAR EIGENVAL IBZKPT PCDAT vasprun.xml REPORT XDATCAR")
        #-----------------------------------------------------
        # Step 3: record the energy and volume from oszicar and outcar
        #-----------------------------------------------------
        E.append(oszicar())
        V.append(outcar())
        #-----------------------------------------------------
        # Step 4: Will not change CONTCAR to POSCAR until acceptance has been
        #         determined
        #-----------------------------------------------------
        copyfile('CONTCAR', 'CONTCAR{}'.format(i))
        #-----------------------------------------------------
    return(E,V)

def build_X(cells):
    """

    cells --> an array that holds each cell

    N --> an array that gives total number of atoms in cells 1 and 2

    Builds the X matrix which tells us how many of species i = 1,...,m are in
    cell j = 1,...,m

    example: x_11 tells me how many of species 1 are in cell 1
                x_21 tells me how many of species 2 are in cell 1, and so on

    """

    m = len(cells)
    types = np.arange(1, m+1)
    X = []
    for atype in types:
        xrow = []
        for acell in cells:
            type_count = []
            tot = len(acell)
            for apos in acell:
                if apos[3] == atype:
                    type_count.append(1)
            xrow.append(sum(type_count)/tot)
        X.append(xrow)

    return(np.array(X))


from scipy.linalg import qr, solve_triangular

def calc_f(X,C):
    
    n = np.loadtxt('data/atoms_per_cell')
    
    X_T = X*n
    C_T = C*n
    fluctuation = 1e-2 # relax the concentration by up to 1/4 of a percent
    # Add regularization to X
    try:
        n = X.shape[0]
        F_list = []
        
        for _ in range(100):
            # Introduce small random fluctuations to C
            C_fluctuated = C + np.random.uniform(-fluctuation, fluctuation, size=C.shape)
            
            # QR decomposition
            Q, R = qr(X)
            
            # Solve R*F = Q^T*C_fluctuated
            QTC = np.dot(Q.T, C_fluctuated)
            F = solve_triangular(R, QTC)
            F_list.append(F)
        
        # Average the results
        f = np.mean(F_list, axis=0)   
    
    except:
        
        f = gmres(X_T,C_T, tol=1e-12)[0]
    
    f = f/np.sum(f) # normalize f
    
    return(f)

def lever_check(X,f,C,step):
    new_c = X@f
    residual = abs(np.linalg.norm(new_c) - np.linalg.norm(C))
    if residual <= 5e-3:
        return(1)
    else:
        np.savetxt("data/failed-residual.txt", [residual])
        np.savetxt("data/failed-concentration.txt", new_c)
        return(0)
