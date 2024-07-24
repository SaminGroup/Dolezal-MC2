import os
import numpy as np
from glob import glob
from ase.io import read
from math import isclose
from scipy.sparse.linalg import gmres


def build_X(unique_atomic_symbols, r=None):
    """

    cells --> an array that holds each cell

    N --> an array that gives total number of atoms in cells 1 and 2

    Builds the X matrix which tells us how many of species i = 1,...,m are in
    cell j = 1,...,m

    example: x_11 tells me how many of species 1 are in cell 1
                x_21 tells me how many of species 2 are in cell 1, and so on

    """


    sim_cells = glob("POSCAR[0-9]")
    sim_cells.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))

    each_row = []
    for acell in sim_cells:

        atoms = read(acell,format='vasp')

        # Calculate the total number of atoms
        total_atoms = len(atoms)

        # Calculate the concentration of each atomic species
        concentrations = []
        for symbol in unique_atomic_symbols:
            count = atoms.get_chemical_symbols().count(symbol)
            concentration = count / total_atoms
            concentrations.append(concentration)
        each_row.append(concentrations)

    X = np.array(each_row).T



    if r is not None:
        acell = "POSCAR"

        atoms = read(acell,format='vasp')

        # Calculate the total number of atoms
        total_atoms = len(atoms)

        # Calculate the concentration of each atomic species
        concentrations = []
        for symbol in unique_atomic_symbols:
            count = atoms.get_chemical_symbols().count(symbol)
            concentration = count / total_atoms
            concentrations.append(concentration)

        concentrations = np.array(concentrations)
        X[:,r] = concentrations

    return(np.array(X))


def mofac_test(F, tol=1e-6):
    
    # Check if all elements of F are within the range [0, 1]
    for f in F:
        if f < 0.0 or f > 1.0:
            return 0
    
    # Check if the sum of F is approximately 1, considering numerical tolerance
    if isclose(sum(F), 1.0, rel_tol=tol, abs_tol=tol):
        return 1
    
    return 0

def calc_f(X,C):
    
    try:
        
        f = np.linalg.solve(X,C)       
    
    except:
        
        f = gmres(X,C, tol=1e-12)[0]
    
    f = f/np.sum(f) # normalize f
    
    return(f)

def potcar(names,opsys):
    """
    names is a list of length m and each index is a string of the element name
    we will write to the POTCAR file in the order that is provided in the names
    list

    names = ['Type 1', 'Type 2',.., 'Type m']
    syntax of names should be captial letter followed by lowercase:

    He, Li, Au, etc
    """
    m = len(names)

    if opsys == "windows":
        beginning = 'type '
        file_name = ' potcars\{}_POTCAR'
    else:
        beginning = 'cat '
        file_name = ' potcars/{}_POTCAR'

    for i in range(m):
        ending = ' > POTCAR{}'.format(i+1)
        file_names = []
        for k in range(m):
            file_names.append(file_name.format(names[k]))
            middle = "".join(file_names)
            phrase = beginning+middle+ending

        os.system(phrase)


def formatter(alist):
    newline = '  '.join(str(num) for num in alist)+"\n"
    return(newline)
