import numpy as np
import os


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
    file_name = ' potcars\{}_POTCAR'

    if opsys == "windows":
        beginning = 'type '
    else:
        beginning = 'cat '

    for i in range(m):
        ending = ' > POTCAR{}'.format(i+1)
        file_names = []
        for k in range(m):
            file_names.append(file_name.format(names[k]))
            middle = "".join(file_names)
            phrase = beginning+middle+ending

        os.system(phrase)


def formatter(alist):
    space1 = " "
    newline = '  '.join(str(num) for num in alist)+"\n"
    return(newline)


def generate_initial_X(C, Natoms):
    # Number of species
    num_species = len(C)

    # Generate a random initial X matrix
    X = np.random.rand(num_species, num_species)

    # Normalize X to get the atomic percents
    X = X / np.sum(X, axis=0)

    # Calculate molar fraction matrix F
    F = np.linalg.lstsq(X, C, rcond=None)[0]

    # Calculate the actual concentration from F
    actual_C = np.dot(X, F)

    # Adjust X to match the user-defined concentration C
    X *= (C / actual_C)

    scaled_X = X*48

    X = np.round(scaled_X).astype(int)

    return X
