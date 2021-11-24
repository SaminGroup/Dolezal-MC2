import ase.io.vasp
import numpy as np
from potcar import potcar
from random import uniform
from scipy.sparse.linalg import gmres
from formatter import formatter


structures = {
'sc'  : 1,
'bcc' : 2,
'fcc' : 4,
'hcp' : 4
}

dims32 = {
'sc'  : [2,4,4],
'bcc' : [2,2,4],
'fcc' : [2,2,2],
'hcp' : [2,2,2]
}

dims48 = {
'sc'  : [3,4,4],
'bcc' : [2,3,4],
'fcc' : [2,2,3],
'hcp' : [2,2,3]
}

def poscar(dims,lattice, Ncells, typecount, names, opsys):
    dir="base-structures/"
    for i in range(Ncells):
        if lattice[i] == "simple cubic":
            cell = ase.io.vasp.read_vasp(dir+"POSCAR_simple_cubic")
        else:
            cell = ase.io.vasp.read_vasp(dir+"POSCAR_{}".format(lattice[i]))
        ase.io.vasp.write_vasp("POSCAR{}".format(i+1),
                               cell*(int(dims[i][0]),int(dims[i][1]),int(dims[i][2])),
                               label = 'Cell {}'.format(i+1),direct=True,sort=True)

        with open('POSCAR{}'.format(i+1)) as f:
            lines = f.readlines()

        potcar(names,opsys)
        lines[5] = (' '.join(names)+'\n')
        lines[6] = formatter(typecount[i])
        with open('POSCAR{}'.format(i+1), 'w') as f:
            f.writelines(lines)



opsys = input("Op. System (windows, linux, unix)? ")
names = input("Species Names (Au Pt Zr Ti etc)? ").split()
m = len(names)

C = np.asarray(input("Total system concentration? ").split(),dtype=float)
choice = input("Do you wish to create your own simcells or use default (y or d)? ")

if choice == 'd':
    N = int(input("How many atoms total per sim cell (32 or 48)? "))
    cells = ['bcc','fcc','hcp']
    lattice = []
    for i in range(m):
        r = int(uniform(0,3))
        lattice.append(cells[r])

else:
    N = int(input("How many atoms total per sim cell? "))
    lattice = input("Lattice type for each cell (sc bcc fcc hcp) ").split()

A = C*np.ones((m,m))
b = C*N*np.ones(m,)
x = gmres(A.T,b)[0]
x = np.array([int(j) for j in x])

typecount = []
for i in range(m):
    typecount.append([])
    for j in range(m):
        typecount[i].append(int(x[j]))
    if sum(typecount[i]) != N:
        add = N-sum(x)
        select = int(uniform(0,m))
        typecount[i][select] += add


if choice == 'd':
    if N == 32:
        pick = dims32
    else:
        pick = dims48

    dims = np.zeros((m,3))
    for i in range(m):
        dims[i] = pick[lattice[i]]

else:
    dims = np.zeros((m,3))
    for i in range(m):
        val = input("Input cell dimensions for Cell {}, {}, which has {} atoms in the unit cell (x y z): ".format(i+1,lattice[i],structures[lattice[i]]))
        val = np.asarray(val.split(),dtype=int)
        dims[i] = val

poscar(dims,lattice,m,typecount,names,opsys)
