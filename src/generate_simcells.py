import ase.io.vasp
import numpy as np
from ase.build import bulk
from random import uniform
import simcell_funcs as fun
from scipy.sparse.linalg import gmres


dims32 = {
'sc'  : (2,4,4),
'bcc' : (2,2,4),
'fcc' : (2,2,2),
'hcp' : (2,2,2)
}

dims48 = {
'sc'  : (3,4,4),
'bcc' : (2,3,4),
'fcc' : (2,2,3),
'hcp' : (2,2,3)
}

def poscar(Natoms, Ncells, typecount, names, opsys, genpot):
    for i in range(Ncells):
        try:
            cell = bulk(names[i],cubic=True)
        except:
            cell = bulk(names[i],orthorhombic=True)

        atomcount = len(cell)

        if atomcount == 1:
            if Natoms == 32:
                dims = dims32['sc']
            else:
                dims = dims48['sc']

        if 1 < atomcount < 4:
            if Natoms == 32:
                dims = dims32['bcc']
            else:
                dims = dims48['bcc']

        else:
            if Natoms == 32:
                dims = dims32['fcc']
            else:
                dims = dims48['fcc']

        ase.io.vasp.write_vasp("POSCAR{}".format(i+1),
                               cell*dims,
                               label = 'Cell {}'.format(i+1),direct=True,sort=True)

        with open('POSCAR{}'.format(i+1)) as f:
            lines = f.readlines()

        if genpot:
            fun.potcar(names,opsys)
        lines[5] = (' '.join(names)+'\n')
        lines[6] = fun.formatter(typecount[i])
        with open('POSCAR{}'.format(i+1), 'w') as f:
            f.writelines(lines)


print("\n --------------------- User Inputs ---------------------")
opsys = input(" 1. Op. System (windows, linux, unix)? ")
names = input(" 2. Species Names (Au Pt Zr Ti etc)? ").split()
m = len(names)

C = np.asarray(input(" 3. Total system concentration (e.g., 0.5, 0.5)? ").split(),dtype=float)

N = int(input(" 5. How many atoms total per sim cell (32 or 48)? "))

A = C*np.ones((m,m))
b = C*N*np.ones(m,)
x = gmres(A.T,b)[0]
x = np.array([int(j) for j in x])

# atom count for each cell is set to indices of X vector
typecount = []
for i in range(m):
    typecount.append([])
    for j in range(m):
        typecount[i].append(int(x[j]))
    if sum(typecount[i]) != N:
        add = N-sum(x)
        select = int(uniform(0,m))
        typecount[i][select] += add


num_swaps = int(4**3)

for _ in range(num_swaps):
    cell1, cell2 = np.random.choice(Ncells, 2, replace=False)
    type1, type2 = np.random.choice(len(names), 2, replace=False)
    
    # Swap one species with another in one cell
    if typecount[cell1][type1] > 0 and typecount[cell2][type2] > 0:
        typecount[cell1][type1] -= 1
        typecount[cell1][type2] += 1
    
    # Reverse the swap in another cell
    if typecount[cell2][type2] > 0:
        typecount[cell2][type2] -= 1
        typecount[cell2][type1] += 1
        

genpot = input(" 6. Generate POTCARs (y or n)? ")

if genpot == "y" or genpot == "yes":
    genpot = True
else:
    genpot = False

print(" --------------------- procedure initialized ---------------------")

poscar(N,m,typecount,names,opsys,genpot)

for i in range(m):
    if genpot:
        print(" ---- Generated POSCAR{} and POTCAR{}".format(i+1,i+1))
    else:
        print(" ---- Generated POSCAR{}".format(i+1))
print(" ----------------------- procedure complete -----------------------")
