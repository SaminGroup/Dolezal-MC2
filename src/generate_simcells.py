import ase.io.vasp
import numpy as np
from ase.build import bulk
from random import uniform
import src.simcell_funcs as fun
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
