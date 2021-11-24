import ase.io.vasp
from potcar import potcar
from formatter import formatter

def poscar(dims, Ncells, typecount, names,opsys):
    dir="base-structures/"
    for i in range(Ncells):
        if lattice[i] == "simple cubic":
            cell = ase.io.vasp.read_vasp(dir+"POSCAR_simple_cubic")
        else:
            cell = ase.io.vasp.read_vasp(dir+"POSCAR_{}".format(lattice[i]))
        ase.io.vasp.write_vasp("POSCAR{}".format(i+1),
                               cell*(dims[i][0],dims[i][1],dims[i][2]),
                               label = 'Cell {}'.format(i+1),direct=True,sort=True)

        with open('POSCAR{}'.format(i+1)) as f:
            lines = f.readlines()
        # uncomment to create the POTCARS
        potcar(names,opsys)
        lines[5] = (' '.join(names)+'\n')
        lines[6] = formatter(typecount[i])
        with open('POSCAR{}'.format(i+1), 'w') as f:
            f.writelines(lines)
"""
opsys = "windows"
names = ["Al","Nb","Ta","Ti","Zr"]
typecount = [[5,7,3,14,19],[2,3,4,14,25],
             [8,10,1,9,20],[7,6,2,15,18],
             [2,10,2,20,14]]

cell_dims = [[2,3,4],[2,2,3],
             [2,2,3],[2,3,4],
             [2,3,4]]

lattice = ["bcc","hcp","fcc","bcc","bcc"]
Ncells = 5
poscar(cell_dims,Ncells,typecount,names)
"""
opsys = "windows"
names = ["Au","Pt"]
typecount = [[16,16],[16,16]]
cell_dims = [[2,2,2],[2,2,2]]
lattice = ["fcc","fcc"]
Ncells = 2
poscar(cell_dims,Ncells,typecount,names,opsys)
