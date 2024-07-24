import ase.io.vasp
import numpy as np
from ase.build import bulk
from scipy.sparse.linalg import gmres
import src.simcell_funcs as fun


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

def poscar(C, names, shapes, Natoms, genpot, opsys):
    
    m = len(names)
    
    # alphabetize the system and maintain proper species concentration
    sorted_lists = sorted(zip(names, C), key=lambda x: x[0])
    names = [x[0] for x in sorted_lists]
    C = np.array([x[1] for x in sorted_lists])
    
    #grab the location of the most prominent species
    max_index = list(C).index(max(C))
    
    A = C*np.ones((m,m))
    
    np.savetxt('concentration', C)
    
    typecounts = []
    for i in range(m):
        
        # most prominent species
        parent = names[max_index]
        
        try:
            
            a = bulk(parent,cubic=True).cell[0,0]
            c = bulk(parent,cubic=True).cell[2,2]
            
            cell = bulk(parent, a=a, crystalstructure=shapes[i], cubic=True)
        
        except:
            
            a = bulk(parent,orthorhombic=True).cell[0,0]
            c = bulk(parent,orthorhombic=True).cell[2,2]
            
            cell = bulk(parent, a=a, c=c, crystalstructure=shapes[i], orthorhombic=True)
            
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
        
        supercell = cell*dims
        
        N = len(supercell.get_chemical_symbols())
        b = C*N*np.ones(m,)
        
        x = gmres(A.T,b)[0]
        x = np.array([int(j) for j in x])
        
        typecount = []
        for j in range(m):
            typecount.append(int(x[j]))
        if sum(typecount) != N:
            add = N-sum(x)
            typecount[i] += add
        
        typecounts.append(typecount)
        
        
        ase.io.vasp.write_vasp(f"POSCAR{i + 1}",
                               supercell,
                               label='Cell {}'.format(i + 1), direct=True, sort=True)

        with open(f"POSCAR{i + 1}") as f:
            lines = f.readlines()

        lines[5] = (' '.join(names) + '\n')
        lines[6] = ' '.join([str(x) for x in typecount]) + "\n"
        with open(f"POSCAR{i + 1}", 'w') as f:
            f.writelines(lines)
    
    # Perform iterative swaps
    # num swaps = a chance to operate on all cell and each species
    # at least one time (1*) increase if want to operate more than
    # 1 time, i.e., 2x = (2*) etc
    num_swaps = int(5*(m*m)+m)
    
    while True:

        for _ in range(num_swaps):
            cell1, cell2 = np.random.choice(m, 2, replace=False)
            type1, type2 = np.random.choice(len(names), 2, replace=False)

            # Swap one species with another in one cell
            if typecounts[cell1][type1] > 0 and typecounts[cell2][type2] > 0:
                typecounts[cell1][type1] -= 1
                typecounts[cell1][type2] += 1

            # Reverse the swap in another cell
            if typecounts[cell2][type2] > 0:
                typecounts[cell2][type2] -= 1
                typecounts[cell2][type1] += 1
        
        # Write POSCAR files
        for i in range(m):

            with open(f"POSCAR{i + 1}") as f:
                lines = f.readlines()

            lines[6] = ' '.join([str(x) for x in typecounts[i]]) + "\n"
            with open(f"POSCAR{i + 1}", 'w') as f:
                f.writelines(lines)
        
        X = fun.build_X(names)
        F = fun.calc_f(X,C)
        
        if fun.mofac_test(F) == 1:
            singular = 0
            break
        print('Avoiding Singular State')
    
        
    # Write POSCAR files
    for i in range(m):

        with open(f"POSCAR{i + 1}") as f:
            lines = f.readlines()

        lines[6] = ' '.join([str(x) for x in typecounts[i]]) + "\n"
        with open(f"POSCAR{i + 1}", 'w') as f:
            f.writelines(lines)
    
    if genpot:
        fun.potcar(names,opsys)
