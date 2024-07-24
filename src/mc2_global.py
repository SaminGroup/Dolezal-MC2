import copy
import numpy as np
import src.myfuncs as fun
from random import uniform

def mc2(step_limit, C_set, T, w, param, intra, supcomp_phrase):
    """
    performs the overall mc2 loop

    step_limit --> max steps
    Ncells --> the number of cells we will initialize
    T --> simulation temperature
    param --> will either be 'begin' or 'continue'
    intra --> True or False, to perform intraswaps
    supcomp_phrase --> how vasp is called on supercomputer
    """
    C = C_set
    m = len(C)
    cells, names, N = fun.create_cells(m)
    

    if param == "begin":

        E_list, V_list = fun.initial_vasp_run(m,names,cells,supcomp_phrase) # returns two lists, each of len() = m
        cells = fun.global_accepted_state_update(cells,m) # update cells after the
        # initial vasp relaxations

        C = fun.concentrations(cells, N) # concentration of each species: constant
        X = fun.build_X(cells) # m x m matrix
        f0 = fun.calc_f(X,C) # m x 1 vector

        E = [] # these lists will only ever be len=2, the initial and final state
        V = []
        E_data = []
        V_data = []
        for i in range(m):
            E.append([])
            V.append([])
            E_data.append([])
            V_data.append([])
        for i in range(m):
            E[i].append(E_list[i])
            E[i].append(0)
            E_data[i].append(E_list[i])
            V[i].append(V_list[i])
            V[i].append(0)
            V_data[i].append(V_list[i])

        xdata = []
        frac_data = []
        for i in range(m):
            frac_data.append([])
            xdata.append([])
        for i in range(m):
            frac_data[i].append(f0[i])
        for i in range(m):
            for j in range(m):
                xdata[i].append([])

        continue_step = 0
        xdata = fun.atomic_percents(X,xdata)
        fun.data_text([frac_data, E_data, V_data], continue_step ,xdata)
    #--------------------------------
    # Data import
    #--------------------------------
    else:

        C = np.loadtxt("data/concentration")
        X = fun.build_X(cells) # m x m matrix
        f0 = fun.calc_f(X,C) # m x 1 vector

        frac_data, E_data, V_data, stepcount, xdata = fun.import_data(m)
        E = [] # these lists will only ever be len=2, the initial and final state
        V = []
        for i in range(m):
            E.append([])
            V.append([])
        for i in range(m):
            E[i].append(E_data[i][-1])
            E[i].append(0)
            V[i].append(V_data[i][-1])
            V[i].append(0)

        continue_step = stepcount # set the step number at the last completed step


    ev_data = [E_data, V_data]
    Flist = [f0, 0]
    states = [cells, 0]

    step = continue_step+1

    while step != step_limit+1:

        initial_state = states[0]

        did_global = 0

        temp_state = copy.deepcopy(initial_state)

        flip = 1
        intra_test = 1

        if intra:

            flip = int(uniform(0,2))

        if flip == 1:

            new_state, r, X_new, f, did_global = fun.flip_and_lever(temp_state,C)
            '''This logic prevents single-flip attempts from remaining
            stuck between singular states -- pushes to global flips instead'''
            if did_global == 1:
                new_state, X_new, f = fun.global_flip(temp_state,C)
                fun.global_potcar_update(m,names)

            else:
                fun.local_potcar_update(m, names, r)

            Flist[1] = f

        else:
            new_state, r, X_new, intra_test = fun.intraswap(temp_state)
            # only update if the move is possible
            if intra_test == 1:
                fun.local_potcar_update(m, names, r)

            f = Flist[0]
            Flist[1] = f # molar frac does not change on intraswap

        if intra_test == 1:

            states[1] = new_state
            #-----------------------------------------------
            # Step 3: run vasp on the new state
            #         and record E,V
            #-----------------------------------------------
            if did_global == 0:
                new_E, new_V = fun.vasp_run(r,supcomp_phrase)
                E[r][1] = new_E # new energy
                V[r][1] = new_V # new volume
                for i in range(m):
                    if i != r:
                        E[i][1] = E[i][0] # final = initial for unselected cells
                        V[i][1] = V[i][0] # final = initial ''
            else:
                new_E, new_V = fun.global_vasp_run(m,supcomp_phrase,step)
                for i in range(m):
                    E[i][1] = new_E[i] # final = initial for unselected cells
                    V[i][1] = new_V[i] # final = initial ''
            #-----------------------------------------------
            # Step 4: calculate the change dH,dV for the acceptance
            #-----------------------------------------------
            dG, cost = fun.enthalpy_and_volume(E, V, X, X_new, Flist, m, N, T)
            accept = fun.acceptance(dG,cost,w)


        else:
            accept = 2



        if accept == 0:
            for i in range(m):
                E_data[i].append(E_data[i][-1])
                V_data[i].append(V_data[i][-1])
                frac_data[i].append(frac_data[i][-1])
                for j in range(len(xdata[i])):
                    xdata[i][j].append(xdata[i][j][-1])

            fun.data_text([frac_data, E_data, V_data], step, xdata)
            step += 1

        elif accept == 1:

            X = X_new
            xdata = fun.atomic_percents(X,xdata)
            Flist[0] = f

            if did_global == 0:
                fun.update_poscar(r,new_state,names) # save contcar as poscar{} file
                E[r][0] = new_E
                V[r][0] = new_V
                newvars = [new_E, new_V]
                dataset = fun.data_update(r, f, frac_data, newvars, ev_data)
                new_state = fun.accepted_state_update(new_state,r,m)

            else:
                fun.global_poscar_update(m,new_state,names)
                for i in range(m):
                    E[i][0] = new_E[i]
                    V[i][0] = new_V[i]
                dataset = fun.global_data_update(f,frac_data,E_data,new_E,V_data,new_V)
                new_state = fun.global_accepted_state_update(new_state,m)

            states[0] = new_state
            fun.data_text(dataset,step,xdata)
            step += 1
