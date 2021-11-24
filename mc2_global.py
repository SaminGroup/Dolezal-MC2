import os
import copy
import numpy as np
import myfuncs as fun
from random import uniform


def mc2(step_limit, Ncells, T, param, intra, supcomp_phrase):
    """
    performs the overall mc2 loop

    step_limit --> max steps
    Ncells --> the number of cells we will initialize
    T --> simulation temperature
    param --> will either be 'begin' or 'continue'
    intra --> True or False, to perform intraswaps
    supcomp_phrase --> how vasp is called on supercomputer
    """

    cells, names, N = fun.create_cells(Ncells)
    m = Ncells

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

    global_step = continue_step

    singular = 0 ; intra_test = 1 ; failed_count = 0 ; did_global = 0
    for step in range(continue_step+1,step_limit+1):
        initial_state = states[0]

        temp_state = copy.deepcopy(initial_state)

        flip = 1
        intra_test = 1
        if intra:
            flip = int(uniform(0,2))

        if flip == 1:


            if did_global == 0:
                new_state, r, X_new = fun.flip_and_lever(temp_state,singular)
                fun.local_potcar_update(m, names, r)
            else:
                new_state, X_new = fun.global_flip(temp_state)
                fun.global_potcar_update(m,names)

        else:
            new_state, r, X_new, intra_test = fun.intraswap(temp_state)
            fun.local_potcar_update(m, names, r)

        if intra_test == 1:


            states[1] = new_state

            f = fun.calc_f(X_new,C)
            Flist[1] = f


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
                new_E, new_V = fun.global_vasp_run(m,supcomp_phrase)
                for i in range(m):
                    E[i][1] = new_E[i] # final = initial for unselected cells
                    V[i][1] = new_V[i] # final = initial ''
            #-----------------------------------------------
            # Step 4: calculate the change dH,dV for the acceptance
            #-----------------------------------------------
            dH,dV,dX = fun.enthalpy_and_volume(E, V, X, X_new, Flist, m, N, T)
            accept = fun.acceptance(dH, dV, dX)



        else:
            accept = 0
            os.system("rm POSCAR")


        if accept == 0:
            for i in range(m):
                E_data[i].append(E_data[i][-1])
                V_data[i].append(V_data[i][-1])
                frac_data[i].append(frac_data[i][-1])
                for j in range(len(xdata[i])):
                    xdata[i][j].append(xdata[i][j][-1])

            fun.data_text([frac_data, E_data, V_data], step, xdata)

        else:

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
                if step == global_step + 10:
                    did_global = 1
            else:
                fun.global_poscar_update(m,new_state,names)
                for i in range(m):
                    E[i][0] = new_E[i]
                    V[i][0] = new_V[i]
                dataset = fun.global_data_update(f,frac_data,E_data,new_E,V_data,new_V)
                new_state = fun.global_accepted_state_update(new_state,m)
                global_step = step
                did_global = 0

            states[0] = new_state
            fun.data_text(dataset,step,xdata)
