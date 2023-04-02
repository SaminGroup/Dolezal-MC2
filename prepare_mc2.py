import os
import numpy as np
import PySimpleGUI as sg
from random import uniform
from scipy.sparse.linalg import gmres
from src.generate_simcells import poscar

sg.theme("Black")
sg.set_options(font=("Arial", 16),border_width=3)
# Define the GUI layout
layout = [
    [sg.Button('Begin New Run')],
    [sg.Button('Continue Existing Run')],
    [sg.Button('Exit')]
]

begin_layout = [
    [sg.Checkbox('Generate Simcells',key='simcells'), sg.Checkbox('Perform Intraswaps', key='intraswap')],
    [sg.Text('MC Step Limit'), sg.Input(size=(5,1)), sg.Text('Temperature (K)'), sg.Input(size=(5,1))],
    [sg.Text('System command to execute VASP (i.e., vasp or mpirun vasp)'), sg.Input(size=(10,1))],
    [sg.Button('Run Script'), sg.Button('Exit')]
]

need_cell_layout = [
        [sg.Text('Please provide number of cells'), sg.Input(size=(5,1))],
        [sg.Button('Save'), sg.Button('Exit')]
]

simcell_layout = [
        [sg.Checkbox('Windows'), sg.Checkbox('Linux/Unix')],
        [sg.Text('System to be Explored (i.e., Au Pt): '), sg.InputText(size=(10,1))],
        [sg.Text('System Concentrations (i.e., 0.5 0.5): '), sg.InputText(size=(15,1))],
        [sg.Checkbox('32 Atom Cell', default=False), sg.Checkbox('48 Atom Cell', default=False), sg.Checkbox('Generate POTCARs', default=False)],
        [sg.Button('Run Script'), sg.Button('Exit')]
]

continue_layout = [
        [sg.Text('Update MC step limit'), sg.InputText(size=(5,1))],
        [sg.Button('Run Script'), sg.Button('Exit')]
]

# Create the GUI window
window = sg.Window('Multi-cell Monte Carlo (MC)2', layout, resizable=True)

while True:
    event, values = window.read()

    if event == sg.WINDOW_CLOSED or event == 'Exit':
        break

    if event == "Begin New Run":

        if not os.path.exists('data'):
            os.makedirs('data')

        begin_window = sg.Window("Begin Run Options", begin_layout)
        while True:
            begin_event, begin_values = begin_window.read()

            if begin_event == sg.WINDOW_CLOSED or begin_event == 'Exit':
                break

            if begin_values['simcells']:
                if not os.path.exists('potcar'):
                    os.makedirs('potcar')
                simcell_window = sg.Window("Simcell Generator", simcell_layout)

                while True:
                    sim_event, sim_values = simcell_window.read()

                    if sim_event == sg.WINDOW_CLOSED or sim_event == 'Exit':
                        break

                    if sim_event == 'Run Script':
                        if sim_values[0] == True:
                            opsys = 'windows'
                        else:
                            opsys = 'linux'

                        names = sim_values[2].split()

                        m = len(names)


                        C = np.asarray(str(sim_values[3]).split(),dtype=float)

                        if sim_values[4] == True:
                            N = 32
                        else:
                            N = 48

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

                        genpot = sim_values[6]

                        poscar(N,m,typecount,names,opsys,genpot)

                        sg.Popup('Simcells are generated!')
                        break

                simcell_window.close()

            else:
                need_window = sg.Window("Cell Number", need_cell_layout)
                while True:
                    need_event, need_values = need_window.read()

                    if need_event == sg.WINDOW_CLOSED or need_event == 'Exit':
                        break

                    if need_event == 'Save':
                        m = need_values[0]

                        break

                need_window.close()

            if begin_event == 'Run Script':
                with open("src/begin_mc2.py") as f:
                    lines = f.readlines()
                    f.close()


                lines[1] = 'step_limit = {}\n'.format(begin_values[0])
                lines[2] = 'Ncells = {}\n'.format(m)
                lines[3] = 'T = {}\n'.format(begin_values[1])
                lines[4] = "param = 'begin'\n"
                lines[5] = "supcomp_phrase = '{}'\n".format(begin_values[2])
                lines[6] = "intra = {}\n".format(begin_values['intraswap'])

                with open('begin.py', 'w') as f:
                    f.writelines(lines)
                    f.close()

                break

        begin_window.close()

    if event == "Continue Existing Run":
        con_window = sg.Window("Begin Run Options", continue_layout)
        while True:
            con_event, con_values = con_window.read()

            if con_event == sg.WINDOW_CLOSED or con_event == 'Exit':
                break

            if con_event == 'Run Script':
                with open('begin.py') as f:
                    lines = f.readlines()
                    f.close()

                lines[1] = 'step_limit = {}\n'.format(con_values[0])
                lines[4] = "param = 'continue'\n"

                with open('begin.py', 'w') as f:
                    f.writelines(lines)
                    f.close()

                break

        con_window.close()

    sg.Popup('MC2 is ready to go!')
    break

window.close()
