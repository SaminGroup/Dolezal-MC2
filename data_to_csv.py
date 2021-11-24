import pandas as pd
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
import seaborn as sns
from plot_aupt import plot_aupt

def data_to_csv():

    """

    19 Nov 2021: changed how xdata is recorded and imported

    """

    species = input('Please list the species in order as found in the POSCAR: ').split()
    Temp = input("Simulation Temperature: ")
    save_title = "".join(species)
    title = " ".join(species)

    dir = "plots_{}_{}/".format(save_title,Temp)

    filenames = [dir+'data/mofac', dir+'data/energy', dir+'data/volume']
    dataset = [[],[],[]]
    ##----------------------------------
    ## Step 1: read in the datasets from the files and save them to a list
    ##----------------------------------
    for i in range(len(filenames)):
        dataset[i] = np.loadtxt(filenames[i],dtype=float)
    dataset = np.array(dataset)
    m = len(dataset[0])
    ##---------------------------------
    ## Energy/Volume per Atom
    ##---------------------------------
    counts = np.array(input("How many atoms per cell? ").split(),dtype=int)
    for k in range(m):
        dataset[1,k] = dataset[1,k]/counts[k]  # energy per atom
        dataset[2,k] = dataset[2,k]/counts[k]  # volume per atom

    deltaE = np.zeros((m,len(dataset[0][0]))) # phase-weight x energy per atom
    deltaV = np.zeros((m,len(dataset[0][0]))) # phase-weight x volume per atom

    for k in range(m):
        deltaE[k,:] = dataset[0,k]*dataset[1,k] # phase-weight x energy per atom
        deltaV[k,:] = dataset[0,k]*dataset[2,k] # phase-weight x volume per atom
    deltaE = np.array([sum(deltaE[:,j]) for j in range(deltaE.shape[1])])
    deltaV = np.array([sum(deltaV[:,j]) for j in range(deltaV.shape[1])])
    TotalE = deltaE
    TotalV = deltaV
    deltaE = (abs(deltaE[0]) - abs(deltaE))
    deltaV = (abs(deltaV[0]) - abs(deltaV))
    ##----------------------------------
    ## Step 2: read in stepcount and At.% data
    ##----------------------------------
    stepcount = np.loadtxt(dir+"data/stepcount",dtype=int)

    xdata = np.loadtxt(dir+"data/xdata").reshape(m,m,stepcount+1)
    ##-----------------------------------
    ## Step 3: create the pandas dataframe and save to csv file
    ##-----------------------------------
    names = []
    for i in range(m):
        names.append([])
    for i in range(m):
        names[i].append('Molar Fraction {}'.format(i+1))
        names[i].append('Energy {}'.format(i+1))
        names[i].append('Volume {}'.format(i+1))
    ##------------------------------------
    ## Step 3a creating the columns of data
    ##------------------------------------
    df = pd.DataFrame({'Steps':np.arange(stepcount+1)})
    for i in range(len(dataset)):
        for j in range(m):
            df[names[j][i]] = pd.Series(dataset[i][j])
    for i in range(m):
        for j in range(len(xdata[i])):
            df['At.% Cell {}, Species {}'.format(i+1,j+1)] = pd.Series(xdata[i][j])
    ##-----------------------------------
    df["Total E"] = pd.Series(TotalE)
    df["Delta E"] = pd.Series(deltaE)
    df["Total V"] = pd.Series(TotalV)
    df['Date'] = pd.Series(pd.to_datetime('today').normalize())
    df.to_csv(dir+'data/dataset.csv',index=False)
    if title == "Au Pt":
        plot_aupt(m,dir,Temp)
    else:
        plot_data(m,dir,species,Temp)
    ##------------------------------------

def plot_data(m,dir,species,Temp):

    title = " ".join(species)

    df = pd.read_csv(dir+'data/dataset.csv')
    mofac = []
    xdata = []
    for i in range(m):
        mofac.append('Molar Fraction {}'.format(i+1))
        xdata.append([])
        df[mofac[i]] = df[mofac[i]]*100
        for j in range(m):
            xdata[i].append('At.% Cell {}, Species {}'.format(i+1,j+1))
            df['At.% Cell {}, Species {}'.format(i+1,j+1)] = df['At.% Cell {}, Species {}'.format(i+1,j+1)]*100
    ##--------------------------------------------------------------
    ## when we are dealing with more than 2 species we will need to sum
    ## the m-1 values and subtract that sum from 100 to produce the mth
    ## species percent
    ##--------------------------------------------------------------
    labels = []
    en_list = []

    for i in range(m):
        labels.append('Cell {}'.format(i+1))
        en_list.append("Energy {}".format(i+1))

    sns.set_theme()
    sns.set_style('ticks')
    plt.rc('font', family='serif')

    fig,ax = plt.subplots()

    df.plot(x = 'Steps', y="Total E", legend=None)
    plt.ylabel('E (eV/atom)')
    plt.grid(which="major",alpha=0.6)
    plt.grid(which="minor",alpha=0.3)
    plt.tight_layout()
    plt.savefig(dir+'total-energy.png', dpi=400,bbox_inches="tight")
    plt.close()

    equil_step = int(input("Equilibrium Begins: "))
    cell_shape = input("Cell Shapes: ").split()

    for i in range(m):
        fac,fbins = np.histogram(df[mofac[i]],25);fbins=fbins[:-1]
        bar = plt.bar(fbins,fac,label="Cell {}".format(i+1),edgecolor="none")
    plt.xticks(np.arange(0,110,10))
    plt.xlabel("Molar Fraction (%)")
    plt.ylabel("Counts")
    plt.legend(frameon=False)
    plt.savefig(dir+"molar-frac-histogram.png",dpi=400,bbox_inches='tight')


    fig,ax = plt.subplots()
    df.plot(x = 'Steps', y= mofac, legend=None)
    plt.ylabel("Molar Fraction (%)")
    plt.axvline(equil_step,color="k")
    plt.legend(cell_shape,frameon=False)
    plt.ylim(0,100)
    plt.grid(which="major",alpha=0.6)
    plt.grid(which="minor",alpha=0.3)
    plt.tight_layout()
    plt.savefig(dir+'molar-frac.png', dpi=400,bbox_inches="tight")
    plt.close()

    fig,ax = plt.subplots()

    df.plot(x = 'Steps', y=en_list,legend=None)
    plt.ylabel('E (eV/atom)')
    plt.legend(cell_shape,frameon=False)
    plt.grid(which="major",alpha=0.6)
    plt.grid(which="minor",alpha=0.3)
    plt.tight_layout()
    plt.savefig(dir+'energy.png', dpi=400,bbox_inches="tight")
    plt.close()

    fig,ax = plt.subplots()

    df.plot(x = 'Steps', y="Total V", legend=None)
    plt.ylabel("V ($\AA^3$/atom)")
    plt.grid(which="major",alpha=0.6)
    plt.grid(which="minor",alpha=0.3)
    plt.tight_layout()
    plt.savefig(dir+'total-volume.png',dpi=400,bbox_inches="tight")
    plt.close()
    ##--------------------------------
    ## Sort the species for the atomic
    ## percents plots
    ##--------------------------------
    atomicpercents = np.loadtxt(dir+"data/xdata").reshape(m,m,df["Steps"].iloc[-1]+1)
    atomicpercents *= 100
    ##----------------------------------
    ## Generate pie chart of surviving
    ## phases and bar plot of final
    ## atomic percent
    ##----------------------------------
    phase_concentrations = np.zeros(m)
    X = []

    for i in range(m):
        phase_concentrations[i] = df[mofac[i]].iloc[-1]
        X.append('Cell {}'.format(i+1))

    final_at_percent = np.zeros((m,m))
    for i in range(m):
        for j in range(m):
            final_at_percent[i,j] = atomicpercents[i][j][-1]

    for i in range(len(phase_concentrations)):
        if abs(phase_concentrations[i]) < 1e-10:
            phase_concentrations[i] = 0.0


    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10,5))

    error = np.zeros((m,))
    for i in range(m):
        error[i] = np.std(df[mofac[i]][equil_step-1:])

    X_axis = np.arange(m)
    w = (1/(m+0.75))

    bar0 = axes[0].bar(X_axis, phase_concentrations, w, edgecolor='black')
    axes[0].set_xticks(X_axis)
    axes[0].set_xticklabels(cell_shape)
    axes[0].set_ylabel('Molar Fraction (%)')
    axes[0].set_ylim(0,phase_concentrations.max() + 10)

    for bar,anerror in zip(axes[0].patches,error):
        bar_value = bar.get_height()
        text = '{:.0f}$\\pm${:.0f}%'.format(bar_value,anerror)
        text_x = bar.get_x() + bar.get_width() / 1.55
        text_y = (bar.get_y() + bar_value)

        axes[0].text(text_x, text_y, text, ha='center', va='bottom',
                        size=8)

    axes[0].set_xticks(X_axis)
    axes[0].set_xticklabels(cell_shape)
    axes[0].set_ylabel('Molar Fraction (%)')
    axes[0].set_ylim(0,phase_concentrations.max() + 10)

    print(final_at_percent)
    for i in range(m):
        bar1 = axes[1].bar(X_axis+((i)*(w)), final_at_percent[:,i], w, edgecolor='black')
        axes[1].set_xticks(X_axis+(i*w/2))


    axes[1].set_xticklabels(cell_shape)
    axes[1].set_ylabel('Final At.%')
    axes[1].set_ylim(0,final_at_percent.max()+20)
    ## label each bar in the barplot
    for bar in axes[1].patches:
        bar_value = bar.get_height()
        text = '{:.0f}%'.format(bar_value)
        text_x = bar.get_x() + bar.get_width() / 2
        text_y = (bar.get_y() + bar_value)
        axes[1].text(text_x, text_y, text, ha='center', va='bottom',
                        size=8)

    axes[1].legend(species,frameon=False, ncol = m, loc='upper right')
    plt.savefig(dir+"phase_percents.png",dpi=400,bbox_inches="tight")
    plt.close()
    ##----------------------------------
    ## Here is where we plot the atomic
    ## percent in each cell over the
    ## entire MC2 run
    ##----------------------------------
    fig, axes = plt.subplots(nrows=m, ncols=1, sharex=True, sharey=True)
    fig.suptitle(title)
    for i in range(m):
        average = 1 - np.average(atomicpercents[i])
        axes[i].stackplot(df['Steps'], atomicpercents[i],
                    labels=species,
                        edgecolor='black',linewidth=1.25)


        textob = X[i]
        ob = offsetbox.AnchoredText(textob, loc=4, frameon=False, prop=dict(fontsize=10,color="white"))
        axes[i].add_artist(ob)
        plt.subplots_adjust(wspace=0.0, hspace=0.0)
        fig.text(0.03, 0.5, 'Atomic Concentration %', va='center', rotation='vertical')


    axes[0].legend(loc='upper left', frameon=False,ncol=m,fontsize=10, labelcolor='w')
    plt.ylim(0,100)
    plt.xlim(0,df['Steps'].max())
    plt.xlabel('Monte Carlo Steps')
    fig.text(0.77,0.89,"T = {} K".format(Temp),fontsize=10)
    plt.savefig(dir+'atomicpercent.png', dpi=400,bbox_inches="tight")
    plt.close()
    ##----------------------------------

data_to_csv()
