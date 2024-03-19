import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from scipy import constants
import scipy as sc
from mpl_toolkits.mplot3d import Axes3D

def kT_eV(kiloton):
    return kiloton*2.611e+31

#Davy Crockett bomb was 0.1 kT

davy_energy = kT_eV(0.01) #energy in eV
davy_heat = 0.35*davy_energy
davy_gamma = 0.065*davy_energy

#initialise the matplotlib.pyplot parameteres
params = {'axes.grid' : True,
          'axes.axisbelow' : True,
          'figure.figsize': (10,7),
          'axes.titlesize' : 16,
          'axes.labelsize' : 12,
          'font.size': 10,}

pylab.rcParams.update(params)

gammas = np.random.uniform(0,0.75,1)

def loop(threshold,lower,upper):
    total_energy = 0
    while total_energy<threshold:
        energy = np.random.uniform(lower,upper,10000)
        total_energy += np.sum(energy)
        np.hstack((gammas,energy))
        print(total_energy)
    return gammas

def get_energy(bomb_energy):
    loop(bomb_energy*0.7,0,0.75)
    loop(bomb_energy*0.1,0.75,2)
    loop(bomb_energy*0.09,2,4.5)
    loop(bomb_energy*0.09,4.5,8)
    loop(bomb_energy*0.02,8,12)
   

def get_atten(energy):
    return (8e-5*(energy**(-0.492)))


def attenuation(n_gamma, x):
    n_transported=0
    for i in range(0,n_gamma):
        n_transported += np.exp(-get_atten(get_energy())*x)
    return(n_transported)

get_energy(1000000)

print(davy_gamma)
    