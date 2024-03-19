import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from scipy import constants
import scipy as sc
from mpl_toolkits.mplot3d import Axes3D

#initialise the matplotlib.pyplot parameteres
params = {'axes.grid' : True,
          'axes.axisbelow' : True,
          'figure.figsize': (10,7),
          'axes.titlesize' : 16,
          'axes.labelsize' : 12,
          'font.size': 10,}

pylab.rcParams.update(params)
