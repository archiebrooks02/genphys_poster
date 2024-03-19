import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = np.genfromtxt("/Users/archiebrooks/Documents/Uni/Group-Project/genphys_poster/attenuation_data.txt", skip_header=1, delimiter=",")

print(data)

energy = data[:,0]
air = data[:,1]
iron = data[:,2]
lead = data[:,3]

def fitting_fn(x,a,b):
    return(a*(x**(b)))

coeff_fit, cov_fit = curve_fit(fitting_fn, energy, air, [8e-5,-0.492])
x = np.linspace(0.4,10,1000000)
y = coeff_fit[0]*(x**coeff_fit[1])

plt.figure()
plt.plot(energy,air,"r.")
plt.plot(x,y)
plt.show()

