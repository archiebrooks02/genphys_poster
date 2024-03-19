import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = np.genfromtxt("/Users/archiebrooks/Documents/Uni/Group-Project/genphys_poster/attenuation_data.txt", skip_header=1, delimiter=",")

print(data)

energy = data[:,0]
air = data[:,1]
iron = data[:,2]
lead = data[:,3]

def fitting_fn(x,a,b,c):
    return((a*(x**b))+c)

air_coeff_fit, air_cov_fit = curve_fit(fitting_fn, energy, air, [8e-5,-0.492,0])
iron_coeff_fit, iron_cov_fit = curve_fit(fitting_fn, energy, iron, [0.1,-2,-0.23])
x = np.linspace(0.4,10,1000000)
y_air = fitting_fn(x,air_coeff_fit[0], air_coeff_fit[1],air_coeff_fit[2])
y_iron = fitting_fn(x,iron_coeff_fit[0],iron_coeff_fit[1],iron_coeff_fit[2])

plt.figure()
plt.plot(energy,air,"r.")
plt.plot(x,y_air)
plt.show()


plt.figure()
plt.plot(energy, iron, "r.")
plt.plot(x,y_iron)
plt.show()

print(air_coeff_fit)
print(iron_coeff_fit)
