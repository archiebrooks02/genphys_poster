import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = np.genfromtxt("/Users/archiebrooks/Documents/Uni/Group-Project/genphys_poster/attenuation_data.txt", skip_header=1, delimiter=",")

energy_1 = data[:4,0]
lead_1 = data[:4,3]

energy_2 = data[3:,0]
lead_2 = data[3:,3]

lead_2 = np.hstack((lead_2,0.628))
energy_2 = np.hstack((energy_2,15))


def fitting_fn(x,a,b,c):
    return((a*(x**b))+c)

def fitting_fn_2(x,a,b,c,d):
    return ((a*(x**3)+b*(x**2)+c*x+d))


lead_coeff_fit_1, lead_cov_fit_1 = curve_fit(fitting_fn, energy_1, lead_1, [0.4,-1.6,0.4])
lead_coeff_fit_2, lead_cov_fit_2 = curve_fit(fitting_fn_2, energy_2, lead_2, [0.0006,-0.0171,0.1573,0.0862])
x_1 = np.linspace(0.4,3,10000)
x_2 = np.linspace(3,15,10000)

y_lead_1 = fitting_fn(x_1,lead_coeff_fit_1[0],lead_coeff_fit_1[1],lead_coeff_fit_1[2])
y_lead_2 = fitting_fn_2(x_2,lead_coeff_fit_2[0],lead_coeff_fit_2[1],lead_coeff_fit_2[2],lead_coeff_fit_2[3])

plt.figure()
plt.plot(energy_1, lead_1, "r.")
plt.plot(energy_2, lead_2, "r.")
plt.plot(x_1,y_lead_1)
plt.plot(x_2,y_lead_2)
plt.show()

print(lead_coeff_fit_1)
print(lead_coeff_fit_2)
