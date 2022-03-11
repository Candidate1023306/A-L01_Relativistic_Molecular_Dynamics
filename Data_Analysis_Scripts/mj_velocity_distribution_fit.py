import numpy
import cmath
import matplotlib
import math
from math import erf
from matplotlib import pyplot
from cmath import pi, sqrt 
from numpy import mgrid, power, exp, asarray, dot, trapz
import scipy
from scipy.optimize.zeros import RootResults
from scipy.special import kn
from scipy.optimize import curve_fit
from numpy.core.function_base import linspace
from numpy.lib.type_check import real



""" 
Constant variables 
"""
n_atoms = 500
boltz = 1.380649*(10**(-23))
timestep_to_time = 10**(-12) #in seconds, so converting from picoseconds.  
timestep = 1e-6 
c = 2.998*(10**8)
c_metal = 2.998e+6
box_size = 23.0
mass = 39.95 #g/mol
mass_to_kg = (10**(-3))/(6.022*10**23)
length_to_m = 10**(-10)
mass_kg = mass*mass_to_kg
rest_mass_energy = mass_kg*(c**2)
box_size_m = length_to_m*box_size
print(str(n_atoms) + " particles of mass: " + str(mass_kg) + " kg in a box of size: " + str(box_size_m) + "m per side")
mean_density = n_atoms*mass_kg/(box_size_m**3)
print("mean system density = " + str(mean_density) + " kg/m^3")


def Maxwell_Juttner(beta, theta, A): #theta is kT/mc^2. A fixes the normalisation that I can't calculate via algebra. 
    gamma = (1-beta**2)**(-0.5)
    return(A*((((gamma**5)*(beta**2))/(theta*kn(2, 1/theta)))*exp(-gamma/theta)))


filler_lines = 9
hist_file_01c = open("500_Atom_0_1c_Argon_Beta_Distribution.txt", "r")
raw_data_01c = []
for line in hist_file_01c:
    raw_data_01c.append(line.split())

hist_file_03c = open("500_Atom_0_3c_Argon_Beta_Distribution.txt", "r")
raw_data_03c = []
for line in hist_file_03c:
    raw_data_03c.append(line.split())

hist_file_05c = open("500_Atom_0_5c_Argon_Beta_Distribution.txt", "r")
raw_data_05c = []
for line in hist_file_05c:
    raw_data_05c.append(line.split())

bins_01c = []
vals_01c = []
bins_03c = []
vals_03c = []
bins_05c = []
vals_05c = []


for i in range(0, len(raw_data_01c)):
    bins_01c.append(float(raw_data_01c[i][0]))  
    vals_01c.append(float(raw_data_01c[i][1]))

for i in range(0, len(raw_data_03c)):
    bins_03c.append(float(raw_data_03c[i][0]))  
    vals_03c.append(float(raw_data_03c[i][1]))

for i in range(0, len(raw_data_05c)):
    bins_05c.append(float(raw_data_05c[i][0]))  
    vals_05c.append(float(raw_data_05c[i][1]))


bins_01c = asarray(bins_01c)
vals_01c= asarray(vals_01c)
bins_03c = asarray(bins_03c)
vals_03c= asarray(vals_03c)
bins_05c = asarray(bins_05c)
vals_05c= asarray(vals_05c)

"""
guessed_juttner = []
for i in range(0, len(bins)):
    guessed_juttner.append(Maxwell_Juttner(bins[i], 0.085, 0.04))
"""

mj_fit_vals_01c, mj_fit_covar_01c = curve_fit(Maxwell_Juttner, bins_01c, vals_01c, p0 = [0.085, 0.04])
fitted_juttner_01c = []

mj_fit_vals_03c, mj_fit_covar_03c = curve_fit(Maxwell_Juttner, bins_03c, vals_03c, p0 = [0.085, 0.04])
fitted_juttner_03c = []

mj_fit_vals_05c, mj_fit_covar_05c = curve_fit(Maxwell_Juttner, bins_05c, vals_05c, p0 = [0.085, 0.04])
fitted_juttner_05c = []

for i in range(0, len(bins_01c)):
    fitted_juttner_01c.append(Maxwell_Juttner(bins_01c[i], mj_fit_vals_01c[0], mj_fit_vals_01c[1]))

for i in range(0, len(bins_03c)):    
    fitted_juttner_03c.append(Maxwell_Juttner(bins_03c[i], mj_fit_vals_03c[0], mj_fit_vals_03c[1]))

for i in range(0, len(bins_05c)):
    fitted_juttner_05c.append(Maxwell_Juttner(bins_05c[i], mj_fit_vals_05c[0], mj_fit_vals_05c[1]))

print("Temperature from the 0.1c Fit = " + str(mj_fit_vals_01c[0]*rest_mass_energy/boltz) + " with error " + str(mj_fit_covar_01c[0][0]*rest_mass_energy/boltz))
print("----")
print("Temperature from the 0.3c Fit = " + str(mj_fit_vals_03c[0]*rest_mass_energy/boltz) + " with error " + str(mj_fit_covar_03c[0][0]*rest_mass_energy/boltz))
print("----")
print("Temperature from the 0.5c Fit = " + str(mj_fit_vals_05c[0]*rest_mass_energy/boltz) + " with error " + str(mj_fit_covar_05c[0][0]*rest_mass_energy/boltz))



pyplot.plot(bins_01c, vals_01c, "k.", label = "Measured 0.1c Distribution")
pyplot.plot(bins_01c, fitted_juttner_01c, "kx", label = "Fitted 0.1c Distribution")

pyplot.plot(bins_03c, vals_03c, "b.", label = "Measured 0.3c Distribution")
pyplot.plot(bins_03c, fitted_juttner_03c, "bx", label = "Fitted 0.3c Distribution")

pyplot.plot(bins_05c, vals_05c, "g.", label = "Measured 0.5c Distribution")
pyplot.plot(bins_05c, fitted_juttner_05c, "gx", label = "Fitted 0.5c Distribution")

pyplot.xlabel("Atom Beta")
pyplot.ylabel("F(beta)")
pyplot.title("Relativistic Beta Distributions")
pyplot.legend(loc = "upper right")
pyplot.show()
