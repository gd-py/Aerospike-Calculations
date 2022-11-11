import numpy as np
import matplotlib.pyplot as plt
from caeroc.formulae.isentropic import Isentropic

## Exponential Atmospheric Pressure Model # Ref: NASA exponential Atmosphere
# Insert height in feet. Returns temperature in deg Celcius.
def tAlt(h):
    return 15.04-0.00649*h

# Insert height in feet. Returns pressure in K-Pa
def pAlt(h):
    h *= 0.3048
    return 101.29*((tAlt(h)+273.1)/288.08)**5.256

## Conversion factor
# (ft.lbf/(lbm.R)) to (kJ/(kg.K)), This can be used for both specific heat conversion and specific gas constant conversion as R = Cp-Cv
c = 0.3048*4.4482/(453.59237*0.555556)

## Propellant Properties, look Ref. 1 for detail (The values here are taken for Chamber Temperature of 1100k)
M = 41.98 # g/mol (or kg/kmol as in Ansys)
Cp_ = 64.87 # J/mol-K
Cp = Cp_/M # kJ/kg-K, multiply this value by 1000 before inserting into Ansys (Ansys accepts value in J/kg-K)
k = 1.0468 # This ratio of gas constants (Cp/Cv) is used for condensed phase. See Ref. 4
print("Cp is: ", Cp*1000)

## Constants
R = 287 # SI units
m = 1.6 # kg, mass of propellent
g = 9.81

## Combustion Chamber Conditions
P_0 = 3447380 # Pa, (= 500 PSI)
T_0 = 1100 # Kelvin
rho_0 = P_0/(R*T_0)

## Ambient Conditions at altitude of 10k feet
h = 10000
P_a = pAlt(h) # Pa, (= 10,000ft pressure altitude)
T_a = 273 + tAlt(h) # Kelvin, (= -4.8 deg Celcius)

gamma = k

## These values are taken from simulation of Ref. 2
A_t = 0.00029532 # m^2, Throat Area
R_e = 0.0237 # m, Exit radius
A_e = np.pi*R_e**2 # m^2, Exit area
Exp_ratio = 6 # Expansion Ratio

isen = Isentropic(gamma=gamma)

V_t = np.sqrt((2*gamma*R*T_0)/(gamma+1)) # Velocity at Throat
Me = [i for i in isen.M(A_Astar=Exp_ratio) if i>1][0] # Exit Mach Number

rhot_rho0 = isen.rho_rho0(M=1) # ratio of density at throat to density at chamber # At throat, M = 1
rhot = rho_0 * rhot_rho0 # density at throat

m_dot = rhot*A_t*V_t # mass flow rate

# Conditions at nozzle exit
Pe_P0 = isen.p_p0(M=Me)
P_e = Pe_P0*P_0 # Pressure at exit

Te_T0 = isen.T_T0(M=Me)
Te = Te_T0*T_0 # Temperature at exit

a_e = np.sqrt(gamma*R*Te) # Sound speed at exit
V_e = Me*a_e

F = m_dot*V_e + (P_e-P_a)*A_e # Thrust given by nozzle

Veq = V_e + (P_e - P_a)*A_e/m_dot # Equivalent velocity, look Ref. 3 for more
I = m*Veq # Total Impulse
Isp = Veq/g

print("Thrust: ", F)
print("Total Impulse: ", I)
print("Specific Impulse: ", Isp)

print(V_t)



## References
# 1. ../../References/KNSU Propellent/Propellent Properties Table.pdf
# 2. ../aerospike-nozzle-design-gui-master
# 3. ../../References/engine/Specific Impulse.pdf
# 4. ../../References/KNSU Propellant/Propellant Properties Calculation.pdf