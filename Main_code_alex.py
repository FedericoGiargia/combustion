import cantera as ct
import matplotlib.pyplot as plt
from pathlib import Path

#########################################################################
# This function computes the average mixture temperature. This is       # 
# based on iteratively computing the average Cp-value of the injected   #
# mixures. The total enthalpy defined by the total mass and Cp will     #
# give the mixture average temperature                                  #
#########################################################################


phi_rich = 5.0
phi_lean = 0.5
mech_file = "gri30.yaml"
fuel_comp_rich = {"H2": 1.0}
T_in = 675 #K
P_in = 12E5 # Pa        This is typically filled in on the spot in the code




def freeflamestuff(T,P,phi, fuel_comp, strain=False,SL=False,H2O=False):
    #########################################################################
    # This function solves the flame profile based on the gri30 mechanics.  # 
    # The model is computed by the Cantera software                         #
    #########################################################################

    To = T
    Po = P

    # Domain width in metres (this is an arbitrary assumed number, not to small to give the mixture the change to convert chemical energy to heat)
    width = 0.12

    gas = ct.Solution(mech_file)
    if strain== False:
        if H2O!=False:
            gas.set_equivalence_ratio(phi, fuel_comp, {"O2": 1.0, "N2": 3.76,"H2O":H2O})
        else:
            gas.set_equivalence_ratio(phi, fuel_comp, {"O2": 1.0, "N2": 3.76})
        gas.TP = To, Po  
        flame = ct.FreeFlame(gas, width=width)
    
    else:
        gas.set_equivalence_ratio(phi, fuel_comp, {"O2": 1.0, "N2": 3.76})
        gas.TP = To, Po 
        flame = ct.CounterflowPremixedFlame(gas, width=width)
        flame.reactants.T = T
        flame.reactants.X = gas.X
        flame.reactants.mdot =  0.5 * gas.density * strain * width 

    flame.transport_model = 'mixture-averaged'
    flame.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)

    flame.solve(loglevel=1, auto=True)
    return flame


# Question Part A:
# Question Part A laminar flame speed:
lst_flame_speeds = []
print("Methane flame speed analysis: ")
for phi in range(7,15):
    flame = freeflamestuff(295,101325,phi/10,{"CH4": 1.0}) #Note! For this the width manually has to be set to 30cm
    # Extract the spatial profiles as a SolutionArray to simplify plotting specific species
    lst_flame_speeds.append(flame.velocity[0])

for i in range(len(lst_flame_speeds)):
    print(f"\nmixture-averaged flame speed for phi= {i/10+0.7} = {lst_flame_speeds[i]} m/s\n")


#Part B A
flame = (freeflamestuff(675,12E5,5,fuel_comp_rich))
# Extract the spatial profiles as a SolutionArray to simplify plotting specific species
profile = flame.to_array()

#Part B B 
print("QUESTION B")
# space variable =water vapour
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(profile("H2O").Y[:,0]/profile("H2O").Y[-1,0], profile.T, ".-")
ax.set_xlabel("Water Vapour c [-]")
ax.set_ylabel("Temperature [K]")
ax.set_title("Flame Temperature Profile VS water vapour progress variable")
ax.grid()

plt.show()


fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(profile("H2O").Y[:,0]/profile("H2O").Y[-1,0], profile("O2").Y, label="O$_2$")
ax.plot(profile("H2O").Y[:,0]/profile("H2O").Y[-1,0], profile("H2").Y, label="H$_2$")
ax.plot(profile("H2O").Y[:,0]/profile("H2O").Y[-1,0], profile("NO").Y, label="NO")
ax.plot(profile("H2O").Y[:,0]/profile("H2O").Y[-1,0], profile("NO2").Y, label="NO$_2$")

ax.legend(loc="best")
ax.set_xlabel("Water Vapour c [-]")
ax.set_ylabel("Mass fraction [-]")
ax.set_title("Species Mass Fraction Profile VS water vapour progress variable")
ax.grid()
plt.show()


fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(profile("H2O").Y[:,0]/profile("H2O").Y[-1,0], profile("NO").Y, label="NO")
ax.plot(profile("H2O").Y[:,0]/profile("H2O").Y[-1,0], profile("NO2").Y, label="NO$_2$")

ax.legend(loc="best")
ax.set_xlabel("Water Vapour c [-]")
ax.set_ylabel("Mass fraction [-]")
ax.set_title("Species Mass Fraction Profile VS water vapour progress variable")
ax.grid()
plt.show()


# space variable =hydrogen
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot((profile("H2").Y[0,0]-profile("H2").Y[:,0])/(profile("H2").Y[0,0]-profile("H2").Y[-1,0]), profile.T, ".-")
ax.set_xlabel("Hydrogen c [-]")
ax.set_ylabel("Temperature [K]")
ax.set_title("Flame Temperature Profile VS Hydrogen progress variable")
ax.grid()

plt.show()


fig, ax = plt.subplots(figsize=(10, 6))
ax.plot((profile("H2").Y[0,0]-profile("H2").Y[:,0])/(profile("H2").Y[0,0]-profile("H2").Y[-1,0]), profile("O2").Y, label="O$_2$")
ax.plot((profile("H2").Y[0,0]-profile("H2").Y[:,0])/(profile("H2").Y[0,0]-profile("H2").Y[-1,0]), profile("H2").Y, label="H$_2$")
ax.plot((profile("H2").Y[0,0]-profile("H2").Y[:,0])/(profile("H2").Y[0,0]-profile("H2").Y[-1,0]), profile("NO").Y, label="NO")
ax.plot((profile("H2").Y[0,0]-profile("H2").Y[:,0])/(profile("H2").Y[0,0]-profile("H2").Y[-1,0]), profile("NO2").Y, label="NO$_2$")

ax.legend(loc="best")
ax.set_xlabel("Hydrogen c [-]")
ax.set_ylabel("Mass fraction [-]")
ax.set_title("Species Mass Fraction Profile VS Hydrogen progress variable")
ax.grid()
plt.show()


fig, ax = plt.subplots(figsize=(10, 6))
ax.plot((profile("H2").Y[0,0]-profile("H2").Y[:,0])/(profile("H2").Y[0,0]-profile("H2").Y[-1,0]), profile("NO").Y, label="NO")
ax.plot((profile("H2").Y[0,0]-profile("H2").Y[:,0])/(profile("H2").Y[0,0]-profile("H2").Y[-1,0]), profile("NO2").Y, label="NO$_2$")

ax.legend(loc="best")
ax.set_xlabel("Hydrogen c [-]")
ax.set_ylabel("Mass fraction [-]")
ax.set_title("Species Mass Fraction Profile VS Hydrogen progress variable")
ax.grid()
plt.show()



#Question C
fuel_comp_lean = {"H2": profile("H2").X[-1,0],"O2": profile("O2").X[-1,0],"H2O": profile("H2O").X[-1,0],"N2": profile("N2").X[-1,0],"NO": profile("NO").X[-1,0],"NO2": profile("NO2").X[-1,0]}
mass_air = profile("H2").Y[-1,0]/2.018*(3.76*28.02+32) # last bit is for to fix the equivalence ratio (done by hand)

def input_TEMPERATURE(T_1,T_expected,T_2, profile_1,air,injected=False,T_3=False,H2O_mass=False):
    #########################################################################
    # This function computes the average mixture temperature. This is       # 
    # based on iteratively computing the average Cp-value of the injected   #
    # mixures. The total enthalpy defined by the total mass and Cp will     #
    # give the mixture average temperature                                  #
    #########################################################################

    gas1= ct.Solution(mech_file)
    gas2= ct.Solution(mech_file)
    gas3= ct.Solution(mech_file)
    gas4= ct.Solution(mech_file)


    gas1.TPX = T_1,12e5, {"H2": profile_1("H2").X[-1,0],"O2": profile_1("O2").X[-1,0],"H2O": profile_1("H2O").X[-1,0],"N2": profile_1("N2").X[-1,0],"NO": profile_1("NO").X[-1,0],"NO2": profile_1("NO2").X[-1,0]}
    gas2.TPX = T_2,12e5,air
    gas3.TPX = T_expected,12e5, {"H2": profile_1("H2").X[-1,0],"O2": profile_1("O2").X[-1,0],"H2O": profile_1("H2O").X[-1,0],"N2": profile_1("N2").X[-1,0],"NO": profile_1("NO").X[-1,0],"NO2": profile_1("NO2").X[-1,0]}
    gas4.TPX = T_expected,12e5,air

    if injected==False: 
        T_lean = ((gas1.h)+(gas2.h)*mass_air)/((gas1.h+gas3.h)/(T_1+T_expected) + mass_air*(gas2.h+gas4.h)/(T_2+T_expected))
    else:
        gas5= ct.Solution(mech_file)
        gas6= ct.Solution(mech_file)
        gas5.TPX = T_3,12e5,injected
        gas6.TPX = T_expected,12e5,injected
        T_lean = ((gas1.h)+(gas2.h)*mass_air-gas5.h*H2O_mass) / ((gas1.h+gas3.h)/(T_1+T_expected) + mass_air*(gas2.h+gas4.h)/(T_2+T_expected) - H2O_mass*(gas5.h+gas6.h)/(T_3+T_expected))


    return T_lean



T_lean= profile.T[-1] # initializing expected lean mixture combustor inlet temperature 
while abs(input_TEMPERATURE(profile.T[-1],T_lean,T_in,profile,{"O2":1,"N2":3.7})-T_lean)>0.1:
    T_lean=(input_TEMPERATURE(profile.T[-1],T_lean,T_in,profile,{"O2":1,"N2":3.7})+T_lean)/2


flame_lean = freeflamestuff(T_lean,12E5,phi_lean,fuel_comp_lean)
profile_lean = flame_lean.to_array()

print("QUESTION C")
print("T lean: ",T_lean)


# space variable = water vapour
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(profile_lean("H2O").Y[:,0]/profile_lean("H2O").Y[-1,0], profile_lean.T, ".-")
ax.set_xlabel("Water Vapour c [-]")
ax.set_ylabel("Temperature [K]")
ax.set_title("Flame Temperature Profile VS Water Vapour progress variable")
ax.grid()
plt.show()

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(profile_lean("H2O").Y[:,0]/profile_lean("H2O").Y[-1,0], profile_lean("O2").Y, label="O$_2$")
ax.plot(profile_lean("H2O").Y[:,0]/profile_lean("H2O").Y[-1,0], profile_lean("H2").Y, label="H$_2$")
ax.plot(profile_lean("H2O").Y[:,0]/profile_lean("H2O").Y[-1,0], profile_lean("NO").Y, label="NO")
ax.plot(profile_lean("H2O").Y[:,0]/profile_lean("H2O").Y[-1,0], profile_lean("NO2").Y, label="NO$_2$")

ax.legend(loc="best")
ax.set_xlabel("Water Vapour c")
ax.set_ylabel("Mass fraction [-]")
ax.set_title("Species Mass Fraction Profile VS Water Vapour progress variable")
ax.grid()
plt.show()

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(profile_lean("H2O").Y[:,0]/profile_lean("H2O").Y[-1,0], profile_lean("NO").Y, label="NO")
ax.plot(profile_lean("H2O").Y[:,0]/profile_lean("H2O").Y[-1,0], profile_lean("NO2").Y, label="NO$_2$")

ax.legend(loc="best")
ax.set_xlabel("Water Vapour c")
ax.set_ylabel("Mass fraction [-]")
ax.set_title("Species Mass Fraction  Profile of NOx gasses VS Water Vapour progress variable")
ax.grid()
plt.show()


# space variable = hydrogen
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(1-profile_lean("H2").Y[:,0]/profile_lean("H2").Y[0,0], profile_lean.T, ".-")
ax.set_xlabel("Hydrogen c [-]")
ax.set_ylabel("Temperature [K]")
ax.set_title("Flame Temperature Profile VS Hydrogen progress variable")
ax.grid()
plt.show()


fig, ax = plt.subplots(figsize=(10, 6))
ax.plot((profile_lean("H2").Y[0,0]-profile_lean("H2").Y[:,0])/(profile_lean("H2").Y[0,0]-profile_lean("H2").Y[-1,0]), profile_lean("O2").Y, label="O$_2$")
ax.plot((profile_lean("H2").Y[0,0]-profile_lean("H2").Y[:,0])/(profile_lean("H2").Y[0,0]-profile_lean("H2").Y[-1,0]), profile_lean("H2").Y, label="H$_2$")
ax.plot((profile_lean("H2").Y[0,0]-profile_lean("H2").Y[:,0])/(profile_lean("H2").Y[0,0]-profile_lean("H2").Y[-1,0]), profile_lean("NO").Y, label="NO")
#ax.plot((profile_lean("H2").Y[0,0]-profile_lean("H2").Y[:,0])/(profile_lean("H2").Y[0,0]-profile_lean("H2").Y[-1,0]), profile_lean("NO2").Y, label="NO$_2$")

ax.legend(loc="best")
ax.set_xlabel("Hydrogen c [-]")
ax.set_ylabel("Mass fraction [-]")
ax.set_title("Species Mass Fraction Profile VS Hydrogen progress variable")
ax.grid()
plt.show()


fig, ax = plt.subplots(figsize=(10, 6))
ax.plot((profile_lean("H2").Y[0,0]-profile_lean("H2").Y[:,0])/(profile_lean("H2").Y[0,0]-profile_lean("H2").Y[-1,0]), profile_lean("NO").Y, label="NO")
ax.plot((profile_lean("H2").Y[0,0]-profile_lean("H2").Y[:,0])/(profile_lean("H2").Y[0,0]-profile_lean("H2").Y[-1,0]), profile_lean("NO2").Y, label="NO$_2$")

ax.legend(loc="best")
ax.set_xlabel("Hydrogen c [-]")
ax.set_ylabel("Mass fraction [-]")
ax.set_title("Species Mass Fraction Profile of NOx gasses VS Hydrogen progress variable")
ax.grid()
plt.show()


#Question D
Y_H2O = profile_lean("H2O").Y[-1,0]
Y_O2 = profile_lean("O2").Y[-1,0]
Y_N2 = 1-Y_H2O-Y_O2

X_H2O = profile_lean("H2O").X[-1,0]
X_O2 = profile_lean("O2").X[-1,0]
X_N2 = 1-X_H2O-X_O2-X_O2 
X_NO = 2*X_O2
Y_NO = X_NO*28.01/(X_H2O*18.018+X_N2*24.02+X_NO*28.01)
print("QUESTION D")
print("Recomputed mass fraction of NO: ",Y_NO,"\nComputed by Cantera mass fraction of NO: ",profile_lean("NO").Y[-1,0],"\nRatio of NO_recomputed/NO_cantera")

# # E 
# strain_rates = [100,500]
# for i in range(1):
#     strained_flame = freeflamestuff(T_lean,12E5,phi_lean,fuel_comp_lean,strain_rates[i],flame_lean.velocity[0])
#     print(f"\nFlame speed at a strain rate of {strain_rates[i]}: {strained_flame.velocity[0]}m/s")


# Question F
T_H2O_injected =  profile_lean.T[-1]
T_mixture_injected = (T_H2O_injected+T_lean)/2

H2O_MAX = profile_lean("H2O").X[-1,0]

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(profile_lean("H2O").Y[:,0]/profile_lean("H2O").Y[-1,0], profile_lean.T, label=f"T - Water vapour{(0)*25}%")
for i in range(4):
    while abs(input_TEMPERATURE(profile.T[-1],T_mixture_injected,T_in,profile,{"O2":1,"N2":3.7},{"H2O":1},T_H2O_injected,H2O_MAX/4*(i+1))-T_mixture_injected)>0.1:
        T_mixture_injected=(input_TEMPERATURE(profile.T[-1],T_mixture_injected,T_in,profile,{"O2":1,"N2":3.7},{"H2O":1},T_H2O_injected,H2O_MAX/4*(i+1))+T_mixture_injected)/2

    phi_lean= phi_lean*(1+3.7)/(H2O_MAX/4*(i+1)+1+3.7)#Calculation of new phi, due to the added water vapour
    flame_water = freeflamestuff(T_mixture_injected,12E5,phi_lean,fuel_comp_lean,False,False,H2O_MAX/4*(i+1)) 
    profile_water = flame_water.to_array()
    
    ax.plot(profile_water("H2O").Y[:,0]/profile_water("H2O").Y[-1,0], profile_water.T, label=f"T - Water vapour{(i+1)*25}%")


ax.legend(loc="best")
ax.set_xlabel("Water Vapour c [-]")
ax.set_ylabel("Temperature [K]")
ax.set_title("Temperature profile VS Water Vapour progress variable")
ax.grid()
print("QUESTION F")
plt.show()



fig, ax = plt.subplots(figsize=(10, 6))


ax.plot(profile_lean("H2O").Y[:,0]/profile_lean("H2O").Y[-1,0], profile_lean("NO").Y, label=f"NO - Water vapour{(0)*25}%")
for i in range(4):
    while abs(input_TEMPERATURE(profile.T[-1],T_mixture_injected,T_in,profile,{"O2":1,"N2":3.7},{"H2O":1},T_H2O_injected,H2O_MAX/4*(i+1))-T_mixture_injected)>0.1:
        T_mixture_injected=(input_TEMPERATURE(profile.T[-1],T_mixture_injected,T_in,profile,{"O2":1,"N2":3.7},{"H2O":1},T_H2O_injected,H2O_MAX/4*(i+1))+T_mixture_injected)/2
    phi_lean= phi_lean*(1+3.7)/(H2O_MAX/4*(i+1)+1+3.7)
    flame_water = freeflamestuff(T_mixture_injected,12E5,phi_lean,fuel_comp_lean,False,False,H2O_MAX/4*(i+1)) 
    profile_water = flame_water.to_array()
    
    ax.plot(profile_water("H2O").Y[:,0]/profile_water("H2O").Y[-1,0], profile_water("NO").Y, label=f"NO - Water vapour{(i+1)*25}%")

ax.legend(loc="best")
ax.set_xlabel("Water Vapour c progress variable [-]")
ax.set_ylabel("Mass fraction [-]")
ax.set_title("Mass fraction profile of NO at verious levels of injected water vapour")
ax.grid()
plt.show()


#How many moles water vapour
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(1-profile_lean("H2").Y[:,0]/profile_lean("H2").Y[0,0], profile_lean.T, label=f"T - Water vapour{(0)*25}%")
for i in range(4):
    while abs(input_TEMPERATURE(profile.T[-1],T_mixture_injected,T_in,profile,{"O2":1,"N2":3.7},{"H2O":1},T_H2O_injected,H2O_MAX/4*(i+1))-T_mixture_injected)>0.1:
        T_mixture_injected=(input_TEMPERATURE(profile.T[-1],T_mixture_injected,T_in,profile,{"O2":1,"N2":3.7},{"H2O":1},T_H2O_injected,H2O_MAX/4*(i+1))+T_mixture_injected)/2
    phi_lean= phi_lean*(1+3.7)/(H2O_MAX/4*(i+1)+1+3.7)
    flame_water = freeflamestuff(T_mixture_injected,12E5,phi_lean,fuel_comp_lean,False,False,H2O_MAX/4*(i+1)) 
    profile_water = flame_water.to_array()

    ax.plot(1-profile_water("H2").Y[:,0]/profile_water("H2").Y[0,0], profile_water.T, label=f"T - Water vapour{(i+1)*25}%")

ax.legend(loc="best")
ax.set_xlabel("Hydrogen c [-]")
ax.set_ylabel("Temperature [K]")
ax.set_title("Temperature profile VS Hydrogen progress variable")
ax.grid()
print("QUESTION F")
plt.show()


fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(1-profile_lean("H2").Y[:,0]/profile_lean("H2").Y[0,0], profile_lean("NO").Y, label=f"NO - Water vapour{(0)*25}%")
for i in range(4):
    while abs(input_TEMPERATURE(profile.T[-1],T_mixture_injected,T_in,profile,{"O2":1,"N2":3.7},{"H2O":1},T_H2O_injected,H2O_MAX/4*(i+1))-T_mixture_injected)>0.1:
        T_mixture_injected=(input_TEMPERATURE(profile.T[-1],T_mixture_injected,T_in,profile,{"O2":1,"N2":3.7},{"H2O":1},T_H2O_injected,H2O_MAX/4*(i+1))+T_mixture_injected)/2
    phi_lean= phi_lean*(1+3.7)/(H2O_MAX/4*(i+1)+1+3.7)
    flame_water = freeflamestuff(T_mixture_injected,12E5,phi_lean,fuel_comp_lean,False,False,H2O_MAX/4*(i+1)) 
    profile_water = flame_water.to_array()
    
    ax.plot(1-profile_water("H2").Y[:,0]/profile_water("H2").Y[0,0], profile_water("NO").Y, label=f"NO - Water vapour{(i+1)*25}%")



ax.legend(loc="best")
ax.set_xlabel("Hydrogen c progress variable [-]")
ax.set_ylabel("Mass fraction [-]")
ax.set_title("Mass fraction profile of NO at verious levels of injected water vapour")
ax.grid()
plt.show()



# distance wise space variable:
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(profile_lean.grid/profile_lean.grid[-1], profile_lean.T, label=f"T - Water vapour{(0)*25}%")
for i in range(4):
    while abs(input_TEMPERATURE(profile.T[-1],T_mixture_injected,T_in,profile,{"O2":1,"N2":3.7},{"H2O":1},T_H2O_injected,H2O_MAX/4*(i+1))-T_mixture_injected)>0.1:
        T_mixture_injected=(input_TEMPERATURE(profile.T[-1],T_mixture_injected,T_in,profile,{"O2":1,"N2":3.7},{"H2O":1},T_H2O_injected,H2O_MAX/4*(i+1))+T_mixture_injected)/2
    phi_lean= phi_lean*(1+3.7)/(H2O_MAX/4*(i+1)+1+3.7)
    flame_water = freeflamestuff(T_mixture_injected,12E5,phi_lean,fuel_comp_lean,False,False,H2O_MAX/4*(i+1)) 
    profile_water = flame_water.to_array()
    
    ax.plot(profile_water.grid/profile_water.grid[-1], profile_water.T, label=f"T - Water vapour{(i+1)*25}%")


ax.legend(loc="best")
ax.set_xlabel("Hydrogen c [-]")
ax.set_ylabel("Temperature [K]")
ax.set_title("Temperature profile VS Hydrogen progress variable")
ax.grid()
plt.show()


fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(profile_lean.grid/profile_lean.grid[-1], profile_lean("NO").Y, label=f"NO - Water vapour{(0)*25}%")
for i in range(4):
    while abs(input_TEMPERATURE(profile.T[-1],T_mixture_injected,T_in,profile,{"O2":1,"N2":3.7},{"H2O":1},T_H2O_injected,H2O_MAX/4*(i+1))-T_mixture_injected)>0.1:
        T_mixture_injected=(input_TEMPERATURE(profile.T[-1],T_mixture_injected,T_in,profile,{"O2":1,"N2":3.7},{"H2O":1},T_H2O_injected,H2O_MAX/4*(i+1))+T_mixture_injected)/2
    phi_lean= phi_lean*(1+3.7)/(H2O_MAX/4*(i+1)+1+3.7)
    flame_water = freeflamestuff(T_mixture_injected,12E5,phi_lean,fuel_comp_lean,False,False,H2O_MAX/4*(i+1)) 
    profile_water = flame_water.to_array()
    
    ax.plot(profile_water.grid/profile_water.grid[-1], profile_water("NO").Y, label=f"NO - Water vapour{(i+1)*25}%")
    ax.plot(profile_lean.grid * 100, profile_lean.T, label=f"T - Water vapour{(i+1)*25}%")


ax.legend(loc="best")
ax.set_xlabel("Distance c progress variable [-]")
ax.set_ylabel("Mass fraction [-]")
ax.set_title("Mass fraction profile of NO at verious levels of injected water vapour")
ax.grid()
plt.show()



