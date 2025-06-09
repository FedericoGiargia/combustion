import cantera as ct
import matplotlib.pyplot as plt
from pathlib import Path
phi_rich = 5.0
phi_lean = 0.5
mech_file = "gri30.yaml"
fuel_comp_rich = {"H2": 1.0}

def freeflamestuff(T,P,phi, fuel_comp):
    # Inlet temperature in kelvin and inlet pressure in pascals
    # In this case we are setting the inlet T and P to room temperature conditions
    To = T
    Po = P#ct.one_atm

    # Domain width in metres
    width = 0.12

    # Create the object representing the gas and set its state to match the inlet conditions
    gas = ct.Solution(mech_file)
    gas.TP = To, Po
    gas.set_equivalence_ratio(phi, fuel_comp, {"O2": 1.0, "N2": 3.76})

    flame = ct.FreeFlame(gas, width=width)
    #flame.transport_model = 'mixture-averaged'
    flame.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)

    flame.solve(loglevel=1, auto=True)
    return flame
flame = (freeflamestuff(675,12E5,5,fuel_comp_rich))
# Extract the spatial profiles as a SolutionArray to simplify plotting specific species
profile = flame.to_array()

# Question Part A:
# print(f"\nmixture-averaged flame speed = {flame.velocity[0]:7f} m/s\n")



# Part B B 


# fig, ax = plt.subplots(figsize=(10, 6))
# ax.plot(profile.grid * 100, profile.T, ".-")
# ax.set_xlabel("Distance [cm]")
# ax.set_ylabel("Temperature [K]")
# ax.set_title("Flame Temperature Profile")
# ax.grid()
# plt.savefig(output_dir / 'Flame_temperature_profile.png', dpi=300)
# plt.show()

# fig, ax = plt.subplots(figsize=(10, 6))

# # ax.plot(profile.grid * 100, profile("CO2").Y, "--", label="CO$_2$")
# # ax.plot(profile.grid * 100, profile("CH4").Y, "--", label="CH$_4$")
# ax.plot(profile.grid * 100, profile("H2O").Y, label="H$_2$O")
# ax.plot(profile.grid * 100, profile("O2").Y, label="O$_2$")
# ax.plot(profile.grid * 100, profile("H2").Y, label="H$_2$")
# ax.plot(profile.grid * 100, profile("NO").Y, label="NO")
# ax.plot(profile.grid * 100, profile("NO2").Y, label="NO$_2$")


# ax.legend(loc="best")
# ax.set_xlabel("Distance [cm]")
# ax.set_ylabel("Mass fraction [-]")
# ax.set_title("Species Mass Fraction Profile")
# ax.grid()
# plt.savefig(output_dir / 'Species_Mass_Fraction_Profile.png', dpi=300)
# plt.show()

# fig, ax = plt.subplots(figsize=(10, 6))

# # ax.plot(profile.grid * 100, profile("CO2").Y, "--", label="CO$_2$")
# # ax.plot(profile.grid * 100, profile("CH4").Y, "--", label="CH$_4$")
# # ax.plot(profile.grid * 100, profile("H2O").Y, label="H$_2$O")
# # ax.plot(profile.grid * 100, profile("O2").Y, label="O$_2$")
# # ax.plot(profile.grid * 100, profile("H2").Y, label="H$_2$")
# ax.plot(profile.grid * 100, profile("NO").Y, label="NO")
# ax.plot(profile.grid * 100, profile("NO2").Y, label="NO$_2$")


# ax.legend(loc="best")
# ax.set_xlabel("Distance [cm]")
# ax.set_ylabel("Mass fraction [-]")
# ax.set_title("Species Mass Fraction Profile")
# ax.grid()
# plt.savefig(output_dir / 'Species_Mass_Fraction_Profile.png', dpi=300)
# plt.show()



#Question C
fuel_comp_lean = {"H2": profile("H2").X[-1,0],"O2": profile("O2").X[-1,0],"H2O": profile("H2O").X[-1,0],"N2": profile("N2").X[-1,0],"NO": profile("NO").X[-1,0],"NO2": profile("NO2").X[-1,0]}

T_Rich = profile.T[-1]
P_Rich = profile.P[-1]

# check if these CP's are all at the right temperature
H2_cp = 14290
O2_cp = 918
H2O_cp = 1900
N2_cp = 1041
NO_cp = 1214.5 #at 1800K
NO2_cp = 1000
Air_cp = 1005

print(profile.T[-1],profile.P[-1])

mass_air = profile("H2").Y[-1,0]/2.018*(3.76*28.02+32) *2*phi_lean#last bit is for to fix the equivalence ratio
T_lean = (675*mass_air*Air_cp + T_Rich*(profile("H2").Y[-1,0]*H2_cp + profile("O2").Y[-1,0]*O2_cp + profile("H2O").Y[-1,0]*H2O_cp + profile("N2").Y[-1,0]*N2_cp + profile("NO").Y[-1,0]*NO_cp + profile("NO2").Y[-1,0]*NO2_cp)) /(profile("H2").Y[-1,0]*H2_cp + profile("O2").Y[-1,0]*O2_cp + profile("H2O").Y[-1,0]*H2O_cp + profile("N2").Y[-1,0]*N2_cp + profile("NO").Y[-1,0]*NO_cp + profile("NO2").Y[-1,0]*NO2_cp+mass_air*Air_cp)

flame_lean = freeflamestuff(T_lean,P_Rich,phi_lean,fuel_comp_lean)
print(f"\nmixture-averaged flame speed = {flame_lean.velocity[0]:7f} m/s\n")

# Create directory for output files if it doesn't exist
output_dir = Path('Laminar_flame_speed')
output_dir.mkdir(exist_ok=True)

# Extract the spatial profiles as a SolutionArray to simplify plotting specific species
profile_lean = flame_lean.to_array()


fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(profile_lean.grid * 100, profile_lean.T, ".-")
ax.set_xlabel("Distance [cm]")
ax.set_ylabel("Temperature [K]")
ax.set_title("Flame Temperature Profile")
ax.grid()
plt.savefig(output_dir / 'Flame_temperature_profile.png', dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(10, 6))

# ax.plot(profile_lean.grid * 100, profile_lean("CO2").Y, "--", label="CO$_2$")
# ax.plot(profile_lean.grid * 100, profile_lean("CH4").Y, "--", label="CH$_4$")
ax.plot(profile_lean.grid * 100, profile_lean("H2O").Y, label="H$_2$O")
ax.plot(profile_lean.grid * 100, profile_lean("O2").Y, label="O$_2$")
ax.plot(profile_lean.grid * 100, profile_lean("H2").Y, label="H$_2$")
ax.plot(profile_lean.grid * 100, profile_lean("NO").Y, label="NO")
ax.plot(profile_lean.grid * 100, profile_lean("NO2").Y, label="NO$_2$")


ax.legend(loc="best")
ax.set_xlabel("Distance [cm]")
ax.set_ylabel("Mass fraction [-]")
ax.set_title("Species Mass Fraction Profile")
ax.grid()
plt.savefig(output_dir / 'Species_Mass_Fraction_Profile.png', dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(10, 6))

# ax.plot(profile_lean.grid * 100, profile_lean("CO2").Y, "--", label="CO$_2$")
# ax.plot(profile_lean.grid * 100, profile_lean("CH4").Y, "--", label="CH$_4$")
# ax.plot(profile_lean.grid * 100, profile_lean("H2O").Y, label="H$_2$O")
# ax.plot(profile_lean.grid * 100, profile_lean("O2").Y, label="O$_2$")
# ax.plot(profile_lean.grid * 100, profile_lean("H2").Y, label="H$_2$")
ax.plot(profile_lean.grid * 100, profile_lean("NO").Y, label="NO")
ax.plot(profile_lean.grid * 100, profile_lean("NO2").Y, label="NO$_2$")


ax.legend(loc="best")
ax.set_xlabel("Distance [cm]")
ax.set_ylabel("Mass fraction [-]")
ax.set_title("Species Mass Fraction Profile")
ax.grid()
plt.savefig(output_dir / 'Species_Mass_Fraction_Profile.png', dpi=300)
plt.show()
