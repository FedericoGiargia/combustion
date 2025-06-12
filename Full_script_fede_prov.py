import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import copy
import time
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

# dctionary for time
_timers = {}

def tic(label='default'):
    _timers[label] = time.perf_counter()
def toc(label='default'):
    if label not in _timers:
        print(f"[toc] Nessun tic('{label}') trovato")
        return
    elapsed = time.perf_counter() - _timers[label]
    print(f"[toc] {label}: {elapsed:.6f} s")
    return elapsed

#%% start by validating the solution on methane gas

"""methane = ct.Solution('gri30.yaml')
phi = 1.4
T0 = 293
P0 = ct.one_atm

fuel_comp = {'CH4' : 1.0}
oxidizer = {'O2':1.0, 'N2':3.76}
methane.set_equivalence_ratio(phi=phi, fuel='CH4', oxidizer={'O2':1.0, 'N2':3.76})
methane.TP = T0, P0

width = 0.12
tic('validation step')


#creating the flame
flame = ct.FreeFlame(methane, width=width)
flame.transport_model = 'mixture-averaged'
flame.set_refine_criteria(ratio=3, slope=0.07, curve=0.12)   #tared with a couple runs

#solving the flame
flame.solve(loglevel=1, auto=True)

#result visualization + comparison with all the other models
print(f'\nmixture-averaged flame speed = {flame.velocity[0]:7f} m/s\n')
profile = flame.to_array()
fig, ax = plt.subplots()
ax.plot(profile.grid*100, profile.T, ".-")
ax.set_xlabel('Distance [cm]')
ax.set_ylabel('Temperature [K]')
plt.show()

fix, ax = plt.subplots()
ax.plot(profile.grid*100, profile('CH4').X, "--", label = 'CH$_4$')
ax.plot(profile.grid*100, profile('CO2').X, "--", label = 'CO$_2$')
ax.plot(profile.grid*100, profile('N2O').X, "--", label = 'N$_2$O')
ax.plot(profile.grid*100, profile('H2O').X, "--", label = 'H$_2$O')
ax.legend(loc= "best")
ax.set_xlabel('Distance [cm]')
ax.set_ylabel('Mole fraction [-]')
plt.show()
toc('validation step')
"""
#%% creation of the class to initialize flames (useful for comparisons)

class PersonalizedFlame:
    def __init__(
        self,
        mech,
        phi=None,
        reactant=None,
        oxidizer=None,
        X_mix=None,
        T0=None,
        P0=None,
        tp_model='multicomponent',
        width=0.03,
        CounterFlowFlame = False,
        reactants_mdot = None,
        products_mdot = None
    ):
        """
        Initialize a 1D freely propagating flame.

        Parameters:
        - mech: kinetic mechanism file
        - phi: equivalence ratio (if X_mix is None)
        - reactant, oxidizer: dicts defining mixture (if X_mix is None)
        - X_mix: dict or list of mole fractions (if provided, phi/reactant/oxidizer ignored)
        - T0: inlet temperature [K]
        - p: pressure [Pa]
        - tp_model: transport model ('mixture-averaged' or 'multicomponent')
        - width: domain width [m]
        """
        self.mech = mech
        self.phi = phi
        self.reactant = reactant
        self.oxidizer = oxidizer
        self.X_mix = X_mix
        self.T0 = T0
        self.P0 = P0
        self.tp_model = tp_model
        self.width = width
        self.CounterFlowFlame = CounterFlowFlame
        self.reactants_mdot = reactants_mdot
        self.products_mdot = products_mdot

        # Create and solve base flame, if we have enough data (to accommodate for point e)
        if CounterFlowFlame == False:
           self.make_and_solve_base_flame()
        else:
            self.make_and_solve_strained_flame()

    def make_and_solve_base_flame(self):
        gas = ct.Solution(self.mech)
        # Determine mixture
        if self.X_mix is not None:
            # Use provided mole fractions and T0, p
            gas.TPX = self.T0, self.P0, self.X_mix
        else:
            # Use phi/reactant/oxidizer
            gas.set_equivalence_ratio(self.phi, self.reactant, self.oxidizer)
            gas.TP = self.T0, self.P0
        # Initialize flame object
        flame = ct.FreeFlame(gas, width=self.width)
        flame.transport_model = self.tp_model
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
        # Solve
        flame.solve(loglevel=0, refine_grid=True, auto=True)
        self.gas0 = gas
        self.flame = flame
    def make_and_solve_strained_flame(self):
        gas = ct.Solution(self.mech)
        # Determine mixture
        if self.X_mix is not None:
            # Use provided mole fractions and T0, p
            gas.TPX = self.T0, self.P0, self.X_mix
        else:
            # Use phi/reactant/oxidizer
            gas.set_equivalence_ratio(self.phi, self.reactant, self.oxidizer)
            gas.TP = self.T0, self.P0
        # Initialize flame object
        flame = ct.CounterflowPremixedFlame(gas, width=self.width)
        flame.reactants.mdot = self.reactants_mdot
        flame.products.mdot = self.products_mdot
        flame.transport_model = self.tp_model
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
        # Solve
        flame.solve(loglevel=0, refine_grid=True, auto=True)
        self.gas0 = gas
        self.flame = flame

    def get_product(self,species_name):
        idx = self.gas0.species_index(species_name) #the process is in two steps, find the index
        return float(self.flame.X[-1, idx])   #use the index in las position

    def get_mass_product(self, species_name):
        idx = self.gas0.species_index(species_name)  # the process is in two steps, find the index
        return float(self.flame.Y[-1, idx])  # to get the mass one

    def flame_speed(self):
        """Return the scalar laminar flame speed (m/s)."""
        # u[0] is the inlet (unburned) velocity = flame speed
        return float(self.flame.velocity[0])

    def compute_fd_sensitivities(self, dA=0.01, top_n=10):
        """
        Compute finite-difference sensitivities of flame speed wrt Arrhenius A factor.
        """
        SL0 = self.flame_speed()
        n_rxn = self.gas0.n_reactions
        sens = np.full(n_rxn, np.nan)

        for i in range(n_rxn):
            # Recreate gas and flame
            gas = ct.Solution(self.mech)
            gas.set_equivalence_ratio(self.phi, self.reactant, self.oxidizer)
            gas.TP = self.T0, self.P0
            flame = ct.FreeFlame(gas, width=self.width)
            flame.transport_model = self.tp_model
            flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)

            # Only perturb Arrhenius rates
            rxn = gas.reaction(i)
            rate = rxn.rate
            if not isinstance(rate, ct.Arrhenius):
                continue
            # Build perturbed rate
            perturbed_rate = ct.Arrhenius(
                rate.pre_exponential_factor * (1.0 + dA),
                rate.temperature_exponent,
                rate.activation_energy
            )
            rxn.rate = perturbed_rate
            gas.modify_reaction(i, rxn)

            # Solve perturbed flame; catch failures
            try:
                flame.solve(loglevel=0, refine_grid=True, auto=True)
                SLi = float(flame.u[0])
                sens[i] = (np.log(SLi) - np.log(SL0)) / np.log(1.0 + dA)
            except Exception:
                continue

        # Sort and return top_n
        idx_sorted = np.argsort(np.nan_to_num(np.abs(sens), nan=-1))[::-1]
        results = []
        for idx in idx_sorted[:top_n]:
            eq = self.gas0.reaction(idx).equation
            results.append({'index': idx, 'sensitivity': sens[idx], 'equation': eq})
        return results
    def plotting_rounine(self, name):
        profile = self.flame.to_array()
        fig, ax = plt.subplots()
        ax.plot(profile.grid * 100, profile.T, '.-')
        ax.set_xlabel('space (cm)')
        ax.set_ylabel('Temperature (K)')
        ax.set_title(f"Temperature profile of the case {mech}")
        plt.figure()
        fix, ax = plt.subplots()
        ax.plot(profile.grid * 100, profile('H2').X, "--", label='H$_2$')
        ax.plot(profile.grid * 100, profile('O2').X, "--", label='O$_2$')
        ax.plot(profile.grid * 100, profile('NO2').X, "--", label='NO$_2$')
        ax.plot(profile.grid * 100, profile('NO').X, "--", label='NO')
        ax.plot(profile.grid * 100, profile('H2O').X, "--", label='H$_2$O')
        ax.legend(loc="best")
        ax.set_xlabel('Distance [cm]')
        ax.set_ylabel('Mole fraction [-]')
        ax.set_title(f"reaction profile of the case {name}")
        plt.show()
        toc('lean mixture')


# %% point (a) with some different variations

tic('comparison step')
hydrogen_gri30_averaged = PersonalizedFlame(mech ='gri30.yaml',
                                            phi =5.0,
                                            reactant = {'H2' : 1.0},
                                            oxidizer={'O2':1.0, 'N2':3.76},
                                            T0 = 675.0,
                                            P0 =12e5,
                                            tp_model='mixture-averaged',
                                            width = 0.3)

#comparison with a smaller phi, to see if the velocities increases (model validation)
hydrogen_gri30_averaged_phi1 = PersonalizedFlame(mech ='gri30.yaml',
                                                 phi =1.0,
                                                 reactant = {'H2' : 1.0},
                                                 oxidizer={'O2':1.0, 'N2':3.76},
                                                 T0 = 675.0,
                                                 P0 =12e5,
                                                 tp_model='mixture-averaged',
                                                 width = 0.3)

#comparison with multicomponent, should be less accurate
hydrogen_gri30_mult = PersonalizedFlame(mech ='gri30.yaml',
                                        phi =5.0,
                                        reactant = {'H2' : 1.0},
                                        oxidizer={'O2':1.0, 'N2':3.76},
                                        T0 = 675.0,
                                        P0 =12e5,
                                        tp_model='multicomponent',
                                        width = 0.3)


hydrogen_keromnes_mult = PersonalizedFlame(mech ='keromnes2013.yaml',
                                           phi=5.0,
                                           reactant={'H2': 1.0},
                                           oxidizer={'O2': 1.0, 'N2': 3.76},
                                           T0=675.0,
                                           P0=12e5,
                                           tp_model='multicomponent',
                                           width=0.3)

hydrogen_konnov_mult = PersonalizedFlame(mech ='h2-konnov.yaml',
                                         phi=5.0,
                                         reactant={'H2': 1.0},
                                         oxidizer={'O2': 1.0, 'N2': 3.76},
                                         T0=675.0,
                                         P0=12e5,
                                         tp_model='multicomponent',
                                         width=0.3)


# Lista di meccanismi da confrontare

mechanisms = [ hydrogen_gri30_averaged ,hydrogen_gri30_averaged_phi1,
              hydrogen_gri30_mult , hydrogen_konnov_mult , hydrogen_keromnes_mult]

for mech in mechanisms: #creates for cycle
    print(f"\n>>> Meccanismo: {mech}") #tells what is what
    #mech.solve() #solves the flame
    mech.flame_speed() #finds the velocity
    """
    sens_list = mech.compute_fd_sensitivities() #tells the 8 most important reactions
    print(f"Top 8 sensibilità per {mech}:")
    for entry in sens_list:
        idx, sens, eq = entry['index'], entry['sensitivity'], entry['equation'] #creates a dictionary

        print(f"  R{idx:03d}: " #3 decimals
              f"S = {sens:+.3f}" #with sign and decimals
              f"  |  {eq}") #uses the dictionary to find what we need
              """
    print(f"the flame velocity for the mechanism : {mech} is : {mech.flame.velocity}")
    profile = mech.flame.to_array()
    fig, ax = plt.subplots()
    ax.plot(profile.grid * 100, profile.T ,'.-')
    ax.set_xlabel('space (cm)')
    ax.set_ylabel('Temperature (K)')
    ax.set_title(f"Temperature profile of the case {mech}")
    plt.figure()
    fix, ax = plt.subplots()
    ax.plot(profile.grid * 100, profile('H2').X, "--", label='H$_2$')
    ax.plot(profile.grid * 100, profile('O2').X, "--", label='O$_2$')
    ax.plot(profile.grid * 100, profile('NO2').X, "--", label='NO$_2$')
    #ax.plot(profile.grid * 100, profile('NO').X, "--", label='NO')
    ax.plot(profile.grid * 100, profile('H2O').X, "--", label='H$_2$O')
    ax.legend(loc="best")
    ax.set_xlabel('Distance [cm]')
    ax.set_ylabel('Mole fraction [-]')
    ax.set_title(f"reaction profile of the case {mech}")
    plt.show()
toc('comparison step')
#%% excercise b

fix, ax = plt.subplots()
ax.plot(profile.grid*100, profile('H2').X, "--", label = 'H$_2$')
ax.plot(profile.grid*100, profile('O2').X, "--", label = 'O$_2$')
ax.plot(profile.grid*100, profile('NO2').X, "--", label = 'NO$_2$')
ax.plot(profile.grid*100, profile('NO').X, "--", label = 'NO')
ax.plot(profile.grid*100, profile('H2O').X, "--", label = 'H$_2$O')
ax.legend(loc= "best")
ax.set_xlabel('Distance [cm]')
ax.set_ylabel('Mole fraction [-]')
plt.show()




#%% excercise c

tic('lean mixture')
"""
[6 pts] Now mix the composition in output from the R part with air in order to obtain an 
equivalence ratio 𝜙 = 0.5. Indicate the amount of air required and assume it is injected at 675 
K and 12 bar. Compute the freely-propagating 1D flat flame for the new condition and plot 
again the variation of NO, NO2, O2 and H2 mass fractions, and temperature, in progress 
variable space using the two definitions given earlier. Are the considerations about progress 
variable the same as before? Comment on the results. 
"""

#we assume that the air entering has the same operating conditions as the I.C
output_T = hydrogen_keromnes_mult.flame.T[-1] #get the output temperature

X_Rout = hydrogen_keromnes_mult.flame.X[-1, :].copy() #get the output components

species = hydrogen_keromnes_mult.gas0.species_names #finds ONLY the species that we have
#extract the components of the previous iteration that we need for combustion, H2 and AIR (O2 + N2)

H2_moles = hydrogen_keromnes_mult.get_product('H2')
O2_moles = hydrogen_keromnes_mult.get_product('O2')
print(O2_moles)  #very important, this is sort of a "validation step". it is almost 0 , negligible
N2_moles = hydrogen_keromnes_mult.get_product('N2')
print(N2_moles)

#the phi = 0.5 is fuel/air divided by stoichiometric. we just add oxigen
n_base = 1.0
H2_product = H2_moles * n_base
O2_needed = H2_product #compel with the formula
air_needed = O2_needed / 0.21
N2_needed = air_needed * 0.79 # need to sum with the previous N2

moles_R = {
    sp: X_Rout[k] * n_base
    for k, sp in enumerate(species)}
# adding step:
moles_mix = moles_R.copy()
moles_mix['O2'] = moles_mix.get('O2', 0.0) + O2_needed #just to account for the small portion, useless
moles_mix['N2'] = moles_mix.get('N2', 0.0) + N2_needed #sum with the one that is there


n_total = sum(moles_mix.values())
#sort of validation, check how much is the phi
X_mix   = { sp: moles_mix[sp]/n_total  for sp in moles_mix }

hydrogen_keromnes_mult_lean = PersonalizedFlame(mech ='keromnes2013.yaml',
                                                X_mix=X_mix,
                                                reactant={'H2': 1.0},
                                                oxidizer={'O2': 1.0, 'N2': 3.76},
                                                T0=675.0,
                                                P0=12e5,
                                                tp_model='multicomponent',
                                                width=0.3)

X_H2_new = hydrogen_keromnes_mult_lean.flame.X[-1,hydrogen_keromnes_mult_lean.gas0.species_index('H2')]
X_O2_new = hydrogen_keromnes_mult_lean.flame.X[-1,hydrogen_keromnes_mult_lean.gas0.species_index('O2')]
phi_check = (X_H2_new / X_O2_new) / 2.0
print(f"Verify φ ≃ {phi_check:.3f}")

hydrogen_keromnes_mult_lean.flame_speed() #finds the velocity

print(f"the flame velocity for the lean mixture is : {hydrogen_keromnes_mult_lean.flame.velocity}")
profile = hydrogen_keromnes_mult_lean.flame.to_array()
fig, ax = plt.subplots()
ax.plot(profile.grid * 100, profile.T ,'.-')
ax.set_xlabel('space (cm)')
ax.set_ylabel('Temperature (K)')
ax.set_title(f"Temperature profile of the case {mech}")
plt.figure()
fix, ax = plt.subplots()
ax.plot(profile.grid * 100, profile('H2').X, "--", label='H$_2$')
ax.plot(profile.grid * 100, profile('O2').X, "--", label='O$_2$')
ax.plot(profile.grid * 100, profile('NO2').X, "--", label='NO$_2$')
ax.plot(profile.grid * 100, profile('NO').X, "--", label='NO')
ax.plot(profile.grid * 100, profile('H2O').X, "--", label='H$_2$O')
ax.legend(loc="best")
ax.set_xlabel('Distance [cm]')
ax.set_ylabel('Mole fraction [-]')
ax.set_title(f"reaction profile of the case {mech}")
plt.show()
toc('lean mixture')



#%% excercise d
""""[6 pts] Recompute the value of NO using equilibrium chemistry and the equation N2+O2 ⇌ 
2NO. Assume that this reaction is happening at the adiabatic flame temperature obtained 
for the L part of the combustor just calculated, and within a mixture approximated to 
consist only of water vapour, oxygen, and nitrogen, where: i) the value of mass fraction 
of water and oxygen is the same as in output from the L part of the combustor; and ii) 
the value of mass fraction of nitrogen is computed as 𝑌 = 1 − 𝑌 − 𝑌 . Show the 
passages used for the calculation. Compare the value of NO with that obtained from 
Cantera in point c) and explain the reasons of the differences (if any)"""

#key assumption : we have the same rate forward and backwards, need to find coefficients and that is it
T_ad = hydrogen_keromnes_mult_lean.flame.T[-1] #find the adiabatic temperature at which we have reaction

#now there is no need to compute R_out: we alrady know the only components of the reaction
Y_H2O =hydrogen_keromnes_mult_lean.get_mass_product('H2O')
Y_O2 =hydrogen_keromnes_mult_lean.get_mass_product('O2')
Y_N2 = 1- Y_H2O - Y_O2

P0 = 12e5
M = hydrogen_keromnes_mult_lean.gas0.molecular_weights       # [kg/kmol]
Ys = np.array([Y_H2O, Y_O2, Y_N2]) #find the index I want the molecular weights of
Ms = np.array([M[hydrogen_keromnes_mult_lean.gas0.species_index(s)] for s in ('H2O','O2','N2')])
Xs = (Ys/Ms) / np.sum(Ys/Ms) #manually find the composition
X_H2O_manual, X_O2_manual, X_N2_manual = Xs

X_H2O = hydrogen_keromnes_mult_lean.get_mass_product('H2O')
X_O2 = hydrogen_keromnes_mult_lean.get_mass_product('O2')
X_N2 = hydrogen_keromnes_mult_lean.get_mass_product('N2')
print (X_H2O_manual, X_H2O, X_O2_manual, X_O2, X_N2_manual, X_N2)

#we have the components, we now need the equilibrium function for the reaction
"""
hydrogen_adiabatic = PersonalizedFlame(mech ='keromnes2013.yaml',
                                       X_mix = {'H2O': X_H2O, 'O2': X_O2, 'N2': X_N2},
                                       T0=T_ad,
                                       P0=P0,
                                       tp_model='multicomponent')
                                       """
hydrogen_adiabatic = ct.Solution('keromnes2013.yaml')
hydrogen_adiabatic.TPX = T_ad, P0, {'H2O' : X_H2O, 'O2' : X_O2, 'N2' : X_N2}
#hydrogen_adiabatic.equilibrate('TP')
#X_NO = hydrogen_adiabatic['NO'].X[0]
#print(f"Molar fraction of NO using chemical equilibrium: {X_NO:.5e}")
#X_NO_actual = hydrogen_keromnes_mult_lean.get_mass_product('NO')
#print(f"actual NO molar fraction with Cantera:  {X_NO_actual}")


#%% point e
"""[6 pts] Now recalculate the flame in point (c) using strained premixed flamelets (using an 
opposed jets configuration) for two increasing values of strain. Compute the consumption speed 
and compare it to that of the freely propagating flame, and explain the reasons of any relevant 
difference.  

hydrogen_low_strained = PersonalizedFlame(mech ='keromnes2013.yaml',
                                                X_mix=X_mix,
                                                #reactant={'H2': 1.0},
                                                #oxidizer={'O2': 1.0, 'N2': 3.76},
                                                T0=675.0,
                                                P0=12e5,
                                                tp_model='multicomponent',
                                                CounterFlowFlame= True,
                                                reactants_mdot= 0.2,
                                                width=0.3)

hydrogen_high_strained = PersonalizedFlame(mech ='keromnes2013.yaml',
                                                X_mix=X_mix,
                                                #reactant={'H2': 1.0},
                                                #oxidizer={'O2': 1.0, 'N2': 3.76},
                                                T0=675.0,
                                                P0=12e5,
                                                tp_model='multicomponent',
                                                CounterFlowFlame= True,
                                                reactants_mdot= 0.6,
                                                width=0.3)


"""


#
L = hydrogen_keromnes_mult_lean.width  # width between the injectors set as the same as the width of the simulation

strain_rates_lst = [1.0, 100.0, 1500.0, 10000.0]  # values of what I want to iterate on[1/s]

#to compute the velocity with the formula U = a*L/2 I need the components of the strain
#we then need the mass flow which will be directly related to the density

rho_u = hydrogen_keromnes_mult_lean.gas0.density  # density, to compute the INLET MASS FLOW

results = []
for a in strain_rates_lst: #set the iteration loop

    U = a * L / 2.0  # calculate the inlet velocity from the formula said previously
    #note, the velocity depends on the strain, meaning that we will have different MASS FLOW

    mdot = rho_u * U
    print(mdot)


    print(f"\n--- Strain rate a = {a:.1f} results into → U = {U:.4f} m/s, mdot = {mdot:.4f} kg/m2/s ---")

    # creating the simulation object with the same parameters as point c except for the model
    #we can't do var_{variable} but we need to create the name and then assign it

    var_name = f"hydrogen_strain_{a}"
    try:
            inst = PersonalizedFlame(
                mech='keromnes2013.yaml',
                X_mix=X_mix,
                T0=675.0,
                P0=12e5,
                tp_model='multicomponent',
                CounterFlowFlame=True,
                reactants_mdot=mdot,
                width=0.3)
    except Exception as e:
            inst = None
            print(f"Error v={a}: {e}") #suggested by chatGPT, to avoid mistakes
            # dinamically assign the var
    globals()[var_name] = inst

    profile = inst.flame.to_array()
    fig, ax = plt.subplots()
    ax.plot(profile.grid * 100, profile.T, '.-')
    ax.set_xlabel('space (cm)')
    ax.set_ylabel('Temperature (K)')
    ax.set_title(f"Temperature profile of the case {mech}")
    plt.figure()
    fix, ax = plt.subplots()
    ax.plot(profile.grid * 100, profile('H2').X, "--", label='H$_2$')
    ax.plot(profile.grid * 100, profile('O2').X, "--", label='O$_2$')
    ax.plot(profile.grid * 100, profile('NO2').X, "--", label='NO$_2$')
    ax.plot(profile.grid * 100, profile('NO').X, "--", label='NO')
    ax.plot(profile.grid * 100, profile('H2O').X, "--", label='H$_2$O')
    ax.legend(loc="best")
    ax.set_xlabel('Distance [cm]')
    ax.set_ylabel('Mole fraction [-]')
    ax.set_title(f"reaction profile of the case {mech}")
    plt.show()
    toc('lean mixture')


    print(f"  Consumption speed determined from the graph =  m/s")

    freespeed = hydrogen_keromnes_mult_lean.flame_speed()  # the lean flame speed

    print(f"  Freely propagating flame speed (lean φ=0.5) = {freespeed:.4e} m/s")

    results.append({
        'a': a,
        'U': U,
        'mdot': mdot,
        'Sc': 0, #at the moment, discuss with others
        'Su_lean': freespeed
    })


#%% point f

"""
[10 pts] Using a freely-propagating flame again, recompute the flame in point c) by adding to 
the injected air an additional amount of water vapour. This could be the case of an exhaust gas 
recirculation (ERG) system.  Assume the water is reinjected at the same temperature of the 
exhaust gases found in point (c). Try different amounts of water (but pay attention the maximum 
you can inject is limited by the amount of water you have at the exhaust in point (c) and evaluate 
the effect on NO and adiabatic flame temperature at the exhaust. Comment in terms of both 
emissions of NO and thermal efficiency. """

#the maximum amount of water used is the one that is present in the system
# take everything from the exhaust
gas0_lean = hydrogen_keromnes_mult_lean.gas0
flame_lean = hydrogen_keromnes_mult_lean.flame
species = gas0_lean.species_names

#number of moles of hydrogen present, and water present, can be found like in point c)

H2_moles = hydrogen_keromnes_mult.get_product('H2') * n_base #already multiply by dummy index
H2O_moles = hydrogen_keromnes_mult.get_product('H2O') * n_base #already multiply by dummy index
O2_output = hydrogen_keromnes_mult_lean.get_product('O2') * n_base #initial value already multiply by dummy index
N2_output = hydrogen_keromnes_mult_lean.get_product('N2') * n_base #initial value already multiply by dummy index
#again, we could recompute the value of oxigen but it is not needed



n_O2_needed = H2_moles / n_base
n_air = n_O2_needed / 0.21
O2_moles = 0.21 * n_air
N2_moles = 0.79 * n_air

#in this case, differently from before, we assume that all the temperature is the same as lean
# Temperatura di ingresso = temperatura di exhaust lean
T_injected = output_T
# P0 =  12e5 as before

# create a list of fraction of total water available
f_values = [0.0, 0.25, 0.5, 0.75, 1.0]
f_names =[0, 25, 50, 75, 100]
results = []
for f in f_values:

    # create the empty mix dictionary (changes at every computation)
    moles_mix = {}

    H2O_added = f * H2O_moles
    moles_mix['H2O'] = moles_mix.get('H2O', 0.0) + H2O_added  # add the water

    # Add the air to mixture (H2 is already in the correct part)
    moles_mix['O2'] = O2_output + O2_moles
    moles_mix['N2'] = N2_output + N2_moles

    # create from mole to X_mix (should be same since 1 mole)
    n_total = sum(moles_mix.values())
    X_mix = { sp: moles_mix[sp] / n_total for sp in moles_mix }


    # solve the flame creating the name iteratively
    a = f_names[f]
    var_name = f"hydrogen_water_percentage_{a}"
    var_name = f"hydrogen_strain_{a}"
    try:
        inst = PersonalizedFlame(
            mech='keromnes2013.yaml',
            X_mix=X_mix,
            T0=675.0,
            P0=12e5,
            tp_model='multicomponent',
            CounterFlowFlame=True,
            reactants_mdot=mdot,
            width=0.3)
    except Exception as e:
        inst = None
        print(f"Error v={a}: {e}")  # suggested by chatGPT, to avoid mistakes
        # dinamically assign the var
    globals()[var_name] = inst

    profile = inst.flame.to_array()
    fig, ax = plt.subplots()
    ax.plot(profile.grid * 100, profile.T, '.-')
    ax.set_xlabel('space (cm)')
    ax.set_ylabel('Temperature (K)')
    ax.set_title(f"Temperature profile of the case {mech}")
    plt.figure()
    fix, ax = plt.subplots()
    ax.plot(profile.grid * 100, profile('H2').X, "--", label='H$_2$')
    ax.plot(profile.grid * 100, profile('O2').X, "--", label='O$_2$')
    ax.plot(profile.grid * 100, profile('NO2').X, "--", label='NO$_2$')
    ax.plot(profile.grid * 100, profile('NO').X, "--", label='NO')
    ax.plot(profile.grid * 100, profile('H2O').X, "--", label='H$_2$O')
    ax.legend(loc="best")
    ax.set_xlabel('Distance [cm]')
    ax.set_ylabel('Mole fraction [-]')
    ax.set_title(f"reaction profile of the case {mech}")
    plt.show()
    toc('lean mixture')

    # adiabatic flame temperature and NO fraction

    T_ad = inst.flame.T[-1]
    Y_NO = inst.get_product('NO')

    print(f"f={f:.2f}, H2O_added={H2O_added:.4f} mol: T_ad={T_ad:.1f} K, Y_NO={Y_NO:.3e}")

    results.append({
        'f': f,
        'H2O_added': H2O_added,
        'T_ad': T_ad,
        'Y_NO': Y_NO
    })


print("\n--- results EGR water addition ---")
for r in results:
    print(f"  f={r['f']:.2f}: T_ad={r['T_ad']:.1f} K, Y_NO={r['Y_NO']:.3e}")


fs = [r['f'] for r in results]
T_ads = [r['T_ad'] for r in results]
Y_NOs = [r['Y_NO'] for r in results]

plt.figure()
plt.plot(fs, T_ads, 'o-')
plt.xlabel('f (water fraction EGR)')
plt.ylabel('Adiabatic flame T [K]')
plt.title('effect of H2O on T_ad')
plt.show()

plt.figure()
plt.plot(fs, Y_NOs, 'o-')
plt.xlabel('f (H2O fraction EGR)')
plt.ylabel('Y_NO as output')
plt.title('H2O effect on NO emissions')
plt.show()


