import cantera as ct
import matplotlib.pyplot as plt
import numpy as np











methane = ct.Solution('gri30.yaml')
phi = 1.4
T0 = 293
P0 = ct.one_atm

fuel_comp = {'CH4' : 1.0}
oxidizer = {'O2':1.0, 'N2':3.76}
methane.set_equivalence_ratio(phi=phi, fuel='H2', oxidizer={'O2':1.0, 'N2':3.76})
methane.TP = T0, P0

width = 0.12

#creating the flame
flame = ct.FreeFlame(methane, width=width)
flame.transport_model = 'mixture-averaged'
flame.set_refine_criteria(ratio=3, slope=0.07, curve=0.12)   #tared with a couple runs

#solving the flame
flame.solve(loglevel=1, auto=True)

#result visualization + comparison with all the other models
print(f'\nmixture-averaged flame speed = {flame.velocity[0]:7f} m/s\n')

#%%



#I wanted to use Miller_Bowman but not available yet and not worth the struggle to code it
hydrogen = ct.Solution('gri30.yaml')  # use keromnes2013.yaml

#other to compare to
hydrogen_comparison = ct.Solution('h2-konnov.yaml')

#comparison with the one suggested by them
hydrogen_comparison2 = ct.Solution('h2o2.yaml')

#setting operative conditions
phi = 5.0
T0  = 675.0           # K
P0  = 12e5            # Pa (12 bar)

#another way of doing it
#gas.TPX = Tin, p, reactants

#setting the equivalence ratio to match the conditions (T,P, equivalence ratio)
fuel_comp = {'H2' : 1.0}
oxidizer = {'O2':1.0, 'N2':3.76}
hydrogen.set_equivalence_ratio(phi=phi, fuel='H2', oxidizer={'O2':1.0, 'N2':3.76})
hydrogen.TP = T0, P0


#setting the values for the domain
width = 0.3

#creating the flame
flame = ct.FreeFlame(hydrogen, width=width)
flame.transport_model = 'multicomponent'
flame.set_refine_criteria(ratio=3, slope=0.07, curve=0.12)   #tared with a couple runs

#solving the flame
flame.solve(loglevel=1, auto=True)

#result visualization + comparison with all the other models
print(f'\nmixture-averaged flame speed = {flame.velocity[0]:7f} m/s\n')



#%%
# comparisons, none of these are necessary
flame2 = ct.FreeFlame(hydrogen_comparison, width=width)
flame2.transport_model = 'mixture-averaged'
flame2.set_refine_criteria(ratio=2, slope=0.03, curve=0.06)
flame2.solve(loglevel=1, auto=True)
print(f'\nmixture-averaged flame speed with mixture = {flame2.velocity[0]:7f} m/s\n')
plt.plot(flame.grid, flame.X[hydrogen.species_index('H2'), :], label='H2 (Multi)')
plt.plot(flame2.grid, flame2.X[hydrogen.species_index('H2'), :], label='H2 (Mix)')
plt.legend()
plt.xlabel('Position [cm]')
plt.ylabel('Mole fraction')
plt.title('Preferential diffusion effect on H2')
plt.show()

flame3 = ct.FreeFlame(hydrogen_comparison, width=width)
flame3.transport_model = 'multicomponent'
flame3.set_refine_criteria(ratio=2, slope=0.03, curve=0.06)
flame3.solve(loglevel=1, auto=True)
print(f'\nmixture-averaged flame speed with konne = {flame3.velocity[0]:7f} m/s\n')
plt.plot(flame3.grid, flame3.X[hydrogen_comparison2.species_index('H2'), :], label='H2 (Komnov)')
plt.plot(flame.grid, flame.X[hydrogen.species_index('H2'), :], label='H2 (Keromnes)')
plt.legend()
plt.xlabel('Position [cm]')
plt.ylabel('Mole fraction')
plt.title('Keromnes vs Komnov')
plt.show()



#plotting of the favourite one
profile = flame.to_array()
fig, ax = plt.subplots()
ax.plot(profile.grid*100, profile.T, ".-")
ax.set_xlabel('Distance [cm]')
ax.set_ylabel('Temperature [K]')
plt.show()

fix, ax = plt.subplots()
ax.plot(profile.grid*100, profile('H2').X, "--", label = 'H$_2$')
ax.plot(profile.grid*100, profile('O2').X, "--", label = 'O$_2$')
ax.plot(profile.grid*100, profile('N2').X, "--", label = 'N$_2$')
ax.plot(profile.grid*100, profile('H2O').X, "--", label = 'H$_2$O')
ax.legend(loc= "best")
ax.set_xlabel('Distance [cm]')
ax.set_ylabel('Mole fraction [-]')
plt.show()




#%%
#excercise b

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



#%%
#excercise c
"""


[6 pts] Now mix the composition in output from the R part with air in order to obtain an 
equivalence ratio 𝜙 = 0.5. Indicate the amount of air required and assume it is injected at 675 
K and 12 bar. Compute the freely-propagating 1D flat flame for the new condition and plot 
again the variation of NO, NO2, O2 and H2 mass fractions, and temperature, in progress 
variable space using the two definitions given earlier. Are the considerations about progress 
variable the same as before? Comment on the results. 




"need to model the mixing of the two, the temperature will be averaged, extract quantities"

composition_R = flame.X[:, -1] #find the composition of it
T_R = flame.T[-1] #find the temperature of it

mix = ct.Quantity(hydrogen)  # R part output, we get how much of it there is
mix.Y = flame.Y[:, -1]
mix.TP = T_R, P0

#intialization of the air to do the same calculation

air = ct.Solution('keromnes2013.yaml')
air.TP = 675.0, 12e5
air.set_equivalence_ratio(phi=0.0, fuel={}, oxidizer={'O2': 1.0, 'N2': 3.76})

"""


T_Rout = flame.T[-1]
X_Rout = flame.X[-1, :].copy()   #  mole fractions (length matches gas.species())
# Y_Rout = flame_rich.Y[-1, :].copy() # you do NOT need this for mixing (we only need X_Rout)

# ====================
# EXERCISE 3: Build φ=0.5 (lean) by mixing R‐outlet + air
# ====================
# 1) Create a new Solution with the SAME mechanism, then assign TPX = T_Rout, P0, X_Rout:
gas_tmp = ct.Solution('gri30.yaml')
gas_tmp.TPX = T_Rout, P0, X_Rout




# -------------------------------------------------
# 2) COSTRUISCO LA MISCELA “LEAN” φ = 0.5
# -------------------------------------------------

# 2.1 index of H2
idx_H2 = gas_tmp.species_index('H2')
# N moles in H2
n_basis     = 1.0                     # prendo 1 kmol di prodotti R per base
n_H2_R      = X_Rout[idx_H2] * n_basis

# 2.2 Poiché φ = 0.5 ⇒ (n_H2 / n_O2) = 1 ⇒ n_O2_aggiunto = n_H2_R
n_O2_needed = n_H2_R

# 2.3 L’aria è 21% O2 e 79% N2 in frazioni molari
n_air       = n_O2_needed / 0.21      # kmol di aria
n_O2_air    = n_air * 0.21
n_N2_air    = n_air * 0.79

# 2.4 Estraggo i moli di ciascuna specie “R” (1 kmol base)
moles_R = {}
for k, sp in enumerate(gas_tmp.species()):
    moles_R[sp.name] = X_Rout[k] * n_basis

# 2.5 Ora costruisco i moli totali “mischiati” = prodotti R + aria
moles_mix = moles_R.copy()
# Aggiungo O2 e N2 provenienti dall’aria
moles_mix['O2'] = moles_mix.get('O2', 0.0) + n_O2_air
moles_mix['N2'] = moles_mix.get('N2', 0.0) + n_N2_air

# 2.6 Calcolo moli totali finali
n_total = sum(moles_mix.values())

# 2.7 Costruisco il dizionario delle frazioni molari finali X_mix
X_mix = { sp: (moles_mix[sp] / n_total) for sp in moles_mix.keys() }

# 2.8 Creo l’oggetto Cantera per la nuova miscela lean
gas_lean = ct.Solution('gri30.yaml')
gas_lean.TPX = 675.0, P0, X_mix   # T=675 K, P=12 bar

# (Facoltativo) Controllo rapido su phi:
# phi = ( X_H2 / X_O2 ) / 2   → a valori intermedi,
# se tutto è andato a buon fine, ottengo phi ≈ 0.5
X_H2_new = gas_lean.X[gas_lean.species_index('H2')]
X_O2_new = gas_lean.X[gas_lean.species_index('O2')]
phi_check = (X_H2_new / X_O2_new) / 2.0
print(f"Verifica φ della miscela lean: φ ≃ {phi_check:.3f} (dovrebbe essere 0.5)")



# -------------------------------------------------
# 3) RISOLVO LA FIAMMA “LEAN” (φ = 0.5)
# -------------------------------------------------
width2 = 0.03   # usiamo di nuovo 3 cm per il dominio 1D
flame_lean = ct.FreeFlame(gas_lean, width=width2)
flame_lean.transport_model = 'mixture-averaged'
flame_lean.set_refine_criteria(ratio=2, slope=0.03, curve=0.06)
flame_lean.set_max_grid_points(500)
flame_lean.solve(loglevel=0, auto=True)

# Profilo di T e Y nel dominio lean
x_lean = flame_lean.grid.copy()     # ascissa [m]
T_lean = flame_lean.T.copy()        # temperatura [K]
Y_lean = flame_lean.Y.copy()        # frazioni in massa [n_punti × n_specie]

# Estrazione indici specie
idx_H2_lean  = gas_lean.species_index('H2')
idx_H2O_lean = gas_lean.species_index('H2O')
idx_NO_lean  = gas_lean.species_index('NO')
idx_NO2_lean = gas_lean.species_index('NO2')
idx_O2_lean  = gas_lean.species_index('O2')

# Vettori “a colonna” per ogni specie di interesse
Y_H2_lean  = Y_lean[:, idx_H2_lean]
Y_H2O_lean = Y_lean[:, idx_H2O_lean]
Y_NO_lean  = Y_lean[:, idx_NO_lean]
Y_NO2_lean = Y_lean[:, idx_NO2_lean]
Y_O2_lean  = Y_lean[:, idx_O2_lean]

# Valori nei reagenti (punto 0) e nei prodotti (ultimo punto)
Y_H2_init_lean   = Y_H2_lean[0]
Y_H2_final_lean  = Y_H2_lean[-1]
Y_H2O_init_lean  = Y_H2O_lean[0]
Y_H2O_final_lean = Y_H2O_lean[-1]

# Definisco le due progress variable (lean):
c_H2_lean  = (Y_H2_init_lean - Y_H2_lean)  / (Y_H2_init_lean - Y_H2_final_lean)
c_H2O_lean = (Y_H2O_lean - Y_H2O_init_lean) / (Y_H2O_final_lean - Y_H2O_init_lean)

# Ordino gli indici in funzione di c crescente
idx_sort_H2_lean  = np.argsort(c_H2_lean)
idx_sort_H2O_lean = np.argsort(c_H2O_lean)


T_Rout = flame.T[-1]
X_Rout = flame.X[-1, :].copy()   #  mole fractions (length matches gas.species())
# Y_Rout = flame_rich.Y[-1, :].copy() # you do NOT need this for mixing (we only need X_Rout)

# ====================
# EXERCISE 3: Build φ=0.5 (lean) by mixing R‐outlet + air
# ====================
# 1) Create a new Solution with the SAME mechanism, then assign TPX = T_Rout, P0, X_Rout:
gas_tmp = ct.Solution('gri30.yaml')
gas_tmp.TPX = T_Rout, P0, X_Rout




# -------------------------------------------------
# 2) COSTRUISCO LA MISCELA “LEAN” φ = 0.5
# -------------------------------------------------

# 2.1 Indice della specie H2 nel meccanismo
idx_H2 = gas_tmp.species_index('H2')
# Numero di moli H2 in 1 kmol di “R‐outlet”
n_basis     = 1.0                     # prendo 1 kmol di prodotti R per base
n_H2_R      = X_Rout[idx_H2] * n_basis

# 2.2 Poiché φ = 0.5 ⇒ (n_H2 / n_O2) = 1 ⇒ n_O2_aggiunto = n_H2_R
n_O2_needed = n_H2_R

# 2.3 L’aria è 21% O2 e 79% N2 in frazioni molari
n_air       = n_O2_needed / 0.21      # kmol di aria
n_O2_air    = n_air * 0.21
n_N2_air    = n_air * 0.79

# 2.4 Estraggo i moli di ciascuna specie “R” (1 kmol base)
moles_R = {}
for k, sp in enumerate(gas_tmp.species()):
    moles_R[sp.name] = X_Rout[k] * n_basis

# 2.5 Ora costruisco i moli totali “mischiati” = prodotti R + aria
moles_mix = moles_R.copy()
# Aggiungo O2 e N2 provenienti dall’aria
moles_mix['O2'] = moles_mix.get('O2', 0.0) + n_O2_air
moles_mix['N2'] = moles_mix.get('N2', 0.0) + n_N2_air

# 2.6 Calcolo moli totali finali
n_total = sum(moles_mix.values())

# 2.7 Costruisco il dizionario delle frazioni molari finali X_mix
X_mix = { sp: (moles_mix[sp] / n_total) for sp in moles_mix.keys() }

# 2.8 Creo l’oggetto Cantera per la nuova miscela lean
gas_lean = ct.Solution('gri30.yaml')
gas_lean.TPX = 675.0, P0, X_mix   # T=675 K, P=12 bar

# (Facoltativo) Controllo rapido su phi:
# phi = ( X_H2 / X_O2 ) / 2   → a valori intermedi,
# se tutto è andato a buon fine, ottengo phi ≈ 0.5
X_H2_new = gas_lean.X[gas_lean.species_index('H2')]
X_O2_new = gas_lean.X[gas_lean.species_index('O2')]
phi_check = (X_H2_new / X_O2_new) / 2.0
print(f"Verifica φ della miscela lean: φ ≃ {phi_check:.3f} (dovrebbe essere 0.5)")



# -------------------------------------------------
# 3) RISOLVO LA FIAMMA “LEAN” (φ = 0.5)
# -------------------------------------------------
width2 = 0.03   # usiamo di nuovo 3 cm per il dominio 1D
flame_lean = ct.FreeFlame(gas_lean, width=width2)
flame_lean.transport_model = 'mixture-averaged'
flame_lean.set_refine_criteria(ratio=2, slope=0.03, curve=0.06)
flame_lean.set_max_grid_points(500)
flame_lean.solve(loglevel=0, auto=True)

# Profilo di T e Y nel dominio lean
x_lean = flame_lean.grid.copy()     # ascissa [m]
T_lean = flame_lean.T.copy()        # temperatura [K]
Y_lean = flame_lean.Y.copy()        # frazioni in massa [n_punti × n_specie]

# Estrazione indici specie
idx_H2_lean  = gas_lean.species_index('H2')
idx_H2O_lean = gas_lean.species_index('H2O')
idx_NO_lean  = gas_lean.species_index('NO')
idx_NO2_lean = gas_lean.species_index('NO2')
idx_O2_lean  = gas_lean.species_index('O2')

# Vettori “a colonna” per ogni specie di interesse
Y_H2_lean  = Y_lean[:, idx_H2_lean]
Y_H2O_lean = Y_lean[:, idx_H2O_lean]
Y_NO_lean  = Y_lean[:, idx_NO_lean]
Y_NO2_lean = Y_lean[:, idx_NO2_lean]
Y_O2_lean  = Y_lean[:, idx_O2_lean]

# Valori nei reagenti (punto 0) e nei prodotti (ultimo punto)
Y_H2_init_lean   = Y_H2_lean[0]
Y_H2_final_lean  = Y_H2_lean[-1]
Y_H2O_init_lean  = Y_H2O_lean[0]
Y_H2O_final_lean = Y_H2O_lean[-1]

# Definisco le due progress variable (lean):
c_H2_lean  = (Y_H2_init_lean - Y_H2_lean)  / (Y_H2_init_lean - Y_H2_final_lean)
c_H2O_lean = (Y_H2O_lean - Y_H2O_init_lean) / (Y_H2O_final_lean - Y_H2O_init_lean)

# Ordino gli indici in funzione di c crescente
idx_sort_H2_lean  = np.argsort(c_H2_lean)
idx_sort_H2O_lean = np.argsort(c_H2O_lean)


