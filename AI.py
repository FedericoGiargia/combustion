import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# 1) Carico il meccanismo (ad esempio GRI-3.0, che contiene H2/O2)
gas = ct.Solution('gri30.yaml')

# 2) Imposto phi = 5.0, T0 = 675 K, P0 = 12 bar
phi = 5.0
T0  = 675.0           # K
P0  = 12e5            # Pa (12 bar)

# 3) Seleziono il carburante (H2) e l’ossidante (aria approssimata come O2 + 3.76 N2)
gas.set_equivalence_ratio(phi,
                          fuel={'H2': 1.0},
                          oxidizer={'O2': 1.0, 'N2': 3.76})
gas.TP = T0, P0

# 4) Definisco la larghezza del dominio (ad es. 3 cm) e creo l’oggetto FreeFlame
width = 0.03    # 3 cm di dominio, tipico punto di partenza
flame_rich = ct.FreeFlame(gas, width=width)
flame_rich.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)
flame_rich.solve(loglevel=1, auto=True)
print(f"\nmixture-averaged flame speed = {flame_rich.velocity[0]:7f} m/s\n")



#excercise 2
# -------------------------------------------------
# 2) ESTRAZIONE DEI PROFILI DI SPECI E TEMPERATURA
# -------------------------------------------------

# 2.1 Coordinate spaziali (non ci serve per il progress variable)
x = flame_rich.grid.copy()             # array delle posizioni [m]
T_profile = flame_rich.T.copy()        # array di temperatura [K]
Y = flame_rich.Y.copy()                # matrice [n_punti × n_specie] di frazioni in massa

# 2.2 Indici delle specie d’interesse
idx_H2  = gas.species_index('H2')
idx_H2O = gas.species_index('H2O')
idx_NO  = gas.species_index('NO')
idx_NO2 = gas.species_index('NO2')
idx_O2  = gas.species_index('O2')

# 2.3 Costruisco vettori “a colonna” delle frazioni in massa
Y_H2  = Y[:, idx_H2]
Y_H2O = Y[:, idx_H2O]
Y_NO  = Y[:, idx_NO]
Y_NO2 = Y[:, idx_NO2]
Y_O2  = Y[:, idx_O2]

# 2.4 Valori iniziali (nello strato dei reagenti) e finali (nei prodotti)
Y_H2_init   = Y_H2[0]
Y_H2_final  = Y_H2[-1]
Y_H2O_init  = Y_H2O[0]
Y_H2O_final = Y_H2O[-1]


# -------------------------------------------------
# 3) DEFINIZIONE DELLE DUE PROGRESS VARIABLE
# -------------------------------------------------

# 3.1 Progress variable basata sul combustibile (H2)
#     c_H2 = ( Y_H2(reactants) - Y_H2(x) ) / ( Y_H2(reactants) - Y_H2(products) )
c_H2 = (Y_H2_init - Y_H2) / (Y_H2_init - Y_H2_final)

# 3.2 Progress variable basata sul prodotto H2O
#     c_H2O = ( Y_H2O(x) - Y_H2O(reactants) ) / ( Y_H2O(products) - Y_H2O(reactants) )
c_H2O = (Y_H2O - Y_H2O_init) / (Y_H2O_final - Y_H2O_init)


# -------------------------------------------------
# 4) ORDINAMENTO PER PROGRESS VALUE CRESCENTE
# -------------------------------------------------

# Poiché vogliamo plottare “Y vs c” in ordine crescente di c,
# ordiniamo gli indici in funzione di ciascuna c.

idx_sort_H2  = np.argsort(c_H2)
idx_sort_H2O = np.argsort(c_H2O)


# -------------------------------------------------
# 5) PLOT DELLE FRAZIONI IN MASSA
# -------------------------------------------------

plt.figure(figsize=(12, 5))

# 5.1 PRIMA COLONNA: vs c_H2
ax1 = plt.subplot(1, 2, 1)
ax1.plot(c_H2[idx_sort_H2], Y_H2[idx_sort_H2],  label='H$_2$')
ax1.plot(c_H2[idx_sort_H2], Y_O2[idx_sort_H2],  label='O$_2$')
ax1.plot(c_H2[idx_sort_H2], Y_NO[idx_sort_H2],  label='NO')
ax1.plot(c_H2[idx_sort_H2], Y_NO2[idx_sort_H2], label='NO$_2$')
ax1.set_xlabel('Progress variable $c$ (basata su H$_2$)')
ax1.set_ylabel('Frazione in massa')
ax1.set_title('Specie in funzione di $c$ (fuel‐based)')
ax1.legend()
ax1.grid(True)

# 5.2 SECONDA COLONNA: vs c_H2O
ax2 = plt.subplot(1, 2, 2)
ax2.plot(c_H2O[idx_sort_H2O], Y_H2[idx_sort_H2O],  label='H$_2$')
ax2.plot(c_H2O[idx_sort_H2O], Y_O2[idx_sort_H2O],  label='O$_2$')
ax2.plot(c_H2O[idx_sort_H2O], Y_NO[idx_sort_H2O],  label='NO')
ax2.plot(c_H2O[idx_sort_H2O], Y_NO2[idx_sort_H2O], label='NO$_2$')
ax2.set_xlabel('Progress variable $c$ (basata su H$_2$O)')
ax2.set_ylabel('Frazione in massa')
ax2.set_title('Specie in funzione di $c$ (product‐based)')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()


# -------------------------------------------------
# 6) PLOT DELLA TEMPERATURA
# -------------------------------------------------

plt.figure(figsize=(12, 5))

ax3 = plt.subplot(1, 2, 1)
ax3.plot(c_H2[idx_sort_H2], T_profile[idx_sort_H2], 'r-')
ax3.set_xlabel('Progress variable $c$ (basata su H$_2$)')
ax3.set_ylabel('Temperatura [K]')
ax3.set_title('Temperatura vs $c$ (fuel‐based)')
ax3.grid(True)

ax4 = plt.subplot(1, 2, 2)
ax4.plot(c_H2O[idx_sort_H2O], T_profile[idx_sort_H2O], 'r-')
ax4.set_xlabel('Progress variable $c$ (basata su H$_2$O)')
ax4.set_ylabel('Temperatura [K]')
ax4.set_title('Temperatura vs $c$ (product‐based)')
ax4.grid(True)

plt.tight_layout()
plt.show()


#%%
#excercise 3
# Estraggo composizione (Y) e temperatura (T) nell’ultimo nodo (prodotti)
# Extract R‐outlet (last node)
T_Rout = flame_rich.T[-1]
X_Rout = flame_rich.X[-1, :].copy()   # 🔑 mole fractions (length matches gas.species())
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



# -------------------------------------------------
# 4) PLOT SPECIE (NO, NO2, O2, H2) E T in funzione di c
# -------------------------------------------------
plt.figure(figsize=(12, 10))

# 4.1 Specie vs c_H2 (lean)
ax1 = plt.subplot(2, 2, 1)
ax1.plot(c_H2_lean[idx_sort_H2_lean], Y_H2_lean[idx_sort_H2_lean],  label='H$_2$',  linewidth=1.5)
ax1.plot(c_H2_lean[idx_sort_H2_lean], Y_O2_lean[idx_sort_H2_lean],  label='O$_2$',  linewidth=1.5)
ax1.plot(c_H2_lean[idx_sort_H2_lean], Y_NO_lean[idx_sort_H2_lean],  label='NO',     linewidth=1.5)
ax1.plot(c_H2_lean[idx_sort_H2_lean], Y_NO2_lean[idx_sort_H2_lean], label='NO$_2$', linewidth=1.5)
ax1.set_xlabel('Progress variable $c$ (fuel‐based)')
ax1.set_ylabel('Frazione in massa')
ax1.set_title('Specie vs $c$ (fuel‐based), φ = 0.5')
ax1.legend()
ax1.grid(True)

# 4.2 Specie vs c_H2O (lean)
ax2 = plt.subplot(2, 2, 2)
ax2.plot(c_H2O_lean[idx_sort_H2O_lean], Y_H2_lean[idx_sort_H2O_lean],  label='H$_2$',  linewidth=1.5)
ax2.plot(c_H2O_lean[idx_sort_H2O_lean], Y_O2_lean[idx_sort_H2O_lean],  label='O$_2$',  linewidth=1.5)
ax2.plot(c_H2O_lean[idx_sort_H2O_lean], Y_NO_lean[idx_sort_H2O_lean],  label='NO',     linewidth=1.5)
ax2.plot(c_H2O_lean[idx_sort_H2O_lean], Y_NO2_lean[idx_sort_H2O_lean], label='NO$_2$', linewidth=1.5)
ax2.set_xlabel('Progress variable $c$ (water‐based)')
ax2.set_ylabel('Frazione in massa')
ax2.set_title('Specie vs $c$ (water‐based), φ = 0.5')
ax2.legend()
ax2.grid(True)

# 4.3 Temperatura vs c_H2 (lean)
ax3 = plt.subplot(2, 2, 3)
ax3.plot(c_H2_lean[idx_sort_H2_lean], T_lean[idx_sort_H2_lean], 'r-', linewidth=1.5)
ax3.set_xlabel('Progress variable $c$ (fuel‐based)')
ax3.set_ylabel('Temperatura [K]')
ax3.set_title('T vs $c$ (fuel‐based), φ = 0.5')
ax3.grid(True)

# 4.4 Temperatura vs c_H2O (lean)
ax4 = plt.subplot(2, 2, 4)
ax4.plot(c_H2O_lean[idx_sort_H2O_lean], T_lean[idx_sort_H2O_lean], 'r-', linewidth=1.5)
ax4.set_xlabel('Progress variable $c$ (water‐based)')
ax4.set_ylabel('Temperatura [K]')
ax4.set_title('T vs $c$ (water‐based), φ = 0.5')
ax4.grid(True)

plt.tight_layout()
plt.show()

# 4.5 Indico la quantità di aria richiesta:
print(f"Per 1 kmol di prodotto R, ho aggiunto n_air = {n_air:.4f} kmol di aria (a 675 K, 12 bar)")


#%%
#excercise 4
# 1) Creo una solution qualsiasi (per esempio uso encore gri30, basta che contenga N2, O2, NO)
gas = ct.Solution('gri30.yaml')

# 2) Per calcolare Kp per N2+O2<=>2NO, utilizziamo la routine equilibrium_constants():
#    equilibria[‘N2:O2:2NO’] corrisponde (N2 + O2 <=> 2 NO).
#    Attenzione: l’ordine delle specie deve rispettare la stechiometria nel file CTI.
Kp_dict = gas.equilibrium_constants              # dizionario di tutte le reazioni
# Trovo l’indice della reazione
rxn_index = gas.reaction_equations().index('N2 + O2 <=> 2 NO')
# oppure in modo piu’ robusto:
for i, eq in enumerate(gas.reaction_equations()):
    if eq.strip() == 'N2 + O2 <=> 2 NO':
        rxn_index = i

# 3) Imposto T e P nel gas (il valore di pressione non conta per Kp, ma serve per far calcolare le proprietà)
T_ad = 3500.0    # esempio: se so che T_ad ~ 3500 K
gas.TP = T_ad, 12e5

# 4) Estraggo Kp(“N2 + O2 <=> 2 NO”) a temperatura T_ad:
Kp_value = Kp_dict[rxn_index]    # questo è Kp in unità di Pa^(-1).
# Se volessi in bar^(-1), divido per (1e5): Kp_bar = Kp_value * 1e5.

#%%
#excercise 5
# 1. Define mechanism and lean mixture (φ = 0.5, T0 = 675 K, P0 = 12 bar)
gas = ct.Solution('gri30.yaml')
phi = 0.5
T0 = 675.0
P0 = 12e5
gas.set_equivalence_ratio(phi,
                          fuel={'H2': 1.0},
                          oxidizer={'O2': 1.0, 'N2': 3.76})
gas.TP = T0, P0

# 2. Freely propagating flat flame (point c)
width = 0.03  # 3 cm
free_flame = ct.FreeFlame(gas, width=width)
free_flame.transport_model = 'mixture-averaged'
free_flame.set_refine_criteria(ratio=2, slope=0.03, curve=0.06)
free_flame.set_max_grid_points(500)
free_flame.solve(loglevel=0, auto=True)
Su_free = free_flame.velocity[0]

# 3. Strained premixed flamelets (opposed jets) for two strain rates
strains = [100.0, 200.0]  # [1/s]
Su_strained = []

for a in strains:
    cf = ct.CounterflowPremixedFlame(gas=gas, width=width)
    cf.transport_model = 'mixture-averaged'
    cf.set_refine_criteria(ratio=2, slope=0.03, curve=0.06)
    cf.strain_rate = a
    cf.solve(loglevel=0, auto=True)
    rho_un = gas.density
    Su = cf.mass_flux / rho_un
    Su_strained.append(Su)

print(f"Freely propagating Su = {Su_free:.6f} m/s")
for a, Su in zip(strains, Su_strained):
    print(f"Strain rate = {a:.1f} 1/s → Su_strained = {Su:.6f} m/s")


#%%
#excercise f
# -----------------------------
# 1) Solve the lean flame (φ = 0.5) to get exhaust composition and temperature
# -----------------------------
gas_lean = ct.Solution('gri30.yaml')
phi_lean = 0.5
T0_lean  = 675.0       # K
P0       = 12e5        # Pa

# Build lean inlet by mixing R‐outlet with air (as in point c)
gas_rich = ct.Solution('gri30.yaml')
phi_rich = 5.0
T0_rich  = 675.0
gas_rich.set_equivalence_ratio(phi_rich,
                               fuel={'H2': 1.0},
                               oxidizer={'O2': 1.0, 'N2': 3.76})
gas_rich.TP = T0_rich, P0

flame_rich = ct.FreeFlame(gas_rich, width=0.03)
flame_rich.transport_model = 'mixture-averaged'
flame_rich.set_refine_criteria(ratio=2, slope=0.03, curve=0.06)
flame_rich.set_max_grid_points(500)
flame_rich.solve(loglevel=0, auto=True)

Y_Rout = flame_rich.Y[-1, :].copy()
T_Rout = flame_rich.T[-1]

gas_tmp = ct.Solution('gri30.yaml')
gas_tmp.TPX = T_Rout, P0, Y_Rout
X_Rout = gas_tmp.X.copy()

idx_H2 = gas_tmp.species_index('H2')
n_H2_R = X_Rout[idx_H2] * 1.0
n_O2_needed = n_H2_R
n_air = n_O2_needed / 0.21
n_O2_air  = 0.21 * n_air
n_N2_air  = 0.79 * n_air

moles_R = {sp.name: X_Rout[k] * 1.0 for k, sp in enumerate(gas_tmp.species())}

moles_lean = moles_R.copy()
moles_lean['O2'] = moles_lean.get('O2', 0.0) + n_O2_air
moles_lean['N2'] = moles_lean.get('N2', 0.0) + n_N2_air

n_total_lean = sum(moles_lean.values())
X_lean_inlet = {sp: moles_lean[sp] / n_total_lean for sp in moles_lean}

gas_lean.TPX = T0_lean, P0, X_lean_inlet

flame_lean = ct.FreeFlame(gas_lean, width=0.03)
flame_lean.transport_model = 'mixture-averaged'
flame_lean.set_refine_criteria(ratio=2, slope=0.03, curve=0.06)
flame_lean.set_max_grid_points(500)
flame_lean.solve(loglevel=0, auto=True)

X_exit_lean = gas_lean.X.copy()
Y_exit_lean = flame_lean.Y[-1, :].copy()
T_exit_lean = flame_lean.T[-1]

idx_H2O = gas_lean.species_index('H2O')
X_H2O_exit = X_exit_lean[idx_H2O]
Y_H2O_exit = Y_exit_lean[idx_H2O]

# -----------------------------
# 2) Build a set of water‐injection cases
# -----------------------------
alpha_list = [0.0, 0.5, 1.0]  # 0%, 50%, 100% of the H2O in 1 kmol of lean exhaust
results = []

for alpha in alpha_list:
    n_H2O_add = alpha * X_H2O_exit

    moles_base = {sp: X_lean_inlet.get(sp, 0.0) * 1.0 for sp in X_lean_inlet}

    moles_mix = moles_base.copy()
    moles_mix['H2O'] = moles_mix.get('H2O', 0.0) + n_H2O_add

    n_total_mix = sum(moles_mix.values())
    X_new = {sp: moles_mix[sp] / n_total_mix for sp in moles_mix}

    gas_egr = ct.Solution('gri30.yaml')
    gas_egr.TPX = T0_lean, P0, X_new

    flame_egr = ct.FreeFlame(gas_egr, width=0.03)
    flame_egr.transport_model = 'mixture-averaged'
    flame_egr.set_refine_criteria(ratio=2, slope=0.03, curve=0.06)
    flame_egr.set_max_grid_points(500)
    flame_egr.solve(loglevel=0, auto=True)

    Y_egr_out = flame_egr.Y[-1, gas_egr.species_index('NO')]
    T_egr_peak = np.max(flame_egr.T)

    results.append((alpha, Y_egr_out, T_egr_peak))

# -----------------------------
# 3) Plot results
# -----------------------------
alphas = [r[0] for r in results]
Y_NO_outs = [r[1] for r in results]
T_peaks = [r[2] for r in results]

plt.figure(figsize=(10, 4))

# Plot NO mass fraction vs alpha
plt.subplot(1, 2, 1)
plt.plot(alphas, Y_NO_outs, marker='o', linewidth=2)
plt.xlabel('α (fraction of exhaust H2O reinjected)')
plt.ylabel('Y_NO at outlet')
plt.title('NO Emissions vs H2O Injection')
plt.grid(True)

# Plot adiabatic flame temperature vs alpha
plt.subplot(1, 2, 2)
plt.plot(alphas, T_peaks, marker='o', linewidth=2)
plt.xlabel('α (fraction of exhaust H2O reinjected)')
plt.ylabel('Adiabatic Flame Temperature [K]')
plt.title('Flame Temperature vs H2O Injection')
plt.grid(True)

plt.tight_layout()
plt.show()

