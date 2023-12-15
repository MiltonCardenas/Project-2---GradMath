import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from auxiliar_codes import *

# Experimental data log-transformed
# T, Tm, Tw, Tmw,
Run1 = np.log10(np.array([
    [32554830, 134173, 26180, 9818],
    [46645200, 481950, 103950, 18900],
    [64240540, 1230460, 309260, 26320],
    [65563680, 9863280, 3000480, 1364580],
    [36366400, 36545600, 10281600, 28806400]
]))

# T, Tm, Tw, Tmw,
Run2= np.log10(np.array([
    [35855330, 158620, 25235, 10815],
    [48652100, 269500, 73500, 9800],
    [62989640, 1081920, 302680, 25760],
    [79088100, 7907900, 3103100, 900900],
    [47349120, 24613680, 15167880, 22069320]
]))

# T, Tm, Tw, Tmw,
Run3 = np.log10(np.array([
    [32597373, 98175, 22908, 6545],
    [52059000, 315000, 110250, 15750],
    [62362300, 847210, 445900, 38220],
    [77218680, 5576620, 2117920, 478240],
    [49714560, 17922240, 17025120, 16138080]
]))

# Parameters bounds
bounds = [
    (-6.0e-2, 6.0e-2),  # p
    (-2.0e-1, 6.0e-2),  # pm
    (-6.0e-2, 6.0e-2),  # pw
    (-2.0e-1, 6.0e-2),  # pmw
    (0, 1.0e-8),        # km
    (0, 1.0e-8),        # kw
    (0, 1.0e-8),        # kr
    (0, 1.0e-8),        # qm
    (0, 1.0e-8)         # qw
]

# Initial conditions
T0 = 32554830
Tm0 = 134173
Tw0 = 26180
Tmw0 = 9818
initial_conditions = [T0, Tm0, Tw0, Tmw0]

# Initial guesses for the parameters
param_paper = [1.5e-2, -2.29e-2, 7.13e-3, 5.68e-4, 1.51e-9, 1.11e-9, 4.36e-10, 4.15e-9, 1.10e-9]

# initial guesses
difference = 0.05
param0 = vary_vector(param_paper, difference)

#param0 = [1e-2, 1e-2, 1e-2, 1e-2,  1e-8, 1e-8, 1e-8, 1e-8, 1e-8]

# Time span and measurement times
t_span = np.linspace(0, 95, 100)
t_measurements = 70 + np.array([0, 24, 45, 69, 93])

# ODE system
def growth_ODE(t, y, p, pm, pw, pmw, km, kw, kr, qm, qw):

    T, Tm, Tw, Tmw = y
    dT_dt = (p - km * Tm - kw * Tw - kr * Tmw) * T
    dTm_dt = (pm + km * T - qm * Tw) * Tm + 0.25 * kr * Tmw * T
    dTw_dt = (pw + kw * T - qw * Tm) * Tw + 0.25 * kr * Tmw * T
    dTmw_dt = (pmw + 0.5 * kr * T) * Tmw + (qm + qw) * Tw * Tm

    return [dT_dt, dTm_dt, dTw_dt, dTmw_dt]

# Objective function to minimize:
def objective_function(params, t_measurements, Run1, Run2, Run3, initial_conditions):

    # ODE solution
    sol = solve_ivp(growth_ODE, [t_measurements[0], t_measurements[-1]], initial_conditions, t_eval=t_measurements, args=params)

    # Errors between the model and the measurements
    residuals1 = np.log10(sol.y.T) - Run1
    residuals2 = np.log10(sol.y.T) - Run2
    residuals3 = np.log10(sol.y.T) - Run3

    # Squared errors
    return np.sum((residuals1+residuals2+residuals3)**2)

# Optimization using minimize function
result = minimize(objective_function, param0, args=(t_measurements, Run1, Run2, Run3, initial_conditions), bounds=bounds, method='L-BFGS-B')

params_opt = result.x
diff = np.round(np.abs((param_paper-params_opt)/param_paper)* 100, 3)
solution = solve_ivp(growth_ODE, [t_measurements[0], t_measurements[-1]], initial_conditions, t_eval=np.linspace(70, t_measurements[-1]), args=params_opt)

t = solution.t
y = solution.y

# ============================================== GRAPHICS =====================================================================

markers = ['o', 's', '^']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']


fig, ax = plt.subplots(figsize=(10, 8))
ax.plot(t, np.log10(y[0, :]), label='T')  # Plot for T
ax.plot(t, np.log10(y[1, :]), label='Tm')  # Plot for Tm
ax.plot(t, np.log10(y[2, :]), label='Tw')  # Plot for Tw
ax.plot(t, np.log10(y[3, :]), label='Tmw')  # Plot for Tmw

for i in [0, 1, 2, 3]:
    if i == 0:  # Add labels only for the first set
        ax.scatter(t_measurements, Run1[:, i], marker=markers[0], color='black', label='Run1')
        ax.scatter(t_measurements, Run2[:, i], marker=markers[1], color='black', label='Run2')
        ax.scatter(t_measurements, Run3[:, i], marker=markers[2], color='black', label='Run3')

    else:
        ax.scatter(t_measurements, Run1[:, i], marker=markers[0], color=colors[i])
        ax.scatter(t_measurements, Run2[:, i], marker=markers[1], color=colors[i])
        ax.scatter(t_measurements, Run3[:, i], marker=markers[2], color=colors[i])

ax.set_xlabel('Time (hours)')
ax.set_ylabel('Logarithm of HIV Population')
ax.set_title('HIV-1 Population Dynamics')
ax.legend()
fig.show()

plt.savefig("Fitted dynamic model.png")