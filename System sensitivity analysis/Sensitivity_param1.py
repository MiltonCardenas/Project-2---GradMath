# -*- coding: utf-8 -*-
"""
@author: milton

"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# model
def growth_ODE(t, y, p, pm, pw, pmw, km, kw, kr, qm, qw):
    T, Tm, Tw, Tmw = y
    dT_dt = (p - km * Tm - kw * Tw - kr * Tmw) * T
    dTm_dt = (pm + km * T - qm * Tw) * Tm + 0.25 * kr * Tmw * T
    dTw_dt = (pw + kw * T - qw * Tm) * Tw + 0.25 * kr * Tmw * T
    dTmw_dt = (pmw + 0.5 * kr * T) * Tmw + (qm + qw) * Tw * Tm
    return [dT_dt, dTm_dt, dTw_dt, dTmw_dt]

# Initial conditions
T0 = 32554830
Tm0 = 134173
Tw0 = 26180
Tmw0 = 9818
initial_conditions = [T0, Tm0, Tw0, Tmw0]

# Parameters from paper
param_paper = [1.5e-2, -2.29e-2, 7.13e-3, 5.68e-4, 1.51e-9, 1.11e-9, 4.36e-10, 4.15e-9, 1.10e-9]


# Time span
t_span = np.linspace(70, 163, 100)

# Sensitivity analysis parameter
diff_max = 0.2
pmax = param_paper[0] * (1 + diff_max)
pmin = param_paper[0] * (1 - diff_max)
p_values = np.linspace(pmin, pmax, 5)


# Parameter span
p_values = np.linspace(pmin, pmax, 5)


fig, axes = plt.subplots(2, 2, figsize=(16, 10))

for p in p_values:
    params = [p] + param_paper[1:]
    solution = solve_ivp(growth_ODE, [t_span[0], t_span[-1]], initial_conditions, t_eval=t_span, args=params)

    # Plot for each dependent variable in the 2x2 grid
    axes[0, 0].plot(solution.t, np.log10(solution.y[0]), label=f'p = {p:.2e}')  # T
    axes[0, 1].plot(solution.t, np.log10(solution.y[1]), label=f'p = {p:.2e}')  # Tm
    axes[1, 0].plot(solution.t, np.log10(solution.y[2]), label=f'p = {p:.2e}')  # Tw
    axes[1, 1].plot(solution.t, np.log10(solution.y[3]), label=f'p = {p:.2e}')  # Tmw


axes[0, 0].set_title('Sensitivity Analysis of Uninfected Cells')
axes[0, 0].set_ylabel('log10(T)')

axes[0, 1].set_title('Sensitivity Analysis of Cells Infected by Mutant Virus')
axes[0, 1].set_ylabel('log10(Tm)')

axes[1, 0].set_title('Sensitivity Analysis of Cells Infected by Wild-type Virus')
axes[1, 0].set_ylabel('log10(Tw)')
axes[1, 0].set_xlabel('Time')

axes[1, 1].set_title('Sensitivity Analysis of Cells Infected by both Virus')
axes[1, 1].set_ylabel('log10(Tmw)')
axes[1, 1].set_xlabel('Time')

# Adding legends and grid
for ax in axes.flat:
    ax.legend()
    ax.grid(True)

plt.tight_layout()
plt.show()

fig.savefig("Sensitivity analysis Param 1.png")

# ===================== Sensitivity Analysis expanded ===========================================

# Nominal trajectory
nominal_solution = solve_ivp(growth_ODE, [t_span[0], t_span[-1]], initial_conditions, t_eval=t_span, args=param_paper)
nominal_trajectories = nominal_solution.y

# 10% Perturbation
perturbed_trajectories = {}
perturbation_size = 0.1  # 10% perturbation

for i, param in enumerate(param_paper):

    perturbed_params = param_paper.copy()
    perturbed_params[i] *= 1.1

    #perturbed parameter
    perturbed_solution = solve_ivp(growth_ODE, [t_span[0], t_span[-1]], initial_conditions, t_eval=t_span, args=perturbed_params)
    perturbed_trajectories[i] = perturbed_solution.y

# To analyze Tm (Cells infected by mutant virus growth)
variable_index = 1  # Index for Tm

# P change
fractional_changes = {}
for i in perturbed_trajectories:
    fractional_change = (perturbed_trajectories[i][variable_index] - nominal_trajectories[variable_index]) / nominal_trajectories[variable_index]
    fractional_changes[i] = fractional_change

# Normalized sensitivity
normalized_sensitivities = {}
for i in fractional_changes:
    normalized_sensitivity = fractional_changes[i] * (param_paper[i] / (param_paper[i] * perturbation_size))
    normalized_sensitivities[i] = normalized_sensitivity

#
parameter_index = 0
normalized_sensitivity_example = normalized_sensitivities[parameter_index]

fig2, ax2 = plt.subplots(figsize=(10, 6))

ax2.plot(t_span, normalized_sensitivity_example)
ax2.set_xlabel('Time')
ax2.set_ylabel(f'Normalized Sensitivity of T changing p')
ax2.set_title('Local Sensitivity Analysis')
ax2.grid(True)

# Save the figure
fig2.savefig("Normalized sensitivity analysis Param 1.png")
