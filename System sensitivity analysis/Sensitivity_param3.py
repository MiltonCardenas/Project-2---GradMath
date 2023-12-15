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

# Sensitivity analysis parameter
diff_max = 0.2
pmax = param_paper[2] * (1 + diff_max)
pmin = param_paper[2] * (1 - diff_max)
pw_values = np.linspace(pmin, pmax, 5)

# Time span
t_span = np.linspace(70, 163, 100)

# Parameter span
pw_values = np.linspace(pmin, pmax, 5)


fig, axes = plt.subplots(2, 2, figsize=(16, 10))

for pw in pw_values:
    params = [param_paper[0],param_paper[1], pw] + param_paper[3:]
    solution = solve_ivp(growth_ODE, [t_span[0], t_span[-1]], initial_conditions, t_eval=t_span, args=params)

    # Plot for each dependent variable in the 2x2 grid
    axes[0, 0].plot(solution.t, np.log10(solution.y[0]), label=f'pw = {pw:.2e}')  # T
    axes[0, 1].plot(solution.t, np.log10(solution.y[1]), label=f'pw = {pw:.2e}')  # Tm
    axes[1, 0].plot(solution.t, np.log10(solution.y[2]), label=f'pw = {pw:.2e}')  # Tw
    axes[1, 1].plot(solution.t, np.log10(solution.y[3]), label=f'pw = {pw:.2e}')  # Tmw


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

fig.savefig("Normalized sensitivity analysis Param 3.png")