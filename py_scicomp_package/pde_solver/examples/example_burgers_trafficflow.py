"""
Burgers equation applied to traffic flow example
================================================

The first figure is the plot of u against miles. This is what the solver
solves. This is then converted to vehicle density (plot 2) and then vehicle
speed (plot 3).
"""
from pde_solver.solver import Domain, Conditions, burgers_eqn_solve
import matplotlib.pyplot as plt
import numpy as np

#   ---------------- PARAMETERS ----------------
vmax = 50  # Maximum speed
rhomax = 38  # Maximum density
rho1 = 19  # Initial 2 lane density
rho2 = 21  # Initial 1 lane density
miles = 10  # Road length


def f(x):
    """Function (for IC). Input must be array."""
    u1 = 1 - 2*rho1/rhomax  # Transform
    u2 = 1 - 2*rho2/rhomax  # Transform
    IC = u1*np.ones(len(x))
    IC[int(len(IC)*0.4):int(len(IC)*0.7)] = u2  # Lane closure
    return IC


#  ---------------- SOLVE ----------------
Grid = Domain(0, miles, 0.001, 10, dx=0.001)  # Input domain specifications
Conds = Conditions(f, "pbc")  # Initial and boundary condition functions
Sol = burgers_eqn_solve(Grid, Conds, method="lfm")  # Solve with method
Sol.plot(6)  # plots

#   ---------------- PLOT ----------------
plt.figure(figsize=(8, 6), dpi=120)
for t in np.linspace(0, Grid.T, 6):
    ts = str(round(t, 3))
    plt.plot(Grid.xvals(), rhomax*(1-Sol(t))/2, label="t = "+ts, linewidth=2.3)
plt.ylabel("Vehicle density", fontsize=16)
plt.xlabel("Mile", fontsize=16)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid()
plt.title("Vehicle Density vs Road Mile", fontsize=16)
plt.legend(fontsize=11)
plt.show()

plt.figure(figsize=(8, 6), dpi=120)
for t in np.linspace(0, Grid.T, 6):
    ts = str(round(t, 3))
    rho = (rhomax*(1-Sol(t))/2)
    plt.plot(Grid.xvals(), vmax*(1-rho/rhomax), label="t = "+ts, linewidth=2.3)
plt.ylabel("Vehicle speed", fontsize=16)
plt.xlabel("Mile", fontsize=16)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid()
plt.title("Vehicle Speed vs Road Mile", fontsize=16)
plt.legend(fontsize=11)
plt.show()
