"""
System test on scheme convergence
=================================

The expected congergence rates are found from the refinement in dx.
"""
import numpy as np
import matplotlib.pyplot as plt
from pde_solver.solver import Domain, Conditions, advection_eqn_solve


def f(x):
    """Periodic function."""
    return np.sin(2*np.pi*x)


#  Parameters
plt.figure(figsize=(8, 6), dpi=115)
Nrange = np.logspace(2, 3, 10)
dt = 0.00005  # Select dt smaller than all dx to be tested
methods = ["luw", "lfm", "lwm"]
c_wave = 1

#  Analysis
for methoduse in methods:

    errs = [0.0]*len(Nrange)
    for i, N in enumerate(Nrange):

        #  Solve
        Grid = Domain(0, 1, dt, .1, Nx=int(N))
        Conds = Conditions(f, "pbc")
        Sol = advection_eqn_solve(Grid, Conds, c=c_wave, method=methoduse)

        #  Find errors
        err = 0
        for t_n in range(Grid.Nt):
            uex = f(Grid.xvals() - c_wave*t_n*Grid.dt)  # Exact solution
            err = (abs(uex - Sol.solndata[:, t_n])).max()
        errs[i] = err

    #  Plot
    errs = np.array(errs)
    plt.loglog(Nrange, errs, label=methoduse)
    plt.scatter(Nrange, errs)

plt.loglog(Nrange, (Nrange)**-1, ':', label="O(dx^1)")
plt.loglog(Nrange, (Nrange)**-2, ':', label="O(dx^2)")
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel("Nx", fontsize=16)
plt.ylabel("error", fontsize=16)
plt.legend(fontsize=11)
plt.grid()
plt.show()
