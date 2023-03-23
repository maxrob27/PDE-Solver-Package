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
T = 1
dtrange = np.array([1/(2**n) for n in range(6, 11)])
methods = ["luw", "lfm", "lwm"]
c_wave = 1

#  Analysis
for methoduse in methods:

    errs = [0.0]*len(dtrange)
    for i, dt in enumerate(dtrange):

        #  Solve
        assert int(1/dt) == 1/dt
        Grid = Domain(0, 1, dt, T, dx=dt/.5)
        Conds = Conditions(f, "pbc")
        Sol = advection_eqn_solve(Grid, Conds, c=c_wave, method=methoduse)

        #  Find errors
        err = 0
        for t_n in range(Grid.Nt):
            uex = f(Grid.xvals() - c_wave*t_n*Grid.dt)  # Exact solution
            err = max(err, (abs(uex - Sol.solndata[:, t_n])).max())

        errs[i] = err

    #  Plot
    errs = np.array(errs)
    plt.loglog(dtrange, errs, label=methoduse)
    plt.scatter(dtrange, errs)

plt.loglog(dtrange, (dtrange)**1, ':', label="O(dt^1)")
plt.loglog(dtrange, (dtrange)**2, ':', label="O(dt^2)")
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel("dt", fontsize=16)
plt.ylabel("error", fontsize=16)
plt.legend(fontsize=11)
plt.grid()
plt.show()
