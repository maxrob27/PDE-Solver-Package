"""
Example linear advection solution (with flux limiter)
=====================================================

The flux limiter may be removed in the function advection_eqn_solve().
Oscillations will then be visible inthe solution.
"""
from pde_solver.solver import Domain, Conditions, advection_eqn_solve
import numpy as np


def f(x):
    """Function (for IC). Input must be array."""
    IC = 0 * x
    IC[int(len(IC)*0.35):int(len(IC)*0.55)] = 1
    return 0.1 * IC


def exactfunc(x, t):
    """Exact solution is shifted initial conditions (f(x))."""
    n = int(t*Grid.Nx % Grid.Nx)
    return np.hstack((f(x)[-n:], f(x)[:-n]))


Grid = Domain(0, 1, 0.001, 1, dx=0.002)  # Input domain specifications
Conds = Conditions(f, "pbc")  # Initial and boundary condition functions
Sol = advection_eqn_solve(Grid, Conds, method="lwm", fluxlim="vanleer", c=1)
Sol.plot(4, exactfunc)  # plots
