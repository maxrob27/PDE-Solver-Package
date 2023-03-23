"""
MMS system test on Burgers equuation
====================================

The method of manufactured solutions (MMS) are applied to investigate the
applied schemes.
"""
from pde_solver.solver import Domain, Conditions, burgers_eqn_solve
import numpy as np


def f(x):
    """Function (for IC). Input must be array."""
    return np.sin(np.pi*x)+np.cos(2*np.pi*(x-0.25))


def exactsoln(x, t):
    """User defined exact solution."""
    ex = np.sin(np.pi*x)*np.cos(a*t) + np.cos(2*np.pi*(x-0.25))*np.exp(-t)
    return ex


def mms_source(x, t):
    """Method of manufactured solution source term."""
    return (-np.exp(-t)*np.cos(2*np.pi*(x-0.25))
            + (np.exp(-t)*np.cos(2*np.pi*(x-0.25))
            + np.cos(a*t)*np.sin(np.pi*x))*(np.pi*np.cos(a*t)*np.cos(np.pi*x)
            - 2*np.exp(-t)*np.pi*np.sin(2*np.pi*(x-0.25)))
            - a*np.sin(a*t)*np.sin(np.pi*x))


a = 1

Grid = Domain(0, 1, 0.0001, 1, dx=0.01)  # Input domain specifications
Conds = Conditions(f, "pbc")  # Initial and boundary condition functions
sol = burgers_eqn_solve(Grid, Conds, method="lwm", mms=mms_source)  # Solve
sol.plot(5, exactsoln)  # plots
