"""
MMS system test on the linear advection equuation
=================================================

The method of manufactured solutions (MMS) are applied to investigate the
applied schemes.
"""
from pde_solver.solver import Domain, Conditions, advection_eqn_solve
import numpy as np


def f(x):
    """Function (for IC). Input must be array."""
    return np.sin(2*np.pi*(x-0))*np.cos(a*0)


def exactsoln(x, t):
    """User defined exact solution."""
    ex = np.sin(2*np.pi*(x-t))*np.cos(a*t)
    return ex


def mms_source(x, t):
    """Method of manufactured solution source term."""
    return -a*np.sin(2*np.pi*(x-t))*np.sin(a*t)


a = 10
c_wave = 1

dx_ = 0.001
Grid = Domain(0, 1, 0.001, .2, dx=dx_)  # Input domain specifications
Conds = Conditions(f, "pbc")  # Initial and boundary condition functions
sol = advection_eqn_solve(Grid, Conds, c=c_wave, method="lwm", fluxlim="",
                          mms=mms_source)  # Solve with specified method
sol.plot(5, exactsoln)  # plots
