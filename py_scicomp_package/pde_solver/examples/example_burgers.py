"""
Example Burgers' equation
=========================

The numerical solution may be compared to the exact solution until roughly
t = 1.581.
"""
from pde_solver.solver import Domain, Conditions, burgers_eqn_solve
import numpy as np
import scipy.optimize as spo


def f(x):
    """Function (for IC). Input must be array."""
    return 0.1*np.sin(2*np.pi*(x+.5))


def exactsoln(x, t):
    """User defined exact solution."""
    def _f(u):
        return 0.1*np.sin(2*np.pi*(x - u*t + 0.5)) - u

    def _fprime(u):
        return (-np.pi*t/5)*np.cos(2*np.pi*(x - u*t + 0.5)) - 1
    return spo.newton(_f, f(x), _fprime, tol=1e-13)


Grid = Domain(0, 1, 0.0001, 1.581, dx=0.0001)  # Input domain specifications
conds = Conditions(f, "pbc")  # Initial and boundary condition functions
sol = burgers_eqn_solve(Grid, conds, method="lwm")  # Solve with method
sol.plot(5, exactsoln)
