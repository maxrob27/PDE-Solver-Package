"""
Unit tests are performed mainly on the solution finding aspect of the solver
============================================================================

Pytest is used to execute this script.
"""
import pytest
from pde_solver.solver import (Domain, Conditions, Fluxlimiter, Scheme,
                               requiredinputs, extendvec, advection_eqn_solve,
                               initialise)
import numpy as np
import random as rn


def test_solver_correctrequiredinputs():
    """All correct inputs given (order irrelevant)."""
    required = ["a", "b", "c"]
    inputs = ["c", "a", "b", "d"]
    requiredinputs(required, inputs)


def test_solver_incorrectrequiredinputs():
    """Insufficient inputs given (order irrelevant)."""
    required = ["a", "b", "c"]
    inputs = ["c", "b", "f"]
    with pytest.raises(NameError):
        requiredinputs(required, inputs)


def test_solver_correctbeta():
    """Correct beta value extracted from flux limiter."""
    beta = round(2*rn.random(), 2)
    txt = "test"
    inp = txt+"_"+str(beta)
    fl = Fluxlimiter(inp)
    assert fl.extractbeta(inp) == (txt, beta)


def test_solver_incorrectbeta():
    """Incorrect flux limiter input."""
    beta = round(2*rn.random(), 2)
    txt = "test"
    inp = txt+"_"+str(beta)+txt
    with pytest.raises(ValueError):
        Fluxlimiter(inp)


def test_solver_betaerror1():
    """No beta value extracted from flux limiter (when none input)."""
    txt = "test"
    inp = txt+"_"
    with pytest.raises(ValueError):
        Fluxlimiter(inp)


def test_solver_betaerror2():
    """No beta value extracted from flux limiter (false input)."""
    txt = "te_st"
    inp = txt
    with pytest.raises(ValueError):
        Fluxlimiter(inp)


def test_solver_extendcorrect():
    """The vector is extended: first -> last, last -> first, length += 2."""
    first = 1
    last = 5
    u = list(range(first, last+1))
    uext = extendvec(u)
    print(u, uext)
    assert (uext[0], uext[-1], len(uext)) == (last, first, len(u) + 2)


def f(x):
    """Initial condition."""
    return np.sin(np.pi*x) + 1


#  Set up problem. All functions and classes in set up separately tested.
Grid = Domain(0, 1, 0.001, 4, dx=0.01)  # Input domain specifications
Conds = Conditions(f, "pbc")  # Initial and boundary condition functions
solndata = initialise(Grid, Conds)[0]


def test_solver_tridiagonal():
    """Tests all schemes create tridiagonal (or bidiagonal) matrices."""
    elementsum = 0
    for ftype in ["luwpbc", "lfmpbc", "lwmpbc"]:  # Tests all schemes (c > 0)
        #  Setup for scheme
        SelectScheme = Scheme(ftype, None, Grid, solndata, c=1, fluxlim="")
        schemefunc = SelectScheme.schemefunc()
        #  Reconstruct "A" matrix
        Nx = Grid.Nx
        solndata[:, 1:Nx+1] = np.identity(Nx)
        Arecon = np.zeros((Nx, Nx))
        #  Apply function
        for i in range(2, Nx+2):
            Arecon[:, i-2] = schemefunc(0, i)
        for j in range(1, Nx-1):
            #  Sum over all non tridiagonal entries, care taken for pbc.
            elementsum += sum(Arecon[j, (j+2) % Nx:(Nx-1+j) % Nx])
    assert elementsum == 0  # Expect 0 for all entries outside diagonal


def test_solver_falseflux():
    """Error when false fluxlimiter input."""
    falsename = "MMSCisaneasycourse"
    with pytest.raises(NotImplementedError):
        advection_eqn_solve(Grid, Conds, c=1, method="lwm", fluxlim=falsename)


#  pytest.main(["unittests_solver.py"])  # Run tests (optional run script)
