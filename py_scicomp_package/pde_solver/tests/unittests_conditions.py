"""
Unit tests are performed mainly on the PDE condition aspects of the solver
==========================================================================

Pytest is used to execute this script.
"""
import pytest
from pde_solver.solver import Domain, Conditions
import numpy as np


def f(x):
    return np.sin(x)


x0 = 0
x1 = 1
Grid = Domain(x0, x1, 0.01, 1, dx=0.01)


def test_conditions_initialchanged():
    """The initial data column is updated to the specified function."""
    Conds = Conditions(f, "pbc")
    A = np.random.rand(Grid.Nx, Grid.Nt)
    Anew = Conds.initial(Grid, A)
    assert sum((Anew[:, 0])) == sum(np.sin(Grid.xvals()))


def test_conditions_initialrestunchanged():
    """The rest of the solution data is unaffected (checks with random)."""
    Conds = Conditions(f, "pbc")
    A = np.random.rand(Grid.Nx, Grid.Nt)
    Anew = Conds.initial(Grid, A)
    assert sum(sum(A[:, 1:])) == sum(sum(Anew[:, 1:]))


def test_conditions_exceededtolerance():
    """Given function f(x) is not sufficiently consistent at boundaries."""
    tol = 0.5*(abs(f(x1) - f(x0)))  # Encapsulates smallest difference
    Conds = Conditions(f, "pbc", tol)
    A = Conds.initial(Grid, Grid.datamat())
    with pytest.raises(RuntimeError):
        Conds.boundary(Grid, A)


def test_conditions_withintolerance():
    """Given function f(x) is sufficiently consistent at boundaries."""
    fvals = f(Grid.xvals())
    tol = abs(max(fvals) - min(fvals))  # Encapsulates largest difference
    Conds = Conditions(f, "pbc", tol)
    A = Conds.initial(Grid, Grid.datamat())
    A = Conds.boundary(Grid, A)[0]
    assert abs(A[0, 0] - A[-1, 0]) <= tol


def test_conditions_assertpbc():
    """Given function f(x) is sufficiently consistent at boundaries."""
    Conds = Conditions(f, "pbc", np.inf)
    A = Conds.initial(Grid, Grid.datamat())
    assert Conds.boundary(Grid, A)[1] == "pbc"


def test_conditions_assertnopbc():
    """Given function f(x) is sufficiently consistent at boundaries."""
    with pytest.raises(NotImplementedError):
        Conds = Conditions(f, "", np.inf)
        Conds.boundary(Grid, Grid.datamat())


#  pytest.main(["unittests_conditions.py"])  # Run tests (optional run script)
