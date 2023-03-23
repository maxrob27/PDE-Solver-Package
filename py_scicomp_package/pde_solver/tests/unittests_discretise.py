"""
Unit tests are performed mainly on the discretisation aspect of the solver
==========================================================================

Pytest is used to execute this script.
"""
import pytest
from pde_solver.solver import Domain


def test_discretise_missinginputs():
    """`dx` and `Nx` are not specified, no discretisation may be made."""
    with pytest.raises(ValueError):
        Domain(0, 1, 0.01, 1)


def test_discretise_inconsistentinputs():
    """`dx` and `Nx` are both specified inconsistently, no discretisation."""
    with pytest.raises(ValueError):
        Domain(0, 1, 0.01, 1, dx=0.01, Nx=50)


def test_discretise_consistentinputs():
    """`dx` and `Nx` are both specified consistently, a discretisation."""
    dx_test = 0.01
    Nx_test = 101
    Grid = Domain(0, 1, 0.01, 1, dx=dx_test, Nx=Nx_test)
    assert [Grid.dx, Grid.Nx] == [dx_test, Nx_test]


def test_discretise_negativeT():
    """T is negative, no discretisation may be made."""
    with pytest.raises(ValueError):
        Domain(0, 1, 0.01, -1, dx=0.01)


def test_discretise_negativedt():
    """dt is negative, no discretisation may be made."""
    with pytest.raises(ValueError):
        Domain(0, 1, -0.01, 1, dx=0.01)


def test_discretise_inconsistentboundaries():
    """x1 < x0, no discretisation may be made."""
    with pytest.raises(ValueError):
        Domain(1, 0, 0.01, 1, dx=0.01)


def test_discretise_xvals():
    """The length of the array containing the x-values should be = Nx."""
    Grid = Domain(0, 1, 0.01, 1, dx=0.01)
    assert len(Grid.xvals()) == Grid.Nx


def test_discretise_tvals():
    """The length of the array containing the x-values should be = Nx."""
    Grid = Domain(0, 1, 0.01, 1, dx=0.01)
    assert len(Grid.tvals()) == Grid.Nt


#  pytest.main(["unittests_discretise.py"])  # Run tests (optional run script)
