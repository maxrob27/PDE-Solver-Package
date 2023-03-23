"""
Defines Conditions class
========================

This class is called when the conditions of the PDE problem are being set up
and initialised.
"""


class Conditions:
    """The conditions on the PDE problem.

    Attributes
    ----------
    initialfunc : function
        The (continuous) function for the initial condition.
    boundaryfunc : function or str
        The (continuous) function for the initial condition or string
        indicating periodic boundary conditions.
    tol : float, optional
        The tolerance used to perform the user input consistency checks. If
        none is entered, a default of `tol=1e-8` is used.

    """
    def __init__(self, initialfunc_, boundaryfunc_, tol_=1e-8):
        self.initialfunc = initialfunc_
        self.boundaryfunc = boundaryfunc_
        self.tol = tol_

    def __str__(self):
        """Presents overview of Conditions."""
        istr = f"Initial conditions: {self.initialfunc.__name__}(x)"
        if type(self.boundaryfunc) != str:
            bstr = f"Boundary conditions: {self.boundaryfunc.__name__}(x)"
        else:
            bstr = "periodic"
        return ("Initial conditions: "+istr+"\nBoundary conditions: " + bstr +
                f"\nTolerance = {self.tol}")

    def initial(self, Domain_, solndata_):
        """Applies initial condition to data storage matrix.

        Parameters
        ----------
        Domain : Class
            The Domain Class for the PDE problem.
        solndata : Matrix
            The solution in a data matrix form.

        Returns
        -------
        Matrix
            The solution in a data matrix form with boundary conditions applied

        """
        solndata_[:, 0] = self.initialfunc(Domain_.xvals())
        return solndata_

    def boundary(self, Domain_, solndata_):
        """Applies boundary condition to data storage matrix.

        Parameters
        ----------
        Domain : Class
            The Domain Class for the PDE problem.
        solndata : Matrix
            The solution in a data matrix form.

        Returns
        -------
        Matrix
            The solution in data matrix form with boundary conditions applied.
        str
            The use of periodic boundary conditions is specified.

        Raises
        ------
        ValueError
            If the user specified tolerance (Conditions class) exceeded. For
            periodic boundary conditions this is tested at the first and last
            node. For non-periodic boundary conditions the initial and boundary
            condition overlap is tested.

        """
        if self.boundaryfunc == "pbc":
            print("Periodic Boundary conditions registered.")
            toltest = abs(self.initialfunc(Domain_.xvals())[0]
                          - self.initialfunc(Domain_.xvals())[-1])
            if self.tol:
                if toltest > self.tol:
                    raise RuntimeError("The initial conditions are not "
                                       f"periodic. A tolerance of {self.tol} "
                                       f" was used (found: tol = {toltest}).")
            return solndata_, "pbc"
        else:
            raise NotImplementedError("Non-periodic boundary conditions have "
                                      "not been implemented. This will be "
                                      "further dissertation work.")
