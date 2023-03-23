"""
Defines Domain class
====================

This class is called when the domain/grid/mesh PDE problem is being set up and
initialised.
"""
import numpy as np


class Domain:
    """The Domain for the PDE equation.

    Attributes
    ----------
    x0 : float
        The start of the domain.
    x1 : float
        The end of the domain (x1 > x0).
    dt : float
        The spacing between time steps.
    T_orig : float
        The final user specified time.
    T : float
        The final time in the domain (equal to or just smaller than T_orig)
    Nt : float
        The total number of time steps.
    dx : float
        The spatial spacing in the domain.
    Nx : float
        The number of spatial nodes in the domain.

    """
    def __init__(self, x0_, x1_, dt_, T_, **kwargs):
        self.x0 = x0_
        self.x1 = x1_
        self.dt = dt_
        self.T_orig = T_

        self.checks()  # Perform input checks

        self.T = (self.tvals())[-1]  # actual T value
        self.Nt = len(self.tvals())

        self.kwargs = kwargs
        self.dx, self.Nx = self.spatial()

    def checks(self):
        """Performs user input checks. Raises errors for infeasible inputs."""
        if self.dt <= 0:
            raise ValueError("User must enter dt > 0.")
        if self.x0 > self.x1:
            raise ValueError("User must enter boundaries where x1 > x0.")
        if self.T_orig < 0 + self.dt:
            raise ValueError("User must specity a final time (T) >= dt.")

    def spatial(self):
        """Returns `dx` and `Nx` under given parameters.

        If dx is given the domain must be discretised such that the spacing is
        exactly `dx`. Conversely if `Nx` is given, there should be exactly `Nx`
        nodes. Both or neither may not be specified. This function returns the
        user-unspecified parameter.

        Parameters
        ----------


        Returns
        -------
        dx : float
            The spatial spacing in the domain.
        Nx : float
            The number of spatial nodes in the domain.

        Raises
        ------
        RuntimeError
            When neither or both values (`dx` and `Nx`) are specified.

        """
        try:
            dxval = self.kwargs["dx"]
        except Exception:
            dxval = None
        try:
            Nxval = self.kwargs["Nx"]
        except Exception:
            Nxval = None

        if not dxval and Nxval:  # dx is not specified, Nx is
            return (self.x1-self.x0)/(Nxval - 1), Nxval
        elif dxval and not Nxval:  # dx is specified, Nx is not
            return dxval, int(np.floor((self.x1-self.x0)/dxval) + 1)

        else:
            if dxval is None or Nxval is None:  # Inconsistent values
                raise ValueError("User must specify either dx or Nx.")
            elif (self.x1-self.x0)/(Nxval - 1) != dxval:
                raise ValueError("User must specify consistent input values.")
            else:  # Consistent values
                return dxval, Nxval

    def xvals(self):
        """Returns array of x values with Nx points from x0 to x1."""
        return np.arange(self.x0, self.x1+self.dx/2, self.dx)

    def tvals(self):
        """Returns array of t values with spacing dt from 0 to T."""
        return np.arange(0, self.T_orig + self.dt/2, self.dt)

    def datamat(self):
        """Returns empty data matrix of dimension: Nx by len(tvals)."""
        return np.zeros((self.Nx, self.Nt))

    def __str__(self):
        """Presents overview of Domain."""
        pts = f"({self.Nx} points)"
        xstr = f"""x-axis: x={self.x0} to x={self.x1} with dx={self.dx} """+pts
        tstr = f"""t-axis: t=0 to t={self.T} with dt={self.dt}."""
        return xstr+"\n"+tstr
