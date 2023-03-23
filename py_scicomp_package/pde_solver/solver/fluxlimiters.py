"""
Defines Fluxlimiter class
=========================

This class is called when a flux limiter is selected.
"""
import numpy as np


class Fluxlimiter:
    """The available flux limiters for the solver.

    Attributes
    ----------
    fluxlimiter : str
        The flux limiter name.
    beta : float
        The value of beta (if required) for the flux limiter.

    """
    def __init__(self, fluxlimiter_):
        self.fluxlimiter, self.beta = self.extractbeta(fluxlimiter_)

    def extractbeta(self, fluxlimtype):
        """Extracts the value of beta from the user input.

        For the sweby and chakosher flux limiter options, a value for beta is
        required. This is input after "_" the flux limiter name. For the other
        limiters this is ignored.

        Parameters
        ----------
        fluxlimtype : str
            The flux limiter including beta parameter.

        Returns
        -------
        fluxlimiter : str
            The flux limiter name.
        beta : float
            The value of beta (if required) for the flux limiter.

        Raises
        ------
        ValueError
            When the user flux limiter input us invalid.

        """
        if "_" in fluxlimtype:
            _loc = 0
            for loc, i in enumerate(fluxlimtype):
                if ord(i) == 95:
                    _loc = loc
            try:
                beta = float(fluxlimtype[_loc+1:])
            except Exception:
                raise ValueError("User flux limiter input invalid.")
            fluxlimtype = fluxlimtype[:_loc]
            return fluxlimtype, beta
        else:
            return fluxlimtype, None

    def __call__(self, inpvec):
        """Implements a flux limiter in the solution.

        This is represented by phi in the report. For the Sweby and
        Chakravarthy-Osher flux limiters extra tests are made to insure all
        inputs are valid (within appropriate range for beta).


        Parameters
        ----------
        inpvec : array
            The solution at the previous timestep.

        Returns
        -------
        array
            Array for phi.

        """
        vlen = len(inpvec)
        phi = np.zeros((vlen))
        for x in range(0, vlen):  # Careful treatment of 0s in fraction
            rnum = (inpvec[x % vlen]-inpvec[(x-1) % vlen])
            if rnum == 0:
                phi[x] = 0
            else:
                rden = (inpvec[(x+1) % vlen]-inpvec[x % vlen])
                if rden == 0:
                    phi[x] = 1
                else:
                    phi[x] = self.fluxfunc(rnum/rden)
        return phi

    def fluxfunc(self, r):
        """Implements a flux limiter in the solution.

        This is represented by phi in the report. For the Sweby and
        Chakravarthy-Osher flux limiters extra tests are made to insure all
        inputs are valid (within appropriate range for beta).


        Parameters
        ----------
        r : float
            Input for fluxlimiter.

        Returns
        -------
        float
            Evaluation of fluxlimiter(r)

        Raises
        ------
        NotImplementedError
            If the user selected scheme is not (yet) implemented.

        """
        #  Check flux limiter exists
        if (hasattr(Fluxlimiter, self.fluxlimiter)
           and callable(eval("Fluxlimiter."+self.fluxlimiter))):
            return eval("self." + self.fluxlimiter + "(r)")
        else:
            raise NotImplementedError(f"Flux limiter: '{self.fluxlimiter}' not"
                                      " implemented.")

    def vanleer(self, r):
        """Van Leer flux limiter."""
        return (r + np.abs(r))/(1 + np.abs(r))

    def superbee(self, r):
        """Superbee flux limiter."""
        return max(0, min(1, 2*r), min(2, r))

    def ospre(self, r):
        """Ospre flux limiter."""
        return 3*r*(r+1)/(2*(r**2+r+1))

    def centeredlimiter(self, r):
        """Monotised central flux limiter."""
        return max(0, min(2*r, 0.5*(2 + r), 2))

    def minmod(self, r):
        """Minmod flux limiter."""
        return min(1, r)

    def sweby(self, r):
        """Sweby flux limiter (0 <= beta <= 1)."""
        beta = self.beta
        if not beta:
            beta = 0.5  # Default value used if beta not specified
        if not 0 <= beta <= 1:
            raise ValueError("beta must be between 0 and 1 for the Sweby"
                             f"limiter (input: {beta}).")
        return max(0, min(beta*r, 1), min(r, beta))

    def chakosher(self, r):
        """Chakravarthy-Osher flux limiter (1 <= beta <= 2)."""
        beta = self.beta
        if not beta:
            beta = 1.5  # Default value used if beta not specified
        if not 1 <= beta <= 2:
            raise ValueError("Beta must be between 1 and 2 for the "
                             f"Chakravarty-Osher limiter (input: {beta}).")
        return max(0, min(r, beta))

    def __str__(self):
        return f"The {self.fluxlimiter} flux limiter is used."
