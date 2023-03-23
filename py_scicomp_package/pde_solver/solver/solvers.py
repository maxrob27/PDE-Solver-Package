"""
Main solver file
================

This file incldes the Soln and Scheme classes. It also includes solving
functions which generally call other classes or functions in order to solve the
PDE. Finally the specific solvers are included for the linear advection
equation and Burgers' equation.

"""
import numpy as np
import matplotlib.pyplot as plt
from PIL import ImageColor
from pde_solver.solver.fluxlimiters import Fluxlimiter


class Soln:
    """The solution data of the equation. Plotting is included.

    Attributes
    ----------
    Domain : Class
        The Domain Class for the PDE problem.
    solndata : Matrix
        The solution in a data matrix form.

    """
    def __init__(self, solndata_, Domain_):
        self.solndata = solndata_
        self.Domain = Domain_

    def plot(self, numplots=3, exactfunc=None):
        """Plots solution between `t=0` and `t=T`.

        The solution is plotted on a 2D graph with the solution at selected
        times presented as different curves. These selected times are evenly
        distributed between `[t,T]`. Nothing is returned but a plot is visually
        shown. The default number of curves is 3 if none is specified.

        Parameters
        ----------
        numplots : float
            The number of different (evenly spaced) times at which to plot.
        exactfunc : function, optional
            If known, this is the exact solution.

        Returns
        -------

        """
        print("plotting solution...")
        plt.figure(figsize=(8, 6), dpi=115)
        ax = plt.gca()
        lw = 2.3
        if exactfunc:
            sch = "scheme: "
        else:
            sch = ""

        for nplot in np.linspace(0, self.Domain.Nt-1, numplots):
            cx = next(ax._get_lines.prop_cycler)['color']
            plt.plot(self.Domain.xvals(), self.solndata[:, int(nplot)],
                     label=sch+f"t = {round(self.Domain.dt*nplot,3)}",
                     color=cx, linewidth=lw)

            if exactfunc:
                # Darken colour for exact solution
                rgbfact = 0.75  # Darken factor
                rgb = ImageColor.getcolor(cx, "RGB")
                rgbnew = tuple([(rgbfact*rgbcolor)/256 for rgbcolor in rgb])

                t_plot = round(self.Domain.dt*nplot, 3)
                plt.plot(self.Domain.xvals(), exactfunc(self.Domain.xvals(),
                         self.Domain.dt*nplot),
                         label=f"exact:     t = {t_plot}",
                         linestyle='--', color=rgbnew, linewidth=1.5*lw)

        plt.xticks(fontsize=13)
        plt.yticks(fontsize=13)
        plt.grid()
        plt.xlabel("x", fontsize=16)
        plt.ylabel("u", fontsize=16)
        plt.legend(fontsize=11)

        plt.show()

    def __call__(self, t):
        """Returns solution vector at time t."""
        t_eval = int(t/self.Domain.dt)
        if t == -1:
            t_eval = t
        elif t < 0:
            raise ValueError("t must be positive (or -1 for last entry)")
        return self.solndata[:, t_eval]


class Scheme:
    """The selected scheme class.

    Attributes
    ----------
    schemename : str
        The selected scheme name for the PDE problem.
    f : function
        The nonlinear function for the conservation law (not for linear).
    Domain : Class
        The Domain Class for the PDE problem.
    solndata : Matrix
        The solution in a data matrix form.
    **c : int
        The wave speed for the linear advection equation.

    """
    def __init__(self, schemename_, fnonlin_, Domain_, solndata_, **kwargs_):
        self.schemename = schemename_
        self.f = fnonlin_
        self.Domain = Domain_
        self.solndata = solndata_
        self.kwargs = kwargs_

    def schemefunc(self):
        """Creates the update function (func) for the finite difference scheme.

        math:: u^{n+1} = func(u^{n})

        Parameters
        ----------

        Returns
        -------
        function
            The update function which is the selected finite difference scheme.

        Raises
        ------
        NotImplementedError
            The specified method and flux limiter (if used) does not exist.

        """
        try:
            return eval("self." + self.schemename + "()")
        except Exception:
            raise NotImplementedError("Finite difference stencil "
                                      f"{self.schemename} not implemented or "
                                      " invalid conditions.")

    def lwmpbc(self):
        """Lax-Wendroff method (periodic)."""
        try:
            fluxlim = self.kwargs['fluxlim']
        except Exception:
            fluxlim = None
        if not self.f:  # Linear
            Cd = self.kwargs["c"]*self.Domain.dt/self.Domain.dx
            if np.abs(Cd) > 1:  # Stability check
                raise RuntimeError("Aborted - The scheme will not be stable. "
                                   "|c*dt/dx| must be smaller than 1, input "
                                   f"was |c*dt/dx| = {np.abs(Cd)}.")

            am1 = 0.5*Cd*(1+Cd)
            a0 = 1 - Cd**2
            ap1 = -0.5*Cd*(1-Cd)

            if not fluxlim:  # No flux limiter
                print("Using Lax Wendroff Method (no flux limiter).")
                A = (np.diagflat(a0*np.ones(self.Domain.Nx))
                     + np.diagflat(am1*np.ones(self.Domain.Nx-1), -1)
                     + np.diagflat(ap1*np.ones(self.Domain.Nx-1), 1))
                A[0, -2] = am1  # Periodic BCs
                A[-1, 1] = ap1  # Periodic BCs

                def func(t, i):
                    return A@self.solndata[:, i-1]

            elif fluxlim:  # Use of flux limiter
                print("Using Lax-Wendroff method with "
                      "" + fluxlim + " flux limiter.")

                def func(t, i):
                    phi = fluxlimiter(self.solndata[:, i-1],
                                      fluxlim)
                    phibar = np.insert(phi[:-1], 0, phi[-1])

                    A = ((1-Cd)*np.identity(self.Domain.Nx)
                         - ap1*np.diag(phi+phibar)
                         + ap1*np.diagflat(phi[:-1], 1)
                         + Cd*np.diagflat(np.ones(self.Domain.Nx-1), -1)
                         + ap1*np.diag(phi[:-1], -1))

                    A[0, -2] = Cd + ap1*phi[-1]
                    A[-1, 1] = ap1*phi[-1]
                    return A@self.solndata[:, i-1]

            else:
                raise NotImplementedError("Fluxlim: '"+fluxlim
                                          + "' is not implemented.")
            return func

        else:  # Nonlinear
            print("Using the (nonlinear) Lax-Wendroff method [Richtmeyer "
                  "two-step method].")

            def func(t, i):
                mu = (self.Domain.dt/self.Domain.dx)
                f = self.f
                uext = extendvec(self.solndata[:, i-1])  # Extended u vector
                up = (uext[1:-1]+uext[2:])/2-mu*(f(uext[2:])-f(uext[1:-1]))/2
                um = (uext[:-2]+uext[1:-1])/2-mu*(f(uext[1:-1])-f(uext[:-2]))/2
                return uext[1:-1] - mu*(f(up) - f(um))  # (12.26) LeVeque

            return func

    def lfmpbc(self):
        """Lax Friedrichs Method (periodic)."""
        if not self.f:  # Linear
            Cd = self.kwargs["c"]*self.Domain.dt/self.Domain.dx
            if np.abs(Cd) > 1:  # Stability check
                raise RuntimeError("Aborted - The scheme will not be stable."
                                   "|c*dt/dx| must be smaller than 1, input "
                                   f"was |c*dt/dx| = {np.abs(Cd)}.")

            print("Using the (linear) Lax-Friedrichs method")

            am1 = 0.5*(1+Cd)
            ap1 = 0.5*(1-Cd)
            A = (np.diagflat(am1*np.ones(self.Domain.Nx-1), -1)
                 + np.diagflat(ap1*np.ones(self.Domain.Nx-1), 1))
            A[0, -2] = am1  # Periodic BCs
            A[-1, 1] = ap1  # Periodic BCs

            def func(t, i):
                return A@self.solndata[:, i-1]
            return func
        else:  # Nonlinear
            print("Using the (nonlinear) Lax-Friedrichs method")

            def func(t, i):
                uext = extendvec(self.solndata[:, i-1])  # Extended u vector
                return (0.5*(uext[:-2]+uext[2:])
                        - (self.Domain.dt/self.Domain.dx)
                        * (self.f(uext[2:])-self.f(uext[:-2])))  # 12.15LeVeque
            return func

    def luwpbc(self):
        """Left upwind method (periodic)."""
        if not self.f:  # Linear
            print("Using the (linear) left upwind method")
            Cd = self.kwargs["c"]*self.Domain.dt/self.Domain.dx
            if Cd < -1:  # Stability check
                raise RuntimeError("Aborted - The scheme will not be stable. "
                                   "c*dt/dx must be greater than -1, input "
                                   f"was c*dt/dx = {Cd}.")
            if self.kwargs["c"] < 0:
                raise RuntimeError("The wavespeed should be positive for left "
                                   "upwind. Input was c = " + self.kwargs["c"])

            A = (np.diagflat(np.ones(self.Domain.Nx))
                 + np.diagflat(-1*np.ones(self.Domain.Nx-1), -1))
            A[0, -2] = -1  # Periodic BCs
            A = (np.identity(self.Domain.Nx)
                 + (self.Domain.dt/self.Domain.dx)*(-self.kwargs["c"])*A)

            def func(t, i):
                return A@self.solndata[:, i-1]
            return func
        else:  # Nonlinear
            def func(t, i):
                uext = extendvec(self.solndata[:, i-1])  # Extended u vector
                return (uext[1:-1] - (self.Domain.dt/self.Domain.dx)
                        * (self.f(uext[1:-1]) - self.f(uext[:-2])))
            return func

    def ruwpbc(self):
        """Right upwind method (periodic)."""
        if not self.f:  # Linear
            print("Using the (linear) right upwind method")
            Cd = self.kwargs["c"]*self.Domain.dt/self.Domain.dx
            if Cd < -1:  # Stability check
                raise RuntimeError("Aborted - The scheme will not be stable. "
                                   "c*dt/dx must be greater than -1, input was"
                                   f" c*dt/dx = {Cd}.")
            if self.kwargs["c"] > 0:
                raise RuntimeError("The wavespeed should be negative for "
                                   "right upwind. Input was c = "
                                   + str({self.kwargs["c"]}))

            A = (np.diagflat(-1*np.ones(self.Domain.Nx))
                 + np.diagflat(np.ones(self.Domain.Nx-1), 1))
            A[-1, 1] = 1  # Periodic BCs
            A = (np.identity(self.Domain.Nx) + (self.Domain.dt/self.Domain.dx)
                 * (-self.kwargs["c"]) * A)

            def func(t, i):
                return A@self.solndata[:, i-1]
            return func
        else:  # Nonlinear
            def func(t, i):
                uext = extendvec(self.solndata[:, i-1])  # Extended u vector
                return (uext[1:-1] - (self.Domain.dt/self.Domain.dx)
                        * (self.f(uext[2:]) - self.f(uext[1:-1])))
            return func


#  Further discretisation functions may be added here  #


def requiredinputs(required, inputs):
    """Performs a check that all required inputs have been provided."""
    for item in required:
        if item not in inputs:  # Required user inputs
            raise NameError(f"No user value set for {item}.")
    for item in inputs:
        if item not in required and item not in ["mms", "fluxlim"]:
            print(f"[Input: '{item}' not used to solve PDE.]")


def initialise(Domain, Conditions):
    """Initialises storage matrix and applies PDE conditions."""
    solndata = Domain.datamat()  # Initialise storage matrix
    solndata = Conditions.initial(Domain, solndata)  # Initial Conditions
    return Conditions.boundary(Domain, solndata)  # Boundary


def marcher(Domain, solndata, rhs, **kwargs):
    """Marches forward in time from `t=0` to `t=T`.

    Parameters
    ----------
    Domain : Class
        The Domain class for the PDE problem.
    solndata : Matrix
        The solution in a data matrix form.
    rhs : function
        The update function which is the selected finite difference scheme.
    **mms : function, optional
        This is the function for the method of manufactured solutions.

    Returns
    -------
    array
        Solution in matrix form with the solution at the next time step.

    """
    try:
        mms = kwargs['mms']
    except Exception:
        mms = None

    t = 0
    dt = Domain.dt
    for i in range(1, Domain.Nt):
        if i % (int(Domain.Nt/10)) == 0:
            print(int(100*i/Domain.Nt + 1), "%")  # Progress report

        t = dt*i
        solndata[:, i] = rhs(t, i)
        if mms:
            solndata[:, i] += dt*mms(Domain.xvals(), t)

    return solndata


def extendvec(inpvec):
    """Extends input vector by one either end peridodically."""
    return np.insert(inpvec, [0, len(inpvec)], [inpvec[-1], inpvec[0]])


def fluxlimiter(inpvec, fluxlimtype):
    """Implements a flux limiter in the solution.

    This is represented by phi in the report. For the Sweby and
    Chakravarthy-Osher flux limiters extra tests are made to insure all inputs
    are valid (within appropriate range for beta).


    Parameters
    ----------
    inpvec : array
        The solution at the previous timestep.
    fluxlimtype : str
        The selected flux limiter (see report/documents for options).

    Returns
    -------
    array
        Array for phi.

    Raises
    ------
    NotImplementedError
        If the user selected scheme is not yet implemented.

    """
    Fluxer = Fluxlimiter(fluxlimtype)
    phi = Fluxer(inpvec)

    return phi


def fdscheme(Domain, solndata, usepbc, **kwargs):
    """Creates the update function (func) for the finite difference scheme.

    math:: u^{n+1} = func(u^{n})

    Parameters
    ----------
    Domain : Class
        The Domain class for the PDE problem.
    solndata : Matrix
        The solution in a data matrix form.
    usepbc : str
        Specifies whether periodic boundary conditions are to be used.
    **method : str
        User specified method.
    **fnonlin : function
        Nonlinear function for conservation law.

    Returns
    -------
    function
        The update function which is the selected finite difference scheme.

    """

    try:
        fnonlin = kwargs['fnonlin']
    except Exception:
        fnonlin = None
    ftype = kwargs["method"]+usepbc

    SelectScheme = Scheme(ftype, fnonlin, Domain, solndata, **kwargs)
    return SelectScheme.schemefunc()


# ------------------------SPECIFIC-EQUATION-SOLVERS----------------------------


def advection_eqn_solve(Domain, Conditions, **kwargs):
    """Solves the linear advection equation.

    Parameters
    ----------
    Domain : Class
        The Domain class for the PDE problem.
    Conditions : Class
        The (initial and boundary) conditions class for the PDE problem.
    **c : int
        The wave speed for the linear advection equation.
    **method : str
        User specified method
    **fluxlimiter : str, optional
        User specified flux limiter

    Returns
    -------
    Class
        The solution class.

    """
    print('\n\n --- Solving Linear Advection Equation ---\n')
    requiredinputs(["c", "method"], kwargs)  # Check inputs provided

    solndata, pbcstr = initialise(Domain, Conditions)  # Storage and Conditions

    rhs_advection_eqn = fdscheme(Domain, solndata, pbcstr, **kwargs)  # Update
    solndata = marcher(Domain, solndata, rhs_advection_eqn, **kwargs)  # March

    return Soln(solndata, Domain)  # Return solution as class


def burgers_eqn_solve(Domain, Conditions, **kwargs):
    """Solves Burgers equation.

    Parameters
    ----------
    Domain : Class
        The Domain class for the PDE problem.
    Conditions : Class
        The (initial and boundary) conditions class for the PDE problem.
    **method : str
        User specified method
    **fluxlimiter : str, optional
        User specified flux limiter

    Returns
    -------
    Class
        The solution class.

    """
    print("""\n\n --- Solving Burgers' Equation ---\n""")
    requiredinputs(["method"], kwargs)  # Check inputs provided

    def funcnonlin(u):  # Burgers equation
        return 0.5*u**2

    solndata, pbcstr = initialise(Domain, Conditions)  # Storage and Conditions

    rhs_burgers_eqn = fdscheme(Domain, solndata, pbcstr, fnonlin=funcnonlin,
                               **kwargs)
    solndata = marcher(Domain, solndata, rhs_burgers_eqn, **kwargs)  # March

    return Soln(solndata, Domain)  # Return solution as class
