# Advection Solver 

This Python file contains a class called Advection which simulates advection of a perturbation with user-defined properties. The solver includes the following features:

The class takes as input the shape of the perturbation ("top_hat" or "gaussian"), the number of cells, the perturbation velocity and the number of periods.

The class calculates the limits, number of cells, and step size for the domain, and defines the initial conditions for the simulation.
The simulation is carried out using one of six different numerical methods: Godunov, centered difference, minmod limiter, MC limiter, Superbee limiter1, and Superbee limiter2.


## How to use the Advection Solver
Import Advection from the Python file:

    import cfdadvection.advection as advec # Importing module


Create an instance of the class with the following arguments:

    shape: string, the shape of the perturbation, either "top_hat" or "gaussian".
    n_cells: integer, the number of cells.
    p_veloc: float, the perturbation velocity.
    num_periods: integer, the number of periods to simulate.

    advection_solver = advec.Advection(shape="top_hat", n_cells=100, p_veloc=1.0, num_periods=2)

Call the advect method of the Advection class to solve the advection equation using one of the six numerical methods:

    advection_solver.advect(method="godunov", C=0.8)

method can be one of "godunov", "unlimited", "minmod", "mc", "superbee1", and "superbee2". C is the Courant number, a numerical parameter that determines the stability of the simulation.

The simulation results are stored in the Advection object, and can be plotted using the plot method that gives us the arrays with the information:

    advection_solver.plot(a)


## Example Usage

    from advection import Advection

    advection_solver = Advection(shape="top_hat", n_cells=100, p_veloc=1.0, num_periods=2)
    a = advection_solver.advect(method="superbee2", C=0.8)
    advection_solver.plot(a)

This will create a simulation of advection with a top hat profile, using the Superbee limiter 2 numerical method. The simulation will run for 2 periods, with a Courant number of 0.8. The results of the simulation can be plotted using matplotlib by the user.


## Analysis 

Based on the obtained results from the top hat profile and the Gaussian profile, it can be inferred that the use of the Godunov approach and the Superbee1 slope limiter tends to produce less accurate results. On the other hand, when considering the top hat profile, the unlimited method i.e centered different method, was observed to have a lower L2 Norm error. However, upon analyzing the advected profile, it is evident that the unlimited method resulted in an undershoot and overshoot of the profile. Finally, the MC slope limiter yielded the best results on both profiles. The MC slope limiter was able to replicate the initial profile and exhibit a lower L2 norm error.


