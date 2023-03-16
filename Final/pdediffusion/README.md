Heat Equation Solver using MPI

This code consits of two parts:
1) diffusionMPI.py

This code solves the heat equation in one dimension using MPI (Message Passing Interface) to parallelize the code. The heat equation is solved using the FFT method, which allows us to calculate the spatial derivative of the temperature function. The heat equation is given by:

u_t = c^2 u_xx

where u is the temperature as a function of position x and time t, and c is the diffusion constant.

The code has been written in Python and requires the following packages to be installed:

    numpy
    matplotlib
    scipy
    mpi4py

Usage

To run the code, open a terminal in the directory where the code is saved and run:

mpiexec -n <number of processes> python diffusionMPI.py

where <number of processes> is the number of processes to be used in parallelization.

The output of the code is saved in the current and outputfolder directory, which is created if it does not exist. The output includes a logfile and a plot of the temperature as a function of temperature and time.

2)diffusionplot.py
 
This python script that helps users to plot the information collected from the logfiles i.e. number of processes vs computing time. These plot are found in the outputfolder directory
