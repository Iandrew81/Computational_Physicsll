#!/usr/bin/env python
# coding: utf-8

# Import libraries

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sint
from mpi4py import MPI
import logging

# get basic information about the MPI communicator
world_comm = MPI.COMM_WORLD
world_size = world_comm.Get_size()
my_rank = world_comm.Get_rank()


# Configure logging
logging.basicConfig(filename='./outputfolder/logfile'+str(world_size)+'.log', filemode = 'w', level=logging.INFO,
                    format='%(message)s')

# Log MPI information
logging.info(f"Process rank {my_rank} out of {world_size}")


# Define constants
c_2 = 1 # difussivity constant

# Length of the box/domain
L = 200

# Discretisation of the domain
N = 2000

# Define the step size
h = L/N

# Define x-axis
x = np.arange(-L/2, +L/2, h)

# Time step
t_step = 0.025
t_max = 400.

# Time discretisation:
t = np.arange(0, t_max, t_step)

# Heat equation, FFT method
# $u_t = c^2 u_{xx}$

# Wavenumbers = spatial frequencies:
k_numbers = 2*np.pi*np.fft.fftfreq(len(x), d = h)

# Initial conditions
u_0 = np.zeros(len(x), dtype= complex)



# Replace zeroes with cos(4*pi*x/L), alpha should be 4*pi/L
u_0[int((L / 2 - L / 8)/h):int((L / 2 + L / 8)/h)]  = np.cos(4*np.pi*x[int((L / 2 - L / 8)/h):int((L / 2 + L / 8)/h)]/L)


# Fourier transform
u_0_fourier = np.fft.fft(u_0)
u_0_fourier_conc = np.concatenate((u_0_fourier.real, u_0_fourier.imag))

# Construct ODE (RHS of ODE)
# Function to get RHS

def RHS_ODE(u_0_fourier_conc, t, k_numbers, c_2):    
    u_tilde = u_0_fourier_conc[:N] + (1j)*u_0_fourier_conc[N:]
    rhs_u_tilde = -(c_2**2)*(k_numbers**2)*u_tilde
    rhs_ode = np.concatenate((rhs_u_tilde.real, rhs_u_tilde.imag))
    return rhs_ode


# k ODEs: solution
start_time1 = MPI.Wtime()
solution = sint.odeint(RHS_ODE, u_0_fourier_conc, t, args = (k_numbers, c_2))
end_time1 = MPI.Wtime()
if my_rank == 0:
    print("k ODEs: solutions: " + str(end_time1 - start_time1))


# Reconstruct Complex solution:
start_time2 = MPI.Wtime()
u_solution = solution[:, :N] + (1j)*solution[:, N:]
end_time2 = MPI.Wtime()
if my_rank == 0:
    print("Reconstruct ODE solutions: " + str(end_time2 - start_time2))

# Inverse Fourier transform of each u_solution
# For loop with k as index
inv_u_solution = np.zeros((u_solution.shape[0]//(world_size),u_solution.shape[1]), dtype = complex)

# Define the workloads
workloads = [ len(t) // world_size for i in range(world_size) ]
for i in range( len(t) % world_size ):
        workloads[i] += 1
my_start = 0
for i in range( my_rank ):
    my_start += workloads[i]
my_end = my_start + workloads[my_rank]


# We want to parallelise the code from here onwards
start_time3 = MPI.Wtime()

for k in range(my_start, my_end):
    inv_u_solution[k-my_start, :] = np.fft.ifft(u_solution[k, :])

end_time3 = MPI.Wtime()

logging.info(f"{world_size} {end_time3 - start_time3}")



if my_rank == 0:
    print("Inverse fourier transform: " + str(end_time3 - start_time3))

if my_rank == 0:
    Full_inv_u_solution = inv_u_solution
    for i in range(1,world_size):
        store = np.empty(inv_u_solution.shape)
        world_comm.Recv(store,source = i,tag = 77)
        Full_inv_u_solution = np.vstack((Full_inv_u_solution,store))
else:
    store = np.array(inv_u_solution.real)
    world_comm.Send(store,dest = 0,tag = 77)

# Plotting the solution:

start_time4 = MPI.Wtime()
# Add colour
R = np.linspace(1, 0, len(t))
B = np.linspace(0, 1, len(t))
G = 0	

if my_rank == 0:

    plt.figure(figsize= (10, 6))
    for j in range(len(t)):
	    plt.plot(x, Full_inv_u_solution[j, :].real, color = [R[j], G, B[j]])
    plt.xlabel("position [m]")
    plt.ylabel("temperature [$\degree$ C]")
    plt.savefig("parallel_output.png")
    plt.close()
end_time4 = MPI.Wtime()
if my_rank == 0:
    print("Plotting the solution: " + str(end_time4 - start_time4))
