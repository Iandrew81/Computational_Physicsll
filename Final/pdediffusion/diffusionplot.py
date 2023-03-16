import glob
import numpy as np
import matplotlib.pyplot as plt


x_values = []
y_values = []

for filename in sorted(glob.glob('./outputfolder/logfile*.log')):
    with open(filename) as f:
        lines = f.readlines()
        x_value = float(lines[1].split()[0])
        y_value = float(lines[1].split()[1])
        x_values.append(x_value)
        y_values.append(y_value)

#Amdahl's law

#First we use the value for just 1 processor to get the time the calculation taakes in serial

s1 = y_values[0] #time in serial

x_amd = np.arange(1, 9,1)

y_amd = lambda x: s1/x_amd


# Plot the data and save the figure
plt.figure(figsize=(10, 8))
ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.plot(x_values, y_values, '-o', label = "Parallel Run")
plt.plot(x_amd, y_amd(x_amd), label = "Amdahls Law ")
plt.xlabel('Number of cores',size = 22)
plt.ylabel('Computing Time', size = 22)
plt.title('MPI Timings',size = 28)
plt.legend( fontsize = 20)

plt.savefig('./outputfolder/Mpi_plot.png')




