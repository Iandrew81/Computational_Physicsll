#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from PIL import Image

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-im", help="Yes or no initial map", action="store_true")
parser.add_argument("-g", help="Yes or no GIF", action="store_true")
parser.add_argument("-RK", help="-RK takes the Runge–Kutta order (2nd or 3rd)", type=int, default=2)
parser.add_argument("-T", help="-T takes a number as the period", type=int, default=5)
parser.add_argument("-e", help="-c takes a as the exentricity", type=float, default=0.01671)

args = parser.parse_args()

class Earth_orbital_motion():
    
    def __init__(self, RK, T, e):
        
        self.T = T
        self.e = e
        self.RK_order = RK

        self.h_1 = 0.01
        
        self.t_axis = np.arange(0, self.T + self.h_1, self.h_1)
        
        # Constants
        self.M_sun = 1.    # mass of Sun in Solar Masses
        self.a = 1. # semi-major axis in AU
        self.G = 4*np.pi**2 # Gravitational constan
        


        # Initial conditions
        self.r_0 = [0, self.a*(1-self.e)]       # initial position of Earth in x and y
        self.v_0 = [-np.sqrt( (self.G*self.M_sun)/(self.a) * (1 + self.e)/(1-self.e) ), 0] # initial velocity of Earth in x and y
        
        self.directory_name = "Period_"+str(self.T) +"-ecc_"+ str(self.e)
    
    
    def initial_map(self):
    
        plt.figure(figsize = (8,6))

        plt.plot(self.r_0[0], self.r_0[1], 'bo', markersize=10, label='Earth')    # Earth
        plt.plot(0, 0, 'yo', c='y', markersize=20, label='Sun')          # Sun


        plt.tick_params(axis='both', which='major', labelsize=14, width=2, length=7)

        plt.xlabel('$x [AU]$', size=15)
        plt.ylabel('$y [AU]$', size=15)

        plt.xlim([-1.6, 1.6])
        plt.ylim([-1.6, 1.6])

        plt.legend(fontsize=16)
        plt.grid()


        plt.savefig("./outputfolder/"+str(self.directory_name)+"/initial_map{:03d}.png".format(0))

        plt.close()

    
    # The slope for the position is dr/dt = v
    def f_odeR(self, t, r_i, h_1):

        return self.h_1 * (-r_i*(self.G*self.M_sun)/(np.linalg.norm(r_i)**3))


    # The slope for the velocity is dv/dt = -(GM/r^3)r
    def f_odeV(self, r_i):

        return -r_i*(self.G*self.M_sun)/np.linalg.norm(r_i)**3

    
    def RK_2(self, t_axis, h_step, r_0, v_0):

        r_sln = np.zeros((len(self.t_axis),2)) # Array to save the positions

        v_sln = np.zeros((len(self.t_axis),2)) # Array to save the velocities

        r_sln[0] = self.r_0 #Initial position

        v_sln[0] = self.v_0 #Initial velocity

        for k in range(len(self.t_axis) - 1):

            #Slopes for RK2 of the positions
            slope_r1 = v_sln[k] #just v(t_i)

            slope_r2 = v_sln[k] + self.f_odeR(self.t_axis[k], r_sln[k], self.h_1)  # v(t_i+h) = v(t_i) + a(t_i)h

            #calculate r(ti+1)
            r_sln[k + 1] = r_sln[k] + self.h_1*(slope_r1 + slope_r2)/2 


            #Slopes for RK2 of the velocities
            slope_v1 = self.f_odeV(r_sln[k]) #a(t_i) = a(r_i)

            slope_v2 = self.f_odeV(r_sln[k + 1]) #a(t_i+h) = a(r_[i+1])

            #calculate r(ti+1)
            v_sln[k + 1] = v_sln[k] + self.h_1*(slope_v1 + slope_v2)/2

        return r_sln, v_sln

    
    def RK_3(self, t_axis, h_step, r_0, v_0):
    
        r_sln = np.zeros((len(self.t_axis),2)) # Array to save the positions

        v_sln = np.zeros((len(self.t_axis),2)) # Array to save the velocities

        r_sln[0] = self.r_0 #Initial position

        v_sln[0] = self.v_0 #Initial velocity

        for k in range(len(self.t_axis) - 1):

            #Slopes for RK2 of the positions
            slope_r1 = v_sln[k] #just v(t_i)

            slope_r2 = v_sln[k] + self.f_odeR(self.t_axis[k], r_sln[k], self.h_1/2)  # v(t_i+h) = v(t_i) + a(t_i)h/2

            slope_r3 = v_sln[k] + self.f_odeR(self.t_axis[k], r_sln[k], self.h_1)  # v(t_i+h) = v(t_i) + a(t_i)h

            #calculate r(ti+1)
            r_sln[k + 1] = r_sln[k] + self.h_1*(slope_r1 + 4*slope_r2 + slope_r3)/6 


            #Slopes for RK2 of the velocities
            slope_v1 = self.f_odeV(r_sln[k]) #a(t_i)

            slope_v2 = self.f_odeV(r_sln[k] + v_sln[k]*h_step/2) #a(t_i+h/2)= a(r_i + v_i*h/2)

            slope_v3 =  self.f_odeV(r_sln[k + 1]) #a(t_i+h)

            #calculate r(t_i+1)
            v_sln[k + 1] = v_sln[k] + self.h_step*(slope_v1 + 4*slope_v2 + slope_v3)/6

        return r_sln, v_sln

    
    def save_history(self):

        # Create two arrays of dimensions 2xlen(t_axis) to store the data for r and v
        if self.RK_order == 2:
        	R, V = self.RK_2(self.t_axis, self.h_1, self.r_0, self.v_0)

        elif self.RK_order == 3:
        	R, V = self.RK_3(self.t_axis, self.h_1, self.r_0, self.v_0)
        # Write the data to a txt file

        

        if not os.path.exists("./outputfolder/"+str(self.directory_name)):
            os.makedirs("./outputfolder/"+str(self.directory_name))

        with open("./outputfolder/"+str(self.directory_name)+"/history.txt", "w") as f:

            f.write("t (years)\t x (AU) \t y (AU) \t  vx (AU/y) \t vy (AU/y) \n")

            for i in range(len(self.t_axis)):
                f.write("{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\n".format(self.t_axis[i], R[i,0],  R[i,1],  V[i,0], V[i,1]))

                

    def plot_map(self, i, x_values, y_values, vx_values, vy_values):
    
        plt.figure(figsize = (8,6))

        plt.plot(x_values[i], y_values[i], 'bo', markersize=10, label='Earth')    # Earth
        plt.plot(0, 0, 'yo', c='y', markersize=20, label='Sun')          # Sun

        plt.plot(x_values[:i], y_values[:i], linestyle=':', c='r', linewidth=0.5, label='Orbit')          # Sun



        plt.tick_params(axis='both', which='major', labelsize=14, width=2, length=7)

        plt.xlabel('$x [AU]$', size=15)
        plt.ylabel('$y [AU]$', size=15)

        plt.xlim([-1.6, 1.6])
        plt.ylim([-1.6, 1.6])

        plt.legend(fontsize=16)
        plt.grid()
        
        
        if not os.path.exists("./outputfolder/"+str(self.directory_name)+"/orbits_images"):
            os.makedirs("./outputfolder/"+str(self.directory_name)+"/orbits_images")

        plt.savefig("./outputfolder/"+str(self.directory_name)+"/orbits_images/orbit{:03d}.png".format(i))

        plt.close()



class animation_class():
    
    def __init__(self, T, e):
        
        self.T = T
        self.e = e
        
        self.directory_name = "Period_"+ str(self.T) +"-ecc_"+ str(self.e)
        
        self.file_path = "./outputfolder/" + self.directory_name + '/history.txt'



    
    def get_data_from_file(self):
    
        data = np.loadtxt(self.file_path, skiprows=1)

        t_values = data[:, 0]
        x_values = data[:, 1]
        y_values = data[:, 2]
        vx_values = data[:, 3]
        vy_values = data[:, 4]

        return t_values, x_values, y_values, vx_values, vy_values
    
    
    def plot_map(self, i, x_values, y_values, vx_values, vy_values):
    
        plt.figure(figsize = (8,6))

        plt.plot(x_values[i], y_values[i], 'bo', markersize=10, label='Earth')    # Earth
        plt.plot(0, 0, 'yo', c='y', markersize=20, label='Sun')          # Sun

        plt.plot(x_values[:i], y_values[:i], linestyle=':', c='r', linewidth=0.5, label='Orbit')          # Sun



        plt.tick_params(axis='both', which='major', labelsize=14, width=2, length=7)

        plt.xlabel('$x [AU]$', size=15)
        plt.ylabel('$y [AU]$', size=15)

        plt.xlim([-1.6, 1.6])
        plt.ylim([-1.6, 1.6])

        plt.legend(fontsize=16)
        plt.grid()

        
        
        
        if not os.path.exists("./outputfolder/"+str(self.directory_name)+"/orbits_images"):
            os.makedirs("./outputfolder/"+str(self.directory_name)+"/orbits_images")

        plt.savefig("./outputfolder/"+str(self.directory_name)+"/orbits_images/orbit{:03d}.png".format(i))

        plt.close()

    
    
    def create_images(self):
    
        t_values, x_values, y_values, vx_values, vy_values = self.get_data_from_file()

        for i in range(len(t_values)):
            if i%5 == 0:
                self.plot_map(i, x_values, y_values, vx_values, vy_values)

    
    
    def create_gif(self):
    
        images_in = "./outputfolder/"+str(self.directory_name)+"/orbits_images/orbit****.png"

        gif_image_out = "./outputfolder/"+str(self.directory_name)+"/earth_orbtit.gif"

        imgs = (Image.open(f) for f in sorted(glob.glob(images_in)))

        img = next(imgs)

        img.save(fp = gif_image_out, format='GIF', append_images=imgs, save_all=True, duration=100, loop=0)




orbit = Earth_orbital_motion(args.RK,args.T, args.e)


if args.im:
    orbit.initial_map()

orbit.save_history()


orbit_gif = animation_class(args.T, args.e)

if args.g:

	orbit_gif.create_images()

	orbit_gif.create_gif()



    
