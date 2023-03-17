import matplotlib.pyplot as plt
import numpy as np


class Advection:

    def __init__(self,shape,n_cells,p_veloc,num_periods):
        # limits:
        if shape == "top_hat":
            self.xmin = 0.
            self.xmax = 1.
        elif shape == "gaussian":
            self.xmin = -2
            self.xmax = 2

        self.shape = shape
        # Number of cells
        self.nx = n_cells

        # Number of ghost cells on each side
        self.ng = 2

        # Lowest and highest indices
        self.ilo = self.ng
        self.ihi = self.ng + self.nx - 1

        # Step size
        self.dx = (self.xmax - self.xmin)/(self.nx)

        # X axis vector -> we need the cell centres:
        self.x = self.xmin + (np.arange(self.nx + 2*self.ng) - self.ng + 0.5)*self.dx 
        
        # Define perturbation velocity
        self.u = p_veloc

        # Define number of periods
        self.num_periods = num_periods

        # Calculate period 
        self.period = (self.xmax - self.xmin)/self.u

        # Maximmum simulation time
        self.tmax = self.num_periods*self.period
        
    def forms(self):
        if self.shape == "top_hat":

            # Empty vector for the solution in double precision:
            a = np.zeros((self.nx + 2*self.ng), dtype = np.float64)
            
            # Add initial conditions;
            a[:] = 0.
            a[np.logical_and(self.x >= 1./3., self.x <= 2./3.)] = 1.
            return a
        
        elif self.shape == "gaussian":
            # Empty vector for the solution in double precision:
             
            a = np.exp(-0.5*(self.x/-.4)**2)
            # Add initial conditions:
            return a


    def fill_bcs(self, a):
    # Fill up the boundary conditions BCs
        for n in range(self.ng):
            # The left BC is
            a[self.ilo - 1 - n] = a[self.ihi - n]
            # The right BC is
            a[self.ihi + 1 + n] = a[self.ilo + n]
        return a

    def time_step(self,C):
        dt = C*self.dx/self.u
        return dt

 

    def states(self,a, dx, dt, method):
        def minmod(a, b):
            if abs(a) <= abs(b) and a*b > 0.0:
                return a
            if abs(b) <= abs(a) and a*b > 0.0:
                return b
            else:
                return 0.0
        def maxmod(a, b):
            if abs(a) > abs(b) and a*b > 0.0:
                return a
            elif abs(b) > abs(a) and a*b > 0.0:
                return b
            else: 
                return 0.0

        # Empty vector for the slope
        slope = np.zeros((self.nx + 2*self.ng), dtype = np.float64)
        r = np.zeros((self.nx + 2*self.ng), dtype = np.float64)
        
        if method == "godunov":
            # 1st approach: Godunov approach
            slope[:] = 0. 
        elif method == "unlimited":
            # 2nd approach: Centred difference method
            for i in range(self.ilo - 1, self.ihi + 2):
                slope[i] = 0.5*(a[i + 1] - a[i-1])/self.dx
        elif method == "minmod":  
            # 3rd approach: minmod limiter
            for i in range(self.ilo - 1, self.ihi + 2):
                slope[i] = minmod((a[i] - a[i-1])/dx, (a[i+1] - a[i])/dx) 
        elif method == "mc":
            # 4th approach: MC limiter
            for i in range(self.ilo - 1, self.ihi + 2):
                slope[i] = minmod(minmod(2*(a[i] - a[i-1])/self.dx, 2*(a[i+1]\
                       - a[i])/self.dx), 0.5*(a[i+1] - a[i-1])/self.dx) 
        elif method == "superbee1":
            #5th approach: Superbee limiter
            for i in range(self.ilo-1, self.ihi +2):
                r[i] = (a[i]-a[i-1])/(a[i+1]-a[i])
                slope[i] = maxmod(minmod(2*r[i],1),minmod(r[i],2)) 
        elif method == "superbee2":
            for i in range(self.ilo-1, self.ihi +2):
                A = minmod( (a[i+1] - a[i])/dx, 2.0*(a[i] - a[i-1])/self.dx )
                B = minmod( (a[i] - a[i-1])/dx, 2.0*(a[i+1] - a[i])/self.dx )
                slope[i] = maxmod(A, B)     
        elif method == "van_leer":
            #6th approach: Van Leer limiter
            for i in range(self.ilo-1, self.ihi +2):
                r = (a[i]-a[i-1])/(a[i+1]-a[i])
                slope[i] = (r+np.abs(r))/(1+np.abs(r))
 
        # Empty vector dor the L and R states
        al = np.zeros((self.nx + 2*self.ng), dtype = np.float64)
        ar = np.zeros((self.nx + 2*self.ng), dtype = np.float64)
    
        # Add derivative calculation
        for i in range(self.ilo, self.ihi + 2):
        
            # Compute L state
            al[i] = a[i - 1] + 0.5*self.dx*(1. - self.u*dt/self.dx)*slope[i-1]
        
            # Compute R state
            ar[i] = a[i] - 0.5*dx*(1. + self.u*dt/self.dx)*slope[i]
        return al, ar

    def riemann(self, al, ar):
        # Computing the flux
        if self.u > 0.:
            return self.u*al
        else: 
            return self.u*ar

    def conservative_update(self,a, dt, flux):
    
        # Empty vector for updated solution
        anew = np.zeros((self.nx + 2*self.ng), dtype = np.float64)
    
        # Update
        anew[self.ilo:self.ihi+1] = a[self.ilo:self.ihi+1] -\
        dt/self.dx*(flux[self.ilo+1:self.ihi+2] - flux[self.ilo:self.ihi+1])
    
        # Return
        return anew

    def plot(self, a):
        
        x1 = self.x[self.ilo:self.ihi+1]
        y1 = a[self.ilo:self.ihi+1]
        
        return x1, y1
    
    def advect(self, method, CFL): 
        a = self.forms()
        t = 0
        while t < self.tmax:
            a1 = self.fill_bcs(a)
            dt = self.time_step(CFL)
            if t + dt > self.tmax:
                dt = self.tmax - t
            al , ar = self.states(a1, self.dx, dt, method)
            flux = self.riemann(al,ar)
            anew1 = self.conservative_update(a , dt, flux)
            a1[:] = anew1[:]
            t+= dt
	    
        return a1



