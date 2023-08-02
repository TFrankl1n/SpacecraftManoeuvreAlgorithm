import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import math as nm

#Additional Features
import planetary_data as my_PD
import tools as my_T


def Default_Settings():                                                        #convert to features, and set central body COEs
    return{
        'centralbody':None,
        'COEs':[],
        'RV':[],
        'J2': False,
        'aero':False,
        'moon_gravity':False,
        'solar_gravity':False
        }                                                                       #add to this, '', 'frame:j200?', 'additional perturbations 

class Orbit_Propagator:
    def __init__(self, state0, tspan, dt, COEs=False, degrees = False, centralbody=my_PD.earth, perturbations = Default_Settings()):
        
        #initialise all relevant values, and replace any which we wish to change
        
        #initialise COESs and the state inputs as pos and state, as we will propagate in state vectors
        print("Running __init__ ...")
        if COEs:                                                               # if passing coes as initial state 
            self.r0, self.v0 = my_T.Convert_COEs_to_RV(state0, degrees=degrees, mu=centralbody['mu'])
        else:                                                                   #Pass straight in if not coes
            self.r0 = state0[:3]                                                #before 3
            self.v0 = state0[3:]                                                #after 3
            
        #self.y0 = np.vstack((self.r0 + self.v0))
        self.y0 = self.r0.tolist() + self.v0.tolist()
        self.tspan = tspan
        self.dt = dt
        self.centralbody = centralbody
        
        #total no. of steps
        self.n_steps = int(np.ceil(self.tspan/self.dt))                        #cieling, rounds any float to the next whole number (2.1 --> 3)

        #initialise variables
        self.t_vals = np.zeros((self.n_steps,1))
        self.y_vals = np.zeros((self.n_steps,6))                               #preallocate known step values, replacing memory rather than making new array every time
        self.t_vals[0] = 0
        self.y_vals[0,:] = self.y0
        self.step = 1                                                          #zero is first, initial condition, whereas this is the first step beyond initial

        #initialise solver
        self.solve = ode(self.differential_Eq)
        self.solve.set_integrator('DOP853')                                     #If multiple orbits being calculated may reduce error with DOPRI5 / DOP853 in parrallel
        self.solve.set_initial_value(self.y0,0)                                #set i.c.'s                                                              

        #define perturbation dictionary
        self.perturbations = perturbations                                     #call in differential EQ
        self.propagate_Orbit()                                                    

    def propagate_Orbit(self):                                                 #From 2 body ode - for initialising, solving at each instance, and propogating orbit 
        print("Propagating Orbit...")
        # Create arrays to store the time and state values
        self.t_vals = np.zeros((self.n_steps, 1))
        self.y_vals = np.zeros((self.n_steps, 6))                              #each row is position and vel         
        self.t_vals[0] = 0
        self.y_vals[0] = self.y0
        #documentation asks for the solver to be in a while loop to keep integrating
        
        #Propagate orbit
        for i in range(1, self.n_steps):                                       #using vectorized operations, calculate all time steps at once and store in array
            t = self.t_vals[i-1] + self.dt
            y = self.solve.integrate(t)
            self.t_vals[i] = t
            self.y_vals[i] = y
      
   
        self.r_vals = self.y_vals[:, :3]                                       # extract first three columns of y_vals = rvals, position vectors at each timestep x, y, z
        self.v_vals = self.y_vals[:, 3:]                                       # extract last three columns of y_vals = vvals, velocity vectors at each timestep 

        #Value Checks        
        #print("R Values: \t ", self.r_vals)
        #print("V Values: \n ", self.v_vals)
        #print("t Values: \t ", self.t_vals)
        #print("y Values: \t ", self.y_vals)

    def differential_Eq(self,t,y):
        #unpack state
        rx,ry,rz,vx,vy,vz = y                                                  #y is the state, in 3 component directions (3D)
        r = np.array([rx,ry,rz])                                               #position array defined as vector for conveience
        v = np.array([vx,vy,vz])
        
        norm_r = np.linalg.norm(r)                                             #norm of 3d radius vector
        a = -r*self.centralbody['mu'] / norm_r**3                              #two body acceleration
        
        #this can be made into for loop over different perturbative method
        if self.perturbations['J2']:                                           #call the perturbation term if set to true
            z2 = r[2]**2                                                        
            r2 = norm_r**2
            tx = r[0] / norm_r*(5*z2 / r2-1) 
            ty = r[1] / norm_r*(5*z2 / r2-1)     
            tz = r[2] / norm_r*(5*z2 / r2-3)
            
            a_j2 = 1.5*self.centralbody['J2']*self.centralbody['mu'] * self.centralbody['radius']**2 / norm_r**4*np.array([tx,ty,tz])
            a+=a_j2                                                             #add elemnt by element 

        return [vx,vy,vz,a[0],a[1],a[2]]
    
#        if self.perturbations['OTHER!']:
#            CALCS
#            
#           ACCELERATION
#           APPEND
#        
#        return


    def calculate_COEs(self, degrees = True):                                   #after propogating:
        print('Calculating Classical Orbital Elements (COEs)...')
    
        self.COEs = np.zeros((self.n_steps,6))                                  #same as we used to initialise state vectors
        
        for n in range(self.n_steps):
            self.COEs[n,:] = my_T.Convert_RV_to_COEs(self.r_vals[n,:], self.v_vals[n,:], mu = self.centralbody['mu'], degrees = degrees) #selecting all values
        
#    def stop condition functions( self ):

#   to check if spacecraft has matched trajectory of target orbit

    def plot_COEs(self, hours=False, days=False, show_plot=False, save_plot=False, title='Orbital Elements Against Time for Orbit', figsize=(18,6)): #PRODUCES GRID PLOTS OF ORBIT VALS 
        #can change hrs to be days or yrs if we would want that???
        print("Plotting COE relationships...")
        
        #   orbital elements differfor each orbit type.
        fig, axs = plt.subplots(nrows=3,ncols=2,figsize=figsize)
        fig.tight_layout(pad=5.0)
        fig.suptitle(title,fontsize=10)
    
        #a, e, i, ta, aop, raan, [year, month, day, hour]]
        #plot sma
        axs[0,0].plot(self.t_vals, self.COEs[:,0])
        axs[0,0].set_title("\nSemi-Major Axis vs. Time")
        axs[0,0].grid(True)
        axs[0,0].set_ylabel('Semi-Major Axis (km)')
        axs[0,0].set_xlabel('Time')
        
        #plot eccentricity
        axs[1,0].plot(self.t_vals, self.COEs[:,1])
        axs[1,0].set_title("\nEccentricity vs. Time")
        axs[1,0].grid(True)
        axs[1,0].set_ylabel('Eccentricity ()')
        axs[1,0].set_xlabel('Time')

         #plot inclination
        axs[2,0].plot(self.t_vals, self.COEs[:,2])
        axs[2,0].set_title("Inclination vs. Time")
        axs[2,0].grid(True)
        axs[2,0].set_ylabel('Inclination Angle (deg)')
        axs[2,0].set_xlabel('Time')

        #plot true anomaly
        axs[0,1].plot(self.t_vals, self.COEs[:,3])
        axs[0,1].set_title("True Anomaly vs. Time")
        axs[0,1].grid(True)
        axs[0,1].set_ylabel('Angle(deg)')
        axs[0,1].set_xlabel('Time')
        
        #plot arg of periapse
        axs[1,1].plot(self.t_vals, self.COEs[:,4])
        axs[1,1].set_title("Arg of Periapse vs. Time")
        axs[1,1].grid(True)
        axs[1,1].set_ylabel('Arg of Periapse (deg)')
        axs[1,1].set_xlabel('Time')
           
        #plot RAAN
        axs[2,1].plot(self.t_vals, self.COEs[:,5])
        axs[2,1].set_title("RAAN vs. Time")
        axs[2,1].grid(True)
        axs[2,1].set_ylabel(' RAAN (hrs, min, sec?)')
        axs[2,1].set_xlabel('Time')
        
        if show_plot:
            plt.show()
            
        if save_plot:
            plt.savefig(title+'.png', dpi=300)
    

    def plot_3d(self, show_plot=False, save_plot = False, title='Single Trajectory Plot'): 
        print("Constructing diff-eq from state vectors, Plotting 3D figure of Trajectory") 
        fig = plt.figure(figsize=(16,8))
        ax = fig.add_subplot(111,projection = '3d')                            
        ax.plot(self.r_vals[:,0],self.r_vals[:,1],self.r_vals[:,2],'k', label='Trajectory')
        ax.plot([self.r_vals[0,0]], [self.r_vals[0,1]], [self.r_vals[0,2]],'bo',label='Initial Position')     #define initial orbit position wrt visual frame
        
        #plot central body                                                     #create sphere on 3d plot
        lat, long = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]                        #azimuth polar latitude/longitude for grid across surface of sphere, and find coords 
        coord_x = self.centralbody['radius']*np.cos(lat)*np.sin(long)
        coord_y = self.centralbody['radius']*np.sin(lat)*np.sin(long)
        coord_z = self.centralbody['radius']*np.cos(long)
        ax.scatter(coord_x, coord_y, coord_z, s=5, facecolor="grey")
        
        max_val = np.max(np.abs(self.r_vals))

        ax.set_xlim([-max_val,max_val])
        ax.set_ylim([-max_val,max_val])
        ax.set_zlim([-max_val,max_val])

        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')


        ax.set_title(title)
        plt.legend()

        if show_plot:
            plt.subplot_tool()
            plt.show()
        if save_plot:
            plt.savefig(title+".png",dpi=300)