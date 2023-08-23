import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import math as nm

#Additional Features
import planetary_data as my_PD
import tools as my_T


def Default_Settings():                                                         
    return{
        'centralbody':None,
        'COEs':[],
        'SV':[],
        'J2': True,
        'aero':False,
        'moon_gravity':False,
        'solar_gravity':False
        #Additional features to be added to dictionary:                         []
        }

class Orbit_Propagator:
    def __init__(self, state0, tspan, dt, withCOEs=None, degrees=None, centralbody=my_PD.earth, perturbations = Default_Settings(), One_Period=None):
        #initialise COESs and the state inputs as pos and state, as we propagate in state vectors
        if withCOEs:                                                           # if passing coes as initial state 
            self.r0, self.v0 = my_T.Convert_COEs_to_SV(state0, degrees=degrees, mu=centralbody['mu'])#, print_results=print_results)

        else:                                                                  #Pass straight in if not COEs
            self.r0 = state0[:3]                                               #before 3
            self.v0 = state0[3:6]                                              #after 3
                       
        self.y0 = np.concatenate((self.r0, self.v0))
        self.tspan = tspan
        self.dt = dt
        self.centralbody = centralbody

        #total no. of steps
        self.n_steps = int(np.ceil(self.tspan/self.dt))                                      #cieling, rounds any float to the next whole number (2.1 --> 3)

        #initialise variables
        self.t_vals = np.zeros((self.n_steps,1))
        self.y_vals = np.zeros((self.n_steps,6))                               #preallocate known step values, replacing memory rather than making new array every time
        self.t_vals[0] = 0
        self.y_vals[0,:] = self.y0
        self.Altitude = my_T.norm(self.y_vals[:, :3 ], axis = 1) - self.centralbody['radius']
        
        self.step = 1                                                          # First step after defined initial state to avoid removing input data
        
        #initialise solver
        self.solve = ode(self.differential_Eq)
        self.solve.set_integrator('DOP853')                                    #rk order 8, (5,3) If multiple orbits being calculated may reduce error with DOPRI5 / DOP853 in parrallel - Default use
        # self.solve.set_integrator('lsoda')                                   #fastest, with higher uncertainty
        # self.solve.set_integrator('dopri5')                                  #rk45 process
        # self.solve.set_integrator('Radau')                                   #Radau method
        self.solve.set_initial_value(self.y0,0)                                # Set initial conditions of orbit                                                              

        #Define perturbation/stop conditions
        self.perturbations = perturbations
        self.stop_propagation_array = []
        self.stop_Propagation()                                                                                                                                     #[not in use with no active stop conditions]
        if One_Period:
            self.propagate_For_One_Period()
        else:
            self.propagate_Orbit()                                                    


    def assign_stop_Propagation(self):                                         # Assign stop conditions for propagation                                              [currently not in use due to stop conditions setting positions to zero after initial state]
                                                                               # List of stop conditions 
            #while step < step:
            #self.stop_propagation_array.append(self.stop_Deorbit)
            #if hohmann:
            #   self.stop_propagation_array.append(self.stop_periapse)
            
            #self.stop_propagation_array.append(self.stop_TargetAlt)
            #self.stop_propagation_array.append(self.stop_TargetIncl)          #append functions in use during operation
        return self.stop_propagation_array
       
    def stop_Periapsis(self, best_state):#=a(1-e)):                            # Check if perigee condition is met                                                 [currently not in use due to stop conditions setting positions to zero after initial state]

        #x = my_T.norm(self.r_vals[self.step, :3]
        #print("this is x:\n",x)
        #print("this is rmin:\n", best_state)
        #if x <= #min from prop)
        #    return False
        return True
    
    def stop_Deorbit(self):                                                    # Check if deorbit altitude is reached                                              [currently not in use due to stop conditions setting positions to zero after initial state]
        if self.Altitude[self.step] < self.centralbody['deorbit_altitude']:     #defined in my_PD
            print('Deorbit Altitude Reached... Halting Propagation')
            return    False                                                     #if this section of class is true, return FALSE
        return     True                                                         # if condition not met return true: continue propagation
    
    def stop_TargetAlt(self, TargetAltitude):                                  # Check if target altitude is reached                                                [currently not in use due to stop conditions setting positions to zero after initial state]
        if self.Altitude[self.step] == 1.0: #my_T.norm(TargetAltitude[:3]):
            return False
        return True
            
    def stop_TargetIncl(self, inclination):                                    # Check if target inclination is reached                                             [currently not in use due to stop conditions setting positions to zero after initial state]
        if self.Inclination[self.step] == 10000.0: #my_T.norm(inclination[:3]):
            return False
        return True

    def stop_Propagation(self): #check stop conditions                         # Check if any stop conditions are                                                   [currently not in use due to stop conditions setting positions to zero after initial state]                                                  
        for stop_propagation in self.stop_propagation_array:               
            if not stop_propagation():                                          # if condition has been met, this is true, - if ANY condition
                return False                                                    # in the above case, then print False
        return True
        
    def propagate_Orbit(self):                                                  #From 2 body ode - for initialising, solving at each instance, and propogating orbit 
        print("Propagating Orbit...")                                           # Propagate the orbit using the differential equation
         
        # Create arrays to store the time and state values
        self.t_vals = np.zeros((self.n_steps, 1))
        self.y_vals = np.zeros((self.n_steps, 6))                               # Each row contains position and vel         
        self.t_vals[0] = 0
        self.y_vals[0] = self.y0
        
        #Propagate orbit
        while self.solve.successful() and self.step < self.n_steps:
            
            t = self.t_vals[self.step-1] + self.dt
            y = self.solve.integrate(t)
            self.t_vals[self.step] = t
            self.y_vals[self.step] = y
            self.Altitude[self.step] =  my_T.norm(y[:3]) - self.centralbody['radius']
   
            self.step +=1                                                       # increase step
            if not self.stop_Propagation():                               # Update state vectors, check stop conditions, and break if necessary               [currently not in use due to stop conditions setting positions to zero after initial state]
                break
        
        self.periapsis_min_index = np.argmin(self.Altitude)                     # Find index with the lowest altitude (perigee step position)
        self.best_hohmann_start_state = self.y_vals[self.periapsis_min_index]   # Find state vectors at lowest altitude (perigee coords)
        self.apoapsis_min_index = np.argmax(self.Altitude)                      # Find index with the highest altitude (apogee step position)
        self.best_hohmann_end_state = self.y_vals[self.apoapsis_min_index]      # Find state vectors at highest altitude (apogee coords)
        
        self.end_of_segment = self.y_vals[-1]
            
        self.r_vals = self.y_vals[:, :3]                                       # extract first three columns of y_vals = rvals, position vectors at each timestep x, y, z
        self.v_vals = self.y_vals[:, 3:]          
        
        
    def propagate_For_One_Period_Reverse(self):                                #propagates at end of one orbit period back to optimal hohmann location              [currently not in use for Plotting Initial orbit due to Position val = 0 error]
        a = np.linalg.norm(self.r0)
        T = 2 * np.pi * np.sqrt(a**3 / self.centralbody['mu'])                 # Propagate orbit for one orbital period
        n_steps_one_period = int(np.ceil(T / self.dt))                         # no. of steps in period 
        
        # Create arrays to store the time and state values for one period
        self.t_vals_one_period = np.zeros((n_steps_one_period, 1))
        self.y_vals_one_period = np.zeros((n_steps_one_period, 6))
        
        self.t_vals_one_period[-1] = 0
        self.y_vals_one_period[-1] = self.y0
        
        # Propagate orbit for one period
        step = n_steps_one_period-1
        while self.solve.successful() and step >= 0:
            t = self.t_vals_one_period[self.step+1] - self.dt                   # Subtract dt to step backwards in direction
            y = self.solve.integrate(t)
            self.t_vals_one_period[step] = t
            self.y_vals_one_period[step] = y
            step -= 1
        
        self.r_vals = self.y_vals_one_period[:, :3]
        self.v_vals = self.y_vals_one_period[:, 3:]

    def differential_Eq(self,t,y):                                             # Define the differential equations for the orbit propagation
        rx,ry,rz,vx,vy,vz = y                                                  #y is the state, in 3 component directions (3D)
        r = np.array([rx,ry,rz])                                               #position array defined as vector for conveience
        v = np.array([vx,vy,vz])
        norm_r = my_T.norm(r)                                                  #norm of 3d radius vector

        a = -r*self.centralbody['mu'] / norm_r**3                              #two body acceleration
         
        #use same format to loop over additional perturbative methods
        if self.perturbations['J2']:                                           #J2 perturbation effect calc at each step
            z2 = r[2]**2                                                         
            r2 = norm_r**2
            tx = r[0] / norm_r*(5*z2 / r2-1) 
            ty = r[1] / norm_r*(5*z2 / r2-1)     
            tz = r[2] / norm_r*(5*z2 / r2-3)
            
            a_j2 = 1.5*self.centralbody['J2']*self.centralbody['mu'] * self.centralbody['radius']**2 / norm_r**4*np.array([tx,ty,tz])
            a+=a_j2                                                             #add elemnt by element 


        # if self.perturbations['aero']:                                       # Values beyond second term set to zero when implemented                                                     [currently not in use for Unstable effect on orbit eq (likely from definition of rho)]
        #     print("AERO WORKING")
        #     Cd = self.centralbody['Cd']                                      # Drag coefficient (adjust as needed for craft's aero properties)
        #     drag_area = self.centralbody['drag_area']                        #Effective cross-sectional drag area (adjust based on the spacecraft's characteristics)
        #     sc_mass = self.centralbody['sc_mass']                            #spacecraft mass
        #     # Calculate relative velocity of spacecraft with respect to the atmosphere
        #     v_rel = my_T.norm(v)/1000
        #     print(norm_r+self.centralbody['radius'])
        #     # rho = nm.log((100/((norm_r+self.centralbody['radius'])-300.0)))/0.05 + 80.0            #incorrect temporary formula from repoted relations in e.g. \  https://doi.org/10.1029/2020EA001321  /          The Atmospheric Neutral Density Experiment (ANDE) And           Modulating Retroreflector In Space (MODRAS): Combined Flight                Experiments For The Space Test Program
        #     rho = 1 / (nm.exp(0.04*(norm_r-400))) -1
        #     # Calculate drag acceleration
        #     D_drag = -0.5 * Cd * rho * v_rel**2 * drag_area
            
        #     a_drag = D_drag / sc_mass
            
        #     a += a_drag
        """                                                                        
                                                                                ###     Additional Perturbations
                                                                                #       Solar Rad Pressure (Spherical)
                                                                                #       Radiation Pressure (Albedo / thermal)
        """
        return [vx,vy,vz,a[0],a[1],a[2]]
    
    def calculate_COEs(self, degrees = True):                                  # Calculate and store the Classical Orbital Elements (COEs)     #after propogating:
        print('Calculating Classical Orbital Elements (COEs)...')
        self.COEs = np.zeros((self.n_steps,6))                                  #same as we used to initialise state vectors
        for n in range(self.n_steps):
            self.COEs[n,:] = my_T.Convert_SV_to_COEs(self.r_vals[n,:], self.v_vals[n,:], mu = self.centralbody['mu'], degrees = degrees) #selecting all values
        

    #Funtions to generate plots of COEs progression at each timestep 
    def plot_COEs(self, hours=False, days=False, show_plot=False, save_plot=False, title='Orbital Elements Against Time for Orbit', 
                  figsize=(18,6), grid_plot=None, full_plot=None): 
        print("Able to plot COEs if required")
        if grid_plot:                                                           #Plot COEs in grid formation
            print("Plotting COE relationships in grid format...")
            fig, axs = plt.subplots(nrows=3,ncols=2,figsize=figsize)
            fig.tight_layout(pad=5.0)
            fig.suptitle(title,fontsize=10)
        
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
    
        if full_plot:                                                           #Plot all COEs together in format similar to STK
            print("Plotting COE relationships together...")
            fig, ax = plt.subplots(figsize=figsize)
            fig.suptitle(title, fontsize=10)
            
            # Plot inclination
            ax.plot(self.t_vals, self.COEs[:, 2], label="Inclination")
            ax.set_ylabel('Inclination Angle (deg)')
            ax.set_title("Orbital Elements vs. Time")
    
            ax.set_xlabel('Time')
            ax.grid(True)
            # Plot true anomaly
            ax.plot(self.t_vals, self.COEs[:, 3], label="True Anomaly")
            ax.set_ylabel('Angle (deg)')
            
            # Plot arg of periapse
            ax.plot(self.t_vals, self.COEs[:, 4], label="Arg of Periapse")
            ax.set_ylabel('Arg of Periapse (deg)')
            
            # Plot RAAN with a secondary y-axis
            ax.plot(self.t_vals, self.COEs[:, 5], label="RAAN")
            ax.set_ylabel('Angle (degrees)')
            
            # Plot semi-major axis
            ax2 = ax.twinx()
            ax2.plot(self.t_vals, self.COEs[:, 0],'k', label="Semi-Major Axis", )
            ax2.set_ylabel('Distance (km)')
            
            # Adjust spacing and layout
            ax.legend(loc="upper left")
            ax2.legend(loc="upper right")
            fig.tight_layout(pad=3.0)
            
            if show_plot:
                plt.show()
            
            if save_plot:
                plt.savefig(title + '.png', dpi=300)