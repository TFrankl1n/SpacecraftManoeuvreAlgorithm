import os
import sys
import numpy as np
from math import sqrt                                                        
import matplotlib.pyplot as plt

# Get the directory path where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
# Add the script directory to the Python path
sys.path.insert(0, script_dir)

#Import files for calling, and assign variable names
import planetary_data as my_PD
import tools as my_T
from orbit_propagator import Default_Settings
centralbody = my_PD.earth

class Satellite:
    time_moved = 0
    
    def __init__(self, initial_elements, tspan, dt, withCOEs=True, add_perts=True, degrees=True):
        self.initial_elements = initial_elements
        self.tspan = tspan
        self.dt = dt
        self.withCOEs = withCOEs
        self.add_perts = add_perts
        self.degrees = degrees
        self.segments = []                                                      # List to store final state values of previous propagation segment

    def propagate(self, degrees, calculate_COEs, full_plot, grid_plot):
        for orbit_elements in self.initial_elements:                            #Check
            if initial_elements[1] <0.0 or initial_elements[1] >1.0:
                raise ValueError("Eccentricity by definition must be between 0 and 1.")
            # Check if input values are in degrees correctly, defining bounds for true anomaly, argument of perigee, and right ascension of the ascending node
            if degrees:
                if initial_elements[2] <0.0 or initial_elements[2] >180.0:
                    raise ValueError("Inclination in degrees must be between 0 and 180.\n (Prograde orbits between 0 and 90 degrees, retrograde beyond this)")
                if  initial_elements[3] <0.0 or initial_elements[3] >360.0:
                    raise ValueError("True Anomaly in degrees must be between 0 and 360.")
                if  initial_elements[4] <0.0 or initial_elements[4] >360.0:
                    raise ValueError("Argument of Perigee in degrees must be between 0 and 360.")
                if  initial_elements[5] <0.0 or initial_elements[5] >360.0:
                    raise ValueError("Right Asc. of the Ascending Node in degrees must be between 0 and 360 .)")
            else:                                                               # Check if any of the initial elements are beyond the bounds for radians             
                if any(np.abs(initial_elements[2:6]) >= 2 * np.pi):
                    raise ValueError("Initial elements may be in radians or one of the initial conditions may be beyond bounds for radians.\nCurrently taking input in radians, set 'degrees=True' for degrees.")        
            if not calculate_COEs:                                              # Check if any plots have been called without recording arrays to plot them with
                if full_plot:
                    raise ValueError("Cannot plot the COEs over propagation duration without 'calculate_COEs=True'.")
                if grid_plot:
                    raise ValueError("Cannot plot the COEs over propagation duration without 'calculate_COEs=True'.")
        try:        
            if self.time_moved > 0:                                             #check if propagation has already occured so we know whether to use [initial conditions + input values] or [previous state and state vecors]
                print("Beginning subsequent propagation...")
                initial_segment = self.segments[-1]
                initial_segment = np.concatenate((initial_segment.r_vals[-1], initial_segment.v_vals[-1]))

                op = my_T.Propagate_Orbit(                                      # Propagate the orbit using last segment elements and settings
                                    initial_segment,
                                    self.tspan,
                                    self.dt,
                                    withCOEs=False,
                                    add_perts=self.add_perts,
                                    degrees=self.degrees,
                                    calculate_COEs=False,
                                    full_plot=False,
                                    grid_plot=False
                )

            else:
                
                print("Beginning first propagation...")
                op = my_T.Propagate_Orbit(                                      # Propagate the orbit using provided initial elements and settings         
                                    self.initial_elements,
                                    self.tspan,
                                    self.dt,
                                    withCOEs=self.withCOEs,
                                    add_perts=self.add_perts,
                                    degrees=self.degrees,
                                    calculate_COEs=False,
                                    full_plot=False,
                                    grid_plot=False
                                )
                self.time_moved +=1                                             # increment time_moved to only use segment values with later propagation
            
            self.segments.append(op)                                            # Append final state of propagation step to segemnts array

        except:
            raise ValueError("Check input orbital elements type. If using classical orbital elements (or TLE input), set 'withCOEs=True', if using state vectors, set 'withCOES=False")
        
    def perform_hohmann_transfer_optimal(self, smachange): 

        a_new = self.initial_elements[0]+smachange
       
        # Propagate to the end of the initial segment
        initial_segment = self.segments[-1]
        print("this is where perf_hohmann begins", self.segments[-1])

        #calculate initial COEs for transfer orbit
        Hohmann_transfer_initial_COEs = my_T.Convert_SV_to_COEs(initial_segment.r_vals[-1], initial_segment.v_vals[-1])
        # Calculate final COEs for transfer orbit
        Hohmann_transfer_final_COEs = [a_new if i == 0 else self.initial_elements[i] if i==1 else x + 180.0 if i == 3 else x for i, x in enumerate(Hohmann_transfer_initial_COEs)]
        #Hohmann_transfer_final_COEs = [a_new, Hohmann_transfer_initial_COEs[1], Hohmann_transfer_initial_COEs[2], Hohmann_transfer_initial_COEs[3], Hohmann_transfer_initial_COEs[4], Hohmann_transfer_initial_COEs[5]]

        print("Transfer initial COEs: ", Hohmann_transfer_initial_COEs)
        print("Transfer final COEs: ", Hohmann_transfer_final_COEs)

        #Propagate transfer with above conditions
        op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer(
                                                                                                self.initial_elements[0], 
                                                                                                a_new,
                                                                                                COEs0=Hohmann_transfer_initial_COEs, 
                                                                                                COEs1=Hohmann_transfer_final_COEs,
                                                                                                dt=100, 
                                                                                                propagate=True, 
                                                                                                add_perts=True
                                                                                            )  
        
        print('DeltaV 0: \t %.3f (km/s)' % DeltaV_Vals[0])                    # Display Calculated DeltaV and time of flight values
        print('DeltaV 1: \t %.3f (km/s)' % DeltaV_Vals[1])
        print('Total DeltaV: ', DeltaV_Vals[1] + DeltaV_Vals[0])
        print('Eccentricity of Transfer: ', e_transfer)
        print('Transfer Time: \t',tof_transfer,' seconds, ',tof_transfer/3600, "hours")

        self.segments.append(op0)                                               # Append initial segment to segments list (one complete orbit)
        self.segments.append(op_transfer)                                       # Append transfer segment to segments list
        self.segments.append(op1)                                               # Append (post second burn) segment to segments list (one complete orbit)
        self.time_moved +=1      
        pass

    def perform_hohmann_transfer(self, smachange):
        a_new = self.initial_elements[0] + smachange
        
        final_segment = self.segments[-1]                                       # Extract the final state of the last propagated segment
        r_final = final_segment.r_vals[-1]
        v_final = final_segment.v_vals[-1]
        
        # Calculate initial COEs for transfer orbit
        Hohmann_transfer_initial_COEs = my_T.Convert_SV_to_COEs(r_final, v_final)
        # Calculate final COEs for transfer orbit
        Hohmann_transfer_final_COEs = [a_new if i == 0 else x + 180.0 if i == 3 else x for i, x in enumerate(Hohmann_transfer_initial_COEs)]
    
    
        # Propagate transfer with above conditions
        op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer(
                                                                                                self.initial_elements[0], a_new,
                                                                                                COEs0=Hohmann_transfer_initial_COEs,
                                                                                                COEs1=Hohmann_transfer_final_COEs,
                                                                                                dt=100,
                                                                                                propagate=True,
                                                                                                add_perts=True
                                                                                            )
        print('13DeltaV 0: \t %.3f (km/s)' % DeltaV_Vals[0])
        print('13DeltaV 1: \t %.3f (km/s)' % DeltaV_Vals[1])
        print('13Total DeltaV: ', DeltaV_Vals[1] + DeltaV_Vals[0])
        print('13Eccentricity of Transfer: ', e_transfer)
        print('13Transfer Time: \t',tof_transfer,' seconds, ',tof_transfer/3600, "hours")

        self.segments.append(op0)                                               # Append initial segment to segments list (one complete orbit)
        self.segments.append(op_transfer)                                       # Append transfer segment to segments list
        self.segments.append(op1)                                               # Append (post second burn) segment to segments list (one complete orbit)
       
        self.time_moved +=1

        pass

    def perform_inclination_change(self, delta_i):
        #TODO                                                   check statements for value inputs - convert to ensure correct calculations below
        
        # Calculate the inclination change in radians
        delta_i_rad = np.radians(delta_i)
        
        # Get the last propagated state from the segments
        final_segment = self.segments[-1]
        r_final = final_segment.r_vals[-1]
        v_final = final_segment.v_vals[-1]
        
        # Calculate the new inclination value
        new_i_rad = np.arccos(np.cos(np.radians(self.initial_elements[2])) * np.cos(delta_i_rad) - np.sin(np.radians(self.initial_elements[2])) * np.sin(delta_i_rad))
        new_i_deg = np.degrees(new_i_rad)
        
        # Modify the inclination in the COEs
        modified_COEs = self.initial_elements.copy()
        modified_COEs[2] = new_i_deg
        
        # Propagate the orbit using the modified COEs
        op = my_T.Propagate_Orbit(
                                    modified_COEs,
                                    self.tspan,
                                    self.dt,
                                    withCOEs=self.withCOEs,
                                    add_perts=self.add_perts,
                                    degrees=True,  # Always use degrees for COEs
                                    calculate_COEs=False,
                                    full_plot=False,
                                    grid_plot=False
                                )
        
        # Append the propagated segment to the segments list
        self.segments.append(op)
        self.time_moved+=1
        pass

    def plot_orbits(self):
        r_vals_list = [segment.r_vals for segment in self.segments]
        labels = [f"Segment {i}" for i in range(1, len(self.segments) + 1)]
        my_T.Plot_N_Orbits(r_vals_list, labels=labels, show_plot=True, title="Orbit Propagation")




###############################################################################
#Enter .txt file location from main below:

filename_ISS = 'Files/TLE_1_ISS.txt'
file_path_ISS_ = os.path.join(script_dir, filename_ISS)
file_path_ISS = file_path_ISS_.replace('\\','/')

filename_Deb = 'Files/TLE_2_Debris.txt'
file_path_Deb_ = os.path.join(script_dir, filename_Deb)
file_path_Deb = file_path_Deb_.replace('\\','/')

filename_SL = 'Files/TLE_3_Starlink-30095.txt'
file_path_SL_ = os.path.join(script_dir, filename_SL)
file_path_SL = file_path_SL_.replace('\\','/')

print(file_path_ISS)                                                            #Check any file path to ensure correctly denoted in software for calling

# script_dir = os.path.dirname(os.path.abspath(__file__))                       #Check directory of python script 
# print("Current script location:", script_dir)

###############################################################################
###::::::::::::::Enter your Desired Propagation segements below::::::::::::::###










############################################################################### Control segments
####Blank initial required variables:
# initial_elements = [    ,   ,   ,   ,   ,   ]                                 # Enter Object initial conditions in either COE or SV form
# dt =                                                                          # Enter time step of propagation (lower values will produce higher fidelity plot, and increase runtime)

#### Add Satellite Object
#    = Satellite(                                                               # Enter Name of Satellite Object
            # initial_elements, 
            # tspan, 
            # dt, 
            # withCOEs=True, 
            # add_perts=True, 
            # degrees=True
            #     )


#### Perform Hohmann transfer
# smachange =                                                                   # Enter desired change in semimajor axis
# sat_1.perform_hohmann_transfer_optimal(smachange) 

#### Perform inclination change
# delta_i =                                                                     # Enter desired change in inclination 
# sat_1.perform_inclination_change(delta_i)

#### Perform orbit propagation
#tspan =                                                                        # Enter desired propagation duration 
# sat_1.propagate(
#     degrees=True, 
#     calculate_COEs=False, 
#     full_plot=False, 
#     grid_plot=False, 
#               )



############################################################################### Example sequence:
    
# Initialize an orbit
initial_elements = [6630, 0.0, 15.0, 15.0, 30.0, 10.0]
dt = 60.0

tspan = 84500.0
sat_1 = Satellite(initial_elements, tspan, dt, withCOEs=True, add_perts=True, degrees=True)

# Perform propagation
tspan = 84500.0
sat_1.propagate(
    degrees=True, 
    calculate_COEs=False, 
    full_plot=False, 
    grid_plot=False, 
    )




# # Perform inclination change
delta_i = 10
sat_1.perform_inclination_change(delta_i)

# Perform Hohmann transfer
smachange = 6500
sat_1.perform_hohmann_transfer_optimal(smachange) #ORIGINAL HOHMANN TRANSFER NEEDS TO BE DEFINED TO TAKE END OF TRANSFER VALUES

#Perform additional propagation if needed
tspan = 44588
sat_1.propagate(
    degrees=True, 
    calculate_COEs=False, 
    full_plot=False, 
    grid_plot=False,
    )

# Plot all segments
sat_1.plot_orbits()


#TODO
    # Perfect Hohmann transfer 
        # BEGIN OP0 PROPAGATION FROM ONE ORBIT BEFORE OPTIMUM POSITION (BLUE CROSS AND WHERE TRANSFER BEGINS)
        # OP1 DOES NOT START WHERE THE PERFECT TRANSFER ENDS, SO MAKE THIS HAPPEN (CAN JUST ENSURE APPENDED COEs ARE THE FINAL OF THE TRANSFER)
    #Standard Hohmann transfer
        # ALLLLL



"""
def propagate_and_plot_multiple_orbits(initial_orbits, tspan, dt, labels=None, withCOEs=True, add_perts=True, degrees=True, calculate_COEs=True, full_plot=True, grid_plot=True, plot_3d=False):
    orbits = []
    #                               ):
  
    for orbit_elements in initial_orbits:
        if initial_elements[1] <0.0 or initial_elements[1] >1.0:
            raise ValueError("Eccentricity by definition must be between 0 and 1.")
        # Check if input values are in degrees correctly, defining bounds for true anomaly, argument of perigee, and right ascension of the ascending node
        if degrees:
            if initial_elements[2] <0.0 or initial_elements[2] >180.0:
                raise ValueError("Inclination in degrees must be between 0 and 180.\n (Prograde orbits between 0 and 90 degrees, retrograde beyond this)")
            if  initial_elements[3] <0.0 or initial_elements[3] >360.0:
                raise ValueError("True Anomaly in degrees must be between 0 and 360.")
            if  initial_elements[4] <0.0 or initial_elements[4] >360.0:
                raise ValueError("Argument of Perigee in degrees must be between 0 and 360.")
            if  initial_elements[5] <0.0 or initial_elements[5] >360.0:
                raise ValueError("Right Asc. of the Ascending Node in degrees must be between 0 and 360 .)")
        # Check if any of the initial elements are beyond the bounds for radians
        else:                                                                       
            if any(np.abs(initial_elements[2:6]) >= 2 * np.pi):
                raise ValueError("Initial elements may be in radians or one of the initial conditions may be beyond bounds for radians.\nCurrently taking input in radians, set 'degrees=True' for degrees.")        
        # Check if any plots have been called without recording arrays to plot them with
        if not calculate_COEs:
            if full_plot:
                raise ValueError("Cannot plot the COEs over propagation duration without 'calculate_COEs=True'.")
            if grid_plot:
                raise ValueError("Cannot plot the COEs over propagation duration without 'calculate_COEs=True'.")
    try:        
        # Propagate the orbit using provided initial elements and settings
        op = my_T.Propagate_Orbit(orbital_elements, tspan, dt, withCOEs=withCOEs, add_perts=add_perts, degrees=degrees)
        orbits.append(op)
        
        if labels is None:
            labels = [f"Orbit {i+1}" for i in range(len(orbits))]
    
         # Plot the 3D orbit trajectory if requested
        if plot_3d:
            my_T.Plot_N_Orbits([op.r_vals for op in orbits], labels=labels, show_plot=True, title="Multiple TLE Plot")
    except:
        raise ValueError("Check input orbital elements. If using classical orbital elements (or TLE input), set 'withCOEs=True', if using state vectors, set 'withCOES=False'. Ensure tspan or dt are nonzero")

if __name__ == '__main__':
    print('Running...')
    
    # Define initial orbits, time span, and time step
    # orbital_elements = [[7000, 0.0, 0.523599 , 0.802851 , 0.0, 3.71681],
    #                     [8736.1, 0.3, 51.6414, 280.3421, 153.847, 5.5858]]
    # orbital_elements = [7000, 0.0, 0.523599 , 0.802851 , 0.0, 3.71681]
    # tspan = 11470
    # dt = 60.0
    """

"""
    Analyzes an orbit using the given parameters.

    Args:
        initial_elements    (list): List of initial elements either:
                                            classical orbital elements [a, e, i, ta, aop, raan].
                                            state vectors [rx, ry, rz, vx, vy, vz].
        tspan               (float): Time span for propagation.
        dt                  (float): Time step for propagation.
        withCOEs            (bool): Whether the initial_elements are classical orbital elements (=True) or state vectors (=False)
        add_perts           (bool): Whether to include J2 perturbative effects. **********************************************************************add options
        degrees             (bool): Whether the input elements are in degrees (=True), radians (=False).
        calculate_COEs      (bool): Whether to calculate COEs at each step of propagation:
            full_plot           (bool): Whether to plot all classical orbital elements on single chart against time
            grid_plot           (bool): Whether to plot all classical orbital elements seperately against time
        plot_3d             (bool): Whether to plot the orbit trajectory in 3D.

    Returns:
        op: The result of orbit propagation.
"""





"""
    #SUBSEQUENT plotting - Single Hohmann from periapsis
    #[       a,           e,        i,      ta,  aop,   raan]
    tspan = 89400.0                                                          #timespan of initial propagation
    #dt = 60
    dt = 100
    c0 = [10736.1, 8.64434e-16, 51.6414, 28.036, 153.847, 5.5858]                  #Initial COEs
    c0 = [28736.1, 0.75, 51.6414, 229.577, 153.847, 5.5858]                  #Initial COEs
      
    #smachange = 6500                                                       #Desired change in sma after transfer
      
    #a_new = c0[0]+smachange
      
    op = my_T.Propagate_Orbit(c0, tspan, dt, withCOEs=True, add_perts=True, degrees=True, calculate_COEs=False, full_plot=False, grid_plot=False)#grid_plot=None, perturbations=Default_Settings())

    print("SECOND PROPAGATION\n")
    #op_ = my_T.Propagate_Orbit(op.best_hohmann_start_state, tspan, dt, withCOEs=False, add_perts=True, calculate_COEs=False, full_plot=False, grid_plot=False) #grid_plot=None, perturbations=Default_Settings())
    op_2 = my_T.Propagate_Orbit(op.end_of_segment, tspan, dt, withCOEs=False, add_perts=True, calculate_COEs=False, full_plot=False, grid_plot=False) #grid_plot=None, perturbations=Default_Settings())

    my_T.Plot_N_Orbits([op.r_vals, op_2.r_vals],labels=['Initial Path', 'Segment continued'],show_plot=False, title = 'Transfer Plot')
      
    #print(my_T.Convert_SV_to_COEs(op1.r_vals[-1], op1.v_vals[-1]))
           
 
       
        # #################################################################CHECK THE PLOT AGAINST THE T_VAL THAT IS OPTIMUM FOR YOU!!!!! - 
    ##### for coes to ssv multiple conversion
#     initial_orbits = [
#     [7000, 0.0, 0.523599 , 0.802851 , 0.0, 3.71681],
#     [8736.1, 0.3, 51.6414, 280.3421, 153.847, 5.5858]
# ]
# COEs_radial_vectors = Convert_COEs_to_SV_multiple(initial_orbits, degrees=True, mu=my_PD.earth['mu'])1
"""
    
def SINGLE(): #ORBIT PATH FROM COEs flags
            """ 
            print perigee from prop
            print apogee from prop
            print perigee from calc
            print apogee from calc
             
            
            print t_vals ?>
                - print to txt??
            
            print v_vals over time 
                - print to txt??
            plot v_vals over time
            print r_vals over time 
                - print to txt??
            plot r_vals over time
            
            """
            return
###############################################################################
   
###############################################################################
        # #MULTIPLE ORBIT PROPAGATION AND PLOT FROM COEs WITH PERTURBATION
        # tspan = 156470
        # c0 = [7000, 0.0, 0.523599 , 0.802851 , 0.0, 3.71681]
        # c1 = [8736.1, 0.3, 51.6414, 280.3421, 153.847, 5.5858] 
        # op0 = my_T.Propagate_Orbit(c0, tspan, dt, withCOEs = True, add_perts=True, degrees=True)#, perturbations=Default_Settings())
        # op1 = my_T.Propagate_Orbit(c1, tspan, dt, withCOEs = True, add_perts=True, degrees=True)#, perturbations=Default_Settings())
        

        # my_T.Plot_N_Orbits([op0.r_vals, op1.r_vals], labels=["Orbit 1", "Orbit 2"], show_plot=True, title = 'Multiple TLE Plot')#' with J2 Perturbation over 2 Days' )
        
        # print(op0.r_vals)
        # print("MIN VALUES HERE:\n\n",np.amin(my_T.norm(op0.r_vals, True)))
        # print(np.argmin(my_T.norm(op0.r_vals, True)))
        
###############################################################################


###############################################################################
        # #SINGLE ORBIT PROP (From TLEs)
        # print("Calculating opop")
        # opop = my_T.Propagate_Orbit(my_T.TLE_to_CEOsPlus('F:/Users/Thomas/UNIVERSITY/Work [] Placements/WORK/5. AIRBUS/1. Work/Manoeuvres/Astro_Dynamics/8_J2/TLE_1_ISS.txt'),
        #                               tspan, dt, withCOEs=True, degrees=True)
        # opop = my_T.Propagate_Orbit(my_T.TLE_to_CEOsPlus(file_path_ISS),
        #                               tspan, dt, withCOEs=True, degrees=True, add_perts=True)
                                      
        # my_T.plot_3d(opop.r_vals, show_plot=(True), title="ISS single TLE Plot")
###############################################################################


###############################################################################
        # #MULTIPLE ORBIT PROP (From TLEs)
        # op0 = my_T.Propagate_Orbit(my_T.TLE_to_CEOsPlus('F:/Users/Thomas/UNIVERSITY/Work [] Placements/WORK/5. AIRBUS/1. Work/Manoeuvres/Astro_Dynamics/8_J2/TLE_1_ISS.txt'), 
        #               tspan, dt, add_perts=True)
        # op1 = my_T.Propagate_Orbit(my_T.TLE_to_CEOsPlus('F:/Users/Thomas/UNIVERSITY/Work [] Placements/WORK/5. AIRBUS/1. Work/Manoeuvres/Astro_Dynamics/8_J2/TLE_2_Debris.txt'), 
        #               tspan, dt, add_perts=True)#, COEs = True, deg=False, add_perts=True)
        # op2 = my_T.Propagate_Orbit(my_T.TLE_to_CEOsPlus('F:/Users/Thomas/UNIVERSITY/Work [] Placements/WORK/5. AIRBUS/1. Work/Manoeuvres/Astro_Dynamics/8_J2/TLE_3_Starlink-30095.txt'), 
        #               tspan, dt, add_perts=True)#, COEs = True, deg=False, add_perts=True)
       
        
        # op0 = my_T.Propagate_Orbit(my_T.TLE_to_CEOsPlus(file_path_ISS), tspan, dt, withCOEs=True, degrees=True, add_perts=True)
        # op1 = my_T.Propagate_Orbit(my_T.TLE_to_CEOsPlus(file_path_Deb), tspan, dt, withCOEs=True, degrees=True, add_perts=True)#, COEs = True, deg=False, add_perts=True)
        # op2 = my_T.Propagate_Orbit(my_T.TLE_to_CEOsPlus(file_path_SL), tspan, dt, withCOEs=True, degrees=True, add_perts=True)
        # #op0 = my_T.Initiate_Propagate(c0, tspan, dt, COEs=True)
        # #my_T.Plot_N_Orbits([op0.r_vals, op1.r_vals], show_plot=True, title = 'Multiple TLE Plot')#' with J2 Perturbation over 2 Days' )
        # tspan_days = tspan/(3600 * 24.0)
        # my_T.Plot_N_Orbits([op0.r_vals,op1.r_vals,op2.r_vals], labels=['ISS','Debris','Starlink-30095'],show_plot=True,
        #                     title = f"Multiple TLE Plot with J2 Perturbation over {tspan_days:.1f} Day(s)")
###############################################################################


###############################################################################        
# #SINGLE HOHMANN TRANSFER - [Value calc (a, v, tof, e)], []
# r0 = 16878.1370050261875804
# r1 = 20678.1370050261875804 
# # r0 = 6736.14
# # r1 = 7000.67
# #[    a,  e,   i,    ta,   aop, raan, [year, month, day, hour]]
# # c0 = [6736.14, 0.00011, 51.6414, 280.3421, 153.8469, 5.5858] 
# # c1 = [6736.14, 0.00011, 51.6414, 280.3421, 153.8469, 5.5858] 
# # c0 = [7000.0, 0.0, 0.523599 , 0.802851 , 0.0, 3.71681]
# # c0 = [r0, 0.5, 51.6414, 100.3421, 153.8469, 5.5858]   
# # c1 = [r1, 0.5, 51.6414, 280.3421, 153.8469, 5.5858]   
# # c0 = [r0, 0.1149999999999999, 62.5, 0.0, 0.0, 10.0]
# # c1 = [r1, 0.1149999999999999, 62.5, 180.0,0.0, 10.0]
# c0 = [8736.1, 0.1, 51.6414, 229.577, 153.847, 5.5858] #[7000, 0.0, 9.2, 0.802851, 0.0, 3.71681]
# c1 = [15239.1, 0.1, 51.6414, 9.577, 153.847, 5.5858] #[7000, 0.0, 9.2, 0.802851, 0.0, 3.71681]

# op, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer(r0=r0, r1=r0, COEs0=c0, COEs1=c1, dt=100, 
#                                                                                     propagate=True, add_perts=(True))  
# #op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer_Circ(r0=r0, r1=r0, COEs0=c0, COEs1=c1, dt=100, 
# #                                                                                          propagate=True, add_perts=(True))  
# print('DeltaV 0: \t %.3f (km/s)' % DeltaV_Vals[0])
# print('DeltaV 1: \t %.3f (km/s)' % DeltaV_Vals[1])
# print('Total DeltaV: ', DeltaV_Vals[1] + DeltaV_Vals[0])
# print('Eccentricity of Transfer: ', e_transfer)
# print('Transfer Time: \t',tof_transfer,' seconds, ',tof_transfer/3600, "hours")

# print("\nrp = \t %.1f (km)", c0[0]*(1-c0[1]))
# my_T.Plot_N_Orbits([op.r_vals, op1.r_vals, op_transfer.r_vals],labels=['Orbit 1 Initial', 'Orbit 1 Final', 'Orbit 1 Transfer'],show_plot=False, title = 'Transfer Plot')

# r2 = 2887
# r3 = 6465
# c2_rot = [centralbody['radius']+r2, 0.02, 30.0,   0.0, 100.0, 20.0]
# c3_rot = [centralbody['radius']+r3, 0.02, 30.0, 180.0, 100.0, 20.0]
###############################################################################


###############################################################################        
        # #SINGLE HOHMANN TRANSFER - [Value calc (a, v, tof, e)], []
        # r0 = 16878.1370050261875804
        # r1 = 20678.1370050261875804 

        # #[    a,  e,   i,    ta,   aop, raan, [year, month, day, hour]]
        # c0 = [r0, 0.314, 62.5, 0.0, 0.0, 10.0]
        # c1 = [r1, 0.314, 62.5, 180.0,0.0, 10.0]
        # op, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer(r0=r0, r1=r0, COEs0=c0, COEs1=c1, dt=100, 
        #                                                                                           propagate=True, add_perts=(True))  
        # #op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer_Circ(r0=r0, r1=r0, COEs0=c0, COEs1=c1, dt=100, 
        # #                                                                                          propagate=True, add_perts=(True))  
        # print('DeltaV 0: \t %.3f (km/s)' % DeltaV_Vals[0])
        # print('DeltaV 1: \t %.3f (km/s)' % DeltaV_Vals[1])
        # print('Total DeltaV: ', DeltaV_Vals[1] + DeltaV_Vals[0])
        # print('Eccentricity of Transfer: ', e_transfer)
        # print('Transfer Time: \t',tof_transfer,' seconds, ',tof_transfer/3600, "hours")
        
        # print("\nrp = \t ", c0[0]*(1-c0[1]), "(km)")
        # my_T.Plot_N_Orbits([op.r_vals, op1.r_vals, op_transfer.r_vals],labels=['Orbit 1 Initial', 'Orbit 1 Final', 'Orbit 1 Transfer'],show_plot=False, title = 'Transfer Plot')
        
        # r2 = 2887
        # r3 = 6465
        # c2_rot = [centralbody['radius']+r2, 0.02, 30.0,   0.0, 100.0, 20.0]
        # c3_rot = [centralbody['radius']+r3, 0.02, 30.0, 180.0, 100.0, 20.0]
##############################################################################


###############################################################################
        # # PLOTTING TWO  HOHMANN TRANSFERS, ONE PLANAR, ONE WITH INCLUNATION
        
        # r0 = 210
        # r1 = 845
        # r2 = 2887
        # r3 = 6465
        # c2 = [centralbody['radius']+r0, 0.0, 0.0, 0.0,  0.0, 0.0]
        # c3 = [centralbody['radius']+r1, 0.0, 0.0, 180.0,0.0, 0.0]
        # c2_rot = [centralbody['radius']+r2, 0.0, 30.0,   0.0, 100.0, 20.0]
        # c3_rot = [centralbody['radius']+r3, 0.0, 30.0, 180.0, 100.0, 20.0]
        
        # #op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer
        # op2, op3, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer(COEs0 = c2, COEs1=c3, propagate = True, add_perts=(True))  
        # # op2, op3, op_transfer, DeltaV_Vals = my_T.Hohmann_Transfer_Circ(propagate = True)  
        
        # op2_rot, op3_rot, op_transfer_rot, DeltaV_Vals_rot, tof_transfer, e_transfer = my_T.Hohmann_Transfer(COEs0 = c2_rot, COEs1=c3_rot, propagate = True)
        # # op2_rot, op3_rot, op_transfer_rot, DeltaV_Vals_rot = my_T.Hohmann_Transfer_Circ(propagate = True)
        
        # print("op2",len(op2.r_vals))
        # print("op3",len(op3.r_vals))
        # print("op2_rot",len(op2_rot.r_vals))
        # print("op3_rot",len(op3_rot.r_vals))
        # my_T.Plot_N_Orbits([op2.r_vals, op3.r_vals, op_transfer.r_vals, op2_rot.r_vals, op3_rot.r_vals, op_transfer_rot.r_vals], 
        #                     labels=['Orbit 1 Initial', 'Orbit 1 Final', 'Orbit 1 Transfer', 'Orbit 2 Initial', 'Orbit 2 Final', 'Orbit 2 Transfer'],
        #                     show_plot=True, title = 'Transfer Plot')
                            
        # print("Here are the DeltaV Values:\n", DeltaV_Vals)
        # print("And the rotated ones why not:\n", DeltaV_Vals_rot)
###############################################################################


###############################################################################
# # PLOTTING INCLINATION TRANSFER
# r0 = 250 #km
# r1 = 250 #km
# #    [                       a,  e,   i,    ta,   aop, raan]
# c0 = [centralbody['radius']+r0, 0.0, 130.0, 0.0, 30.0, 10.0]
# c1 = [centralbody['radius']+r1, 0.0, 120.0, 0.0, 30.0, 10.0]
# op0, op1, deg_change, DeltaV_Val = my_T.Inclination_Transfer(COEs0=c0, COEs1=c1, propagate=True)

# my_T.Plot_N_Orbits([op0.r_vals, op1.r_vals],
#                     labels=['Orbit 1 Initial', 'Orbit 1 Final', 'Orbit 1 Transfer', 'Orbit 2 Initial', 'Orbit 2 Final', 'Orbit 2 Transfer'],
#                     show_plot=True, title = 'Transfer Plot')
                    
# print("Here is the DeltaV of the transfer:\n", DeltaV_Val)
# ###############################################################################


###############################################################################
        #SUBSEQUENT plotting - periapsis detection - with j2
        #[       a,           e,        i,      ta,  aop,   raan]
# tspan = 1456
# dt=100

# c0 = [8736.1, 0.3, 51.6414, 280.3421, 153.847, 5.5858]  
# #propagate orbit from initial COEs:
# op = my_T.Propagate_Orbit(c0, tspan, dt, withCOEs=True, add_perts=True, degrees=True, calculate_COEs=False, full_plot=False, grid_plot=False)#grid_plot=None, perturbations=Default_Settings())
# #propagate orbit from periapsis
# op_ = my_T.Propagate_Orbit(op.best_hohmann_start_state, tspan, dt, withCOEs=False, add_perts=True, calculate_COEs=False, full_plot=False, grid_plot=False)#grid_plot=None, perturbations=Default_Settings())

# my_T.plot_3d(op.r_vals, show_plot=(True), title="Single COE Plot Initial")
# my_T.plot_3d(op_.r_vals, show_plot=(True), title="Single COE Plot Final")
# my_T.Plot_N_Orbits([op.r_vals, op_.r_vals], labels=['Input Coes', 'RVs from periapse (Best Value)'],show_plot=True, title = 'Transfer Plot')

# print(op_.r_vals)
# print("smallest value:\t",np.amin(my_T.norm(op_.r_vals, True)))
# print("smallest value pos:\t",np.argmin(my_T.norm(op_.r_vals, True)))
# ###############################################################################


###############################################################################
        # SUBSEQUENT plotting - Single Hohmann from periapsis
        # [       a,           e,        i,      ta,  aop,   raan]
        #             """
        #             tspan = 434000.0                                                          #timespan of initial propagation
        #             #dt = 60
        #             dt = 100
        #             c0 = [10736.1, 8.64434e-16, 51.6414, 28.036, 153.847, 5.5858]                  #Initial COEs
        #             c0 = [28736.1, 0.1, 51.6414, 229.577, 153.847, 5.5858]                  #Initial COEs
                  
        #             smachange = 6500                                                       #Desired change in sma after transfer
                  
        #             a_new = c0[0]+smachange
        #             """
        # c0 = [18736.1, 0.0, 51.6414, 280.342, 153.847, 5.5858]   
        # propagate orbit from initial COEs:
        # perigee_rad = c0[0]*(1-c0[1]) 
        # cbperigee_rad = c0[0]*(1-c0[1]) - centralbody['radius']
        # apogee_rad = c0[0]*(1+c0[1])
        # cbapogee_rad = c0[0]*(1+c0[1]) - centralbody['radius']
        # print("rp",perigee_rad)#, cbperigee_rad)
        # print("rp",apogee_rad)#, cbapogee_rad)
      
        # sma_new_rp = perigee_rad / (1-c0[1])
        # sma_new_ra = apogee_rad / (1+c0[1])
      
      
        # user enrty for new rp or ra:
        # rp_new = 18241.65999999999
        # smanew = rp_new / (1-c0[1])
        # print(smanew)
      
        # not strictly right, as we want to raise ra first from rp:
        # rp_change = 7000
        # perigee_rad_new = perigee_rad + rp_change
        # smanew2 = (perigee_rad+rp_change) / (1-c0[1])
        # print(smanew2)
      
        # ra_change = 7000
        # apogee_rad_new = apogee_rad+ra_change
        # smanew3 = (apogee_rad+ra_change) / (1-c0[1])
        # print(smanew3)
      
        # print("new a rp:", sma_new_rp)
        # print("new a ra:", sma_new_ra)
      
      
        # generate target COEs
        # print(Hohmann_transfer_initial_COEs)
        # periapsis = my_T.norm(op.best_hohmann_start_state)
        # apoapsis = my_T.norm(op.best_hohmann_end_state+smachange/2.0)
        # print("THIS IS THE PERIAPSIS",periapsis)
        # propagate orbit from periapsis
        # print("SECOND PROPAGATION\n")
        # op_ = my_T.Propagate_Orbit(op.best_hohmann_start_state, tspan, dt, withCOEs=False, add_perts=True, calculate_COEs=False, full_plot=False, grid_plot=False) #grid_plot=None, perturbations=Default_Settings())
        #         """
        #         op = my_T.Propagate_Orbit(c0, tspan, dt, withCOEs=True, add_perts=True, degrees=True, calculate_COEs=False, full_plot=False, grid_plot=False)#grid_plot=None, perturbations=Default_Settings())
              
        #         # print("SUPPOSED BEST PERIAPSIS",my_T.norm(op.best_hohmann_start_state[:3]))
        #         print("13\n\nTHIS IS THE ONE: ",op.best_hohmann_end_state)
              
        #         Hohmann_transfer_initial_COEs = my_T.Convert_SV_to_COEs(op.best_hohmann_start_state[:3], op.best_hohmann_start_state[3:])
        #         Hohmann_transfer_final_COEs = [a_new if i == 0 else x + 180.0 if i == 3 else x for i, x in enumerate(Hohmann_transfer_initial_COEs)]
        #         print(Hohmann_transfer_initial_COEs)
        #         print(Hohmann_transfer_final_COEs)
                
        #         my_T.plot_3d(op.r_vals, show_plot=(True), title="Single COE Plot")
        #         """
        # #print(op.r_vals)
        # # my_T.plot_3d(op_.r_vals, show_plot=(True), title="Single COE Plot")
        # print("13\n\ninitial sma:", c0[0])
        # print("13Final sma:", a_new, "\n\n")
        # op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer(c0[0], a_new,
        #                                                                                       COEs0=Hohmann_transfer_initial_COEs, COEs1=Hohmann_transfer_final_COEs,
        #                                                                                       dt=100, 
        #                                                                                       propagate=True, add_perts=(True))  
        # print('13DeltaV 0: \t %.3f (km/s)' % DeltaV_Vals[0])
        # print('13DeltaV 1: \t %.3f (km/s)' % DeltaV_Vals[1])
        # print('13Total DeltaV: ', DeltaV_Vals[1] + DeltaV_Vals[0])
        # print('13Eccentricity of Transfer: ', e_transfer)
        # print('13Transfer Time: \t',tof_transfer,' seconds, ',tof_transfer/3600, "hours")
      
        # # print("\nrp = \t %.1f (km)", c0[0]*(1-c0[1]))
        # # my_T.plot_3d(op.r_vals, show_plot=(True), title="Single COE Plot")
        # my_T.Plot_N_Orbits([op.r_vals, op0.r_vals, op1.r_vals, op_transfer.r_vals],labels=['Initial Path','Initial Orbit ', 'Final Orbit', 'Transfer Orbit'],show_plot=False, title = 'Transfer Plot')
      
        # print(my_T.Convert_SV_to_COEs(op1.r_vals[-1], op1.v_vals[-1]))
   
       
       
       # #################################################################CHECK THE PLOT AGAINST THE T_VAL THAT IS OPTIMUM FOR YOU!!!!! - 
       
       
       #  print("13DeltaV of Manoeuvre:\n", DeltaV_Vals)
       
        
       #  print(op0.t_vals[65])
       
       #  #print(op_.r_vals)
       #  #print("smallest value:\t",np.amin(my_T.norm(op_.r_vals, True)))
       #  #print("smallest value pos:\t",np.argmin(my_T.norm(op_.r_vals, True)))
###############################################################################


###############################################################################
        #SUBSEQUENT plotting - Single Hohmann - identifying optimal place for transfer:
        #[       a,           e,        i,      ta,  aop,   raan]
        
        
        # tspan = 65000.0                                                          #timespan of initial propagation
        # dt = 100
        # c0 = [28736.1, 0.10048, 51.6415, 233.164, 153.861, 5.57665]                  #Initial COEs
        
        # smachange = 6500                                                       #Desired change in sma after transfer
        
        # a_new = c0[0]+smachange
    
        # op = my_T.Propagate_Orbit(c0, tspan, dt, withCOEs=True, add_perts=True, degrees=True, calculate_COEs=False, full_plot=False, grid_plot=False)#grid_plot=None, perturbations=Default_Settings())
        
        
        # #loop over time for state vectors at each timestep
        # print(op.y_vals)
        # print(len(op.y_vals))
        # #print("transfer times:", tof_transfer)
        # #print(op.tof_transfer)
        
        # # add these to complete array
        # print
        # #array to store DeltaV_Values
        # Each_DeltaV_Value = []
        # Each_tof_Value = []
        # # for i in range(1, 10, 1):
        # #     print(i)
        
        
        # for y_states in op.y_vals:
        #     print("y_state value: ", y_states)
        #     COEs_initial = my_T.Convert_SV_to_COEs(y_states[:3], y_states[3:])
        #     print(COEs_initial)
        #     op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer(
        #         r0=y_states[0],
        #         r1=a_new,
        #         COEs0=COEs_initial,
        #         COEs1=[a_new if i == 0 else x + 180.0 if i == 3 else x for i, x in enumerate(COEs_initial)],
        #         dt=100,
        #         propagate=True,
        #         add_perts=True
        #     )
        #     Total_DeltaV_sum = sum(DeltaV_Vals)
        #     Each_DeltaV_Value.append([Total_DeltaV_sum])
        #     Each_tof_Value.append(tof_transfer)
        
        # # # Print collected DeltaV_Vals
        # for i, Total_DeltaV_sum in enumerate(Each_DeltaV_Value):
        #     print(f"Run {i+1} - DeltaV 0: {DeltaV_Vals[0]:.3f} (km/s), DeltaV 1: {DeltaV_Vals[1]:.3f} (km/s)")
        # #for i, tof_transfer in enumerate(Each_tof_Value):
        #     #print(f"Run {i+1} - DeltaV 0: {tof_transfer[0]:.3f} (km/s), DeltaV 1: {DeltaV_Vals[1]:.3f} (km/s)")

        
        # print(len(Each_DeltaV_Value))
        # Each_DeltaV_Value = np.array(Each_DeltaV_Value)
        
        


        #plot delta v over what time the transfer occurs
        # plt.plot(op.t_vals, Each_DeltaV_Value, label='Summed DeltaV')
        # plt.xlabel('Time')
        # plt.ylabel('Summed DeltaV')
        # plt.title('Summed DeltaV vs Time')
        # plt.legend()
        # plt.grid(True)
        # plt.show()
        
        #plot time of flight of transfer to cost of DeltaV
        # plt.plot(Each_tof_Value, Each_DeltaV_Value, label='Summed DeltaV')
        # plt.xlabel('Transfer Time of Flight (seconds)')
        # plt.ylabel('Summed DeltaV (km/s)')
        # plt.title('Summed DeltaV vs Time')
        # plt.legend()
        # plt.grid(True)
        # plt.show()
        
        # # Choose the desired time range [start, stop]
        # selected_plot_timespan = [25000, 60000]
        # # Find the indices that fall within the desired time range
        # selected_plot_timespan_indices = np.where((op.t_vals >= selected_plot_timespan[0]) & (op.t_vals <= selected_plot_timespan[1]))
        
        # # Extract the corresponding time values and summed DeltaV values
        # selected_t_vals = op.t_vals[selected_plot_timespan_indices]
        # selected_delta_v_vals = np.array(Each_DeltaV_Value)[selected_plot_timespan_indices]
        # plt.plot(selected_t_vals, selected_delta_v_vals, label='Summed DeltaV')
        # plt.xlabel('Time')
        # plt.ylabel('Summed DeltaV')
        # plt.title('Summed DeltaV vs Time')
        # plt.legend()
        # plt.grid(True)
        # plt.show()
        
        # print(np.amin(selected_delta_v_vals))
        # print(np.argmin(selected_delta_v_vals))
        # print(op.t_vals[184])
        
        # print(np.amin(Each_DeltaV_Value))
        # print(np.argmin(Each_DeltaV_Value))
        
                # """
                # # print("SUPPOSED BEST PERIAPSIS",my_T.norm(op.best_hohmann_start_state[:3]))
                # print("\n\nTHIS IS THE ONE: ",op.best_hohmann_end_state)
                
                # Hohmann_transfer_initial_COEs = my_T.Convert_SV_to_COEs(op.best_hohmann_start_state[:3], op.best_hohmann_start_state[3:])
                # Hohmann_transfer_final_COEs = [a_new if i == 0 else x + 180.0 if i == 3 else x for i, x in enumerate(Hohmann_transfer_initial_COEs)]
                # print(Hohmann_transfer_initial_COEs)
                # print(Hohmann_transfer_final_COEs)
                
                # #my_T.plot_3d(op.r_vals, show_plot=(True), title="Single COE Plot")
                # #print(op.r_vals)
                # # my_T.plot_3d(op_.r_vals, show_plot=(True), title="Single COE Plot")
                # print("13\n\ninitial sma:", c0[0])
                # print("13Final sma:", a_new, "\n\n")
                # op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer(c0[0], a_new,
                #                                                                                       COEs0=Hohmann_transfer_initial_COEs, COEs1=Hohmann_transfer_final_COEs,
                #                                                                                       dt=100, 
                #                                                                                       propagate=True, add_perts=(True))  
                # print('13DeltaV 0: \t %.3f (km/s)' % DeltaV_Vals[0])
                # print('13DeltaV 1: \t %.3f (km/s)' % DeltaV_Vals[1])
                # print('13Total DeltaV: ', DeltaV_Vals[1] + DeltaV_Vals[0])
                # print('13Eccentricity of Transfer: ', e_transfer)
                # print('13Transfer Time: \t',tof_transfer,' seconds, ',tof_transfer/3600, "hours")
                
                # # print("\nrp = \t %.1f (km)", c0[0]*(1-c0[1]))
                # # my_T.plot_3d(op.r_vals, show_plot=(True), title="Single COE Plot")
                # my_T.Plot_N_Orbits([op.r_vals, op0.r_vals, op1.r_vals, op_transfer.r_vals],labels=['Initial Path','Initial Orbit ', 'Final Orbit', 'Transfer Orbit'],show_plot=False, title = 'Transfer Plot')
                
                # print(my_T.Convert_SV_to_COEs(op1.r_vals[-1], op1.v_vals[-1]))
                
                # print("13DeltaV of Manoeuvre:\n", DeltaV_Vals)
                
                # """
        #print(op_.r_vals)
        #print("smallest value:\t",np.amin(my_T.norm(op_.r_vals, True)))
        #print("smallest value pos:\t",np.argmin(my_T.norm(op_.r_vals, True)))
###############################################################################


###############################################################################

#broken******************
#SUBSEQUENT MANOEUVRES - INCLINATION CHANGE AND PROPAGATION WITH J2 PERTERBATION
#INCLINATION TRANSFER -----**requires correction for transfer to take J2 position inclination impulse 
# r0 = 250 #km
# r1 = 250 #km
# #    [                       a,  e,   i,    ta,   aop, raan]
# c0 = [centralbody['radius']+r0, 0.0, 0.0, 0.0, 30.0, 10.0]
# c1 = [centralbody['radius']+r1, 0.0, 25.0, 0.0, 30.0, 10.0]
# op0, op1, deg_change, DeltaV_Val = my_T.Inclination_Transfer(COEs0=c0, COEs1=c1, propagate=True, add_perts=(True))

# print("Here is the DeltaV of the transfer:\n", DeltaV_Val)
# #PROPAGATE ORBIT OVER 2 DAYS
# tspan = 3600.0*24.0*6.0
# dt=100
# M2 = my_T.Convert_SV_to_COEs(op0.r_vals[-1], op0.v_vals[-1], print_results=False, degrees=(True))
# op = my_T.Propagate_Orbit(M2, tspan, dt, withCOEs=True, add_perts=True)#, perturbations=Default_Settings())

# my_T.Plot_N_Orbits([op0.r_vals, op1.r_vals],
#                     labels=['Orbit Initial', 'Orbit Post Transfer'],
#                     show_plot=True, title = 'Transfer Plot')
###############################################################################


###############################################################################
        
        # #PLOT POSITION AND VELOCITY VECTORS OVER ORBIT (initial)

        # #Plot of v and r for propagated states throughout one complete orbital period 

        # op_Position_ = my_T.norm(op.r_vals, axis=True)
        # op_Velocity_ = my_T.norm(op.v_vals, axis=True)
        
        # #Position and Velocity Plots
        # fig, axs = plt.subplots(nrows=2,ncols=2, figsize = (18, 10))
        # fig.tight_layout(pad=5.0)
        # fig.suptitle("Position and Velocity Propagation [Radau]",fontsize=10)
        
        # axs[0,0].plot(op.t_vals/(3600), op.r_vals)
        # axs[0,0].set_title("Initial Orbit Position Vectors vs. Time")
        # axs[0,0].grid(True)
        # axs[0,0].set_ylabel('Position Vector (km/s)')
        # axs[0,0].set_xlabel('Time (hrs)')
        
        # axs[1,0].plot(op.t_vals/(3600), op_Position_)
        # axs[1,0].set_title("Initial Orbit Position vs. Time")
        # axs[1,0].grid(True)
        # axs[1,0].set_ylabel('Position Vector (km/s)')
        # axs[1,0].set_xlabel('Time (hrs)')
        
        # axs[0,1].plot(op.t_vals/(3600), op.v_vals)
        # axs[0,1].set_title("Initial Orbit Velocity vs. Time")
        # axs[0,1].grid(True)
        # axs[0,1].set_ylabel('Velocity Vector (km/s)')
        # axs[0,1].set_xlabel('Time (hrs)')
        
        # axs[1,1].plot(op.t_vals/(3600), op_Velocity_)
        # axs[1,1].set_title("Initial Orbit Velocity vs. Time")
        # axs[1,1].grid(True)
        # axs[1,1].set_ylabel('Velocity Vector (km/s)')
        # axs[1,1].set_xlabel('Time (hrs)')
        
        # plt.show()



        # #PLOT POSITION AND VELOCITY VECTORS OVER ONE COMPLETE ORBIT (final and initial)
        # r0 = 6570
        # r1 = 42160
        # #    [ a,   e,    i,  ta,  aop, raan]
        # c0 = [r0, 0.0, 16.0, 0.0, 30.0, 10.0]
        # c1 = [r1, 0.0, 16.0, 180.0,30.0, 10.0]
        # op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer(r0=r0, r1=r0, COEs0=c0, COEs1=c1, altitude=False, dt=100, propagate=True)  

        # #Plot of v and r for propagated states throughout one complete orbital period 
        # op0_Position_ = my_T.norm(op0.r_vals, axis=True)
        # op1_Position_ = my_T.norm(op1.r_vals, axis=True)
        # op0_Velocity_ = my_T.norm(op0.v_vals, axis=True)
        # op1_Velocity_ = my_T.norm(op1.v_vals, axis=True)
        # #op_Position_ = my_T.norm(op.r_vals, axis=True)
        # #op_Velocity_ = my_T.norm(op.v_vals, axis=True)
        
        # op0_Position = np.sqrt(np.sum(op0.r_vals**2, axis=1))
        # op1_Position = np.sqrt(np.sum(op1.r_vals**2, axis=1))
        # op0_Velocity = np.sqrt(np.sum(op0.v_vals**2, axis=1))
        # op1_Velocity = np.sqrt(np.sum(op1.v_vals**2, axis=1))

        # #Position and Velocity Plots
        # fig, axs = plt.subplots(nrows=2,ncols=4, figsize = (18, 10))
        # fig.tight_layout(pad=5.0)
        # fig.suptitle("Position and Velocity Propagation",fontsize=10)
    
        # axs[0,0].plot(op0.t_vals/(3600), op0.r_vals)
        # axs[0,0].set_title("Initial Orbit Position Vectors vs. Time")
        # axs[0,0].grid(True)
        # axs[0,0].set_ylabel('Position Vector (km/s)')
        # axs[0,0].set_xlabel('Time (hrs)')
        
        # axs[1,0].plot(op0.t_vals/(3600), op0_Position_)
        # axs[1,0].set_title("Initial Orbit Position vs. Time")
        # axs[1,0].grid(True)
        # axs[1,0].set_ylabel('Position Vector (km/s)')
        # axs[1,0].set_xlabel('Time (hrs)')

        # axs[0,1].plot(op0.t_vals/(3600), op0.v_vals)
        # axs[0,1].set_title("Initial Orbit Velocity vs. Time")
        # axs[0,1].grid(True)
        # axs[0,1].set_ylabel('Velocity Vector (km/s)')
        # axs[0,1].set_xlabel('Time (hrs)')

        # axs[1,1].plot(op0.t_vals/(3600), op0_Velocity)
        # axs[1,1].set_title("Initial Orbit Velocity vs. Time")
        # axs[1,1].grid(True)
        # axs[1,1].set_ylabel('Velocity Vector (km/s)')
        # axs[1,1].set_xlabel('Time (hrs)')
        
        # axs[0,2].plot(op1.t_vals/(3600), op1.r_vals)
        # axs[0,2].set_title("Final Orbit Position Vectors vs. Time")
        # axs[0,2].grid(True)
        # axs[0,2].set_ylabel('Position Vector (km/s)')
        # axs[0,2].set_xlabel('Time (hrs)')
        
        # axs[1,2].plot(op1.t_vals/(3600), op1_Position_)
        # axs[1,2].set_title("Final Orbit Position vs. Time")
        # axs[1,2].grid(True)
        # axs[1,2].set_ylabel('Position Vector (km/s)')
        # axs[1,2].set_xlabel('Time (hrs)')

        # axs[0,3].plot(op1.t_vals/(3600), op1.v_vals)
        # axs[0,3].set_title("Final Orbit Velocity vs. Time")
        # axs[0,3].grid(True)
        # axs[0,3].set_ylabel('Velocity Vector (km/s)')
        # axs[0,3].set_xlabel('Time (hrs)')

        # axs[1,3].plot(op1.t_vals/(3600), op1_Velocity)
        # axs[1,3].set_title("Final Orbit Velocity vs. Time")
        # axs[1,3].grid(True)
        # axs[1,3].set_ylabel('Velocity Vector (km/s)')
        # axs[1,3].set_xlabel('Time (hrs)')
        
        # plt.show()
###############################################################################
        




#CALCULATING RESULTS

# L1 = [ 1975.092981,     5755.050166,     6994.405950,
#  1555.918924,     5709.751529,     6988.984694,
#  12553.604192,     2423.263859,     1844.133241]

# L2 = [ 1975.0914033,   5755.04984318,  6994.40666112,
# 1555.85979119,  5709.74427148,  6988.98367143,
# 12656.14403594 , 2294.81356652 , 1672.47615788]

# L3 = []

# print(L2[0])

# for i in range(9):
#     print(i)
#     L3.append(L1[i] - L2[i])
    
# print(L3)
