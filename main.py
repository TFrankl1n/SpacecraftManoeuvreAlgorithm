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
        self.DeltaV_Sat_Object = []

    def propagate(self, degrees, calculate_COEs, full_plot, grid_plot):
        for orbit_elements in self.initial_elements:                            #Check inputs
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
        print("End of previous segment", self.segments[-1])

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
        last_segment_r_final = final_segment.r_vals[-1]
        last_segment_v_final = final_segment.v_vals[-1]
        #print("Transfer start coords:", last_segment_r_final)
        # Calculate initial COEs for transfer orbit
        Hohmann_transfer_initial_COEs = my_T.Convert_SV_to_COEs(last_segment_r_final, last_segment_v_final)
        # print("[main] Hohmann Transfer start array:\n",Hohmann_transfer_initial_COEs)
        # print(last_segment_r_final, last_segment_v_final)
        # Calculate final COEs for transfer orbit
        Hohmann_transfer_final_COEs = [a_new if i == 0 else x + 180.0 if i == 3 else x for i, x in enumerate(Hohmann_transfer_initial_COEs)]
        # print("[main] Hohmann Transfer end array:\n",Hohmann_transfer_final_COEs)
        # print(my_T.Convert_COEs_to_SV(Hohmann_transfer_final_COEs))
        # Propagate transfer with above conditions
        op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer(
                                                                                                Hohmann_transfer_initial_COEs[0], a_new,
                                                                                                COEs0=Hohmann_transfer_initial_COEs,
                                                                                                COEs1=Hohmann_transfer_final_COEs,
                                                                                                dt=100,
                                                                                                propagate=True,
                                                                                                add_perts=True
                                                                                            )                                      
        # print('13DeltaV 0: \t %.3f (km/s)' % DeltaV_Vals[0])
        # print('13DeltaV 1: \t %.3f (km/s)' % DeltaV_Vals[1])
        # print('13Total DeltaV: ', DeltaV_Vals[1] + DeltaV_Vals[0])
        # print('13Eccentricity of Transfer: ', e_transfer)
        # print('13Transfer Time: \t',tof_transfer,' seconds, ',tof_transfer/3600, "hours")
        self.segments.append(op0)                                              # Append initial segment to segments list (one complete orbit)
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

# print(file_path_ISS)                                                            #Check any file path to ensure correctly denoted in software for calling

# script_dir = os.path.dirname(os.path.abspath(__file__))                         #Check directory of python script 
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
initial_elements = [10563, 0.2, 25.0, 50.0, 30.0, 10.0]
dt = 60.0

tspan = 2000.0
sat_1 = Satellite(initial_elements, 
                  tspan, 
                  dt, 
                  withCOEs=True,
                  add_perts=True, 
                  degrees=True
                  )

# Perform propagation
sat_1.propagate(
                degrees=True, 
                calculate_COEs=False, 
                full_plot=False, 
                grid_plot=False, 
                )
#Immediate Hohmann Transfer
smachange = 8754
sat_1.perform_hohmann_transfer(smachange) 

initial_elements = [10563, 0.2, 25.0, 50.0, 30.0, 10.0]
dt = 60.0

tspan = 9900.0
sat_2 = Satellite(initial_elements, 
                  tspan, 
                  dt, 
                  withCOEs=True,
                  add_perts=True, 
                  degrees=True
                  )
# Perform propagation
sat_2.propagate(
                degrees=True, 
                calculate_COEs=False, 
                full_plot=False, 
                grid_plot=False, 
                )
#Immediate Hohmann Transfer
smachange = 8754
sat_2.perform_hohmann_transfer(smachange) 

# Plot all segments
sat_1.plot_orbits()
sat_2.plot_orbits()
