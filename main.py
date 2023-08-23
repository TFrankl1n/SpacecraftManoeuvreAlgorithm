import os
import sys
import csv
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
    
    def __init__(self, initial_elements, dt, withCOEs=None, add_perts=True, degrees=True):
        # Initialise the Satellite class with provided parameters
        self.initial_elements = initial_elements
        self.dt = dt
        self.withCOEs = withCOEs
        self.add_perts = add_perts
        self.degrees = degrees
        self.segments = []                                                     # List to store final state values of previous propagation segment
        self.DeltaV_Sat_Object = []                                            # List to store DeltaV values of manoeuvre series

    def propagate(self, tspan, calculate_COEs, full_plot, grid_plot):
        self.tspan = tspan
        for orbit_elements in self.initial_elements:                           # Check inputs
            if initial_elements[1] <0.0 or initial_elements[1] >1.0:
                raise ValueError("Eccentricity by definition must be between 0 and 1.")
            if self.degrees:                                                   # Check if input values are in degrees correctly, defining bounds for true anomaly, argument of perigee, and right asc of the ascending node
                if initial_elements[2] <0.0 or initial_elements[2] >180.0:
                    raise ValueError("Inclination in degrees must be between 0 and 180.\n (Prograde orbits between 0 and 90 degrees, retrograde beyond this)")
                if  initial_elements[3] <0.0 or initial_elements[3] >360.0:
                    raise ValueError("True Anomaly in degrees must be between 0 and 360.")
                if  initial_elements[4] <0.0 or initial_elements[4] >360.0:
                    raise ValueError("Argument of Perigee in degrees must be between 0 and 360.")
                if  initial_elements[5] <0.0 or initial_elements[5] >360.0:
                    raise ValueError("Right Asc. of the Ascending Node in degrees must be between 0 and 360 .)")
            else:                                                              # Check if any of the initial elements are beyond the bounds for radians             
                if any(np.abs(initial_elements[2:6]) >= 2 * np.pi):
                    raise ValueError("Initial elements may be in radians or one of the initial conditions may be beyond bounds for radians.\nCurrently taking input in radians, set 'degrees=True' for degrees.")        
            if not calculate_COEs:                                             # Check if any plots have been called without recording arrays to plot them with
                if full_plot:
                    raise ValueError("Cannot plot the COEs over propagation duration without 'calculate_COEs=True'.")
                if grid_plot:
                    raise ValueError("Cannot plot the COEs over propagation duration without 'calculate_COEs=True'.")
        try:        
            if self.time_moved > 0:                                            # Check if propagation has already occured so we know whether to use [initial conditions + input values] or [previous state and state vecors]
                print("Beginning next propagation...")
                initial_segment = self.segments[-1]
                initial_segment = np.concatenate((initial_segment.r_vals[-1], initial_segment.v_vals[-1]))

                op = my_T.Propagate_Orbit(                                                  # Propagate the orbit using last segment elements and settings
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
                self.time_moved +=1                                                         # increment time_moved to only use segment values with later propagation

            else:
                
                print("Beginning first propagation...")
                op = my_T.Propagate_Orbit(                                                  # Propagate the orbit using provided initial elements and settings         
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
                self.time_moved +=1                                                         # increment time_moved to only use segment values with later propagation        
            self.segments.append(op)                                                        # Append final state of propagation step to segemnts array

        except:
            raise ValueError("Check input orbital elements type. If using classical orbital elements (or TLE input), set 'withCOEs=True', if using state vectors, set 'withCOES=False")
        
    def perform_hohmann_transfer_optimal(self, smachange, dt): 
        
        if self.time_moved > 0:                                                #check if propagation has already occured, to use [initial conditions + input values] or [previous state and state vecors]
            print("Beginning Hohmann Transfer propagation...")
            a_new = self.initial_elements[0]+smachange
            
            # Propagate to the end of the initial segment
            initial_segment = self.segments[-1]
            
            Hohmann_transfer_initial_COEs = my_T.Convert_SV_to_COEs(initial_segment.r_vals[-1], initial_segment.v_vals[-1])                                                                  # Calculate initial COEs for transfer orbit
            Hohmann_transfer_final_COEs = [a_new if i == 0 else self.initial_elements[i] if i==1 else x + 180.0 if i == 3 else x for i, x in enumerate(Hohmann_transfer_initial_COEs)]       # Calculate final COEs for transfer orbit
            
            #Propagate transfer with above conditions
            op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer_ideal(
                                                                                                    self.initial_elements[0], 
                                                                                                    a_new,
                                                                                                    COEs0=Hohmann_transfer_initial_COEs, 
                                                                                                    COEs1=Hohmann_transfer_final_COEs,
                                                                                                    dt=100, 
                                                                                                    propagate=True, 
                                                                                                    add_perts=True
                                                                                                )
            self.time_moved +=1                                                # increment time_moved to only use segment values with later propagation

        else:                                                                  # No prior prop completed 
            print("Beginning Hohmann Transfer as first propagation...")
            a_new = self.initial_elements[0]+smachange
            
            #calculate initial COEs for transfer orbit
            Hohmann_transfer_initial_COEs = my_T.Convert_COEs_deg_to_rad(self.initial_elements)

            # Calculate final COEs for transfer orbit
            Hohmann_transfer_final_COEs = [a_new if i == 0 else self.initial_elements[i] if i==1 else x + 180.0 if i == 3 else x for i, x in enumerate(Hohmann_transfer_initial_COEs)]
            # Hohmann_transfer_final_COEs = my_T.Convert_COEs_deg_to_rad(Hohmann_transfer_final_COEs)
            
            #Propagate transfer with above conditions
            op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer_ideal(
                                                                                                    self.initial_elements[0], 
                                                                                                    a_new,
                                                                                                    COEs0=Hohmann_transfer_initial_COEs, 
                                                                                                    COEs1=Hohmann_transfer_final_COEs,
                                                                                                    dt=100, 
                                                                                                    propagate=True, 
                                                                                                    add_perts=True
                                                                                                )  
            self.time_moved +=1                                                 # increment time_moved to only use segment values with later propagation
        
        print('DeltaV Burn 1: \t %.3f (km/s)' % DeltaV_Vals[0])                 # Display Calculated DeltaV and time of flight values
        print('DeltaV Burn 2: \t %.3f (km/s)' % DeltaV_Vals[1])
        print('Total DeltaV: ', DeltaV_Vals[1] + DeltaV_Vals[0], "(km/s)")
        print('Eccentricity of Transfer: ', e_transfer)
        print('Transfer Time: \t',tof_transfer,' seconds,\nTransfer Time: \t ',tof_transfer/3600, "hours")
        
        self.segments.append(op0)                                               # Append initial segment to segments list (one complete orbit)
        self.segments.append(op_transfer)                                       # Append transfer segment to segments list
        self.segments.append(op1)                                               # Append (post second burn) segment to segments list (one complete orbit)
        self.time_moved +=1      
        pass

    def perform_Phasing(self, smachange, dt):                                                                                                           #[Further DeltaV calculation investigation required ]
        
        
        final_segment = self.segments[-1]                                       # Extract the final state of the last propagated segment
        last_segment_r_final = final_segment.r_vals[-1]
        last_segment_v_final = final_segment.v_vals[-1]
        # Calculate initial COEs for transfer orbit
        Phasing_transfer_initial_COEs = my_T.Convert_SV_to_COEs(last_segment_r_final, last_segment_v_final)
        a_new = self.initial_elements[0] + smachange

        # Propagate transfer with above conditions
        op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer(
                                                                                                Phasing_transfer_initial_COEs[0], a_new,
                                                                                                COEs0=Phasing_transfer_initial_COEs,
                                                                                                COEs1=Phasing_transfer_initial_COEs,
                                                                                                dt=dt,
                                                                                                propagate=True,
                                                                                                add_perts=True
                                                                                            )                                      

        print('Eccentricity of Transfer: ', e_transfer)
        print('Transfer Time: \t',tof_transfer,' seconds, ',tof_transfer/3600, "hours")
        self.segments.append(op_transfer)                                       # Append transfer segment to segments list
        self.time_moved +=1
    
        pass
    def perform_inclination_change(self, delta_i):
        
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
        
        # Calculate and print DeltaV
        DeltaV_Val = np.linalg.norm(op.v_vals[-1] - v_final)
        print("DeltaV:", DeltaV_Val, "(km/s)")
        pass


    def plot_orbits(self, show_plot=True, title="My Custom Plot Title", dark_theme=False):
        if show_plot:
            r_vals_list = [segment.r_vals for segment in self.segments]
            labels = [f"Segment {i}" for i in range(1, len(self.segments) + 1)]
            my_T.Plot_N_Orbits(r_vals_list, labels=labels, show_plot=True, title="Hohmann Transfer Propagated Orbit", dark_theme=dark_theme)


    def save_to_csv(self, filename):
        try:
            with open(filename, mode='w', newline='') as file_csv:
               write2csv = csv.writer(file_csv)
               write2csv.writerow(["Time (s)", "x (km)", "y (km)", "z (km)", "vx (km/sec)", "vy (km/sec)", "vz (km/sec)"])  # Write export file headers
    
               for segment in self.segments:
                   for i in range(len(segment.t_vals)):
                       t_vals = segment.t_vals[i]
                       r_vals = segment.r_vals[i]
                       v_vals = segment.v_vals[i]
                       write2csv.writerow([t_vals[0], r_vals[0], r_vals[1], r_vals[2], v_vals[0], v_vals[1], v_vals[2]])
        except Exception as excep:
            print(f"Error saving data to {filename}: {excep}")

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
####Blank initial required variables array:
# initial_elements = [    ,   ,   ,   ,   ,   ]                                 # Enter Object initial conditions in either COE or SV form
# dt =                                                                          # Enter time step of propagation (lower values will produce higher fidelity plot, and increase runtime)

####Inport TLE data
# initial_elements = my_T.TLE_to_CEOsPlus(         ),                           # Enter File_Path_Name to propagate from 
# dt =   
#### Add Satellite Object
#    = Satellite(                                                               # Enter Name of Satellite Object
            # initial_elements, 
            # dt, 
            # withCOEs=True, 
            # add_perts=True, 
            # degrees=True
            #     )
            

################################### Manoeuvres:

#### Perform Hohmann transfer
# smachange =                                                                   # Enter desired change in semimajor axis increase or decrease
# [object name] .perform_hohmann_transfer_optimal(smachange, dt) 

#### Perform inclination change
# delta_i =                                                                     # Enter desired change in inclination 
# [object name] .perform_inclination_change(delta_i)

#### Perform Phasing Manoeuvre
# sma_change =                                                                  # Enter difference in sma for Phasing Manoeuvre                              
# [Object name] .perform_Phasing(sma_change, dt)

#### Perform orbit propagation
#tspan =                                                                        # Enter desired propagation duration 
# [object name] .propagate(
#     tspan = tspan
#     calculate_COEs=False, 
#     full_plot=False, 
#     grid_plot=False, 
#               )




################################### Tools:

    
#### Plot Orbit
# [Object name] .plot_orbits()

### Save the data to a CSV file
# filename = "orbit_data.csv"
# file_path = os.path.join(script_dir, filename)
# [Object name ] .save_to_csv(file_path)

### Print save directory
# print(f"Script Directory: {script_dir}")
# print(f"CSV File Path: {file_path}")
# print(f"Data saved to {filename}")
        


"""
                                                                                                Analyzes an orbit using the given parameters.
                                                                                            
                                                                                                Args:
                                                                                                    initial_elements    (list): List of initial elements either:
                                                                                                                                        classical orbital elements [a, e, i, ta, aop, raan].
                                                                                                                                        state vectors [rx, ry, rz, vx, vy, vz].
                                                                                                    tspan               (float): Time span for propagation.
                                                                                                    dt                  (float): Time step for propagation.
                                                                                                    withCOEs            (bool): Whether the initial_elements are classical orbital elements (=True) or state vectors (=False)
                                                                                                    add_perts           (bool): Whether to include J2 perturbative effects.                                                                             **********************************************************************add options
                                                                                                    degrees             (bool): Whether the input elements are in degrees (=True), radians (=False).
                                                                                                    calculate_COEs      (bool): Whether to calculate COEs at each step of propagation:
                                                                                                        full_plot           (bool): Whether to plot all classical orbital elements on single chart against time
                                                                                                        grid_plot           (bool): Whether to plot all classical orbital elements seperately against time
                                                                                                    plot_3d             (bool): Whether to plot the orbit trajectory in 3D.
                                                                                            
                                                                                                Returns:
                                                                                                    op: The result of orbit propagation.
"""
############################################################################### Example sequence:



initial_elements = [15500.0, 0.1, 15.0, 15.0, 30.0, 10.0]                                 # Enter Object initial conditions in either COE or SV form
dt =  60.0                                                                        # Enter time step of propagation (lower values will produce higher fidelity plot, and increase runtime)


### Add Satellite Object
first_Sat = Satellite(                                                               # Enter Name of Satellite Object
            my_T.TLE_to_CEOsPlus(file_path_ISS),                   # Enter File_Path_Name to propagate from   , 
            dt, 
            withCOEs=True, 
            add_perts=True, 
            degrees=True
                )


### Perform orbit propagation
tspan = 3600.0 * 24.0 * 5.0                                                        # Enter desired propagation duration 
first_Sat.propagate(
    tspan=tspan,
    calculate_COEs=False, 
    full_plot=False, 
    grid_plot=False, 
              )

### Perform Hohmann transfer
smachange =  16500.0                                                              # Enter desired change in semimajor axis increase or decrease
first_Sat.perform_hohmann_transfer_optimal(smachange, dt) 


### Perform inclination change
delta_i =  12.5                                                                   # Enter desired change in inclination 
first_Sat.perform_inclination_change(delta_i)

tspan =  3600.0 * 24.0 * 5.0                                                                  # Enter desired propagation duration 
first_Sat.propagate(
    tspan=tspan,
    calculate_COEs=False, 
    full_plot=False, 
    grid_plot=False, 
              )

### Perform Phasing Manoeuvre
sma_change =  8000.0                                                               # Enter difference in sma for Phasing Manoeuvre                              
first_Sat.perform_Phasing(sma_change, dt)

### Perform orbit propagation
tspan =  3600.0 * 24.0 * 2.0                                                                  # Enter desired propagation duration 
first_Sat.propagate(
    tspan=tspan,
    calculate_COEs=False, 
    full_plot=False, 
    grid_plot=False, 
              )
    
### Plot Orbit
first_Sat.plot_orbits()


























        
    



# ____________________________________________________________________________________________________________________________________________________________________________ #
############################################################################### DEV ENVIRONMENT: ###############################################################################

#Immediate Hohmann Transfer:
    # def perform_hohmann_transfer(self, smachange, dt):
    #     if smachange>0:
    #         a_new = self.initial_elements[0] + smachange
    #     elif smachange<0:
    #         a_new = self.initial_elements[0] - smachange
    #     else:
    #         raise ValueError("If performing altitude change, the semi major axis requires a non-zero value.")
        
    #     final_segment = self.segments[-1]                                       # Extract the final state of the last propagated segment
    #     last_segment_r_final = final_segment.r_vals[-1]
    #     last_segment_v_final = final_segment.v_vals[-1]
    #     # Calculate initial COEs for transfer orbit
    #     Hohmann_transfer_initial_COEs = my_T.Convert_SV_to_COEs(last_segment_r_final, last_segment_v_final)
    #     # Calculate final COEs for transfer orbit
    #     Hohmann_transfer_final_COEs = [a_new if i == 0 else x + 180.0 if i == 3 else x for i, x in enumerate(Hohmann_transfer_initial_COEs)]
    #     # Propagate transfer with above conditions
    #     op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer(
    #                                                                                             Hohmann_transfer_initial_COEs[0], a_new,
    #                                                                                             COEs0=Hohmann_transfer_initial_COEs,
    #                                                                                             COEs1=Hohmann_transfer_final_COEs,
    #                                                                                             dt=dt,
    #                                                                                             propagate=True,
    #                                                                                             add_perts=True#,
    #                                                                                             #degrees = True
    #                                                                                         )                                      
    #     print('DeltaV 0: \t %.3f (km/s)' % DeltaV_Vals[0])
    #     print('DeltaV 1: \t %.3f (km/s)' % DeltaV_Vals[1])
    #     print('Total DeltaV: ', DeltaV_Vals[1] + DeltaV_Vals[0])
    #     print('Eccentricity of Transfer: ', e_transfer)
    #     print('Transfer Time: \t',tof_transfer,' seconds, ',tof_transfer/3600, "hours")
    #     self.segments.append(op0)                                              # Append initial segment to segments list (one complete orbit)
    #     self.segments.append(op_transfer)                                      # Append transfer segment to segments list
    #     self.segments.append(op1)                                              # Append (post second burn) segment to segments list (one complete orbit)
    #     self.time_moved +=1

    #     pass

                                                                               ################Errors:
                                                                                   #Calculation of total DeltaV due to nested array storage of velocity componenets
                                                                                   # matching final(target) orbit and end of transfer path




###############################################################################
        
   
       
       



########################################################################################################################## [Run without above class]
###################################################################################################Function to loop through initial states, and at each position complete a transfer, storing the DeltaV values to observe how tof, DV, etc. changes at different transfer start times for one initial state
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








##############################################################################################                     [To be implemented as flags in Satellite class to print velocity and position vectors over time]

#PLOT POSITION AND VELOCITY VECTORS OVER ORBIT (initial)

# #Plot of v and r for propagated states throughout one complete orbital period 

# op_Position_ = my_T.norm(op0.r_vals, axis=True)
# op_Velocity_ = my_T.norm(op0.v_vals, axis=True)

# #Position and Velocity Plots
# fig, axs = plt.subplots(nrows=2,ncols=2, figsize = (18, 10))
# fig.tight_layout(pad=5.0)
# fig.suptitle("Position and Velocity Propagation [Radau]",fontsize=10)

# axs[0,0].plot(op0.t_vals/(3600), op0.r_vals)
# axs[0,0].set_title("Initial Orbit Position Vectors vs. Time")
# axs[0,0].grid(True)
# axs[0,0].set_ylabel('Position Vector (km/s)')
# axs[0,0].set_xlabel('Time (hrs)')

# axs[1,0].plot(op0.t_vals/(3600), op_Position_)
# axs[1,0].set_title("Initial Orbit Position vs. Time")
# axs[1,0].grid(True)
# axs[1,0].set_ylabel('Position Vector (km/s)')
# axs[1,0].set_xlabel('Time (hrs)')

# axs[0,1].plot(op0.t_vals/(3600), op0.v_vals)
# axs[0,1].set_title("Initial Orbit Velocity vs. Time")
# axs[0,1].grid(True)
# axs[0,1].set_ylabel('Velocity Vector (km/s)')
# axs[0,1].set_xlabel('Time (hrs)')

# axs[1,1].plot(op0.t_vals/(3600), op_Velocity_)
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











#Code additional ideas:
    
    #TODO
        # Perfect Hohmann transfer 
            # Identify perfect Periapsis from orbit prop
        
        # Inclination change manoeuvre
            # Plotting manoeuvre from optimal position (at zaxis=0 (Detect AN an DN))
            #subsequent plotting of manoeuvre from prev. conditions from AN or DN
            # subsequent manoeuvre from any position.
            # track or input desired RAAN change to induce inclination change
        # Eccentricity change version of Phasing manoeuvre
        # Fast Hohmann (with immediate burn after new altitude reached. - using stop conditions to identify when altitude = target)
        

        #USER MANUAL:
            #customising graph outputs????
            #graphs of position and velocity vectors automated
            #loop function for plotting delta v vs tof. (or vs when manoeuvre begins?)
            #TODO***Fix ideal inclination change
            #TODO**Sort immediate inclination change!!
            #HOW  to turn on and off stop conditions
            # How to sleect a certain perturbative effect [adding to perts=True] init CHANGE PROPAGATOR CENTRAL BODY PERTURBATIONS -
            
        #DESCRIPTIONS:
            #Diff eq:
                # PERTURBATION INPUTS AT TOP LEVEL
                # AERO DRAG CALC + INPUTS
                # RADIATION PRESSURE + INPUTS