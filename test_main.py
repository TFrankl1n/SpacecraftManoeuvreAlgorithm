import numpy as np
from math import sqrt                                                        
import matplotlib.pyplot as plt

import os
import sys     
current_directory = os.getcwd()
files_folder = os.path.join(current_directory, 'Files')
filename_1 = 'TLE_1_ISS.txt'
file_path = os.path.join(files_folder, filename_1)

import planetary_data as my_PD
import tools as my_T
from orbit_propagator import Default_Settings

centralbody = my_PD.earth

tspan = 3600 * 24.0 *1.0
dt = 100.0

if __name__ == '__main__':    
    print('Running...')
    Initiaite_Propagation = True
    
    if Initiaite_Propagation:
        print("Main Running Successfully:")
           
###############################################################################
        # #SINGLE ORBIT PATH FROM COEs
        # c0 = [centralbody['radius']+r0, 0.0, 28.0, 15.0, 7.0, 6.0]     
        # op = my_T.Propagation_Orbit(c0, tspan, dt, add_perts=True)#, perturbations=Default_Settings())

        # my_T.plot_3d(op.r_vals, show_plot=(True), title="Single COE Plot")
###############################################################################


###############################################################################
        # #MULTIPLE ORBIT PROPAGATION AND PLOT FROM COEs WITH PERTURBATION
        # op0 = my_T.Propagation_Orbit(c0, tspan, dt, add_perts=True)#, perturbations=Default_Settings())
        # op1 = my_T.Propagation_Orbit(c1, tspan, dt, add_perts=True)#, perturbations=Default_Settings())

        # my_T.Plot_N_Orbits([op0.r_vals, op1.r_vals], labels=["Orbit 1", "Orbit 2"], show_plot=True, title = 'Multiple TLE Plot')#' with J2 Perturbation over 2 Days' )
###############################################################################


###############################################################################
        # #SINGLE ORBIT PROP (From TLEs)
        # print("Calculating opop")
        # opop = my_T.Propagation_Orbit(my_T.TLE_to_CEOsPlus('F:/Users/Thomas/UNIVERSITY/Work [] Placements/WORK/5. AIRBUS/1. Work/Manoeuvres/Astro_Dynamics/8_J2/TLE_1_ISS.txt'),
        #                               tspan, dt)#, COEs = True, deg=False, add_perts=True)
        
        # my_T.plot_3d(opop.r_vals, show_plot=(True), title="ISS single TLE Plot")
        
        # #IMPORT WITH file_path
        
        # opop1 = my_T.Propagation_Orbit(my_T.TLE_to_CEOsPlus(file_path),tspan, dt)#, COEs = True, deg=False, add_perts=True)
        
        # my_T.plot_3d(opop1.r_vals, show_plot=(True), title="ISS single TLE Plot")
###############################################################################


###############################################################################
        #  #MULTIPLE ORBIT PROP (From TLEs)
        # op0 = my_T.Propagation_Orbit(my_T.TLE_to_CEOsPlus('F:/Users/Thomas/UNIVERSITY/Work [] Placements/WORK/5. AIRBUS/1. Work/Manoeuvres/Astro_Dynamics/8_J2/TLE_1_ISS.txt'), 
        #              tspan, dt, add_perts=True)
        # op1 = my_T.Propagation_Orbit(my_T.TLE_to_CEOsPlus('F:/Users/Thomas/UNIVERSITY/Work [] Placements/WORK/5. AIRBUS/1. Work/Manoeuvres/Astro_Dynamics/8_J2/TLE_2_Debris.txt'), 
        #              tspan, dt, add_perts=True)#, COEs = True, deg=False, add_perts=True)
        # op2 = my_T.Propagation_Orbit(my_T.TLE_to_CEOsPlus('F:/Users/Thomas/UNIVERSITY/Work [] Placements/WORK/5. AIRBUS/1. Work/Manoeuvres/Astro_Dynamics/8_J2/TLE_3_Starlink-30095.txt'), 
        #              tspan, dt, add_perts=True)#, COEs = True, deg=False, add_perts=True)
       
        # #op0 = my_T.Initiate_Propagation(c0, tspan, dt, COEs=True)
        # #my_T.Plot_N_Orbits([op0.r_vals, op1.r_vals], show_plot=True, title = 'Multiple TLE Plot')#' with J2 Perturbation over 2 Days' )
        # tspan_days = tspan/(3600 * 24.0)
        # my_T.Plot_N_Orbits([op0.r_vals,op1.r_vals,op2.r_vals], labels=['ISS','Debris','Starlink-30095'],show_plot=True,
        #                    title = f"Multiple TLE Plot with J2 Perturbation over {tspan_days:.1f} Day(s)")
###############################################################################


###############################################################################        
        # #SINGLE HOHMANN TRANSFER - [Value calc (a, v, tof, e)], []
        # r0 = 6570
        # r1 = 42160
        # #[    a,  e,   i,    ta,   aop, raan, [year, month, day, hour]]
        # c0 = [r0, 0.0, 16.0, 0.0, 30.0, 10.0]
        # c1 = [r1, 0.0, 16.0, 180.0,30.0, 10.0]
        # op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer_Circ(r0=r0, r1=r0, COEs0=c0, COEs1=c1, altitude=False, dt=100, propagate=True)  
        # print('DeltaV 0: \t %.3f (km/s)' % DeltaV_Vals[0])
        # print('DeltaV 1: \t %.3f (km/s)' % DeltaV_Vals[1])
        # print('Total DeltaV: ', DeltaV_Vals[1] + DeltaV_Vals[0])
        # print('Eccentricity of Transfer: ', e_transfer)
        # print('Transfer Time: \t',tof_transfer,' seconds, ',tof_transfer/3600, "hours")
        
        # my_T.Plot_N_Orbits([op0.r_vals, op1.r_vals, op_transfer.r_vals],labels=['Orbit 1 Initial', 'Orbit 1 Final', 'Orbit 1 Transfer'],show_plot=False, title = 'Transfer Plot')
##############################################################################


###############################################################################
        # #PLOT POSITION AND VELOCITY VECTORS OVER ONE COMPLETE ORBIT
        # r0 = 6570
        # r1 = 42160
        # #    [ a,   e,    i,  ta,  aop, raan]
        # c0 = [r0, 0.0, 16.0, 0.0, 30.0, 10.0]
        # c1 = [r1, 0.0, 16.0, 180.0,30.0, 10.0]
        # op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer_Circ(r0=r0, r1=r0, COEs0=c0, COEs1=c1, altitude=False, dt=100, propagate=True)  

        # #Plot of v and r for propagated states throughout one complete orbital period 
        # op0_Position_ = my_T.norm(op0.r_vals, axis=True)
        # op1_Position_ = my_T.norm(op1.r_vals, axis=True)
        # op0_Velocity_ = my_T.norm(op0.v_vals, axis=True)
        # op1_Velocity_ = my_T.norm(op1.v_vals, axis=True)
        
        # op0_Position = np.sqrt(np.sum(op0.r_vals**2, axis=1))
        # op1_Position = np.sqrt(np.sum(op1.r_vals**2, axis=1))
        # op0_Velocity = np.sqrt(np.sum(op0.v_vals**2, axis=1))
        # op1_Velocity = np.sqrt(np.sum(op1.v_vals**2, axis=1))

        # #Position and Velocity Plots
        # fig, axs = plt.subplots(nrows=2,ncols=4, figsize = (18, 10))
        # fig.tight_layout(pad=5.0)
        # fig.suptitle("Position and Velocity Propagation [Radau]",fontsize=10)
    
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
        

###############################################################################
        # PLOTTING TWO  HOHMANN TRANSFERS, ONE PLANAR, ONE WITH INCLUNATION
        
        # r0 = 210
        # r1 = 845
        # r2 = 2887
        # r3 = 6465
        # c2 = [centralbody['radius']+r0, 0.0, 0.0, 0.0,  0.0, 0.0]
        # c3 = [centralbody['radius']+r1, 0.0, 0.0, 180.0,0.0, 0.0]
        # c2_rot = [centralbody['radius']+r2, 0.0, 30.0,   0.0, 100.0, 20.0]
        # c3_rot = [centralbody['radius']+r3, 0.0, 30.0, 180.0, 100.0, 20.0]
        
        # #op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer
        # op2, op3, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer_Circ(COEs0 = c2, COEs1=c3, propagate = True, add_perts=(True))  
        # # op2, op3, op_transfer, DeltaV_Vals = my_T.Hohmann_Transfer_Circ(propagate = True)  
        
        # op2_rot, op3_rot, op_transfer_rot, DeltaV_Vals_rot, tof_transfer, e_transfer = my_T.Hohmann_Transfer_Circ(COEs0 = c2_rot, COEs1=c3_rot, propagate = True)
        # # op2_rot, op3_rot, op_transfer_rot, DeltaV_Vals_rot = my_T.Hohmann_Transfer_Circ(propagate = True)
        
        # print("op2",len(op2.r_vals))
        # print("op3",len(op3.r_vals))
        # print("op2_rot",len(op2_rot.r_vals))
        # print("op3_rot",len(op3_rot.r_vals))
        # my_T.Plot_N_Orbits([op2.r_vals, op3.r_vals, op_transfer.r_vals, op2_rot.r_vals, op3_rot.r_vals, op_transfer_rot.r_vals], 
        #                    labels=['Orbit 1 Initial', 'Orbit 1 Final', 'Orbit 1 Transfer', 'Orbit 2 Initial', 'Orbit 2 Final', 'Orbit 2 Transfer'],
        #                     show_plot=True, title = 'Transfer Plot')
                            
        # print("Here are the DeltaV Values:\n", DeltaV_Vals)
        # print("And the rotated ones why not:\n", DeltaV_Vals_rot)
###############################################################################


###############################################################################
        # PLOTTING INCLINATION TRANSFER
        # op0, op1, deg_change, DeltaV_Val = my_T.Inclination_Transfer(COEs0=c0, COEs1=c1, propagate=True)
        
        # my_T.Plot_N_Orbits([op0.r_vals, op1.r_vals],
        #                     labels=['Orbit 1 Initial', 'Orbit 1 Final', 'Orbit 1 Transfer', 'Orbit 2 Initial', 'Orbit 2 Final', 'Orbit 2 Transfer'],
        #                     show_plot=True, title = 'Transfer Plot')
                            
        # print("Here is the DeltaV of the transfer:\n", DeltaV_Val)
###############################################################################


###############################################################################
        # #SUBSEQUENT MANOEUVRES - INCLINATION CHANGE AND PROPAGATION WITH J2 PERTERBATION
        #INCLINATION TRANSFER -----**requires correction for transfer to take J2 position inclination impulse 
        r0 = 250 #km
        r1 = 250 #km
        #    [                       a,  e,   i,    ta,   aop, raan]
        c0 = [centralbody['radius']+r0, 0.0, 130.0, 0.0, 30.0, 10.0]
        c1 = [centralbody['radius']+r1, 0.0, 120.0, 0.0, 30.0, 10.0]
        op0, op1, deg_change, DeltaV_Val = my_T.Inclination_Transfer(COEs0=c0, COEs1=c1, propagate=True, add_perts=(True))

        print("Here is the DeltaV of the transfer:\n", DeltaV_Val)
        #PROPAGATE ORBIT OVER 2 DAYS
        tspan = 3600.0*24.0*2.0
        dt=100
        M2 = my_T.Convert_RV_to_COEs(op0.r_vals[-1], op0.v_vals[-1], print_results=False, degrees=(False))
        op = my_T.Propagation_Orbit(M2, tspan, dt, add_perts=True)#, perturbations=Default_Settings())
        
        my_T.Plot_N_Orbits([op0.r_vals, op1.r_vals, op.r_vals],
                            labels=['Orbit Initial', 'Orbit Post Transfer', 'Orbit Propagation Post Transfer'],
                            show_plot=True, title = 'Transfer Plot')
###############################################################################