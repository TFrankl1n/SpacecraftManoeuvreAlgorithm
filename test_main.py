import numpy as np
from math import sqrt                                                        
import matplotlib.pyplot as plt
from scipy.integrate import ode
from scipy.integrate import quad
from mpl_toolkits.mplot3d import Axes3D


import os
import sys     
sys.path.append('F:/Users/Thomas/UNIVERSITY/Work [] Placements/WORK/5. AIRBUS/1. Work/Manoeuvres/Manoeuvre_Alg/')
#sys.path.append('../Resources/')
# Get the directory from sys.path, identifying filename and constructing path
#directory = sys.path[-1]
#filename_1 = 'TLE_1_ISS.txt'
#file_path = os.path.join(directory, filename_1)

from orbit_propagator import Orbit_Propagator as my_OP
from orbit_propagator import Default_Settings
import planetary_data as my_PD
import tools as my_T

centralbody = my_PD.earth

tspan = 3600 * 24.0 *5.0
dt = 1000.0

if __name__ == '__main__':    
    print('Running...')
    perturbations = Default_Settings()
    perturbations['J2'] = True
    r0 = 6570
    r1 = 42160
    
    #PLOTTING ONE  HOHMANN TRANSFERS, ONE PLANAR, - ERROR
    #a, e, i, ta, aop, raan, [year, month, day, hour]]
    c2 = [centralbody['radius']+r0, 0.0, 16.0, 0.0, 30.0, 10.0]
    c3 = [centralbody['radius']+r1, 0.0, 16.0, 180.0,30.0, 10.0]
    c2 = [r0, 0.0, 16.0, 0.0, 30.0, 10.0]
    c3 = [r1, 0.0, 16.0, 180.0,30.0, 10.0]

#   op0 = my_OP(c0, tspan, dt, COEs=True, deg=False, perturbations=perturbations)   

    
    op2, op3, op_transfer, DeltaV_Vals, tof_transfer, e_transfer = my_T.Hohmann_Transfer(COEs0 = c2, COEs1=c3, propagate = True)  
        
    print("op2",len(op2.r_vals))
    print("op3",len(op3.r_vals))
    print(DeltaV_Vals)
    my_T.Plot_N_Orbits([op2.r_vals, op3.r_vals, op_transfer.r_vals],labels=['Orbit 1 Initial', 'Orbit 1 Final', 'Orbit 1 Transfer'],show_plot=False, title = 'Transfer Plot')

    print('DeltaV 0: \t %.3f (km/s)' % DeltaV_Vals[0])
    print('DeltaV 1: \t %.3f (km/s)' % DeltaV_Vals[1])
    print('Total DeltaV: ', DeltaV_Vals[1] + DeltaV_Vals[0])
    print('Eccentricity of Transfer: ', e_transfer)
    print('Transfer Time: \t',tof_transfer,' seconds, ',tof_transfer/3600, "hours")   
    
    print("\n\nTo transfer from these two altitudes, the engines must perform a total velocity change of approx.", round(DeltaV_Vals[1] + DeltaV_Vals[0], 4),"(km/s)")
    print("This manoeuvre will take approx.", round(tof_transfer/3600, 3), "hours")
    
    #op2.plot_3d(show_plot=(False))
    #op2.calculate_COEs()
    #op3.plot_COEs(show_plot=(False))
    op2.plot_vs(show_plot=True, hours=True)
    #op2 = my_OP(my_T.TLE_to_CEOsPlus(file_path), tspan, dt, COEs = True, deg=False, perturbations=perturbations)

    
    #PLOTTING TWO  HOHMANN TRANSFERS, ONE PLANAR, ONE WITH INCLUNATION
    # c2 = [centralbody['radius']+r0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # c3 = [centralbody['radius']+r1, 0.0, 0.0, 180.0,0.0, 0.0]
    # c2_rot = [centralbody['radius']+r0, 0.0, 30.0, 0.0, 100.0, 20.0]
    # c3_rot = [centralbody['radius']+r1, 0.0, 30.0, 180.0, 100.0, 20.0]
    
    
    # op2, op3, op_transfer, DeltaV_Vals = my_T.Hohmann_Transfer(COEs0 = c2, COEs1=c3, propagate = True)  
    # # op2, op3, op_transfer, DeltaV_Vals = my_T.Hohmann_Transfer(propagate = True)  
    
    # op2_rot, op3_rot, op_transfer_rot, DeltaV_Vals_rot = my_T.Hohmann_Transfer(COEs0 = c2_rot, COEs1=c3_rot, propagate = True)
    # # op2_rot, op3_rot, op_transfer_rot, DeltaV_Vals_rot = my_T.Hohmann_Transfer(propagate = True)
    
    # print("op2",len(op2.r_vals))
    # print("op3",len(op3.r_vals))
    # print("op2_rot",len(op2_rot.r_vals))
    # print("op3_rot",len(op3_rot.r_vals))
    # my_T.Plot_N_Orbits([op2.r_vals, op3.r_vals, op_transfer.r_vals, op2_rot.r_vals, op3_rot.r_vals, op_transfer_rot.r_vals],
    #                     labels=['Orbit 1 Initial', 'Orbit 1 Final', 'Orbit 1 Transfer', 'Orbit 2 Initial', 'Orbit 2 Final', 'Orbit 2 Transfer'],
    #                     show_plot=True, title = 'Transfer Plot')
                        #az = -45.0, el = 0.0, axes = 9500,
                        
    #print("Here are the DeltaV Values:\n", DeltaV_Vals)
    #print("ANd the rotated ones why not:\n", DeltaV_Vals_rot)

    

    
