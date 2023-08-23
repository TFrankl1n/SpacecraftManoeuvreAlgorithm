import math as nm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import datetime

import planetary_data as my_PD
from orbit_propagator import Orbit_Propagator as my_OP
from orbit_propagator import Default_Settings

deg_to_rad = np.pi/180.0
rad_to_deg = 180.0/np.pi
    
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

def norm(v, axis=None):
    if axis is not None:
        return np.linalg.norm(v, axis=1)
    else:
        return np.linalg.norm(v)

def plot_3d(r_vals, show_plot=False, save_plot = False, title='Single Trajectory Plot', centralbody=my_PD.earth): #THIS IS UNECESSARY???? is for this one clearly
    print("Constructing diff-eq from state vectors, Plotting 3D figure of Trajectory") 
    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(111,projection = '3d')                                             #add to screen (plot space) [one plot one one fig, first col, row, value]

    ax.plot(r_vals[:,0],r_vals[:,1],r_vals[:,2],'k', label='Trajectory')                    # plot trajectory         
    ax.plot([r_vals[0,0]], [r_vals[0,1]], [r_vals[0,2]],'ko',label='Initial Position')      # define initial orbit position - visual of whats going on, where begin + prientation wrt visual frame
    
    #plot central body - assuming orbiting perfect earth for now                            - create sphere on 3d plot
    lat, long = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]                                         #azimuth polar latitude/longitude for grid across surface of sphere, and find coords 
    coord_x = centralbody['radius']*np.cos(lat)*np.sin(long)
    coord_y = centralbody['radius']*np.sin(lat)*np.sin(long)
    coord_z = centralbody['radius']*np.cos(long)                                             #calculating height pos  
    ax.scatter(coord_x, coord_y, coord_z, s=5, facecolor="grey")

    pos_Periapse = norm(r_vals, True).argmin()
    val_Periapse = r_vals[pos_Periapse]
    print("Here is the periapse value - after prop, should match value and position seen in propagator", pos_Periapse)
    ax.plot(val_Periapse[0], val_Periapse[1], val_Periapse[2], "bx", label = 'Periapse')    #Identify periapsis of orbit with blue x
    
    pos_Apoapse = np.argmax(norm(r_vals, True))
    val_Apoapse = r_vals[pos_Apoapse]
    print("Here is the apoapse value- after prop, should match value and position seen in propagator", pos_Apoapse)
    ax.plot(val_Apoapse[0], val_Apoapse[1], val_Apoapse[2], "rx", label = 'Apoapse')        #Identify apoapsis of orbit with red x
    
    pos_final = r_vals[-1]
    ax.plot(pos_final[0], pos_final[1], pos_final[2], 'wo', label='Final Position')          #Identify end of orbit with white spot   
    
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)
    plot_radius = 0.5*max([x_range, y_range, z_range])                                      # The plot bounding box - sphere - call half the max range the plot radius.
    
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_title(title)
    plt.legend(loc=3)                                                                       #if not labelled prev., use: plt.legend(['Trajectory','Initial Pos.'])
    if show_plot:
        plt.subplot_tool()
        plt.show()
    if save_plot:
        plt.savefig(title+".png",dpi=300)                                  

def Plot_N_Orbits(r_vals, labels=['Trajectory','Initial Pos.'], centralbody = my_PD.earth, show_plot = False, save_plot = False, title='Multiple Orbit Plotting',
                  dark_theme=False):                                                         # algorithm to take COE to r and v vectors in inertial frame to calculate
    print("Plotting Multiple Orbits...")
    
    if dark_theme:
        plt.style.use('dark_background')
        plt.rcParams['text.color'] = 'white'                                                # Set text color to white
        plt.rcParams['xtick.color'] = 'white'                                               # Set x-axis tick color to white
        plt.rcParams['ytick.color'] = 'white'                                               # Set y-axis tick color to white
        plt.rcParams['axes.labelcolor'] = 'white'                                           # Set axes label color to white
    else:
        plt.style.use('default')
        plt.rcParams['text.color'] = 'black'                                                # Reset text color to black
        plt.rcParams['xtick.color'] = 'black'                                               # Reset x-axis tick color to black
        plt.rcParams['ytick.color'] = 'black'                                               # Reset y-axis tick color to black
        plt.rcParams['axes.labelcolor'] = 'black'                                           # Reset axes label color to black

        
    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(111,projection='3d')                                               # Add to screen (plot space and 3D axis)
    #plot trajectory from initial position
    initial_position_label = False                                                          #check if label exists
    n=0                                                                                     # For label iteration!
    for r in r_vals:
        ax.plot(r[:,0].reshape(-1),r[:,1].reshape(-1),r[:,2].reshape(-1), label=labels[n])  # Plot whole trajectory 
        i=0
        while i != np.shape(r_vals)[0]:
            if not initial_position_label:
                ax.plot([r[0,0]], [r[0,1]], [r[0,2]], 'ko', label='Initial Position')       # Define initial orbit position, black spot
                initial_position_label = True            
            if initial_position_label:
                ax.plot([r[0,0]], [r[0,1]], [r[0,2]], 'ko')
            i+=1
        n+=1
        
    #Plot final position
    final_position_label = False                                                            #check if label exists                                 
    for r in r_vals:
        pos_final = r_vals[-1][-1]
        i=0
        while i != np.shape(r_vals)[0]:
            if not final_position_label:
                ax.plot([r[-1,0]], [r[-1,1]], [r[-1,2]], 'wo', label='F Position', markersize=2) 
                final_position_label = True                                                 #Label now exists
            if final_position_label: 
                ax.plot([r[-1,0]], [r[-1,1]], [r[-1,2]], 'wo', markersize=2) 
            i+=1
            
    #plot central body (sphere on 3d plot)
    lat, long = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]                            #azimuth polar latitude/longitude for grid across surface of sphere, and find coords 
    coord_x = centralbody['radius']*np.cos(lat)*np.sin(long)
    coord_y = centralbody['radius']*np.sin(lat)*np.sin(long)
    coord_z = centralbody['radius']*np.cos(long)                               #calculating height pos  
    wireframe=False                                                           #check if label exists
    if wireframe:
        ax.plot_wireframe(coord_x, coord_y, coord_z, color = "grey")
    else:
        if dark_theme:
            ax.scatter(coord_x, coord_y, coord_z, s=5, facecolor="white")
        else:
            ax.scatter(coord_x, coord_y, coord_z, s=5, facecolor="grey")
                                                                                                                                                                        
    #Plot Orbit's Nodes (Ascending and Descending - visually equal)                                                                                                                 #[Some node positions missed occasionally without evident correlation]
    Node_label = False                                                         #check if label exists
    for r in r_vals:                    # Get the coordinates of the point
        i=0
        while i != np.shape(r_vals)[0]:
            idx_for_nodes = np.argmin(np.abs(r[:, 2]))                              # Find index of point closest to zero vertically
            pos_node = r[idx_for_nodes] 
            if not Node_label:
                ax.plot(pos_node[0], pos_node[1], pos_node[2], 'g^', markersize=6, label='Orbit Node')  # Use a green triangle marker
                Node_label = True                                                           #Label now exists
            if Node_label:
                ax.plot(pos_node[0], pos_node[1], pos_node[2], "g^")
            i+=1
            
    #Plot Perigee and Apogee crosses on each orbit prop.
    i=0
    Periapse_label = False                                                                  #check if label exists
    Apoapse_label = False                                                                   #check if label exists
    while i != np.shape(r_vals)[0]:
        pos_Periapsis = np.argmin(norm(r_vals[i], True))
        val_Periapsis = r_vals[i][pos_Periapsis]
        if not Periapse_label:
            ax.plot(val_Periapsis[0], val_Periapsis[1], val_Periapsis[2], "bx", label = 'Periapse')
            Periapse_label = True                                                           #Label now exists
        pos_Apoapsis = np.argmax(norm(r_vals[i], True))
        val_Apoapsis = r_vals[i][pos_Apoapsis]
        if not Apoapse_label:
            ax.plot(val_Apoapsis[0], val_Apoapsis[1], val_Apoapsis[2], "rx", label = 'Apoapse') 
            Apoapse_label = True                                                            #Label now exists
        if Periapse_label:
            ax.plot(val_Periapsis[0], val_Periapsis[1], val_Periapsis[2], "bx")
        if Apoapse_label:
            ax.plot(val_Apoapsis[0], val_Apoapsis[1], val_Apoapsis[2], "rx") 
        i+=1

    #Adjust plot limits to fit graphing
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)
    plot_radius = 0.5*max([x_range, y_range, z_range])                                      # The plot bounding box - sphere - call half the max range the plot radius.
    
    #Plot within bounds
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    
    #Set labels and title
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_title(title)
    if dark_theme:
        ax.xaxis.pane.set_edgecolor('black')
        ax.yaxis.pane.set_edgecolor('black')
        ax.zaxis.pane.set_edgecolor('black')
        ax.grid(True, color=(0.121, 0.121, 0.121, 1.0))
    plt.legend(loc=3)
    plt.tight_layout()
    
    if dark_theme:
        ax.w_xaxis.set_pane_color((0.121, 0.121, 0.121, 1.0))                               # Adjust the RGB values as needed
        ax.w_yaxis.set_pane_color((0.121, 0.121, 0.121, 1.0))
        ax.w_zaxis.set_pane_color((0.121, 0.121, 0.121, 1.0))
    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(title+".png",dpi=300)   
    return r_vals
       


def Convert_COEs_to_SV(COEs=[],degrees=None,mu=my_PD.earth['mu'], print_results=False):     #convert classical orb. to r and v vectors  #pass in coes, defining degrees (convert to radians, pass in mu)
    
    if len(COEs)==7:                                                                        #Pulling COEs from TLE and saving to array
        a, e, i, raan, aop, ta = [(val) for val in COEs[:-1]]
        if print_results:
            print('a =\t %.1f \ne =\t %.2f \ni =\t %.3f \nta =\t %.4f \naop =\t %.6f \nraan =\t %.6f' % (a, e, i, ta, aop, raan))
    elif degrees:                                                                           #Converting COEs to radians (for calculation)
        a = COEs[0]
        e = COEs[1]
        i, ta, aop, raan = [np.radians(val) for val in COEs[2:6]]
        if print_results:
            print('a =\t %.1f \ne =\t %.2f \ni =\t %.3f \nta =\t %.4f \naop =\t %.6f \nraan =\t %.6f' % (a, e, i, ta, aop, raan))
    else:                                                                                   #pass in COEs as already in radians
        a, e, i, ta, aop, raan = COEs
        if print_results:
            print('a =\t %.1f \ne =\t %.2f \ni =\t %.3f \nta =\t %.4f \naop =\t %.6f \nraan =\t %.6f' % (a, e, i, ta, aop, raan))

    Ecc_Anom = Calc_Eccentric_Anomaly([ta,e],'tae')                                         #Use True Anomoly method to calculate Ecc Anom
    r_norm = a*(1-e**2)/(1+e*nm.cos(ta))                                                    

    #calc r and v vectors in perifocal frame                                                #no motion in z axis as perifocal in in plane
    r_perif = r_norm*np.array([nm.cos(ta),nm.sin(ta),0])
    v_perif = nm.sqrt(mu*a)/r_norm*np.array([-nm.sin(Ecc_Anom),nm.cos(Ecc_Anom)*nm.sqrt(1-e**2),0])

    perif_to_eci = (PQW_to_EarthCentralInertial(raan,aop,i))
    
    #calculate r and v vectors in inertial frame
    r = np.dot(perif_to_eci,r_perif)                                            
    v = np.dot(perif_to_eci,v_perif)
    return r,v

def Convert_COEs_to_SV_multiple(COEs_list, degrees=None, mu=my_PD.earth['mu'], print_results=None):
    if isinstance(COEs_list[0], list):                                                      # If COEs_list contains multiple sets of COEs        
        result = []
        for COEs in COEs_list:
            r, v = Convert_COEs_to_SV(COEs, degrees=degrees, mu=mu, print_results=print_results)
            result.append((r, v))
        return result
    else:                                                                                   # If COEs_list contains a single set of COEs
        return [Convert_COEs_to_SV(COEs_list, degrees=degrees, mu=mu, print_results=print_results)]

def Convert_SV_to_COEs(r, v, mu=my_PD.earth['mu'],degrees=None, print_results=None):        # Converting Radial Vectors to COEs
    r_norm = norm(r)
    h = np.cross(r,v)                                                                       #specific angular momentum                                                    
    h_norm=norm(h)

    i = nm.acos(h[2]/h_norm)                                                                #inclination
    e = ((norm(v)**2 - mu/r_norm) * r-np.dot(r,v) * v)/mu                                   #eccentricity vector
    e_norm = norm(e)                                                                        #eccentricity scalar
    if e_norm>1:
        print("Eccentricity value greater than expected, setting to max")
        e_norm=1.0
    if e_norm<0:
        print("Eccentricity value smaller than expected, setting to min")
        enorm=0.0
    
    N = np.cross([0,0,1],h)                                                                 #node line
    N_norm = norm(N)
    raan = nm.acos(N[0]/N_norm)                                                             #raan
    if N[1]<0:
        raan = 2*np.pi - raan                                                               #quadrant check (as cos functions have limits of 180 to -80) 
    aop = nm.acos(np.dot(N,e)/N_norm/e_norm)                                                #argument of perigee
    if e[2]<0:
        aop = 2*np.pi - aop
    ta = nm.acos(np.dot(e,r)/e_norm/r_norm)                                                 #true anomaly
    if np.dot(r,v)<0:
        ta = 2*np.pi - ta
    a = r_norm * (1+e_norm * nm.cos(ta)) / (1-e_norm**2)                                    #SMA                    
    if print_results:
        print("Printing results\na =", a, "\ne =", e_norm, "\ni =", i * rad_to_deg, "\nraan =", raan * rad_to_deg, "\naop =", aop * rad_to_deg, "\nta =", ta * rad_to_deg)
    if degrees: 
        return [a, e_norm, i*rad_to_deg, ta*rad_to_deg, aop*rad_to_deg, raan*rad_to_deg]
    else:
        return [a, e_norm, i, ta, aop, raan]

def Convert_COEs_deg_to_rad(COEs_deg):                                                      #says on the tin!!! (degrees to radians)
    a_deg, e, i_deg, ta_deg, aop_deg, raan_deg = COEs_deg
    i_rad = np.radians(i_deg)
    ta_rad = np.radians(ta_deg)
    aop_rad = np.radians(aop_deg)
    raan_rad = np.radians(raan_deg)
    return [a_deg, e, i_rad, ta_rad, aop_rad, raan_rad]

def TLE_to_CEOsPlus(tle_filename, mu=my_PD.earth['mu']):                                    #convert TLE info to COE + additional useful parameters
    print("Reading TLE file...")
    with open(tle_filename,'r') as f:
        lines = f.readlines()
        
    #seperate into three lines
    line0 = lines[0].strip()                                                    #name of sat
    line1 = lines[1].strip().split()                                            #dataline 1      strip removes newline characters, split removes and splits string by spaces  
    line2 = lines[2].strip().split()                                            #dataline 2                     ""                             ""

    #epoch (yr and day)
    epoch = line1[3]                                                            
    year, month, day, hour = Calc_Epoch(epoch)
    #collect coes:
    i = float(line2[2])*deg_to_rad                                                          #inclination, rad                                                           (#third term of second line)
    raan = float(line2[3])*deg_to_rad                                                       #right ascention of ascending node, rad
    e_string = line2[4]                                                                     #eccentricity - corrected calc below to place after 0.
    e = float('0.'+e_string)
    aop = float(line2[5])*deg_to_rad                                                        #argument of perigee, rad
    MeanAnom = float(line2[6])*deg_to_rad                                                   #mean anomaly, rad
    mean_motion = float(line2[7])                                                           #mean motion, revs/day
    T = 1/mean_motion*24*3600                                                               #period, seconds
    a = (T**2*mu/4.0/np.pi**2)**(1/3.0)                                                     #SMA
    Ecc_Anom = Calc_Eccentric_Anomaly([MeanAnom,e],'newton')                                 #calc eccentric anomaly (using iterative method for stable eqs)
    ta = Calc_True_Anomaly([Ecc_Anom,e])                                                    #calc True anomaly
    return [a, e, i, raan, aop, ta, [year, month, day, hour]]

def TLE_to_SV(tle_filename):                                                                # read in .txt TLE file and convert to State Vectors
    return Convert_COEs_to_SV(TLE_to_CEOsPlus(tle_filename))

def PQW_to_EarthCentralInertial(raan,aop,i):                                                #Convert to ICRF [if need to convert out of ICRF, use this function transposed]
    row0 = [nm.cos(raan)*nm.cos(aop) - nm.sin(raan)*nm.sin(aop)*nm.cos(i) , -nm.cos(raan)*nm.sin(aop) - nm.sin(raan)*nm.cos(aop)*nm.cos(i) ,  nm.sin(raan)*nm.sin(i)]
    row1 = [nm.sin(raan)*nm.cos(aop) + nm.cos(raan)*nm.sin(aop)*nm.cos(i) , -nm.sin(raan)*nm.sin(aop) + nm.cos(raan)*nm.cos(aop)*nm.cos(i) , -nm.cos(raan)*nm.sin(i)]
    row2 = [nm.sin(aop)*nm.sin(i) ,                                          nm.cos(aop)*nm.sin(i) ,                                           nm.cos(i)            ] 
    return np.array([row0, row1, row2])

def Calc_Eccentric_Anomaly(arr,method,tol=1e-8):                                            # calculate ecentric anomaly  (for TLEs)
    if method =='newton':                                                                   # Newton's method - performs random guesses of x, iteratively finding Eccentric anomaly (good for converging values)
        Me,e=arr
        if Me < np.pi/2.0: 
            Ecc_Anom0 = Me + e / 2.0
        else:
            Ecc_Anom0 = Me - e
        for i in range(200):                                                                #arbitrary no. of steps:
            ratio = (Ecc_Anom0-e*np.sin(Ecc_Anom0)-Me)/(1-e*np.cos(Ecc_Anom0));
            if abs(ratio)<tol:
                if i==0:
                    return Ecc_Anom0                                                        #did converge    
                else:
                    return Ecc_Anom0-ratio                                                  #did converge                                                        
            else:
                Ecc_Anom1 = Ecc_Anom0-ratio                                 
                Ecc_Anom0 = Ecc_Anom1
                return False                                                                #did not converge
    elif method == 'tae':                                                                   #Direct calculation method!
        ta,e = arr
        return 2.0*nm.atan(nm.sqrt((1-e)/(1+e))*nm.tan(ta/2.0))
    else:
        print('Invalid method for eccentric anomaly')
                
def Calc_Epoch(epoch):                                                                      #called by TLE reader for dates (for TLEs)
    #epoch year
    year = int('20'+epoch[:2])
    epoch = epoch[2:].split('.')
    day_of_year = int(epoch[0])-1
    hour_of_day = float('0.'+epoch[1])*24.0
    
    #combine for year-month-day
    date = datetime.date(year,1,1)+datetime.timedelta(day_of_year)                          # counts in leap year and proper dayes and dates 1, 1, gives first day of year, add no of days into year NOW
    
    #extract month and day 
    month = float(date.month)
    day = float(date.day)
    
    return year, month, day, hour_of_day

def Calc_True_Anomaly(arr):                                                                 #Calculate True Anomaly from Ecc Anom and eccentricity (for TLEs)
    Ecc_Anom, e =arr
    if e >= 1:
        raise ValueError("Eccentricity must be between 0 and 1")
    return 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(Ecc_Anom/2.0))



############################################################################################ Orbit propagation and Transfer Setup:

    

def Propagate_Orbit(COEs, tspan, dt, centralbody=my_PD.earth, withCOEs=None, add_perts=True, perturbations=Default_Settings(), calculate_COEs=None, 
                      degrees=None, print_results=None, grid_plot=None, full_plot=None):
    if add_perts:
        perturbations['J2'] = True   
    else:
        perturbations['J2'] = False
    op = my_OP(COEs, tspan, dt, withCOEs=withCOEs, degrees=degrees, perturbations=perturbations)
    if calculate_COEs:
        op.calculate_COEs()
        op.plot_COEs(show_plot=(True),hours = True, grid_plot=grid_plot, full_plot=full_plot)
    return op


def Hohmann_Transfer(r0=0, r1=0, centralbody=my_PD.earth, COEs0=[], COEs1=[], altitude=True, #Being used as the Phasing transfer function called by the Satellite Class
                          dt=100, names=['Initial', 'Final', 'Transfer'], propagate=True,
                          add_perts=True, perturbations=Default_Settings()):
    if COEs0:                                                                               # Check if COEs are passed in correctly
        a0 = COEs0[0]
        # e0 = COEs0[1]
        COEs1[3] = COEs0[3]+180  # *****NOT NECESSARY?
    if propagate:                                                                           # Check if we are propagating orbit correctly
        if add_perts:
            perturbations['J2'] = True   
        else:
            perturbations['J2'] = False
        if not COEs0: 
            print("Classical Orbital Elements have not been provided.")    
        
        T_0 = 2 * np.pi * (a0**3 / centralbody['mu'])**0.5
        op0 = my_OP(COEs0, T_0, dt, withCOEs=True, degrees=False, perturbations=perturbations)
        op0_COEs = Convert_SV_to_COEs(op0.r_vals[0], op0.v_vals[0])
        
        COEs1 = op0_COEs
        COEs1[0] = op0_COEs[0] + r1
        COEs1[3] = COEs0[3] + 180
        T_1 = 2 * np.pi * (COEs1[0]**3 / centralbody['mu'])**0.5
        op1 = my_OP(COEs1, T_1, dt, withCOEs=True, degrees=False, perturbations=perturbations)

        #apogee and perigee calcs
        r_p0 = norm(op0.r_vals[-1])

        #op1 is the perfect transfer, we want to find the 
        r_a1 = norm(op1.r_vals[0])
        
        a_transfer = (r_p0 + r_a1) / 2.0                                        # Calculate the transfer radius and velocities 
        v_circ_init = nm.sqrt(centralbody['mu'] / r_p0)                         # Calculate the transfer initial velocities 
        V_circ_final = nm.sqrt(centralbody['mu'] / r_a1)                        # Calculate the transfer final velocities 

        v0_transfer = nm.sqrt(centralbody['mu'] * (2 / r_p0 - 1 / a_transfer))
        v1_transfer = nm.sqrt(centralbody['mu'] * (2 / r_a1 - 1 / a_transfer))

        
        # Given velocity vector
        velocity_direction0 = op0.v_vals[0]                                         #np.array([, , , ])
        direction_norm0 = norm(op0.v_vals[0])                                    # Calculate the norm of the direction vector        
        normalized_direction0 = velocity_direction0 / direction_norm0              # Normalize the direction vector
        scaled_velocity_components0 = normalized_direction0 * v0_transfer         # Calculate the scaled velocity components
        
        velocity_direction1 = op1.v_vals[0]                                     #np.array([, , ])
        direction_norm1 = norm(op1.v_vals[0])                                    # Calculate the norm of the direction vector        
        normalized_direction1 = velocity_direction1 / direction_norm1              # Normalize the direction vector
        scaled_velocity_components1 = normalized_direction1 * v1_transfer         # Calculate the scaled velocity components

        SV_transfer_initial = np.concatenate((op0.r_vals[0], scaled_velocity_components0))
        tof_transfer = 2 * np.pi * nm.sqrt((a_transfer**3) / centralbody['mu'])     # Calculate time of flight and delta-v values
        DeltaV_Vals = [0.5*abs(scaled_velocity_components0 - v_circ_init), 0.5*abs(V_circ_final - scaled_velocity_components1)]
        e_transfer = 1 - r_p0 / a_transfer 
        
        if altitude: # Use central body radius if not
            r_p0 += centralbody['radius']
            r_a1 += centralbody['radius']
        
        op_transfer = my_OP(SV_transfer_initial, tof_transfer, dt, withCOEs=False)

        
        COEs_transfer_final = Convert_SV_to_COEs(op_transfer.r_vals[-1], op_transfer.v_vals[-1])
        
        return op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer
    else:
        print("We didn't propagate anything")
        return DeltaV_Vals, tof_transfer, e_transfer
    

def Hohmann_Transfer_ideal(r0=0, r1=0, centralbody=my_PD.earth, COEs0=[], COEs1=[], altitude=True, 
                          dt=100, names=['Initial', 'Final', 'transfer'], propagate=True,
                          add_perts=True, perturbations=Default_Settings()):
    # Check if COEs are passed in correctly
    if COEs0:
        a0 = COEs0[0]
        e0 = COEs0[1]
        a1 = COEs1[0]
        e1 = COEs1[1]
        #apogee and perigee calcs
        rp0 = a0*(1-e0)
        ra0 = a0*(1+e0)
        rp1 = a1*(1-e1)
        ra1 = a1*(1+e1)
    elif altitude: # Use central body radius if not
        rp0 += centralbody['radius']
        ra0 += centralbody['radius']
        rp1 += centralbody['radius']
        ra1 += centralbody['radius']
        
    # Calculate the transfer radius and velocities from initial orbit properties
    r_transfer = (rp0 + ra1) / 2.0
    v_circ_init = nm.sqrt(centralbody['mu'] / rp0)
    v_circ_final = nm.sqrt(centralbody['mu'] / ra1)
    v0_transfer = nm.sqrt(centralbody['mu'] * (2 / rp0 - 1 / r_transfer))
    v1_transfer = nm.sqrt(centralbody['mu'] * (2 / ra1 - 1 / r_transfer))
    
    # Calculate time of flight and delta-v values
    tof_transfer = np.pi * nm.sqrt((r_transfer**3) / centralbody['mu'])
    DeltaV_Vals = [0.5*abs(v0_transfer - v_circ_init), 0.5*abs(v_circ_final - v1_transfer)]
    e_transfer = 1 - rp0 / r_transfer 
    
    # Propagate orbit for DeltaV
    if propagate:
        if add_perts:
            perturbations['J2'] = True   
        else:
            perturbations['J2'] = False
        if not COEs0: 
            print("Classical Orbital Elements have not been provided.")            
        COEs_transfer_initial = [r_transfer, e_transfer, COEs0[2], 0.0, COEs0[4], COEs0[5]] #Construct initial state of transfer at optimal position
        
        op_transfer = my_OP(COEs_transfer_initial, tof_transfer, dt, withCOEs=True)         #Propagate optimal transfer orbit
        
        COEs_transfer_final = Convert_SV_to_COEs(op_transfer.r_vals[-1], op_transfer.v_vals[-1]) #Take final position after tof transfer complete (half orbit) and define final state in COEs
        
        COEs1[3] = COEs_transfer_final[3]                                                   #Set final True anomaly as found
        
        #Time Period for one complete orbit for initial and final orbit states (will use to propagate with)
        T_0 = 2 * np.pi * (a0**3 / centralbody['mu'])**0.5                                  
        T_1 = 2 * np.pi * (a1**3 / centralbody['mu'])**0.5
        
        op0 = my_OP(COEs0, T_0, dt, withCOEs=True, degrees=False, perturbations=perturbations)#                                                    [,One_Period=True) removed due to incorrect values being produced after structure change - did employ stepping time back method to ensure initial oribit final = transfer orbit initial states]
        op1 = my_OP(COEs1, T_1, dt, withCOEs=True, degrees=False, perturbations=perturbations)
        return op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer
    else:
        print("We didn't propagate anything")
    return DeltaV_Vals, tof_transfer, e_transfer

def Inclination_Transfer(i0 = 0, i1 = 0, centralbody=my_PD.earth, COEs0=[], COEs1=[], altitude=True, 
                      dt = 100, names=['Initial', 'Final', 'transfer'], propagate = True,
                      add_perts=True, perturbations=Default_Settings()):    
    #check COEs are passed in correctly
    if COEs0:
        r0 = COEs0[0]
        r1 = COEs0[0]
        i0 = COEs0[2]

    if altitude: #use central body rad if not
        #print("setting Altitude")    
        r0 += centralbody['radius']
        r1 += centralbody['radius']
        
    if propagate:
        if add_perts:
            perturbations['J2'] = True   
        else:
            perturbations['J2'] = False
        if not COEs0:                                                           # if altitudes passed in = equitorial plot
            #a, e, i, ta, aop, raan
            COEs0 = [r0, 0.0, 0.0, 0.0, 0.0, 0.0]
            COEs1 = [r1, 0.0, 0.0, 0.0,0.0, 0.0]
                
        #initial and final orbit period
        T_0 = 2 * np.pi * (r0**3 / centralbody['mu'])**0.5                      # to plot one complete orbit
        T_1 = 2 * np.pi * (r1**3 / centralbody['mu'])**0.5 
        
        #propagate orbit for DeltaV
        op0 = my_OP(COEs0, T_0, dt, withCOEs=True, degrees=True, perturbations=perturbations)                                             #change these to tspan**** and include tspan in initialising steps
        op1 = my_OP(COEs1, T_1, dt, withCOEs=True, degrees=True, perturbations=perturbations)   
        
        r0, v0 = Convert_COEs_to_SV(COEs=COEs0, degrees=False, mu=my_PD.earth['mu'], print_results = None)
        r1, v1 = Convert_COEs_to_SV(COEs=COEs1, degrees=False, mu=my_PD.earth['mu'], print_results = None)
        
        i0_rad = np.radians(i0)
        i1_rad = np.radians(i1)
        
        #Calculate the inclination change (plane change) in radians
        plane_change_rad = np.arccos(np.sin(i0_rad) * np.sin(i1_rad) + np.cos(i0_rad) * np.cos(i1_rad))         #
        
        op0_v_vals0_norm = norm(op0.v_vals[0])
        op1_v_vals0_norm = norm(op1.v_vals[0])        

        DeltaV_Val = 2*op0_v_vals0_norm*nm.sin(plane_change_rad/2)

        print("DeltaV_Vals", DeltaV_Val)
        return op0, op1, plane_change_rad, DeltaV_Val
    else:
        return plane_change_rad, DeltaV_Val





############################################################################################ Instant Hohmann Transfer - end of transfer orbit does not correlate with target orbit - Solution may be to apply apse line rotation which accounts for rotation beyond periapsis
# def Hohmann_Transfer(r0=0, r1=0, centralbody=my_PD.earth, COEs0=[], COEs1=[], altitude=True, 
#                           dt=100, names=['Initial', 'Final', 'transfer'], propagate=True,
#                           add_perts=True, perturbations=Default_Settings()):
#     print("Hohmann_Transfer Function Running...")
#     # Check if COEs are passed in correctly
#     if COEs0:
#         a0 = COEs0[0]
#         e0 = COEs0[1]
#         a1 = COEs1[0]
#         e1 = COEs1[1]
#         #apogee and perigee calcs
#         rp0 = a0*(1-e0)
#         ra0 = a0*(1+e0)
#         rp1 = a1*(1-e1)
#         ra1 = a1*(1+e1)
#     elif altitude: # Use central body radius if not
#         rp0 += centralbody['radius']
#         ra0 += centralbody['radius']
#         rp1 += centralbody['radius']
#         ra1 += centralbody['radius']
        
#     # Calculate the transfer radius and velocities
#     r_transfer = (rp0 + ra1) / 2.0
#     v_circ_init = nm.sqrt(centralbody['mu'] / rp0)
#     V_circ_final = nm.sqrt(centralbody['mu'] / ra1)
#     v0_transfer = nm.sqrt(centralbody['mu'] * (2 / rp0 - 1 / r_transfer))
#     v1_transfer = nm.sqrt(centralbody['mu'] * (2 / ra1 - 1 / r_transfer))
    
#     # Calculate time of flight and delta-v values
#     tof_transfer = np.pi * nm.sqrt((r_transfer**3) / centralbody['mu'])
#     DeltaV_Vals = [abs(v0_transfer - v_circ_init), abs(V_circ_final - v1_transfer)]
#     e_transfer = 1 - rp0 / r_transfer 
    
#     # Propagate orbit for DeltaV
#     if propagate:
#         if add_perts:
#             perturbations['J2'] = True   
#         else:
#             perturbations['J2'] = False
#         if not COEs0: 
#             print("Classical Orbital Elements have not been provided.")    
        
#         COEs_transfer = [r_transfer, e_transfer, COEs0[2], 0.0, COEs0[4], COEs0[5]] 
#         T_0 = 2 * np.pi * (a0**3 / centralbody['mu'])**0.5
#         T_1 = 2 * np.pi * (a1**3 / centralbody['mu'])**0.5
        
#         print("THING 1\n", COEs0)
#         print("THING 2\n", COEs1)
#         print("T_0\t", T_0)
#         print("T_1\t", T_1)
        
#         op0 = my_OP(COEs0, T_0, dt, withCOEs=True, degrees=False, perturbations=perturbations)
#         op1 = my_OP(COEs1, T_1, dt, withCOEs=True, degrees=False, perturbations=perturbations)
#         op_transfer = my_OP(COEs_transfer, tof_transfer, dt, withCOEs=True)
        
#         print("Delta vals: \t", DeltaV_Vals)
#         return op0, op1, op_transfer, DeltaV_Vals, tof_transfer, e_transfer
#     else:
#         print("We didn't propagate anything")
#     return DeltaV_Vals, tof_transfer, e_transfer