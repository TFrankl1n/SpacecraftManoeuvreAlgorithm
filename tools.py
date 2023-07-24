# -*- coding: utf-8 -*-

import math as nm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import datetime

# plt.style.use('ggplot')

import planetary_data as my_PD
from orbit_propagator import Orbit_Propagator as my_OP

deg_to_rad = np.pi/180.0
rad_to_deg = 180.0/np.pi

def norm(v):
    return np.linalg.norm(v)

r_count = 1

# algorithm to take COE to r and v vectors in inertial frame to worht with 
def Plot_N_Orbits(r_vals, labels=['Trajectory','Initial Pos.'], centralbody = my_PD.earth, show_plot = False, save_plot = False, title='Multiple Orbit Plotting - From TLE Files'):
    print("Plotting Multiple Orbits...")
    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(111,projection = '3d')                                # add to screen (plot space) [one plot one one fig, first col, row, value]

    #plot trajectory - how many position arrays there are in r_vals variable = use for loop
    #print("r_vals:\n", r_vals)
    n = 0
                                                                  
    for r in r_vals:
        ax.plot(r[:,0],r[:,1],r[:,2],  label=labels[n])                        # Plot whole trajectory 
        ax.plot([r[0,0]], [r[0,1]], [r[0,2]], 'ko')                            #define initial orbit position (state- where orbit begin + orientation wrt visual frame)
        n+=1
        
 
    #plot central body - assuming orbiting perfect earth for now               - create sphere on 3d plot
    lat, long = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]                            #azimuth polar latitude/longitude for grid across surface of sphere, and find coords 
    coord_x = centralbody['radius']*np.cos(lat)*np.sin(long)
    coord_y = centralbody['radius']*np.sin(lat)*np.sin(long)
    coord_z = centralbody['radius']*np.cos(long)                                    #calculating height pos  
    #ax.plot_wireframe(coord_x, coord_y, coord_z, color = "grey")
    ax.scatter(coord_x, coord_y, coord_z, s=5, facecolor="grey")
 
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)
    # The plot bounding box - sphere - call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])
    
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_title(title)
    plt.legend(loc=3)                                                       
 
    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(title+".png",dpi=300)                                         
    
    print("rvals",len(r_vals))
    print("r_count??",r_count)
    #plt.plot(t_vals, r_vals)
    #r_count = r_count+1
    
    return r_vals
        
#convert classical orb. to r and v vectors
def Convert_COEs_to_RV(COEs=[],deg=False,mu=my_PD.earth['mu'], print_results = False):      #pass in coes, defining degrees (convert to radians, pass in mu)
    print("Converting COEs to Radial Vectors...")
    print("COES values:", COEs)
    print("length of COEs:", len(COEs))
    
    if deg:
        print("Converting COEs to radians...")
        a, e, i, ta, aop, raan = [np.radians(val) for val in COEs]
        if print_results:
            print('a =\t %.1f \ne =\t %.2f \ni =\t %.3f \nraan =\t %.4f \naop =\t %.6f \nta =\t %.6f' % (a, e, i, raan, aop, ta))
    else:
        print("COEs already in radians...")
        a, e, i, ta, aop, raan = COEs
        if print_results:
            print('a =\t %.1f \ne =\t %.2f \ni =\t %.3f \nraan =\t %.4f \naop =\t %.6f \nta =\t %.6f' % (a, e, i, raan, aop, ta))

    Ecc_Anom = Calc_Eccentric_Anomaly([ta,e],'tae')
    
    r_norm = a*(1-e**2)/(1+e*np.cos(ta))
    
    #calc r and v vectors in perifocal frame []
    r_perif = r_norm*np.array([nm.cos(ta),nm.sin(ta),0])                                                    # no motion in z axis
    v_perif = nm.sqrt(mu*a)/r_norm*np.array([-nm.sin(Ecc_Anom),nm.cos(Ecc_Anom)*nm.sqrt(1-e**2),0])
    
    #rotation matrix from perifocal to Earth Centred Inertial
    perif_to_eci = np.transpose(EarthCentralInertial_to_PeriFocal(raan,aop,i))
    
    #calculate r and v vectors in inertial frame
    r = np.dot(perif_to_eci,r_perif)                                            
    v = np.dot(perif_to_eci,v_perif)
    #print("Calculated State Vectors:\n r: \t %.1f \n v: %.2f\t" % (r, v))
    # print("r:\n", r)
    # print("v:\n", v)
    return r,v


def Convert_RV_to_COEs(r, v, mu=my_PD.earth['mu'],degrees = False, print_results = False):
    #print("Converting Radial Vectors to COEs ...")
    #norm of position vector
    r_norm = norm(r)
    #specific angular momentum                                                    
    h = np.cross(r,v)                                                           #EQ[2.18]
    h_norm=norm(h)
    #inclination
    i = nm.acos(h[2]/h_norm)
    #eccentricity vector
    e = ((norm(v)**2 - mu/r_norm) * r-np.dot(r,v) * v)/mu
    #eccentricity scalar
    e_norm = norm(e)
    #node line
    N = np.cross([0,0,1],h)                                                     
    N_norm = norm(N)
    #raan
    raan = nm.acos(N[0]/N_norm)
    if N[1]<0:
        raan = 2*np.pi - raan                                                   #quadrant check (as cos functions have limits of 180 to -80) EQ[6]
        
    #argument of perigee
    aop = nm.acos(np.dot(N,e)/N_norm/e_norm)
    if e[2]<0:
        aop = 2*np.pi - aop
    
    #true anomaly
    ta = nm.acos(np.dot(e,r)/e_norm/r_norm)
    if np.dot(r,v)<0:
        ta = 2*np.pi - ta
    
    #SMA
    a = r_norm * (1+e_norm * nm.cos(ta)) / (1-e_norm**2)
    
    if print_results:
        print('a = ',a)
        print('e = ',e_norm)
        print('i = ',i*rad_to_deg)
        print('raan = ', raan*rad_to_deg)
        print('aop = ', aop*rad_to_deg)
        print('ta = ', ta*rad_to_deg)
        
    if degrees: 
        return [a, e_norm, i*rad_to_deg, ta*rad_to_deg, aop*rad_to_deg, raan*rad_to_deg]
    else:
        return [a, e_norm, i, ta, aop, raan]

def TLE_to_CEOsPlus(tle_filename, mu=my_PD.earth['mu']):                        #convert TLE info to COE + additional useful parameters
    print("Reading TLE file")
 #read TLE file
    with open(tle_filename,'r') as f:
        lines = f.readlines()
        
    #seperate into three lines
    line0 = lines[0].strip()                                                    #name of sat
    line1 = lines[1].strip().split()                                            #dataline 1      strip removes newline characters, split removes and splits string by spaces  
    line2 = lines[2].strip().split()                                            #dataline 2                     ""                             ""
    
    #epoch (yr and day)
    epoch = line1[3]                                                            
    year, month, day, hour = Calc_Epoch(epoch)
    
    #collect coes
    
    #inclination, rad                                                           #third term of second line
    i = float(line2[2])*deg_to_rad
    #right ascention of ascending node, rad
    raan = float(line2[3])*deg_to_rad
    #eccentricity - corrected calc below with 0.
    e_string = line2[4]
    e = float('0.'+e_string)
    #argument of perigee, rad
    aop = float(line2[5])*deg_to_rad
    #mean anomaly, rad
    Me = float(line2[6])*deg_to_rad
    #mean motion, revs/day
    mean_motion = float(line2[7]) 
    #period, seconds
    T = 1/mean_motion*24*3600
    #SMA
    a = (T**2*mu/4.0/np.pi**2)**(1/3.0)
    #calc eccentric anomaly                                                    #Newton's method to solve for ecc anomaly
    Ecc_Anom = Calc_Eccentric_Anomaly([Me,e],'newton')    
    #calc rue anomaly
    ta = Calc_True_Anomaly([Ecc_Anom,e])
    return [a, e, i, ta, aop, raan, [year, month, day, hour]]

def TLE_to_RadVects(tle_filename):
    print("Converting TLE COEs to Radial Vectors")
    return Convert_COEs_to_RV(TLE_to_CEOsPlus(tle_filename))

#inertial to perifocal rotation matrix
def EarthCentralInertial_to_PeriFocal(raan,aop,i):
    print('Frame of reference fixing...? - central inertial to perifocal')
    row0 = [-nm.sin(raan)*nm.cos(i)*nm.sin(aop) + nm.cos(raan)*nm.cos(aop), nm.cos(raan)*nm.cos(i)*nm.sin(aop) + nm.sin(raan)*nm.cos(aop), nm.sin(i)*nm.sin(aop)]
    row1 = [-nm.sin(raan)*nm.cos(i)*nm.cos(aop) - nm.cos(raan)*nm.sin(aop), nm.cos(raan)*nm.cos(i)*nm.cos(aop) - nm.sin(raan)*nm.sin(aop), nm.sin(i)*nm.cos(aop)]
    row2 = [nm.sin(raan)*nm.sin(i), -nm.cos(raan)*nm.sin(i),nm.cos(i)]
    return np.array([row0, row1, row2])
   
    
def Calc_Eccentric_Anomaly(arr,method,tol=1e-8):                               #calculate ecentric anomaly - required for other calcs
    print("Calculating Eccentric Anomaly...")
    #we have true anomaly and eccentricity which can be subbed in straight from the coes 
    if method =='newton':                                                       #random guesses of x, find slope, mostly for orbits WILL converge so is good numerical method
        #newton's method for iteratively finding Eccentric anomaly (Ecc_Anom))  #OTHER METHODS ARE AVAILABLE SO COMPARE OPTIONS FOR THIS******
        Me,e=arr
        if Me < np.pi/2.0: 
            Ecc_Anom0 = Me + e / 2.0
        else:
            Ecc_Anom0 = Me - e
        for i in range(200): #arbitrary no. of stes:
            ratio = (Ecc_Anom0-e*np.sin(Ecc_Anom0)-Me)/(1-e*np.cos(Ecc_Anom0));
            if abs(ratio)<tol:
                if i==0:
                    return Ecc_Anom0
                else:
                    return "There may be an error here, Calc_Eccentric_Anomaly value for i"
            else:
                Ecc_Anom1 = Ecc_Anom0-ratio
                Ecc_Anom0 = Ecc_Anom1
                #did not converge
        return False
    elif method == 'tae':
        ta,e = arr
        return 2.0*nm.atan(nm.sqrt((1-e)/(1+e))*nm.tan(ta/2.0))
    else:
        print('Invalid method for eccentric anomaly')
                

def Calc_Epoch(epoch):
    print("Calculating Epoch...")
    #epoch year
    year = int('20'+epoch[:2])
    epoch = epoch[2:].split('.')
    day_of_year = int(epoch[0])-1
    hour_of_day = float('0.'+epoch[1])*24.0
    
    #combine for year-month-day
    date = datetime.date(year,1,1)+datetime.timedelta(day_of_year)             # counts in leap year and proper dayes and dates 1, 1, gives first day of year, add no of days into year NOW
    
    #extract month and day 
    month = float(date.month)
    day = float(date.day)
    
    return year, month, day, hour_of_day


def Calc_True_Anomaly(arr):
    print("Calculating True Anomaly")
    Ecc_Anom, e =arr
    if e >= 1:
        raise ValueError("Eccentricity must be between 0 and 1")
    return 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(Ecc_Anom/2.0))

def Calc_RV_to_Perifocal():
    return "Calc_RV_to_Perifocal??"




def Hohmann_Transfer_Values(r0, r1, altitude = True, centralbody = my_PD.earth):
    print("Hohmann_Transfer_Values Function Running")
    #if using altitude instead of sma
    if altitude:
        r0 += centralbody['radius']
        r1 += centralbody['radius']
    #assuming e=0
    a_transfer = (r0+r1)/2.0                                                    #calc sma of transfer
    v_circ_init = nm.sqrt(centralbody['mu']  / r0)                                       #calc circular orbit velocity
    V_circ_final = nm.sqrt(centralbody['mu'] / r1)
    v0_transfer = nm.sqrt(centralbody['mu'] * (2/r0 - 1/a_transfer))                     #calc initial and final transfer orbit vel - use in prop to change in
    v1_transfer = nm.sqrt(centralbody['mu'] * (2/r1 - 1/a_transfer))
    
    tof_transfer = np.pi * nm.sqrt (a_transfer**3 / centralbody['mu'])
    DeltaV_vals = [abs(v0_transfer - v_circ_init), abs(V_circ_final - v1_transfer)]
    e_transfer = 1 - r0/a_transfer
    
    return DeltaV_vals, tof_transfer, e_transfer                                             #DeltaV values for each burn.



def Hohmann_Transfer(r0 = 0, r1 = 0, centralbody=my_PD.earth, COEs0=[], COEs1=[], altitude=True, 
                     dt = 100, names=['Initial', 'Final', 'transfer'], propagate = True, ):
    #check COEs are passed in correctly
    
    print("Hohmann_Transfer Function Running")

    if COEs0:
        r0 = COEs0[0]
        r1 = COEs1[0]

    elif altitude: #use central body rad if not
        r0 += centralbody['radius']
        r1 += centralbody['radius']

    #assuming e=0 [circular orbit]
    a_transfer = (r0+r1)/2.0                                                    #calc sma of transfer
    v_circ_init = nm.sqrt(centralbody['mu']  / r0)                                       #calc circular orbit velocity
    V_circ_final = nm.sqrt(centralbody['mu'] / r1)
    v0_transfer = nm.sqrt(centralbody['mu'] * (2/r0 - 1/a_transfer))                     #calc initial and final transfer orbit vel
    v1_transfer = nm.sqrt(centralbody['mu'] * (2/r1 - 1/a_transfer))
    
    t_transfer = np.pi * nm.sqrt (a_transfer**3 / centralbody['mu'])
    DeltaV_Vals = [abs(v0_transfer - v_circ_init), abs(V_circ_final - v1_transfer)]
    e_transfer = 1 - r0/a_transfer 
    
    #propagate orbit for DeltaV
    if propagate:
        if not COEs0:                                                           # if altitudes passed in = equitorial plot
            COEs0 = [r0, 0.0, 0.0, 0.0, 0.0, 0.0]
            COEs1 = [r1, 0.0, 0.0, 180.0,0.0, 0.0]
                                                 # eccentricity of transfer orbit [as ta=0 at prigee]
        #transfer orbit COEs
        COEs_transfer = [a_transfer, e_transfer, COEs0[2], 0.0, COEs0[4], COEs0[5]]
        
        #initial and final orbit period
        T_0 = 2 * np.pi * (r0**3 / centralbody['mu'])**0.5                      # to plot one complete orbit
        T_1 = 2 * np.pi * (r1**3 / centralbody['mu'])**0.5
    

        op0 = my_OP(COEs0, T_0, dt, COEs=True )
        op1 = my_OP(COEs1, T_1, dt, COEs=True )                                 #change these to tspan****
        op_transfer = my_OP(COEs_transfer, t_transfer, dt, COEs=True )          # all passed in must be in the same frame, always start with ta=0, rotating through 180 deg, keep other vals the same
        
        
        print("Delta vals: \t",DeltaV_Vals)
        return op0, op1, op_transfer, DeltaV_Vals, t_transfer, e_transfer
    else:
        print("We didn't propagate anything")
    return DeltaV_Vals, t_transfer, e_transfer     

                    