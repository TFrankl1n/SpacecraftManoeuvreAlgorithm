G_meters = 6.67408e-11

G = G_meters*10**-9

earth = {
    'name':'Earth',
    'mass':5.974e24,
    'mu':5.974e24*G,
    'radius':6378.0, #km
    'J2':1.08263e-3,
    'deorbit_altitude':400, #km
    'Cd': 2.2,                                                                 # Drag coefficient (adjust as needed for craft's aero properties)
    'drag_area': 20.0,                                                         #Effective cross-sectional drag area (adjust based on the spacecraft's characteristics) m^2
    'sc_mass': 500.0,                                                          #Spacecraft mass in kg (data from default Astrogator prop values)         
    #sma # km 
    #SOI??
    #
}

moon = {
    'name':'Moon',
    'mass':5,
    'mu':5*G,
    'radius': 1737.4, #km
    'J2':202.7e-6
}

sun = {
    'name':'Sun',
    'mass':1.989e30,
    'mu':1.32712e11,
    'radius':695700.0, #km              696340.0?
    'J2' : 0
}
