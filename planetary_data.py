# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 13:27:18 2023

@author: Thomas
"""

G_meters = 6.67408e-11

G = G_meters*10**-9

earth = {
    'name':'Earth',
    'mass':5.974e24,
    'mu':5.974e24*G,
    'radius':6378.0,
    'J2':1.08263e-3
}

moon = {
    'name':'Moon',
    'mass':5,
    'mu':5*G,
    'radius': 'small',
    'J2':202.7e-6
}

sun = {
    'name':'Sun',
    'mass':1.989e30,
    'mu':1.32712e11,
    'radius':695700.0, 
    'J2' : 0
}
    
