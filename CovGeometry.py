# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 13:19:30 2024

@author: sina
"""

# Physical geometrical coverage constraints
# Spherical Earth assumption


# Parameter definition:
    # h = Altitude [km]
    # theta = Coverage HALF-angle [rad]
    # beta = Field Of View HALF-angle [rad]
    # el = Elevation angle [rad]
    # Re = Mean radius of Earth at Equator [km]

from math import pi,cos,sin,asin

# Following function converts the available satellite altitude to the half-coverage angle
# (should also provide either minimum elevation angle 'e' or the field of view half angle 'b')
# returns the coverage half-angle, fov half-angle, and the minimum elevation angle
def h2theta(h, C, char, Re = 6378.14):
    r = h+Re
    max_beta = asin(Re/r)
    max_theta = pi/2-max_beta
    if char == 'e':
        el = C
        if (el>pi/2)or(el<0):
            print("Incorrect minimum elevation angle input")
            return
        beta = asin(Re*sin(el+pi/2)/r)
        theta = pi/2-(beta+el)
        return theta,beta,el
    elif char == 'b':
        beta = C
        if beta>max_beta:
            return max_theta, max_beta, 0
        el = pi/2 - asin(r*sin(beta)/Re) # since asin output is in the [-pi/2,pi/2] range
        theta = pi/2-(beta+el)
        return theta,beta,el

# Following function converts the coverage half-angle to the required satellite altitude
# (should also provide either minimum elevation angle 'e' or the field of view half angle 'b')
# returns the altitude, fov half-angle, and the minimum elevation angle
def theta2h(theta, C, char, Re = 6378.14):
    max_beta = pi/2 - theta
    min_h = Re/cos(theta) - Re
    if char == 'e':
        el = C
        if (el>pi/2)or(el<0):
            print("Incorrect minimum elevation angle input")
            return
        beta = pi/2-(theta+el)
        r = Re*sin(pi/2+el)/sin(beta)
        return r-Re,beta,el
    elif char == 'b':
        beta = C
        if beta>max_beta:
            return min_h,max_beta,0
        el = pi/2 - (theta+beta)
        r = Re*sin(pi/2+el)/sin(beta)
        return r-Re,beta,el