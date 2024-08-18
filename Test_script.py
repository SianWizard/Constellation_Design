# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:31:39 2024

@author: sina
"""

# Test script
from math import degrees,radians,ceil,pi
import CovGeometry
import ConstConcept

### Introduction

# This script finds a fast analytical LEO constellation design for EO purposes
# It satisfies at least 2 obsevations per 24 hours with minimum number of sats
# for a desired latitude band. This program assumes 2BP problem with spherical 
# Earth. The constellation contains symmetrical circular orbits only.
# Further analysis should be done to visualize the constellation.
# Orbital phasing are independant of the phasing in the neighboring planes


### Inputs

i = radians(90) # Desired satellite inclination [rad]
h = 600 # Maximum satellite altitude [km]
latitude_band = (radians(0),radians(90)) # Desired latitude band between 0 and 90 [deg]
max_gap_time = (0.5)*3600 # Maximum gap time in target coverage [s]
### Setting either target minimum elevation angle or the satellite sensor half-angle
min_el = radians(40)
#beta = radians(45)






### Preliminary calculations

theta,beta,el = CovGeometry.h2theta(h, min_el, 'e')
#theta,beta,el = CovGeometry.h2theta(h, beta, 'b')



### For Equatorial orbits
if abs(pi/2-i)>(pi/2-theta):
    miu = 398600.44
    Re = 6378.14
    Tsat = 2*pi*((h+Re)**3/miu)**0.5
    min_num_sat = ceil(Tsat/max_gap_time)
    raans = [0]
    P = 1

else:
    ### Minimum number of satellites per plane to cover a minimum latitude

    min_num_sat,delta_t_gap = ConstConcept.minSatperPlane(h, theta , i, latitude=latitude_band[0], max_allow_gap=max_gap_time)
    #print("Number of Satellites per plane needed for minimum 2 revisits per 24 hrs = ",min_num_sat)
    #print("The minimum inclination for the orbits for optimal coverage on the desired latitude band is = ",ceil(degrees(latitude_band[-1]))," [deg]")
    if abs(pi/2-i)>abs(pi/2-latitude_band[1]):
        print("\nInput inclination is lower than the required minimum inclination angle!!!")

    ### Minimum number of Planes needed to satisfy maximum Gap time
    P,raans,DELTA = ConstConcept.NumPlanesNeeded(i,latitude_band[0],theta, max_gap_time, delta_t_gap)
    #print("Number of orbital planes needed with ", degrees(i), " [deg] inclination angle, is = ",P)




    ########################## Checking different latitudes ########################
    lats = [(latitude_band[0]+latitude_band[1])/2,latitude_band[1]]
    for lat in lats:
        P_temp,raans_temp,DELTA_temp = ConstConcept.NumPlanesNeeded(i,lat,theta, max_gap_time, delta_t_gap)
        if P_temp>P:
            P,raans,DELTA = P_temp,raans_temp,DELTA_temp
        elif P_temp==P:
            if DELTA_temp<DELTA:
                P,raans,DELTA = P_temp,raans_temp,DELTA_temp
            
            
    ## Checking
    P1,raans1,DELTA1 = ConstConcept.NumPlanesNeeded(i,latitude_band[0],theta, max_gap_time, delta_t_gap)
    P2,raans2,DELTA2 = ConstConcept.NumPlanesNeeded(i,(latitude_band[0]+latitude_band[1])/2,theta, max_gap_time, delta_t_gap)
    P3,raans3,DELTA3 = ConstConcept.NumPlanesNeeded(i,latitude_band[1],theta, max_gap_time, delta_t_gap)
    #raans = raans1

raans_deg = [0] * len(raans)

for u in range(len(raans)):
    raans_deg[u] = degrees(raans[u])

### Final Output
print("\n\n***************************************************")
print("************  Suggested Constellation  ************")
print("***************************************************\n")

print("Selected inputs:\n")
print("Orbits altitude = ", h, " [km]")
print("Orbits inclination = ", degrees(i), " [deg]"," (Required minimum Inc = ", ceil(degrees(latitude_band[-1]))," [deg])")
print("Minimum target elevation angle = ", degrees(el), " [deg]")
print("Satellite FOV half angle = ", degrees(beta), " [deg]")
print("Maximum visibility gap time = ", (max_gap_time/3600), " [hr]")

print("\nSuggested constellation:\n")
print("Total number of satellites = ", P*min_num_sat)
print("Number of satellites per plane= ", min_num_sat)
print("Number of orbital planes = ", P)
print("RAAN distribution = ", *raans_deg)