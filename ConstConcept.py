# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 15:11:37 2024

@author: sina
"""
import numpy as np
from math import pi,ceil,sin,cos,acos,asin

# minimum number of satellites required per plane to have visibility at least once per day
# spherical earth and circular orbit assumption (2BP)
# Orbital phasing are independant of the phasing in the neighboring planes
def minSatperPlane(h,theta, i=pi/2,latitude = 0 , max_allow_gap=24*3600,miu = 398600.44,Re = 6378.14,omega_e = 7.292115e-5):

    if abs(pi/2-i)>abs(pi/2-latitude):
        print("Can't find the optimal nominal number of sats per plane!!")
    r = h+Re
    Tsat = 2*pi*(r**3/miu)**0.5

    delta_landa = omega_e * Tsat
    
    
    if (sin(theta)/cos(latitude))>1:
        min_no_SatPerPlane = 1
    else:
        theta_lat = asin(sin(theta)/cos(latitude)) 
        min_no_SatPerPlane = 1+ceil((delta_landa-2*theta_lat)/(2*theta_lat)) # This gives the min number of sats for 1 visit per strip pass
        if min_no_SatPerPlane<1:
            min_no_SatPerPlane=1
    # Until now minimum number of sats per plane for at least 1 visit per day has been calculated
    # Now the maximum allowable gap time constraint        
    
    delta_t_gap = Tsat/min_no_SatPerPlane
    if delta_t_gap>max_allow_gap:
        min_no_SatPerPlane = ceil(Tsat/max_allow_gap)
        delta_t_gap = Tsat/min_no_SatPerPlane
    
    return min_no_SatPerPlane,delta_t_gap





# Calculating the orbit variation from 180 deg at a specific latitde
def NumPlanesNeeded(i,lat,theta, max_gap, delta_t_gap ,omega_e = 7.292115e-5):
    # calculations for any single plane at a specific latitude   
    # From here on out, calculations are done over a latitude plane and angles
    Om_lat_lat=OmegaLat(i,lat)
    
    
    ######
    # if (sin(theta)/cos(lat))>1:
    #     return 1,[0],0
    # theta_lat = asin(sin(theta)/cos(lat))
    # print(degrees(2*theta_lat))
    #######
    
    if (sin(theta)/cos(lat))>1:
        theta_lat = pi
    else:
        theta_lat = asin(sin(theta)/cos(lat))
    
    
    
    max_gap_angle = max_gap * omega_e
    sep = pi-2*Om_lat_lat
    
    #first_el_prev = 0
    #last_el_prev = sep
    
    #first_el_new = first_el_prev + 0
    #last_el_new = last_el_prev + 0
    
    P=1
    raans = [0]*360
    ind = 0
    #################### Snippet 1
    # DELTA = max_gap_angle-2*theta_lat
    
    # while (last_el_prev-first_el_new+2*theta_lat)>max_gap_angle:
    #     P=P+1
    #     first_el_new = first_el_new + DELTA
    #     last_el_new = first_el_new+sep
        
    #     ind = ind + 1
    #     raans[ind] = raans[ind-1] + DELTA
    
    # while (2*pi-last_el_new+2*theta_lat)>max_gap_angle:
    #     P=P+1
    #     first_el_new = last_el_new + DELTA
    #     last_el_new = first_el_new+sep
        
    #     ind = ind + 1
    #     raans[ind] = raans[ind-1] + sep + DELTA
    
    # del raans[ind+1:]
    #####################
    
    ##################### Snippet 2
    # DELTA = max_gap_angle
    
    # while (last_el_prev-first_el_new)>max_gap_angle:
    #     P=P+1
    #     first_el_new = first_el_new + DELTA
    #     last_el_new = first_el_new+sep
        
    #     ind = ind + 1
    #     raans[ind] = raans[ind-1] + DELTA
    
    # while (2*pi-last_el_new)>max_gap_angle:
    #     P=P+1
    #     first_el_new = last_el_new + DELTA
    #     last_el_new = first_el_new+sep
        
    #     ind = ind + 1
    #     raans[ind] = raans[ind-1] + sep + DELTA
    
    # del raans[ind+1:]
    #####################
    
    ##################### Snippet 3
    # delta_t_gap_angle = delta_t_gap*omega_e
    # DELTA = max_gap_angle-2*delta_t_gap_angle+2*theta_lat
    
    
    # #inner gap
    # gap = (last_el_prev-first_el_new)-2*theta_lat+2*delta_t_gap_angle
    # #print("sep: ",sep,", DELTA: ",DELTA,", gap: ",gap,", max_gap_angle: ",max_gap_angle,", ind: ",ind)
    # while gap>max_gap_angle:
    #     P=P+1
    #     first_el_new = first_el_new + DELTA
    #     ind = ind + 1
    #     raans[ind] = first_el_new
    #     gap = (last_el_prev-first_el_new)-2*theta_lat+2*delta_t_gap_angle
    
    # #outer gap
    # last_el_new = first_el_new + sep
    # gap = (2*pi - last_el_new)-2*theta_lat+2*delta_t_gap_angle
    # #print("sep: ",sep,", DELTA: ",DELTA,", gap: ",gap,", max_gap_angle: ",max_gap_angle)
    # while gap>max_gap_angle:
    #     P=P+1
    #     first_el_new = last_el_new + DELTA
    #     last_el_new = first_el_new+sep
    #     ind = ind + 1
    #     raans[ind] = first_el_new
    #     gap = (2*pi - last_el_new)-2*theta_lat+2*delta_t_gap_angle
    #     #print(gap)
    
    # del raans[ind+1:]
    #####################
    
    ##################### Snippet 4
    delta_t_gap_angle = delta_t_gap*omega_e
    DELTA = max_gap_angle-2*delta_t_gap_angle+2*theta_lat
    
    ind = -1
    P = 0
    raan_start = 0
    raan_end = -DELTA
    gap_out = (2*pi - raan_end)-2*theta_lat+2*delta_t_gap_angle
    while gap_out>max_gap_angle:
        P=P+1
        ind = ind + 1
        raan_start = raan_end + DELTA
        raans[ind]=raan_start
        raan_end = raan_start + sep
        raan_fix = raan_end
        gap_in = (raan_fix-raan_start)-2*theta_lat+2*delta_t_gap_angle
        while gap_in>max_gap_angle:
            P=P+1
            ind = ind+1
            raan_start = raan_start + DELTA
            raans[ind] = raan_start
            raan_end = raan_start + sep
            gap_in = (raan_fix-raan_start)-2*theta_lat+2*delta_t_gap_angle
            gap_break = (2*pi - raan_start)-2*theta_lat+2*delta_t_gap_angle
            if gap_break<max_gap_angle:
                break
        gap_out = (2*pi - raan_end)-2*theta_lat+2*delta_t_gap_angle
        gap_break = (2*pi - raan_start)-2*theta_lat+2*delta_t_gap_angle
        if gap_break<max_gap_angle:
            break
    del raans[ind+1:]
    #####################
    
    # Prograde vs Retrograde
    # if i>pi/2:
    #     for m in range(len(raans)):
    #         if m == 0:
    #             continue
    #         raans[m]=raans[m]-sep
    
    return P,raans,DELTA

# calculating the deviation of the orbit with respect to ascension point on specific latitude
def OmegaLat(i,lat):
    if i>pi/2:
        i = pi-i
    O = np.array([0,0,0])
    A = np.array([1,0,0])
    A_lat = np.array([cos(lat),0,sin(lat)])
    e = np.array([0,-sin(i),cos(i)])
    O_lat = np.array([0,0,sin(lat)])
    
    if (sin(lat)/sin(i)>1)and(abs(sin(lat)/sin(i)-1)<1e-5):
        landa = pi/2
    else:
        landa = asin(sin(lat)/sin(i))
    #print("gamma in degrees: ", degrees(landa))
    
    ri = A-O
    rf = cos(landa)*ri + sin(landa)*np.cross(e, ri)+(1-cos(landa))*(np.dot(e,ri))*e # Euler axis rotation
    
    if abs(rf[-1] - sin(lat))>1e-4:
        print("Something wrong with the rotation")
        
    r1 = A_lat-O_lat
    r2 = rf - O_lat
    Om_lat_lat = acos(np.dot(r1, r2)/(np.linalg.norm(r1)*np.linalg.norm(r2)))
    return Om_lat_lat