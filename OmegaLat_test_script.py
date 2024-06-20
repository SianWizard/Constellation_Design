# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 13:28:02 2024

@author: sina
"""

# OmegaLat test script

from ConstConcept import OmegaLat
from math import degrees,radians,pi
from numpy import  arange


i = radians(130)
max_lat = pi/2 - abs(pi/2-i)
lats = arange(0,max_lat+radians(1),radians(5))
omega_lats = [0]*len(lats)
index = range(len(lats))
for m in index:
    omega_lats[m] = degrees(OmegaLat(i, lats[m]))