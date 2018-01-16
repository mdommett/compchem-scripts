#!/usr/bin/env python2

import numpy as np
from itertools import chain
from sys import argv
import sys
from periodic import element
import math

coords = open(argv[1],"r")
coords_lines = coords.read().splitlines()
natoms = int(coords_lines[0])
symbol,x,y,z,q,atomic_no = [],[],[],[],[],[]

# Get atomic coordinates from the .xyz file of the dimer

for line in coords_lines[2:natoms+2]:
    coord = line.split()
    symbol.append(coord[0])
    x.append(float(coord[1]))
    y.append(float(coord[2]))
    z.append(float(coord[3]))
    
# Calculate the Centre of Mass (CoM)

Xnumerator,Ynumerator,Znumerator,Xdenominator,Ydenominator,Zdenominator=0,0,0,0,0,0

for i in range(natoms):
    Xnumerator += element(symbol[i]).mass*x[i]
    Ynumerator += element(symbol[i]).mass*y[i]
    Znumerator += element(symbol[i]).mass*z[i]
    Xdenominator += element(symbol[i]).mass
    Ydenominator += element(symbol[i]).mass
    Zdenominator += element(symbol[i]).mass

Xicm = Xnumerator/Xdenominator
Yicm = Ynumerator/Ydenominator
Zicm = Znumerator/Zdenominator

print Xicm,Yicm,Zicm
