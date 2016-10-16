#!/usr/bin/env python

import numpy as np
from sys import argv

xyz = open(argv[1],'r').read().splitlines()
unit = raw_input("Would you like output in Bohr (B) or Angstroms (A)?\n") 
xc =[]
yc =[]
zc =[]

for line in xyz[2:]:
    ln = line.split()
    xc.append(float(ln[1]))
    yc.append(float(ln[2]))
    zc.append(float(ln[3]))

xmin = min(xc)
xmax = max(xc)
ymin = min(yc)
ymax = max(yc)
zmin = min(zc)
zmax = max(zc)
b2a = float(0.52917721067)
xdist = (xmax-xmin)
ydist = (ymax-ymin)
zdist = (zmax-zmin)

if unit == 'B':

	print "x distance : {0}".format(xdist)
	print "y distance : {0}".format(ydist)
	print "z distance : {0}".format(zdist)
	print "All distances are in Bohr"
elif unit == 'A':

	print "x distance : {0}".format(xdist*b2a)
        print "y distance : {0}".format(ydist*b2a)
        print "z distance : {0}".format(zdist*b2a)
        print "All distances are in Angstroms"
else: 
	print "I did not understand your units, enter B or A"
