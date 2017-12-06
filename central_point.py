#!/usr/bin/env python

import numpy as np
from sys import argv

xyz = open(argv[1],'r').read().splitlines()

xc =[]
yc =[]
zc =[]

for line in xyz[2:]:
    ln = line.split()
    xc.append(float(ln[1]))
    yc.append(float(ln[2]))
    zc.append(float(ln[3]))

x = np.average(xc)
y = np.average(yc)
z =np.average(zc)
print x,y,z
