#!/usr/bin/env python2

import numpy as np
from itertools import chain
from sys import argv
import sys
from periodic import element

NBO_in = open(argv[1],"r").read().splitlines()

start = []
stop = []
for i,j in enumerate(NBO_in):

    if j == "    Atom  No    Charge         Core      Valence    Rydberg      Total":
        start.append(i)
    if j == " =======================================================================":
        stop.append(i)

symbol,number,natural_charge,core,valence,rydberg,total = [],[],[],[],[],[],[]

if len(start)==len(stop):
    for i in range(len(start)):
        for line in NBO_in[start[i]+2:stop[i]]:
            cols = line.split()
            symbol.append(cols[0])
            number.append(int(cols[1]))
            natural_charge.append(float(cols[2]))
            core.append(float(cols[3]))
            valence.append(float(cols[4]))
            rydberg.append(float(cols[5]))
            total.append(float(cols[6]))
    

for i in range(len(symbol)):
    print "{0:<2}  {1:>3}  {2:>8.5f}".format(symbol[i], number[i], natural_charge[i])
