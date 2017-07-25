#!/usr/bin/env python

import numpy as np
from itertools import chain
from sys import argv,exit
from periodic import element
import matplotlib.pyplot as plt

NBO_in = open(argv[1],"r").read().splitlines()

start = []
stop = []

count=0
for i,j in enumerate(NBO_in):

    if j == "    Atom  No    Charge         Core      Valence    Rydberg      Total":

        start.append(i)
    if j == " =======================================================================":
        stop.append(i)
alphasymbol,alphanumber,alphanatural_charge = [],[],[]
betasymbol,betanumber,betanatural_charge = [],[],[]

if len(start)==len(stop):
    for line in NBO_in[start[1]+2:stop[1]]:
        cols = line.split()
        alphasymbol.append(cols[0])
        alphanumber.append(int(cols[1]))
        alphanatural_charge.append(float(cols[2]))
    for line in NBO_in[start[2]+2:stop[2]]:
        cols = line.split()
        betasymbol.append(cols[0])
        betanumber.append(int(cols[1]))
        betanatural_charge.append(float(cols[2]))
else:
    sys.exit("Error reading charges")
density=[]
if alphasymbol==betasymbol and alphanumber==betanumber:
    for i in range(len(alphasymbol)):
        print "{0:<2}  {1:>3}  {2:>8.5f}".format(alphasymbol[i], alphanumber[i], betanatural_charge[i]-alphanatural_charge[i])
        density.append(betanatural_charge[i]-alphanatural_charge[i])
else:
    sys.exit("Error")

plt.plot(alphanumber,density)
plt.show()
