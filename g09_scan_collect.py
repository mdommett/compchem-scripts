#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from sys import argv

infile = open(argv[1],"r").read().splitlines()
scf = []
for line in infile:
    split_line = line.split()
    if len(split_line) == 9:
        if split_line[0] == 'SCF':
            scf.append(float(split_line[4]))
    if len(split_line) == 5:
        if split_line[0] == 'Total':
            scf.append(float(split_line[4]))
    if len(split_line) == 2:
        if split_line[0] == 'Optimization':
            scf.append(str(split_line[0]))


count = 0
s0_au = []
s1_au = []

for i in scf:    
    if isinstance(i, basestring):
        s0_au.append(scf[count-2])
        s1_au.append(scf[count-1])
    count +=1
    
x = range(1,len(s0_au)+1)

fig,ax = plt.subplots()
plt.plot(x,s0_au,'o-b',linewidth=2)
plt.plot(x,s1_au,'o-r',linewidth=2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_tick_params(labelsize=12, width=1.5)
ax.yaxis.set_tick_params(labelsize=12, width=1.5)
plt.xlabel('X', fontsize=12)
plt.ylabel('Energy (au)', fontsize=12)
plt.show()
