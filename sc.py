#!/usr/bin/env python

from sys import argv

infile = open(argv[1],"r").read().splitlines()
scf = []
for line in infile:
    split_line = line.split()
    if len(split_line) == 9:
        if split_line[0] == 'SCF':
            scf.append(float(split_line[4]))
    if len(split_line) == 2:
        if split_line[0] == 'Optimization':
            scf.append(str(split_line[0]))

count = 0
for i in scf:
    
    if isinstance(i, basestring):
        print scf[count-1]
    count +=1