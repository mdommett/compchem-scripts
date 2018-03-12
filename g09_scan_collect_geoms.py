#!/usr/bin/env python2

from sys import argv,exit
from periodic import element



finp = open(argv[1],"r").read().splitlines()
fout = open(argv[1]+"-opt.xyz","w")

lines = filter(None, (line.rstrip() for line in finp))
hold_geom=[]
for i,j in enumerate(lines):
    if j.split()[0]=="NAtoms=":
        natoms=int(j.split()[3])
        break
for i,j in enumerate(lines):
    if j == " Number     Number       Type             X           Y           Z":
        hold_geom=[]
        count = int(i)
        for line in lines[count+2:count+natoms+2]:
            splt = line.split()
            coords=(element(splt[1]).symbol,float(splt[3]),float(splt[4]),float(splt[5]))
            hold_geom.append(coords)
    if j==" Optimization completed.":
        fout.write("{}\n\n".format(natoms))
        for line in hold_geom:
            symbol,x,y,z=line
            fout.write("{0:<2} {1:>13.9f} {2:>13.9f} {3:>13.9f}\n".format(symbol,x,y,z))
