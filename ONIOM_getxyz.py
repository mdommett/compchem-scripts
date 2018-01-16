#!/usr/bin/env python2
from sys import argv
from periodic import element



finp = open(argv[1],"r").read().splitlines()
fout = open(argv[1]+".xyz","w")
#natoms = int(raw_input("How many atoms are in high layer? "))
lines = filter(None, (line.rstrip() for line in finp))
for i,j in enumerate(lines):
    if j.split()[0]=="NAtoms=":
	       natoms=int(j.split()[3])
for i,j in enumerate(lines):
    if j == " Number     Number       Type             X           Y           Z":
        count = int(i)
        fout.write("{}\n\n".format(natoms))
        for line in lines[count+2:count+natoms+2]:
            splt = line.split()
            symbol=(element(splt[1]).symbol)
            x=float(splt[3])
            y=float(splt[4])
            z=float(splt[5])
            fout.write("{0:<2} {1:>13.9f} {2:>13.9f} {3:>13.9f}\n".format(symbol,x,y,z))
