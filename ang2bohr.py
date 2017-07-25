#!/usr/bin/env python
from datetime import datetime
startTime = datetime.now()
import numpy as np
from sys import argv
ang2bohr=1.88973
infile=open(argv[1],"r").read().splitlines()
outfile=open(argv[1]+"-bohr","w")
for i in infile:
    ix=i.split()
    outfile.write("{0:>13.9f}   {1:>13.9f}   {2:>13.9f}   {3:>13.9}\n".format(float(ix[0])*ang2bohr,float(ix[1])*ang2bohr,float(ix[2])*ang2bohr,ix[3]))
print datetime.now() - startTime
