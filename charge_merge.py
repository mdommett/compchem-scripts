#!/usr/bin/env python
import numpy as np
from sys import argv

charges=open(argv[1],"r").read().splitlines()
coords=open(argv[2],"r").read().splitlines()


x = [i.split()[1] for i in coords[2:]]
y = [i.split()[2] for i in coords[2:]]
z = [i.split()[3] for i in coords[2:]]

with open(argv[3],"w") as o:
	for i in range(len(charges)):
		o.write("{0:>13.9f} {1:>13.9f} {2:>13.9f} {3:>9.6f}\n".format(float(x[i]),float(y[i]),float(z[i]),float(charges[i])))
