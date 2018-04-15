#!/usr/bin/env python

import sys
with open(sys.argv[1],"r") as f:
    for i,line in enumerate(f):
        if " [GEOMETRIES] (XYZ)" in line:
            natoms=int(f[i+1].split()[0])
        break
print(i)
print(natoms)
