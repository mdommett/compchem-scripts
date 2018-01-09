#!/usr/bin/env python

import read_file as rf
import edit_file as ef
import handle_atoms as ha
import numpy as np
import sys

cluster=rf.read_xyz(sys.argv[1])[-1]
top=rf.read_xyz(sys.argv[2])[-1]
newcluster=[]
for i,j in enumerate(cluster):
    if j in top:
        newcluster.insert(0,j)
    else:
        newcluster.append(j)
ef.write_xyz("new-cluster.xyz",newcluster)

