#!/usr/bin/env python

from sys import argv,exit
from collections import defaultdict
import numpy as np

template = open(argv[1],"r").read().splitlines()
split_template = (template[0]).split()
natoms_template = int(split_template[0])
atm1 = int(argv[2])
lin1 = template[atm1]

cols = lin1.split()
colslen = len(cols)

cons = []
for con in cols[6:colslen]:
    
    con = int(con)
    cons.append(con)
    
for con in cons:
    cols = template[con].split()
    colslen = len(cols)
    ne
