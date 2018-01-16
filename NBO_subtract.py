#!/usr/bin/env python2

import numpy as np
from itertools import chain
from sys import argv
import sys
from periodic import element

set_1 = open(argv[1],"r").read().splitlines()
set_2 = open(argv[2],"r").read().splitlines()

if len(set_1)==len(set_2):
    
    for i in range(len(set_1)):
        
        print float(set_1[i].split()[2])-float(set_2[i].split()[2])
                  
