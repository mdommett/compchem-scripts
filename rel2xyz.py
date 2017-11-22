#!/usr/bin/env python
import sys
import read_file as rf
import edit_file as ef
import handle_atoms as ha
import numpy as np

QE_f = sys.argv[1]
atoms=rf.read_qe(QE_f)
ef.write_xyz(QE_f[:-4]+"-opt.xyz",atoms)

vectors=np.loadtxt("vectors")
mega=ha.make_mega_cell(atoms,2,2,2,vectors,alt_multi=True)
ef.write_xyz(QE_f[:-4]+"-mega.xyz",mega)
