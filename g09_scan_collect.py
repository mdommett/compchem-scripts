#!/usr/bin/env python

from sys import argv,exit
from periodic import element
import matplotlib.pyplot as plt
from datetime import datetime
startTime = datetime.now()

finp = open(argv[1],"r").read().splitlines()
fout = open(argv[1]+"-energies","w")

lines = filter(None, (line.rstrip() for line in finp))
hold_scf=[]
hold_s1=[]

for j in lines:

    if  " SCF Done:" in j:
        scf=float(j.split()[4])

    if " Total Energy" in j:
        s1=float(j.split()[4])

    if j==" Optimization completed.":
        hold_scf.append(scf)
        hold_s1.append(s1)
        del scf,s1

scf_final=[(i-hold_scf[0])*27.2114 for i in hold_scf]
s1_final=[(i-hold_scf[0])*27.2114 for i in hold_s1]

if len(scf_final)==len(s1_final):
        for i in range(len(scf_final)):
            fout.write("{0:>2.3f}  {1:>2.3f}\n".format(scf_final[i],s1_final[i]))
else:
    sys.exit("\n\nCheck the G09 output file, there is a problem\n\n")
print datetime.now() - startTime
plotgraph = raw_input("\nDo you want to plot a graph? Yes/No : ")
x=range(len(s1_final))
if plotgraph in ("yes","YES","Yes","y","Y"):
    fig,ax = plt.subplots()
    ax.plot(x,scf_final,'o-b',linewidth=2)
    ax.plot(x,s1_final,'o-r',linewidth=2)
    plt.xlabel('X', fontsize=12)
    plt.ylabel('Energy (eV)', fontsize=12)
    plt.show()
elif plotgraph in ("no","NO","No","n","N"):
    print "\nOkay, I will just print the energies\n"
else:
    print "\nI did not understand your answer. Please enter 'Yes' or 'No'\n"
