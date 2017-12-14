#!/usr/bin/env python
from os import system,remove
iterlog = open("iter.log").read().splitlines()

s0 = []
s1 = []

for line in iterlog[1:]:
    splt = line.split()
    s1.append(float(splt[6])*27.2114)
    s0.append(float(splt[7])*27.2114)


with open("data","w") as d:
    for i in range(len(s1)):
        d.write("{0}  {1} {2}\n".format(i+1,s0[i]-s0[0],s1[i]-s0[0]))
    d.close()
        
with open("plot","w") as p:
    p.write("set xlabel \"Step\"\nset ylabel \"Energy (eV)\"\n")
    p.write("plot 'data' using 1:2 title 'S0' with lines linecolor black linewidth 2,\\\n")
    p.write("     'data' using 1:3 title 'S1' with lines linecolor 'red' linewidth 2")
    p.close()

system("gnuplot -persist plot")
system("rm -f plot data")
