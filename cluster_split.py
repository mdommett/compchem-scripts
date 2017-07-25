#!/usr/bin/env python

from sys import argv

infile =  open(argv[1],"r").read().splitlines()
out  = str(argv[2])
natoms_total = int(infile[0].split()[0])
natoms_permolecule = int(raw_input("How many atoms per molecule?\n"))
nmolecules = int(natoms_total/natoms_permolecule)
print "If there are {0} atoms in the molecule, then there are {1} molecules in total".format(natoms_permolecule, nmolecules)
answer = str(raw_input("Would you like the output to separate xyz files? [Y/N]\n"))


if answer == "Y" or answer == "y" or answer == "yes":

    print "Printing {} xyz files... \n".format(nmolecules)

    atomcount = 1
    for i in range(nmolecules):
        outfile = open(out+"."+str(i+1)+".xyz","w")
        for line in infile[2+i*natoms_permolecule:2+i*natoms_permolecule+natoms_permolecule]:
            if atomcount == 1:
                outfile.write( "{}\n\n".format(natoms_permolecule))
    
            if atomcount < natoms_permolecule:
                outfile.write("{}\n".format(line))
                atomcount +=1
        
            else: 
                outfile.write("{}\n".format(line))
                atomcount = 1
        print "Molecule {} printed".format(i+1)
        outfile.close()

elif answer == "N" or answer == "n" or answer == "no":
    outfile = open(out+".xyz","w")
    atomcount = 1
    for line in infile[2:]:
        if atomcount == 1:
            outfile.write( "{}\n\n".format(natoms_permolecule))
    
        if atomcount < natoms_permolecule:
            outfile.write("{}\n".format(line))
            atomcount +=1
        
        else: 
            outfile.write("{}\n".format(line))
            atomcount = 1

    outfile.close()
