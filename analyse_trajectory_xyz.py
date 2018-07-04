import numpy as np
from sys import argv,exit
import argparse

if __name__ == "__main__":
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--dihedral", help="Dihedral angle between atoms A,B,C,D", nargs=4))
    parser.add_argument("-a","--angle",help="Angle between atoms A,B,C",nargs=3)
    parser.add_argument("-p","--pyramidal",help="pyramidal angle between bond AB plane CD",
    nargs=4)
