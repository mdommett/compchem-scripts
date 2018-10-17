#!/usr/bin/env python
import numpy as np
from sys import argv,exit
from sklearn.metrics.pairwise import pairwise_distances

def SVD(X):
    """
    Singular value decomposition method.
    Source: https://gist.github.com/lambdalisue/7201028
    """
    # Find the average of points (centroid) along the columns
    C = np.average(X, axis=0)
    # Create CX vector (centroid to point) matrix
    CX = X - C
    # Singular value decomposition
    U, S, Vt = np.linalg.svd(CX)
    # The last row of V matrix indicate the eigenvectors of
    # smallest eigenvalues (singular values).
    N = Vt[-1]

    # Extract a, b, c, d coefficients.
    x0, y0, z0 = C
    a, b, c = N
    d = -(a * x0 + b * y0 + c * z0)
    return a, b, c, d

def file_to_matrix(infile):
    xyz=np.zeros((len(infile),3))
    for i in range(len(infile)):
        xyz[i,0]=infile[i].split()[1]
        xyz[i,1]=infile[i].split()[2]
        xyz[i,2]=infile[i].split()[3]
    return xyz

def project_to_plane(a,b,c,centroid,point):
    x,y,z=point[0],point[1],point[2]
    d,e,f=centroid[0],centroid[1],centroid[2]
    t=((a*d)-(a*x)+(b*e)-(b*y)+(c*f)-(c*z))/((a**2)+(b**2)+(c**2))
    projection=np.array([x+(t*a),y+(t*b),z+(t*c)])
    return projection

infile=open(argv[1],'r').readlines()
xyz=file_to_matrix(infile[2:])
ref_no=int(argv[2])-1
basal_no=[int(i)-1 for i in argv[3:7]]
ref=xyz[ref_no,:]
basal=xyz[basal_no,:]
a,b,c,d=SVD(basal)
cent=np.average(basal, axis=0)
p=project_to_plane(a,b,c,cent,ref)
dist=np.linalg.norm(ref-p)
print("Distance from best fitting basal plane\n={:.2f}".format(dist))
