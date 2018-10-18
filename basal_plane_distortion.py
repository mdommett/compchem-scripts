#!/usr/bin/env python
import numpy as np
from sys import argv,exit
from sklearn.metrics.pairwise import pairwise_distances
import argparse

def SVD(X):
    """
    Singular value decomposition method.
    Source: https://gist.github.com/lambdalisue/7201028
    The fitting of a plane to a set of points can be solved by using the SVD. The
    best fitting plane can be found for minimising the mean orthogonal distance
    between each point and the plane. First the centroid of the points is taken, and
    the coordinates moved to that point. The SVD of the matrix is taken, where the final
    column of the right singluar values (V) are the a,b,c values (the normal vector) of the
    plane

    Parameters
    ----------
    X: Nx3 numpy array
    Returns
    ----------
    a,b,c,d : floats for the plane equation aX + bY +cZ + d = 0
    """
    # Find the average of points (centroid) along the columns
    C = np.average(X, axis=0)
    # Create CX vector (centroid to point) matrix
    CX = X - C
    # Singular value decomposition
    U, S, Vt = np.linalg.svd(CX)
    # The last row of Vt matrix indicate the eigenvectors of
    # smallest eigenvalues (singular values).
    N = Vt[-1]

    # Extract a, b, c, d coefficients.
    x0, y0, z0 = C
    a, b, c = N
    d = -(a * x0 + b * y0 + c * z0)
    return a, b, c, d

def file_to_matrix(infile):
    """
    Converts an xyz coordinate file to an Nx3 matrix

    Parameters
    ----------
    infile: list of N lists
    Returns
    ----------
    xyz: Nx3 matrix
    """
    xyz=np.zeros((len(infile),3))
    for i in range(len(infile)):
        xyz[i,0]=infile[i].split()[1]
        xyz[i,1]=infile[i].split()[2]
        xyz[i,2]=infile[i].split()[3]
    return xyz

def project_to_plane(a,b,c,d,point):
    """
    Projects the coordinates of a point onto a plane described by a,b,c,d

    Parameters
    ----------
    a,b,c,d: floats
    point: 1x3 np array
    Returns
    ----------
    projection: 1x3 np array
    """
    m,n,o=point[0],point[1],point[2]
    t=(-d-a*m-b*n-c*o)/((a**2)+(b**2)+(c**2))
    projection=np.array([m+(t*a),n+(t*b),o+(t*c)])
    return projection

def squared_distance_to_plane(a,b,c,d,point):
    """
    Find the mean orthogonal squared distance from a set of points (or point) to the plane

    Parameters
    ----------
    a,b,c,d: floats
    point: Nx3 np array
    Returns
    ----------
    distance: float of the mean squared othrogonal distance
    """
    x,y,z=point.T
    N=np.array([a,b,c])
    distance=np.mean(((a*x+b*y+c*z+d)**2/np.linalg.norm(N)))
    return distance

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input",help="XYZ file",type=str)
    parser.add_argument("--ref","-r",help="Reference atom",type=int)
    parser.add_argument("--basal","-b",help="Atoms to fit plane to", type=int, nargs='*')
    user_input = argv[1:]
    args = parser.parse_args(user_input)

    infile=open(args.input,'r').readlines()
    xyz=file_to_matrix(infile[2:])
    ref_no=args.ref-1
    basal_no=[int(i)-1 for i in args.basal]
    ref=xyz[ref_no,:]
    basal=xyz[basal_no,:]
    a,b,c,d=SVD(basal)
    mean_squared_error_of_plane=squared_distance_to_plane(a,b,c,d,basal)
    p=project_to_plane(a,b,c,d,ref)
    squared_distance=squared_distance_to_plane(a,b,c,d,ref)
    distance=np.sqrt(squared_distance)
    print("MSE of fit:\n{:.2f}".format(mean_squared_error_of_plane))
    print("Distance from best fitting basal plane:\n{:.2f}".format(distance))
