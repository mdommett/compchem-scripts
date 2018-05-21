#!/usr/bin/env python
import read_file as rf
from sys import argv
import numpy as np
import edit_file as ef
from sklearn.metrics.pairwise import pairwise_distances

def coordinate_matrix(atoms):
    coords=np.zeros((len(atoms),3))
    for i in range(len(atoms)):
        coords[i,0]=atoms[i].x
        coords[i,1]=atoms[i].y
        coords[i,2]=atoms[i].z
    return coords

def long_axis(distance_matrix):
    # Gets the maximum argument, which is a flattened index, and gives it back as a 2d index
    return np.array(np.unravel_index(np.argmax(distance_matrix, axis=None), distance_matrix.shape))

def rotation_matrix(axis, theta):

    """Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians."""

    axis = np.asarray(axis)
    axis = axis/np.sqrt(np.dot(axis, axis))

    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def rotate_x(v,theta):
    rot=np.array([[1,0,0],
                [0,np.cos(theta),-np.sin(theta)],
                [0,np.sin(theta),np.cos(theta)]])
    return np.dot(rot,v)

def rotate(axis,theta):
    axis = axis/np.sqrt(np.dot(axis, axis))
    a11=np.cos(theta)+(axis[0]**2*(1-np.cos(theta)))
    a12=axis[0]*axis[1]*(1-np.cos(theta))-axis[2]*np.sin(theta)
    a13=axis[0]*axis[2]*(1-np.cos(theta))+axis[1]*np.sin(theta)
    a21=axis[1]*axis[0]*(1-np.cos(theta))+axis[2]*np.sin(theta)
    a22=np.cos(theta)+(axis[1]**2*(1-np.cos(theta)))
    a23=(axis[1]*axis[2]*(1-np.cos(theta)))-axis[0]*np.sin(theta)
    a31=(axis[2]*axis[0]*(1-np.cos(theta)))-axis[1]*np.sin(theta)
    a32=(axis[2]*axis[1]*(1-np.cos(theta)))+axis[0]*np.sin(theta)
    a33=np.cos(theta)+(axis[2]**2*(1-np.cos(theta)))
    R=np.array([[a11,a12,a13],
                [a21,a22,a23],
                [a31,a32,a33]])
    return R

def align_z(coordinates,point1,point2):
    """first need to move a point to the origin
    then rotate into XZ plane
    """
    # translate chosen point to origin
    translated=coordinates-point1
    # long axis vector
    vector=point2-point1
    # unit_vector  of vector
    unit_vector=vector/np.sqrt(np.dot(vector, vector))
    # z-axis vector
    z=np.array([0,0,1])
    # cos(theta) of angle to between current vector and z-axis
    rot_angle=(np.arccos(np.dot(z,vector)/(np.linalg.norm(z)*np.linalg.norm(vector))))%(2*np.pi)# modulo covers for negative angle
    #print(np.degrees(rot_angle))
    #  rotation axis that connects vector and z-axis
    rot_axis=np.cross(z,vector)/np.linalg.norm(np.cross(z,vector))
    x,y,z=rot_axis
    c=np.cos(rot_angle)
    #print(c)
    C=1-c

    s=np.sqrt(1-(C*C))
    # https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
    # rotation matrix :
    rmat = np.array([[x*x*C+c,x*y*C-z*s,x*z*C+y*s],
                    [y*x*C+z*s,y*y*C+c,y*z*C-x*s],
                    [z*x*C-y*s,z*y*C+x*s,z*z*C+c]])
    new_coords=np.array([np.dot(i,rmat) for i in translated])
    return new_coords

def translate(coordinates,point):
    return coordinates-point

def magnitude(vector):
    return np.linalg.norm(vector)

def costheta(vector1,vector2):
    return np.dot(vector1,vector2)/((magnitude(vector1)*(magnitude(vector2))))

def distance(point1,point2):
    return np.sqrt(np.sum([(b - a)**2 for (a, b) in zip(point1,point2)]))

def dihedral_angle(p0,p1,p2,p3):
    # q0 = p1-p0 vector connecting p0 to p1
    q0 = p1 - p0
    q1 = p2 - p1
    q2 = p3 - p2

    q0q1=np.cross(q0,q1) # vector normal to plane p0p1p2
    q1q2=np.cross(q1,q2) # normal to plane p1p2p3

    n0=q0q1/magnitude(q0q1) # unit normal
    n1=q1q2/magnitude(q1q2)

    u0=n1 # orthogonal unit vectors
    u2=q1/magnitude(q1)
    u1=np.cross(u0,u2)

    cos_theta=np.dot(n0,u0)
    sin_theta=np.dot(n0,u1)
    theta=-np.arctan2(sin_theta,cos_theta)
    return np.degrees(theta)

def get_snap_vector((mon_0_0,mon_0_1),(mon_1_0,mon_1_1),mon_1):
    """
    Pass two coordinates per molecule (the extreme coordinates for the long axis)
    to generate the best vector to overlap the molecules
    the molecule with points 3 and 4 is translated
    """

    # intermolecular vectors between the extreme points of each molecule
    extreme_vectors=np.array([mon_1_0-mon_0_0,mon_1_0-mon_0_1,
    mon_1_1-mon_0_0,mon_1_1-mon_0_1])
    # translate each extreme point in mon_1 to each extreme point in molecule 0
    mon_1_translations=np.array([np.array((mon_1_0,mon_1_1)-vector) for vector in extreme_vectors])
    # new intermolecular extreme vectors
    new_extreme_vectors=np.array(mon_1_translations-(mon_0_0,mon_0_1))
    # best translation is the one that minimises the distances between the intermolecular extremes
    snap_vector=extreme_vectors[np.array([abs(magnitude(vector)) for vector in new_extreme_vectors]).argmin()]
    return snap_vector

if __name__=='__main__':

    dimer=rf.read_pos(argv[1])
    natoms_dimer=len(dimer)
    natoms_monomer=natoms_dimer/2
    dimer_array=np.zeros((natoms_dimer,3))

    mon_0=coordinate_matrix(dimer[0:natoms_monomer])
    mon_1=coordinate_matrix(dimer[natoms_monomer:])
    dim=np.concatenate((mon_0,mon_1),axis=0)
    mon_0_distances=pairwise_distances(mon_0)
    mon_1_distances=pairwise_distances(mon_1)
    mon_0_extreme_atoms=long_axis(mon_0_distances)
    mon_1_extreme_atoms=long_axis(mon_1_distances)#+natoms_monomer # back to index of original dimer
    mon_0_extreme_coords=np.array([mon_0[mon_0_extreme_atoms[0]],mon_0[mon_0_extreme_atoms[1]]])
    mon_1_extreme_coords=np.array([mon_1[mon_1_extreme_atoms[0]],mon_1[mon_1_extreme_atoms[1]]])

    mon_0_long_axis=mon_0_extreme_coords[1]-mon_0_extreme_coords[0]
    mon_1_long_axis=mon_1_extreme_coords[1]-mon_1_extreme_coords[0]
    #print(dihedral_angle(mon_0[mon_0_extremes[0]],mon_0[mon_0_extremes[1]],mon_1[mon_1_extremes[0]],mon_1[mon_1_extremes[1]]))

    # vectors between the extreme points of each molecule
    # in matrix form to get every combination
    snap_vector=get_snap_vector(mon_0_extreme_coords,mon_1_extreme_coords,mon_1)
    magnitude_snap_vector=magnitude(snap_vector)
    mon_1_snapped=mon_1-snap_vector
    #print(mon_1_snapped[mon_1_extreme_atoms[0]])
    snapped_dim=np.concatenate((mon_0,mon_1_snapped))
    for i ,j in enumerate(dimer):
        j.x,j.y,j.z=snapped_dim[i]
    ef.write_xyz("{}_translated.xyz".format(argv[1][:-4]),dimer)
    """
    for no,vector in enumerate(extreme_vectors):
        #print(magnitude(vector))
        #theta = np.radians(angle)
        #rotated=np.array([np.dot(rotate(np.array([0,0,1]),theta), i) for i in new_coords])
        mon_1_snapped=mon_1-vector
        print(np.mean(pairwise_distances(mon_0,mon_1_snapped)))
        snapped_dim=np.concatenate((mon_0,mon_1_snapped))
        #print(distance(snapped_dim[19],snapped_dim[7]))
        #translated_dim=align_z(snapped_dim,mon_0_extreme_coords[0],mon_0_extreme_coords[1])
        for i,j in enumerate(dimer):
            j.x,j.y,j.z=snapped_dim[i]
        ef.write_xyz("{}.xyz".format(no),dimer)"""
