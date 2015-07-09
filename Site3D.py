import sys
import scipy.spatial as sp
import math
import pickle
import multiprocessing
import optparse
import pylab
from numpy import *
import glob
import os

# Helper Functions
def eval_distance(mapped_state_distances, cutoff):
    # assumes mapped states and state_distances
    mapped_cutoff_states=[]
    for state in xrange(len(mapped_state_distances)):
        # if mindist to protein > cutoff, add to unbound frames 
        if mapped_state_distances[state] > cutoff:
            mapped_cutoff_states.append(state)
    mapped_cutoff_states=array([int(i) for i in mapped_cutoff_states])
    return mapped_cutoff_states



def cart2sph(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r = math.sqrt(XsqPlusYsq + z**2)               # r
    elev = math.acos((z/math.sqrt(XsqPlusYsq)))       # theta
    az = math.atan((y/x))                           # phi
    return r,  math.degrees(az), math.degrees(elev)


def sph2cart(r, az, elev):
    elev=math.radians(elev)
    az=math.radians(az)
    x=r*math.cos(az)*sin(elev)
    y=r*math.sin(az)*sin(elev)
    z=r*math.cos(elev)
    return x,y,z
    
def pocket_sphere_coors(radius, center, resolution):
    sphere_cart_coors=[]
    for dx in arange(-radius,radius, resolution):
        for dy in arange(-radius,radius, resolution):
            for dz in arange(-radius,radius, resolution):
                sphere_cart_coors.append((center[0]+dx, center[1]+dy, center[2]+dz))
    return sphere_cart_coors     
            

def parse_pocket_file(pocket_file, resolution=0.5):
    pocket_data=dict()
    fhandle=open(pocket_file)
    n=0
    for line in fhandle.readlines():
        if len(line) < 2:
            break
        pocket_data[n]=dict()
        atom=line[0:3]
        atomnum=line[6:10]
        atomname=line[12:15]
        resname=line[16:20]
        chain=line[21]
        resnum=line[22:25]
        xcoor=float(line[30:37])
        ycoor=float(line[38:45])
        zcoor=float(line[46:53])
        # official
        occupancy=line[56:61]
        beta=line[62:67]
        print "sphere %s radius %s" % (n, beta)
        # may want this flex
        #occupancy=line[67:70]
        #beta=line[54:59]
        radius=float(beta)
        center=(xcoor, ycoor, zcoor)
        pocket_data[n]['center']=center
        pocket_data[n]['radius']=radius
        n+=1
    return pocket_data


def get_pocket_minmax(pocket_data, allcoor, pad=3.0, resolution=0.5):
    reduced_coors=dict()
    xmin=100000
    xmax=0
    ymin=100000
    ymax=0
    zmin=100000
    zmax=0
    mins=[xmin, ymin, zmin]
    maxes=[xmax, ymax, zmax]
    for sphere in sorted(pocket_data.keys()):
        for frame in xrange(allcoor.shape[0]):
            for n in range(0,3):
                new=allcoor[frame][abs(allcoor[frame][:,n]-pocket_data[sphere]['center'][n]) < pocket_data[sphere]['radius']]    
                if min(new[:,n]) < mins[n]:
                    mins[n]=min(new[:,n])
                elif max(new[:,n]) > maxes[n]:
                    maxes[n]=max(new[:,n])
    lengths=dict()
    for n in range(0,3):
        maxes[n]=maxes[n]+pad
        mins[n]=mins[n]-pad
        lengths[n]=int(((round(maxes[n]))-round(mins[n])))
    box_volume=lengths[0]*lengths[1]*lengths[2]
    print "pocket box volume %s angstroms^3" % box_volume
    total=max(lengths.values())
    ranges=dict()
    for n in range(0,3):
        ranges[n]=arange(int(round(mins[n])), int(round(maxes[n])), resolution)
    reduced_coors=dict()
    for frame in xrange(allcoor.shape[0]):
        matches=where((allcoor[frame][:,0]>=mins[0])&(allcoor[frame][:,0]<maxes[0])&(allcoor[frame][:,1] >= mins[1])&(allcoor[frame][:,1] < maxes[1])&(allcoor[frame][:,2] >= mins[2])&(allcoor[frame][:,2] < maxes[2]))
        reduced_coors[frame]=allcoor[frame][matches]
    return reduced_coors, ranges[0], ranges[1], ranges[2], box_volume




# Class for 3D Grid
class Site3D:
    def __init__(self, resolution=0.5, pops=None, xaxis=None, yaxis=None, zaxis=None, grid=None):
        self.pops=pops
        self.dx =resolution
        self.dy =resolution
        self.dz=resolution
        self.xaxis = array(xaxis)
        self.yaxis = array(yaxis)
        self.zaxis = array(zaxis)
        self.tol=self.dx/2.0
        X,Y,Z=meshgrid(self.xaxis,self.yaxis,self.zaxis)
        self.pocketgrid=vstack((X.ravel(), Y.ravel(), Z.ravel())).T

 
    def map_sphere_occupancy_grid(self, pocketdata, reduced_coor, cutoff=1.4, pad=None):
        # pocketdata has spheres n with center and radii
        # all coor is all protein coors
        pocketoccup=zeros((len(self.xaxis), len(self.yaxis), len(self.zaxis)))
        # loop over all grid points
        for frame in xrange(len(reduced_coor.keys())):
            frameoccup=ones((len(self.xaxis), len(self.yaxis), len(self.zaxis)))
            distances=sp.distance.cdist(self.pocketgrid, reduced_coor[frame])
            # array of gridpoint index, protein coor
            occupied=where(distances< cutoff)
            if occupied[0].size:
                for index in occupied[0]:
                    coor=self.pocketgrid[index]
                    testarray=array([38.0, 47.5, 25.0])
                    i=where(self.xaxis==coor[0])[0]
                    j=where(self.yaxis==coor[1])[0]
                    k=where(self.zaxis==coor[2])[0]
                    frameoccup[i,j,k]=0
            pocketoccup+=frameoccup
        return pocketoccup


    def write_pdb(self, dir, matrix, frequency):
        count=0
        atomname='DUM'
        resid=1
        resname='DUM'
        occupancy=0.00
        beta=0.00
        testarray=array([38.0, 47.5, 25.0])
        ohandle=open('%s/pocketgrid_open%0.1f.pdb' % (dir, frequency), 'w')
        for i in xrange(len(self.xaxis)):    
            for j in xrange(len(self.yaxis)):    
                for k in xrange(len(self.zaxis)):    
                    #if i==18 and j==19 and k==4:
                    if matrix[i,j,k]==frequency:
                        xcoor=self.xaxis[i]
                        ycoor=self.yaxis[j]
                        zcoor=self.zaxis[k]
                        atomnum=count+1
                        line='ATOM{0: >7}{1: >4} {2:>4} X{3:>4}    {4: >8.3f}{5: >8.3f}{6: >8.3f}{7: >6.2f}{8: >6.2f} \n'.format(atomnum, atomname, resname, resid, xcoor, ycoor, zcoor, occupancy, beta)
                        ohandle.write(line)
                        count+=1
                    else:
                        count+=1
        ohandle.close()
        return

    def write_dx(self, GD, dir, filename):
        newfile=open('%s/%s_sitefrequency.dx' % (dir, filename), 'w')
        newfile.write('# Data calculated Pocket open frequency\n')
        newfile.write('object 1 class gridpositions counts %s %s %s\n' % (GD.shape[0], GD.shape[1], GD.shape[2]))
        newfile.write('origin %s %s %s\n' % (self.xaxis[0], self.yaxis[0], self.zaxis[0]))
        newfile.write('delta %s 0 0\n' % self.dx)
        newfile.write('delta 0 %s 0\n' % self.dy)
        newfile.write('delta 0 0 %s\n' % self.dz)
        newfile.write('object 2 class gridconnections counts %s %s %s\n' % (GD.shape[0], GD.shape[1], GD.shape[2]))
        newfile.write('object 3 class array type double rank 0 items %s data follows\n' % (GD.shape[0]*GD.shape[1]*GD.shape[2]))
        intergrid=zeros((GD.shape[0], GD.shape[1], GD.shape[2]))
        count=0
        for i in range(0, GD.shape[0]):
            for j in range(0, GD.shape[1]):
                for k in range(0, GD.shape[2]):
                    if count==2:
                        if GD[i][j][k]==0:
                            newfile.write('%s\n' % int(GD[i][j][k]))
                        else:
                            newfile.write('%s\n' % GD[i][j][k])
                        count=0
                    else:
                        if GD[i][j][k]==0:
                            newfile.write('%s\t' % int(GD[i][j][k]))
                        else:
                            newfile.write('%s\t' % GD[i][j][k])
                        count+=1
        newfile.write('\nobject "ligand free energy" class field')
        newfile.close()

    def pmfvolume(self, GD):
        sum=0
        oldx=self.xaxis[0]
        oldval=0
        for i in range(0, len(self.xaxis)):
            oldy=self.yaxis[0]
            for j in range(0, len(self.yaxis)):
                oldz=self.zaxis[0]
                for k in range(0, len(self.zaxis)):
                    if abs(self.zaxis[k]-oldz) !=0:
                        if GD[i,j, k] != 0:
                            sum+=GD[i,j, k]*abs(self.xaxis[i]-oldx)*abs(self.yaxis[j]-oldy)*abs(self.zaxis[k]-oldz)
                    oldz=self.zaxis[k]
                oldy=self.yaxis[j]
            oldx=self.xaxis[i]
        return float(sum)
 


