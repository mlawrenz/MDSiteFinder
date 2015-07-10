import sys
import scipy.spatial as sp
import math
import pickle
import multiprocessing
import optparse
import pylab
import numpy
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
    mapped_cutoff_states=numpy.array([int(i) for i in mapped_cutoff_states])
    return mapped_cutoff_states

def format_pdb_line(atomnum, atomname, resname, resnum, xcoor, ycoor, zcoor, occupancy, beta):
    line='ATOM{0: >7}{1: >4} {2:>4} X{3:>4}    {4: >8.3f}{5: >8.3f}{6: >8.3f}{7: >6.2f}{8: >6.2f} \n'.format(atomnum, atomname, resname, resnum, xcoor, ycoor, zcoor, occupancy, beta)
    return line


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
    for dx in numpy.arange(-radius,radius, resolution):
        for dy in numpy.arange(-radius,radius, resolution):
            for dz in numpy.arange(-radius,radius, resolution):
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
    extra=3.0 #protein search for edge effects
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
        ranges[n]=numpy.arange(int(round(mins[n])), int(round(maxes[n])), resolution)
    reduced_coors=dict()
    for frame in xrange(allcoor.shape[0]):
        # include protein just outside of box
        #matches=numpy.where((allcoor[frame][:,0]>=mins[0])&(allcoor[frame][:,0]< maxes[0])&(allcoor[frame][:,1] >= mins[1])&(allcoor[frame][:,1] < maxes[1])&(allcoor[frame][:,2] >= mins[2])&(allcoor[frame][:,2] < maxes[2]))
        matches=numpy.where((allcoor[frame][:,0]>=(mins[0]-extra))&(allcoor[frame][:,0]< (maxes[0]+extra))&(allcoor[frame][:,1] >= (mins[1]-extra))&(allcoor[frame][:,1] < (maxes[1]+extra))&(allcoor[frame][:,2] >= (mins[2]-extra))&(allcoor[frame][:,2] < (maxes[2]+extra)))
        reduced_coors[frame]=allcoor[frame][matches]
    return reduced_coors, ranges[0], ranges[1], ranges[2], box_volume




# Class for 3D Grid
class Site3D:
    def __init__(self, total_frames, resolution=0.5, pops=None, xaxis=None, yaxis=None, zaxis=None, reduced_coors=None):
        self.total_frames=total_frames
        self.pops=pops
        self.dx =resolution
        self.dy =resolution
        self.dz=resolution
        self.xaxis = numpy.array(xaxis)
        self.yaxis = numpy.array(yaxis)
        self.zaxis = numpy.array(zaxis)
        self.tol=self.dx/2.0
        X,Y,Z=numpy.meshgrid(self.xaxis,self.yaxis,self.zaxis)
        # gives a grid that is indexed by looping over y,x,z and shape of total
        # gridpoints
        self.pocketgrid=numpy.vstack((X.ravel(), Y.ravel(), Z.ravel())).T
        # initialize tally with zero, will reshape it later
        self.pocketoccup=numpy.zeros((self.pocketgrid.shape[0]))
        self.reduced_coors=reduced_coors

    def map_sphere_occupancy_grid(self, pocketdata, cutoff=3.0, pad=None):
        # pocketdata has spheres n with center and radii
        # all coor is all protein coors
        # loop over all grid points
        # save frameoccupancies
        init=0
        for frame in xrange(len(self.reduced_coors.keys())):
            frameoccup=numpy.ones((self.pocketgrid.shape[0]))
            distances=sp.distance.cdist(self.pocketgrid, self.reduced_coors[frame])
            # array of gridpoint index, protein coor
            # shape of occupied is gridpoint, protein atom close
            # count gridpoints with at least 1 protein atom close
            occupied=numpy.where(distances< cutoff)
            unique_occupied=numpy.unique(occupied[0]) 
            frameoccup[unique_occupied]=0
            check_solvated=numpy.where(distances<5.0)
            ref=numpy.arange(0, distances.shape[0])
            unique_solvated=numpy.unique(check_solvated[0])
            solvated_inds=numpy.setdiff1d(ref, unique_solvated) 
            frameoccup[solvated_inds]=0
            self.pocketoccup+=frameoccup
            if init==0:
                framelog=frameoccup
                init+=1
            else:
                framelog=numpy.vstack((framelog, frameoccup))
            # turn off check solvated search, can do visualization post
            # processing
        return framelog


    def write_pdb(self, dir, outname, frequency_indices, freq_val, buffer=1.0):
        count=0
        atomname='DUM'
        resid=1
        resid_alpha='ABCDEFGHIJKLMOPQRSTUVWXYZ'
        resid_alpha=resid_alpha*100
        resid_count=0
        resname='DUM'
        occupancy=0.00
        beta=0.00
        ohandle=open('%s/%s_open%0.1f.pdb' % (dir, outname, freq_val), 'w')
        for index in frequency_indices:
            xcoor=self.pocketgrid[index][0]
            ycoor=self.pocketgrid[index][1]
            zcoor=self.pocketgrid[index][2]
            atomnum=count+1
            # keep spheres as diff residues, but can't be over 4 digits (add
            # alpha)
            if len(str(resid))>3:
                resid=1
                resid_count+=1
            resnum='%s%s' % (resid_alpha[resid_count], resid)
            print resnum
            resid+=1
            line=format_pdb_line(atomnum, atomname, resname, resnum, xcoor, ycoor, zcoor, occupancy, beta)
            ohandle.write(line)
            count+=1
        ohandle.close()
        return

    def write_dx(self, dir, filename):
        # reshape freq due to ravel in order to format for OpenDX
        reshape_freq=numpy.zeros((len(self.xaxis), len(self.yaxis),len(self.zaxis)))
        count=0
        for j in range(0, len(self.yaxis)):
            for i in range(0, len(self.zaxis)):
                for k in range(0, len(self.zaxis)):
                    reshape_freq[i,j,k]=self.pocketoccup[count]/self.total_frames
                    count+=1
        newfile=open('%s/%s_sitefrequency.dx' % (dir, filename), 'w')
        newfile.write('# Data calculated Pocket open frequency\n')
        newfile.write('object 1 class gridpositions counts %s %s %s\n' % (reshape_freq.shape[0], reshape_freq.shape[1], reshape_freq.shape[2]))
        newfile.write('origin %s %s %s\n' % (self.xaxis[0], self.yaxis[0], self.zaxis[0]))
        newfile.write('delta %s 0 0\n' % self.dx)
        newfile.write('delta 0 %s 0\n' % self.dy)
        newfile.write('delta 0 0 %s\n' % self.dz)
        newfile.write('object 2 class gridconnections counts %s %s %s\n' % (reshape_freq.shape[0], reshape_freq.shape[1], reshape_freq.shape[2]))
        newfile.write('object 3 class array type double rank 0 items %s data follows\n' % (reshape_freq.shape[0]*reshape_freq.shape[1]*reshape_freq.shape[2]))
        intergrid=numpy.zeros((reshape_freq.shape[0], reshape_freq.shape[1], reshape_freq.shape[2]))
        count=0
        for i in range(0, reshape_freq.shape[0]):
            for j in range(0, reshape_freq.shape[1]):
                for k in range(0, reshape_freq.shape[2]):
                    if count==2:
                        if reshape_freq[i][j][k]==0:
                            newfile.write('%s\n' % int(reshape_freq[i][j][k]))
                        else:
                            newfile.write('%s\n' % reshape_freq[i][j][k])
                        count=0
                    else:
                        if reshape_freq[i][j][k]==0:
                            newfile.write('%s\t' % int(reshape_freq[i][j][k]))
                        else:
                            newfile.write('%s\t' % reshape_freq[i][j][k])
                        count+=1
        newfile.write('\nobject "ligand free energy" class field')
        newfile.close()

 


