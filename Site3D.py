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
            

def parse_all_pocket_files(pocketdir, resolution=0.5):
    framedata=dict()
    index=0
    print  sorted(glob.glob('%s/*pdb' % pocketdir))
    for file in sorted(glob.glob('%s/*pdb' % pocketdir)):
        framedata[index]=dict()
        sites=numpy.loadtxt(file, usecols=(3,))
        xcoor=numpy.loadtxt(file, usecols=(5,))
        ycoor=numpy.loadtxt(file, usecols=(6,))
        zcoor=numpy.loadtxt(file, usecols=(7,))
        radii=numpy.loadtxt(file, usecols=(9,))
        centers=numpy.dstack((xcoor, ycoor, zcoor))
        centers=centers.reshape(centers.shape[1],  centers.shape[2])
        framedata[index]['centers']=centers
        framedata[index]['radii']=radii
        index+=1
        #pocket_coors=[]
        #for i in xrange(len(centers)):
        #    coors=pocket_sphere_coors(radii[i], centers[i], resolution)
        #    if i==0:
        #        pocketcoors=coors
        #        i+=1
        #    else:
        #        pocketcoors=numpy.vstack((pocketcoors, coors))
        #framedata[index]=pocketcoors
        #index+=1
        #for site in set(sites):
        #    framedata[index][site]=dict()
        #    frames=numpy.where(sites==site)[0]            
        #    centers=numpy.dstack((xcoor[frames], ycoor[frames], zcoor[frames]))
        #    framedata[index][site]['centers']=centers.reshape(centers.shape[1], centers.shape[2])
        #    framedata[index][site]['radii']=radii[frames]
    return framedata


def protein_grid(allcoor, resolution=0.5):
    extra=3.0 #protein search for edge effects
    xmin=100000
    xmax=0
    ymin=100000
    ymax=0
    zmin=100000
    zmax=0
    mins=[xmin, ymin, zmin]
    maxes=[xmax, ymax, zmax]
    for frame in xrange(allcoor.shape[0]):
        for n in range(0,3):
            mintest=numpy.min(allcoor[:,:,n].flatten())
            maxtest=numpy.max(allcoor[:,:,n].flatten())
            if mintest < mins[n]:
                mins[n]=mintest
            if maxtest > maxes[n]:
                maxes[n]=maxtest
    lengths=dict()
    for n in range(0,3):
        maxes[n]=maxes[n]+resolution
        mins[n]=mins[n]-resolution
        lengths[n]=int(((round(maxes[n]))-round(mins[n])))
    box_volume=lengths[0]*lengths[1]*lengths[2]
    print "protein box volume %s angstroms^3" % box_volume
    total=max(lengths.values())
    ranges=dict()
    for n in range(0,3):
        ranges[n]=numpy.arange(int(round(mins[n])), int(round(maxes[n])), resolution)
    return ranges[0], ranges[1], ranges[2], box_volume


# Class for 3D Grid
class Site3D:
    def __init__(self, total_frames, resolution=0.5, xaxis=None, yaxis=None, zaxis=None):
        self.total_frames=total_frames
        self.dx =resolution
        self.dy =resolution
        self.dz=resolution
        self.xaxis = numpy.array(xaxis)
        self.yaxis = numpy.array(yaxis)
        self.zaxis = numpy.array(zaxis)
        
        # gives a grid that is indexed by looping over y,x,z and shape of total
        X,Y,Z=numpy.meshgrid(self.xaxis,self.yaxis,self.zaxis)
        self.pocketgrid=numpy.vstack((X.ravel(), Y.ravel(), Z.ravel())).T
        self.pocketoccup=numpy.zeros((self.pocketgrid.shape[0]))
        #self.reduced_coors=reduced_coors

    def map_sphere_occupancy_grid(self, pocketdata, cutoff=3.0, pad=None):
        # pocketdata has spheres n with center and radii
        # all coor is all protein coors
        # loop over all grid points
        # save frameoccupancies
        init=0
        for frame in range(0, self.total_frames):
            frameoccup=numpy.zeros((self.pocketgrid.shape[0]))
            distances=sp.distance.cdist(self.pocketgrid, pocketdata[frame]['centers'])
            for index in xrange(pocketdata[frame]['centers'].shape[0]):
                frames=numpy.where(distances[:,index] < pocketdata[frame]['radii'][index])[0]
                frameoccup[frames]=1
            self.pocketoccup+=frameoccup
            if init==0:
                framelog=frameoccup
                init+=1
            else:
                framelog=numpy.vstack((framelog, frameoccup))
        return framelog

    def write_pdb(self, outfile, frequency_indices):
        count=0
        atomname='DUM'
        resid=1
        resid_alpha='ABCDEFGHIJKLMOPQRSTUVWXYZ'
        resid_alpha=resid_alpha*100
        resid_count=0
        resname='DUM'
        occupancy=0.00
        beta=0.00
        ohandle=open(outfile, 'w')
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
            resid+=1
            line=format_pdb_line(atomnum, atomname, resname, resnum, xcoor,
ycoor, zcoor, occupancy, beta)
            ohandle.write(line)
            count+=1
        ohandle.close()
        return


    def write_dx(self, outfile):
        # reshape freq due to ravel in order to format for OpenDX
        reshape_freq=numpy.zeros((len(self.xaxis), len(self.yaxis),len(self.zaxis)))
        count=0
        for j in range(0, len(self.yaxis)):
            for i in range(0, len(self.xaxis)):
                for k in range(0, len(self.zaxis)):
                    reshape_freq[i,j,k]=self.pocketoccup[count]/self.total_frames
                    count+=1
        newfile=open(outfile, 'w')
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

