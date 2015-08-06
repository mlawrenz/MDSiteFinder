

import scipy.spatial as sp
import operator
import Site3D
import time
import mdtraj
import optparse
import pylab
import numpy
import glob
import os
import sys
import pickle
from scipy import interpolate, integrate

"""
MDSiteFinder
===========

** To be used on MOE SiteFinder results  **

Program requires as input (see options with -h flag):
*
*

Program output is a OpenDX file SiteFrequency.dx (Occupancy for passed in site
at default 0.5 Angstrom gridspace.
This file can be loaded into VMD:
* you can run the command: vmd -pdb reference.pdb -dx SiteFrequency.dx
* make a new representation in Isosurface, and scale the isovalues ranging from
* % opened, from 0 to 1.
If your 3D surface is not centered on your protein, check that your trajectory
is aligned to the reference structure.

"""

def main(pocketfile, trajfile, topo, outname, writedx=True):
    xcoors=numpy.loadtxt(pocketfile, usecols=(5,), ndmin=1)
    ycoors=numpy.loadtxt(pocketfile, usecols=(6,), ndmin=1)
    zcoors=numpy.loadtxt(pocketfile, usecols=(7,), ndmin=1)
    xcoors=xcoors.reshape(-1,1)
    ycoors=ycoors.reshape(-1,1)
    zcoors=zcoors.reshape(-1,1)
    tmp=numpy.hstack((xcoors, ycoors))
    opencoors=numpy.hstack((tmp, zcoors))
    pad=8.0
    # load traj
    dir=os.path.dirname(trajfile)
    print "REMOVING HYDROGENS FROM CALC (REFLECTS CUTOFF)"
    # check if you already ran the calc, load framelog if so
    start=float(time.time())
    traj=mdtraj.load(trajfile, top=topo)
    indices=[]
    for i in traj.topology.atoms:
        if 'H' not in i.name:
            indices.append(i.index)
    newcoors=10*traj.xyz[:,indices,:] #multiply by 10 bc mdtraj scalers
    cutoff=3.0
    open_frames=dict()
    count=0
    for frame in xrange(newcoors.shape[0]):
        distances=sp.distance.cdist(opencoors, newcoors[frame])
        open=numpy.where(distances<3.0)[0]
        if open.size < 20:
            print frame, open.size
        if not open.size:
            count+=1
    tally=float(count)/newcoors.shape[0]
    print "open for %s" % tally
            

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-o', '--outname', dest='outname',
                      help='output dx name')
    parser.add_option('-y', '--topo', dest='topo',
                      help='topology')
    parser.add_option('-t', '--trajfile', dest='trajfile',
                      help='MD traj file')
    parser.add_option('-p', '--pocketfile', dest='pocketfile',
                      help='protein pocket file output from MOE SiteFinder')
    #parser.add_option('-f', action="store_true", dest="writefree", help='perform free energy calc with P-L COM distances')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    #if options.writefree==True:
    main(pocketfile=options.pocketfile, trajfile=options.trajfile, topo=options.topo, outname=options.outname)


