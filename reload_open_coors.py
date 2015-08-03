test
test=numpy.ones((2,2,2))
test
test[0,1]=0
test
test+=test
test
test[0,0]=0
test
test+=test
test
test=numpy.ones((2,2,2))
test
test[0,1]=0
test
test=numpy.ones((2,2,2))
orig=numpy.ones((2,2,2))
new=numpy.ones((2,2,2))
new[0,1]=0
new
test
orig=orig+new
orig
new=numpy.ones((2,2,2))
new[0,1]=0
new[0,0]=0
new
orig=orig+new
orig
import mdtraj
mdtraj??
traj=mdtraj.load('test/aln-holo.pdb')
traj.??
traj.?
traj.n_atoms
traj.xyz[0][1327]
maxgrid=20
import numpy
import scipy.spatial as sp
x=numpy.linspace(0,20,1)
print x
x=numpy.linspace(0,20,20)
print x
x=numpy.linspace(0,20)
x
x=numpy.linspace(0,20,20)
x
x=numpy.linspace(0,20,19)
print x
x
x=numpy.linspace(0,20,20)
print x
y=numpy.linspace(0,20,20)
z=numpy.linspace(0,20,20)
center=numpy.array(4,4,4)
cutoff=2
x=x[abs(x-4)<cutoff]
y=x[abs(y-4)<cutoff]
y=y[abs(y-4)<cutoff]
z=z[abs(z-4)<cutoff]
print x
X,Y,Z=numpy.meshgrid(x,y,z)
print X
data=numpy.vstack((X.ravel(), Y.ravel(), Z.ravel())).T
print data.shape
data
distance=sp.distance.dist(data, center, reshape(0,-1)).ravel()
distance=sp.distance.cdist(data, center, reshape(0,-1)).ravel()
distance=sp.distance.cdist(data, center.reshape(0,-1)).ravel()
print center
center=numpy.array(4,4,4)
center=numpy.array([4,4,4])
distance=sp.distance.cdist(data, center.reshape(0,-1)).ravel()
distance=sp.distance.cdist(data, center.reshape(1,-1)).ravel()
points_in_sphere=data[distance < cutoff]
print points_in_sphere
import readlinle
import readline
readline.write_history_file('plot3D.py')
import numpy
numpy.ones((100,100,1000))
test=numpy.ones((100,100,1000))
print test[0]
print test[0].shape
print test[1].shape
print test[2].shape
print test[0]=1
print test[0,:,:]=1
print test[0,0:1000]=1
print test[0]
print test[0].shape
test[0]
import mdtraj
test=mdtraj.load('test/aln-apo.pdb')
test.xyz.shape
test.xyz[0]
test.xyz[0].shape
test.xyz[0][0]
test.xyz[0][:,]
test.xyz[0]
test.xyz[0][0]
test.xyz[0][0][0]
test.xyz[0][:][0]
test.xyz[0][0:-1][0]
test.xyz[0][0:-1]
test.xyz[:,:,0]
test.xyz[0][0:5]
numpy.zeros((5,5,5))
import numpy
numpy.zeros((5,5,5))
test=numpy.zeros((5,5,5))
test[0,2,2]=1.4
test[0,1,1]=1.8
test[0,0,0]=1.0
test
test[0,1,1]=1.85
test[0,2,2]=1.42
test[0,0,0]=1.01
test
round(test,1)
numpy.round(test, decimal=1)
numpy.round(test, decimals=1)
import numpy, pylab
shape_names=numpy.loadtxt('emolecule_shape_overlap.txt', usecols=(1,))
shape_names=numpy.loadtxt('emolecule_shape_overlap.txt', usecols=(0,))
shape_scores=numpy.loadtxt('emolecule_shape_overlap.txt', usecols=(1,))
!head emolecule_amgen_correspond.txt
emolec_names=numpy.loadtxt('emolecule_amgen_correspond.txt', usecols=(0,), skiprows=1)
amgen_names=numpy.loadtxt('emolecule_amgen_correspond.txt', usecols=(1,), skiprows=1)
activity_names=numpy.loadtxt('C6emolecules_results.txt', usecols=(0,))
activity_scores=numpy.loadtxt('C6emolecules_results.txt', usecols=(1,))
activity_scores=numpy.loadtxt('C6emolecules_results.txt', usecols=(1,), dtype=str)
convert=dict()
data=dict()
for (emol, amgen) in zip(emolec_names, amgen_names):
    convert[emol]=amgen
convert
convert=dict()
for (emol, amgen) in zip(emolec_names, amgen_names):
    key=int(emol)
    convert[key]=int(amgen)
convert'
convert
data['activity']=dict()
data['shape']=dict()
for (name, activity) in zip(activity_names, activity_scores):
    data['activity'][int(name)]=activity
for (name, score) in zip(shape_names, shape_scores):
    data['shape'][int(name)]=int(score)
data
score
activity
data['activity']
data['score']
data['shape']
shape_scores
data['shape']=dict()
for (name, score) in zip(shape_names, shape_scores):
    data['shape'][int(name)]=score
data['shape
data['shape']
import readline
readline.write_history_file('plot_shape_vs_activity.py')
num=0.129
(num+0.5)/0.5
import picklet
import pickle
test=dict()
test[0]=1
ohandle=open('test', 'wb')
ohandle=open('test.pickle', 'wb')
picklet.dump(test, ohandle)
pickle.dump(test, ohandle)
ohandle.close()
test
test=222
fhandle=open('test.pickle', 'rb')
pickle.open(fhandle)
pickle.?
pickle??
pickle.load(fhandle)
fhandle.close()
1-1
0.679-1
abs(0.679-1)
import numpy
test=numpy.ones(4,4,4)
test=numpy.ones((4,4,4))
test[2,2]=0.55
test[2,0]=0.23
test[0,0]=0.77
test
test=numpy.ones((2,2,2))
test[0,0]=0.77
test[0,1]=0.23
test[1,1]=0.98
test
test-1
abs(test-1)
test
import numpy, pylab
file=open('sites_allcmd.txt')
line=file.readline()
line
data=dict()
for label in line.split():
    data[label]=[]
data
label.rstrip('"')
label
line.split()
line.split(',')
for label in line.split(','):
    data[label]=[]
label
label.rstrip('"')
label.rstrip().rstrip('"')
label.rstrip().rstrip('"').rstrip('"')
label.rstrip().rstrip('"')
label.rstrip().rstrip('"').split('"')[1]
for label in line.split(','):
    label=label.rstrip().rstrip('"').split('"')[1]
    data[label]=[]
data
line=fhandle.readline()
line=file.readline()
line
line.split(',')
len(line.split(','))
len(data.keys())
line
line.split(',')
labels
data.keys()
data=dict()
for label in line.split(','):
    label=label.rstrip().rstrip('"').split('"')[1]
    data[label]=[]
file=open('sites_allcmd.txt')
line=file.readline()
data=dict()
for label in line.split(','):
    label=label.rstrip().rstrip('"').split('"')[1]
    data[label]=[]
len(data.keys())
line
order=[i.rstrip().rstrip('"').split('"')[1] for i in line.split(',')]
print order
start=0
column=0
file=open('sites_allcmd.txt')
for line in file.readlines():
    if start==0:
        start=1
        continue
    else:
        for entry in line.split(','):
            key=order[column]
            data[key].append(entry.rstrip().rstrip('"').split('"')[1])
            column+=1
data
import readline
readline.write_history_file('scrape_mdb.py')
import glob
import os
ls frame10000.pdb
ls format_frame10000.pdb
import pdb
import readline
readline.write_history_file('change_name.py')
test1=numpy.array([1,2,3,4,5])
import numpy
test1=numpy.array([4,5,6,7,8])
test1=numpy.array([4,5,6,7,8,9,10])
test2=numpy.array([4,5,6,7,8,9,10])
test1=numpy.array([4,5,6,7,8])
test1=numpy.array([1,2,3,4,5])
numpy.logical_and(numpy.logical_and(test1 !=0, test2 !=0), test1 == test2 )
test2=numpy.array([4,5,6,7,8])
numpy.logical_and(numpy.logical_and(test1 !=0, test2 !=0), test1 == test2 )
T=test1-test2
print T
test1
if test1 == test2.any():
    print True
if test1.any() == test2.any():
    print True
numpy.allclose(test1, test2)
numpy.intersect1d(test1, test2)
test2=numpy.array([4,5,6,7,8,9,10])
numpy.intersect1d(test1, test2)
import scipy as sp
test1=(0.00, 0.00, 0.00)
test2=(2.00, 2.00, 2.00)
sp.distance.cdist(test1, test2)
!vi Site3D.py
import scipy.spatial as sp
sp.distance.cdist(test1, test2)
test.reshape(1,-1)
test2.reshape(1,-1)
test2=numpy.array([2.00, 2.00, 2.00])
import numpy
test2=numpy.array([2.00, 2.00, 2.00])
test1=numpy.array([0.00, 0.00, 0.00])
sp.distance.cdist(test1, test2)
sp.distance.cdist(test1, test2.reshape(1,-1))
sp.distance.cdist(test1.reshape(1,-1), test2.reshape(1,-1))
test1
test2
test1.reshape(1,-1)
test1.reshape(1,-1).shape
test1=numpy.array([0.00, 0.00, 0.00], ndmin=1)
test1.shape
sp.distance.cdist(test1.reshape(1,-1), test2.reshape(1,-1))
numpy.sqrt(2**2+2**2)
numpy.sqrt(2**2+2**2+2**2)
import numpy
test=numpy.meshgrid((3,3,3))
import numpy
range=numpy.arange(0,3,0.5)
test=numpy.meshgrid(range, range, range)
test.shape
test
test[0].shape
test[0,:]
test=numpy.meshgrid(range, range, range)
test[0].shape
test[0]
test=numpy.grid(range, range, range)
test=numpy.mgrid(range, range, range)
numpy.hstack((range, range))
test=numpy.vstack((range, range))
test=numpy.vstack((test, range))
test.shape
test
test[0]
!vi sphere_distance_3D.py
X,Y,Z=numpy.meshgrid(x,y,z)
ax.scatter(X,Y,Z, c='r', alpha=0.5)
data=numpy.vstack((X.ravel(), Y.ravel(), Z.ravel())).T
X,Y,Z=numpy.meshgrid(x,y,z)
X,Y,Z=numpy.meshgrid(range,range,range)
data=numpy.vstack((X.ravel(), Y.ravel(), Z.ravel())).T
data.shape
data
data.shape
data[0]
data[1]
data[2]
import time
start=time.time()
end=time.time()
print start-end
pritn end-start
print end-start
line=['A']*20
line
line=['A'*20]
line
line=[' '*20]
line
line[15:28]='aaa'
line
line='ATOM CA  10   10   10'
line='ATOM CA  10   10   {0 >3}'.format(120)
line='ATOM CA  10   10   {0: >3}'.format(120)
print line
line='ATOM CA  10   10   10'
newline='ATOM CA  10   10   {0: >3}'.format(120)
print line
print newline
newline='ATOM CA  10   10  {0: >3}'.format(120)
print newline
print line
newline='ATOM CA  10   10  {0}: > 3}'.format(120
)
newline='ATOM CA  10   10  {0: > 3}'.format(120
)
print newline
prniot line
print line
print newline
newline='ATOM CA  10   10 {0: > 3}'.format(120
)
print newline
print line
newline='ATOM CA  10   10 {0: > 3}'.format(10)
print newline
print line
line
newline='ATOM CA  10   10  {0: >3}'.format(120)
newline
line
newline='ATOM CA  10   10  {0: >3}'.format(10)
newline
line
import numpy
test1=numpy.array([1,2,3])
test2=numpy.array([1,2,3])
numpy.array_equal(test1, test2)
import scipy.spatial
import numpy
test1=numpy.array([38.0, 47.5, 25.0])
test2=numpy.array([38.0, 47.5, 25.0])
test2=numpy.array([38.33499, 47.451, 25.496])
print scipy.spatial.distance(test1, test2)
import scipy.spatial as sp
print sp.cdist(test1, test2)
!grep cdist AnaSites.py
!grep dist AnaSites.py
!grep dist Site3D.py
print sp.distance.cdist(test1, test2)
print sp.distance.cdist(test1.reshape(1,-1), test2.reshape(1,-1))
print sp.distance.cdist(test1.reshape(1,-1), test2.reshape(1,-1))[0]
import mdtraj as md
test=mdtraj.load('../SiteTest/aln-holo.pdb')
test=md.load('../SiteTest/aln-holo.pdb')
test.restrict_atoms??
test.restrict_atoms?
test.atom_slice?
test.topology.atoms_by_name
test.topology.atoms_by_name?
for i in test.topology.atoms_by_name:
    if 'H' not in i:
        print i
test.topology.atoms
for i in test.topology.atoms:
    if 'H' not in i:
        print i
print i
test.topology.
test.topology.atoms
for i in test.topology.atom:
    if 'H' not in i:
        print i
names=[]
for i in test.topology.atoms:
    names.append(i)
names
len(names)
traj.xyz.shape
test.xyz.shape
for (n,name) in enumerate(names):
    if 'H' not in name.split('-')[1]:
        print n, name.split('-')[1]:
for (n,name) in enumerate(names):
    if 'H' not in name.split('-')[1]:
        print n, name.split('-')[1]
names
names[0]
len(names[0])
name.name
for i in test.topology.atoms:
    if 'H' not in i.name:
        priont i.index
for i in test.topology.atoms:
    if 'H' not in i.name:
        print i.index
print i.name
indices=[]
for i in test.topology.atoms:
    if 'H' not in i.name:
        indices.append(i.index)
len(indices)
l
test.shape
test.xyz.shape
test.xyz[0][indices,:]
reduced_coors=test.xyz[0][indices,:]
reduced_coors.shape
indices[0]
indices[200]
test.xyz[0][200]
reduced_coors[200]
test.xyz[0][423]
import readline
readline.write_history_file('exclude hydrogens.py')
import numpy
test=numpy.zeros((2,2))
test
test=numpy.zeros((3, 2, 2,2))
test
test[0]
test[0][0][0]
test[0][0][0][0]
test.shape
test[0]
test[1]
test[2]
test=numpy.zeros((5, 3))
test
test[0]
test[0][1]=1
test[1][2]=5
test
test[3][0]=2
test
test.shape
test.ravel()
test=numpy.zeros((3,3,3))
test
test.shape
count=0
for i in range(0,3):
    for j in range(0,3):
        for k in range(0,3):
            if count==2 or count==6:
                test[i,j,k]=1
                count=+1
            else:
                count+=1
test
test[0,2,2]
count
test
test.shape
count=0
test=numpy.zeros((5,4,3))
for i in range(0,5):
    for j in range(0,4):
        for k in range(0,3):
            if count==2 or count==6:
                test[i,j,k]=1
                count+=1
            else:
                count+=1
test
test.shape
new=test.evale()
new=test.ravel()
new
new.shape
new
test
test[0]
test[0,0,1]
new[0]
test=numpy.meshgrid((range(0,5), range(0,4), range(0,3))
)
test=numpy.meshgrid((range(0,5), range(0,4), range(0,3)))
)
test=numpy.meshgrid(range(0,5), range(0,4), range(0,3))
)
test=numpy.meshgrid(range(0,5), range(0,4), range(0,3))
test.shape
test[0]
new=test.ravel().T
new=test.ravel()
!vi Site3D.py
new=numpy.vstack((test[0].ravel(), test[1].ravel(), test[2].ravel())
)
new
new.shape
new=numpy.vstack((test[0].ravel(), test[1].ravel(), test[2].ravel())
).T
new.shape
test[0]
test[0].shape
test
test=numpy.meshgrid(range(0,5), range(0,4), range(0,3))
X,Y,Z=numpy.meshgrid(range(0,5), range(0,4), range(0,3))
new=numpy.vstack((X.ravel(), Y.ravel(), Z.ravel())
).T
new.shape
new[0]
new[1]
new[2]
new[3]
test[0]
test[0][0]
test[0][0][0]
test[0].shape
new.shapenew[3]
new
new[0]
new[1]
new[2]
new[3]
new[4]
new[5]
new[6]
new[7]
new[8]
new[9]
new[10]
new[11]
new[12]
new[13]
new[14]
new[15]
new=numpy.vstack((X.ravel(), Y.ravel(), Z.ravel())
)
new.shape
new.shape[0]
new=new.reshape((60,3)
)
new.shape
new[0]
new[1]
new[2]
new[3]
new[4]
new[5]
new[6]
new[7]
new[8]
new[9]
new[55]
X.ravel()
new=numpy.vstack((X.ravel(), Y.ravel(), Z.ravel())
).T
new=numpy.vstack((X.ravel(), Y.ravel())).T
new[0]
new[1]
new[2]
new[3]
new[4]
Y.ravel()
test[0]
test[0][0]
X.ravel()
Y.ravel()
Z.ravel()
orig.shape
orig=numpy.zeros((5,4,3))
orig.shape
X,Y,Z=numpy.meshgrid(range(0,5), range(0,4), range(0,3))
new=numpy.vstack((X.ravel(), Y.ravel(), Z.ravel())).T
new.shape
count=0
for j in range(0,4):
    for x in range(0,5):
        for z in range(0,3):
            print new[count]
            count+=1
new[60]
new[59]
new.shape
X,Y,Z=numpy.meshgrid(range(0,1), range(0,3), range(0,1))
new=numpy.vstack((X.ravel(), Y.ravel(), Z.ravel())).T
count=0
for j in range(0,3):
    for x in range(0,1):
        for z in range(0,1):
            print new[count]
            count+=1
new.shape
new
new.shape[0]
test='ABCD'
test[0]
test=400
len(test)
test.size
len(str{test))
len(str(test))
test='ABCD'
test*2
numpy.loadtxt('frame1_sites.pdb', usecols=(5,6,7))
import numpy
test=numpy.loadtxt('frame1_sites.pdb', usecols=(5,6,7))
test=numpy.loadtxt('frame1_sites.pdb', usecols=(5))
test=numpy.loadtxt('frame1_sites.pdb', usecols=(5), skiprows=3)
test=numpy.loadtxt('frame1_sites.pdb', usecols=(5,), skiprows=3)
import numpy
test1=numpy.array([1,2,3,4])
test2=numpy.array([1,2,3,4])
test3=numpy.array([1,2,3,4])
numpy.hstack((test1, test2, test3))
numpy.vstack((test1, test2, test3))
numpy.concatenate((test1, test2, test3))
numpy.concatenate((test1, test2, test3), axis=1)
numpy.concatenate((test1, test2, test3), axis=2)
numpy.concatenate((test1, test2, test3), axis=3)
numpy.vstack((test1, test2, test3))
numpy.vstack((test1.T, test2.T, test3.T))
test1.T
test1
test1.T()
numpy.vstack((test1.T, test2.T, test3.T))
test1
test1.T
numpy.dstack((test1, test2, test3))
import numpy
!vi Site3D.py
import scipy.spatial as sp
test1=numpy.zeros((5,3))
test1
test2=numpy.ones((2,10,3))
test2[0].shape
test2[0][0]
result=sp.distances.cdist(test1, test2)
result=sp.distance.cdist(test1, test2)
test1
test2
result=sp.distance.cdist(test2, test1)
import os
os.path.dirname('../MoeSiteTest/cmd_1_reduced_frames/format_frame01000.pdb')
data=numpy.loadtxt('gyr-NIST-25.concat.colvars.traj', usecols=(1,))
import numpy
data=numpy.loadtxt('gyr-NIST-25.concat.colvars.traj', usecols=(1,))
data.shape
numpy.min(data)
print '04i' % 4
print '%04i' % 4
print '%04i' % 40
import numpy, pylab
data=numpy.loadtxt('rgyr-NIST-25.dat', usecols=(1,), skiprows=1)
pylab.plot(dataa)
pylab.plot(data)
pylab.show()
import numpy, pylab
data=numpy.loadtxt('rgyr-BSA-15.dat', usecols=(1,), skiprows=1)
pylab.plot(data)
pylab.show()
import os
for file in glob.glob('min*pdb'):
    num=int(file.split('rmsd-')[1].split('.pdb'))
    print 'min-long_cmd-holo-rmsd-%05i.pdb' % num
import glob
for file in glob.glob('min*pdb'):
    num=int(file.split('rmsd-')[1].split('.pdb'))
    print 'min-long_cmd-holo-rmsd-%05i.pdb' % num
files=glob.glob('min*pdb')
for file in files:
    num=int(file.split('rmsd-')[1].split('.pdb')[0])
    print 'min-long_cmd-holo-rmsd-%05i.pdb' % num
ls
num
for file in files:
    num=int(file.split('rmsd-')[1].split('.pdb')[0])
    print 'min-long_cmd-holo-rmsd-%05i.pdb' % num
num
print '%05i' % num
file
mv min-long_cmd-holo-rmsd.pdb ..
for file in files:
    num=int(file.split('rmsd-')[1].split('.pdb')[0])
    print 'min-long_cmd-holo-rmsd-%05i.pdb' % num
file
!mv min-long_cmd-holo-rmsd.pdb ..
files=glob.glob('min*pdb')
for file in files:
    num=int(file.split('rmsd-')[1].split('.pdb')[0])
    print 'min-long_cmd-holo-rmsd-%05i.pdb' % num
files
for file in files:
    num=int(file.split('rmsd-')[1].split('.pdb')[0])
    os.system('mv  %s min-long_cmd-holo-rmsd-%05i.pdb' % (file, num))
pwd
cd ../amd_frames
files=glob.glob('min*pdb')
for file in files:
    num=int(file.split('rmsd-')[1].split('.pdb')[0])
    os.system('mv  %s min-amd-holo-rmsd-%05i.pdb' % (file, num))
file
! mv min-amd-dih-holo-rmsd.pdb ..
for file in files:
    num=int(file.split('rmsd-')[1].split('.pdb')[0])
    os.system('mv  %s min-amd-holo-rmsd-%05i.pdb' % (file, num))
file
files=glob.glob('min*pdb')
for file in files:
    num=int(file.split('rmsd-')[1].split('.pdb')[0])
    os.system('mv  %s min-amd-holo-rmsd-%05i.pdb' % (file, num))
file
!mv min-amd-holo-rmsd.pdb ..
files=glob.glob('min*pdb')
for file in files:
    num=int(file.split('rmsd-')[1].split('.pdb')[0])
    os.system('mv  %s min-amd-holo-rmsd-%05i.pdb' % (file, num))
import glob, os
files=glob.glob('aln*pdb')
for file in files:
    num=int(file.split('cmd-')[1].split('.pdb')[0])
    os.system('mv  %s aln-strip-all-1.2micros-cmd-%05i.pdb' % (file, num))
fhandle=open('Compounds_pyrazole.smi')
pyrazole=[]
for line in fhandle.readlines():
    pyrazole.append(line)
pyrazole
pyridine=[]
fhandle=open('Compounds_pyridine.smi')
for line in fhandle.readlines():
    pyridine.append(line)
pyridine
fhandle=open('Compounds.smi')
all=[]
for line in fhandle.readlines():
    all.append(line)
for i in pyridine:
    if i in pyrazole:
        print i
for i in pyrazole:
    if i in pyridine:
        print i
for i in all:
    if i not in pyridine:
        if i not in pyrazole:
            print i
ohandle=open('Compounds_missed.smi', 'w')
for i in all:
    if i not in pyridine:
        if i not in pyrazole:
            ohandle.write(i)
ohandle.close()
!more Compounds_missed.smi
import readline
readline.write_history_file('check.py')
import numpy, pylab
data=numpy.loadtxt('rgyr-BSA-15.dat')
high_data=numpy.loadtxt('rgyr-high-BSA-15.dat')
pylab.plot(data, 'low')
data
ls
!ls *gyr*dat
! head rgyr-BSA-15.dat
data=numpy.loadtxt('rgyr-BSA-15.dat', usecols=(1,))
high_data=numpy.loadtxt('rgyr-high-BSA-15.dat', usecols=(1,))
pylab.plot(data, label='low')
pylab.plot(high_data, label='hi')
pylab.legend()
pylab.show()
import numpy, pylab
data=numpy.loadtxt('gyr-NIST-25-all.pmf', usecols=(1,))
coor=numpy.loadtxt('gyr-NIST-25-all.pmf', usecols=(0,))
pylab.scatter(coor, data)
pylab.show90
pylab.show()
import numpy, pylab
data=numpy.loadtxt('gyr-BSA-15-all.pmf', usecols=(1,))
coor=numpy.loadtxt('gyr-BSA-15-all.pmf', usecols=(0,))
pylab.scatter(coor, data)
pylab.show()
pwd
cd ../../../namd-gas-simulations/
cd BSA-
ls
ls *dat
data=numpy.loadtxt('rgyr-BSA-15.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
import numpy, pylab
data=numpy.loadtxt('rgyr-BSA-15.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
import numpy, pylab
data=numpy.loadtxt('rgyr-high-BSA-15.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
import numpy, pylab
data=numpy.loadtxt('rgyr-high-BSA-15-6.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
import numpy, pylab
data=numpy.loadtxt('rgyr-BSA-15.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
import numpy, pylab
data=numpy.loadtxt('rgyr-BSA-15-9.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
import numpy, pylab
data=numpy.loadtxt('rgyr-NIST-25.dat', usecols=(1,))
coor=numpy.loadtxt('rgyr-NIST-25.dat', usecols=(0,))
pylab.plot(data)
pylab.show()
import numpy, pylab
coor=numpy.loadtxt('rgyr-NIST-25-4.dat', usecols=(0,))
data=numpy.loadtxt('rgyr-NIST-25-5.dat', usecols=(1,))
!more ptraj-NIST-25.in
pylab.plot(data)
pylab.show()
coor=numpy.loadtxt('rgyr-NIST-25.dat', usecols=(0,))
data=numpy.loadtxt('rgyr-NIST-25.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
data=numpy.loadtxt('rgyr-NIST-25-5.dat', usecols=(1,))
import numpy, pylab'
import numpy, pylab
data=numpy.loadtxt('rgyr-NIST-25-5.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
import numpy, pylab'
import numpy, pylab
data=numpy.loadtxt('rgyr-CytC-7.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
data=numpy.loadtxt('rgyr-Ub-5.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
import numpy, pylab
data=numpy.loadtxt('rgyr-CytC-7.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
import numpy, pylab
data=numpy.loadtxt('rgyr-equil-NIST-25-5.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
import numpy, pylab
data=numpy.loadtxt('rgyr-equil-BSA-15-9.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
systems=['CytC-7', 'BSA-15', 'NIST-25']
types=['anderson-highcutoff', 'langevin-highcutoff', 'langevin-lowcutoff']
data=dict()
for type in types:
    data[type]=dict()
    for sys in systems:
        files=glob.glob('%s/*%s*' % (type, sys))
        if files:
            tmp=numpy.loadtxt(file, usecols=(1,))
            data[type][sys]=tmp
import glob
data=dict()
for type in types:
    data[type]=dict()
    for sys in systems:
        files=glob.glob('%s/*%s*' % (type, sys))
        if files:
            tmp=numpy.loadtxt(file, usecols=(1,))
            data[type][sys]=tmp
import numpy, pylab
data=dict()
for type in types:
    data[type]=dict()
    for sys in systems:
        files=glob.glob('%s/*%s*' % (type, sys))
        if files:
            tmp=numpy.loadtxt(file, usecols=(1,))
            data[type][sys]=tmp
files
data=dict()
for type in types:
    data[type]=dict()
    for sys in systems:
        files=glob.glob('%s/*%s*' % (type, sys))
        if files:
            for file in files:
                data[type][sys]=numpy.loadtxt(file, usecols=(1,))
data
for sys in systems:
    pylab.figure()
    for key in data.keys():
        if sys in data[key].keys():
            pylab.plot(data[key][sys], label=type)
            pylab.title(sys)
pylab.show()
for sys in systems:
    pylab.figure()
    for key in data.keys():
        if sys in data[key].keys():
            pylab.plot(data[key][sys], label=type)
            pylab.title(sys)
            pylab.legend()
pylab.show()
for sys in systems:
    pylab.figure()
    for key in data.keys():
        if sys in data[key].keys():
            pylab.plot(data[key][sys], label=key)
            pylab.title(sys)
            pylab.legend()
pylab.show()
import readline
readline.write_history_file('make_plots.py')
import readline
readline.write_history_file('make_plots.py')
import readline
readline.write_history_file('tmp')
import numpy, pylab
import glob
files=glob.glob('*rmsd*dat')
data=dict()
files=glob.glob('*holo*rmsd*dat')
files=glob.glob('*rmsd*dat')
head holo-kif18-pocket-rmsd.dat
!head holo-kif18-pocket-rmsd.dat
import mdtraj
traj=mdtraj.load('aln-strip-amd-step10.dcd', topo='../aln-apo.pdb')
traj=mdtraj.load('aln-strip-amd-step10.dcd', top='../aln-apo.pdb')
traj.xyz
traj.xyz.shape
traj=mdtraj.load('aln-strip-all-1.2micros-cmd.dcd', top='../aln-apo.pdb')
traj.xyz.shape
import numpy
array1=numpy.array([2,2,2])
array1=numpy.array([1,1,1])
array2=numpy.array([2,2,2])
ref=numpy.array([4,4,4])
matrix=numpy.vstack(array1, array2)
matrix=numpy.vstack((array1, array2))
matrix.shape
matrix
!grep cipy Site3D.py
import scipy.spatial as sp
!grep sp Site3D.py
sp.distance.cdist(ref, matrix)
ref
mesh=numpy.meshgrid(ref, ref)
mesh.shape
sp.distance.cdist(mesh, matrix)
mesh[0]
mesh[1]
matrix.shape
matrix.shape[0]
numpy.vstack((ref, ref))
test=numpy.vstack((ref, ref))
sp.distance.cdist(test, matrix)
matrix
test
import numpy
import scipy.spatial as sp
array1=([2,2,2])
array1=numpy.array([2,2,2])
array1=numpy.array([1,1,1])
array2=numpy.array([2,2,2])
matrix=numpy.zeros((1,2,3))
matrx
matrix
matrix[0]=array1
matrix[2]=array2
matrix[1]=array2
matrix[0]=array1
matrix
matrix[0][0]=array1
matrix[0][1]=array2
matrix
matrix.shape
ref=numpy.array([4,4,4])
sp.distance.cdist(ref, array1)
ref
array1
sp.distance.cdist((ref, array1))
sp.distance.cdist(ref, matrix[0])
matrix.shape
ref.shape
ref
arrat
sp.distance.cdist(ref, matrix[0][:])
matrix[0].shape
len(matrix[0])
len(ref)
ref
sp.distance.cdist(matrix[0][:], ref)
matrix
ref
ref_matrix=numpy.vstack((ref, ref))
ref_matrix.shape
ref_matrix
sp.distance.cdist(ref, matrix[0])
sp.distance.cdist(ref, matrix)
sp.distance.cdist(ref, array1)
sp.distance.cdist(ref, matrix0])
sp.distance.cdist(ref_matrix, matrix[0])
ref_matrix
matrix[0]
!vi Site3D.py
ref=numpy.array([0,0,0])
ref2=numpy.array([4,4,4])
ref_matrix=numpy.vstack((ref, ref2))
ref_amtrix
ref_matrix
sp.distance.cdist(ref_matrix, matrix[0])
matrix[0]
matrix
matrix.shape
new_matrix=numpy.zeros((2,2,3))
new_matrix[0][0]=array1
new_matrix[0][1]=array2
new_matrix[1][0]=array1
new_matrix[1][1]=array2
new_matrix
new_matrix.shape
ref
ref2
ref_matrix=numpy.zeros((2,2,3))
ref_matrix[0][0]=ref
ref_matrix[0][1]=ref2
ref_matrix[1][0]=ref
ref_matrix[1][1]=ref2
ref_matrix
ref_matrix.shape
sp.distance.cdist(ref_matrix, new_matrix)
ref_matrix.shape
new_matrix.shape
import readline
readline.write_history_file('check.py')
import check
test, ref=check.make()
test
ref
import scipy.spatial as sp
sp.distance.cdist(ref, test)
ref.shape
sp.distance.cdist(ref[0], test)
ref[0].shape
sp.distance.cdist(ref[0], test[0])
sp.distance.cdist(ref[0], test[0], 'euclidean')
ref.shape
ref.T
ref.T.shape
sp.distance.cdist(ref.T, test.T, 'euclidean')
sp.distance.cdist(ref[:], test.T, 'euclidean')
sp.distance.cdist(ref[:], test[:], 'euclidean')
test[:]
test[:].shape
test
ref
new=627.5*(data[file][type]-min)
min=-506.2375
data=-506.2355
new=627.5*(data[file][type]-min)
new=627.5*(data-min)
new
import numpy
data=numpy.loadtxt('rgyr-equil-BSA-15-9.dat', usecols=(1,))
import pylab
pylab.plot(data)
pylab.show()
import numpy, pylab
data=numpy.loadtxt('rgyr-high-BSA-15.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
!ls -l *dat
data=numpy.loadtxt('rgyr-high-BSA-15-6.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
import numpy, pylab
data=numpy.loadtxt('rgyr-equil-NIST-25-5.dat', usecols=(1,))
pylab.plot(data)
pylab.show()
import numpy, pylab
data=numpy.loadtxt('rgyr-BSA-15-all.pmf')
pwd
ls *all.pmf
data=numpy.loadtxt('rgyr-BSA-15-all.pmf')
data=numpy.loadtxt('gyr-BSA-15-all.pmf')
data.shape
pylab.plot(data[0], data[1])
pylab.show()
data[:,0]
data[:,1]
pylab.plot(data[:,0], data[:,1])
pylab.show()
import numpy, pylab
data=numpy.loadtxt('gyr-NIST-25-all.pmf')
pylab.plot(data[:,0], data[:,1])
pylab.show()
import glob
files=glob.glob('../Mcl-1/analysis/trajs/open*gridopt.pdb')
files
fhandle=open(files[0])
fhandle.readline()
line=fhandle.readline()
line.split()
line.split()[5]
line.split()[6]
xcoors=numpy.loadtxt(files[0], usecols=(5,))
import numpy
xcoors=numpy.loadtxt(files[0], usecols=(5,))
xcoors=numpy.loadtxt(files[0], usecols=(6,))
xcooers
xcoors
ycoors=numpy.loadtxt(files[0], usecols=(7,))
zcoors=numpy.loadtxt(files[0], usecols=(8,))
numpy.vstack((xcoors, ycoors))
test=numpy.hstack((xcoors, ycoors))
all=numpy.hstack((test, zcoors))
all
all.shape
test=numpy.hstack((xcoors.T, ycoors.T))
test
test.shape
test=numpy.hstack((xcoors, ycoors))
test
test[0]
test.shape
xcoors=numpy.loadtxt(files[0], usecols=(6,), ndmin=1)
ycoors=numpy.loadtxt(files[0], usecols=(7,), ndmin=1)
zcoors=numpy.loadtxt(files[0], usecols=(8,), ndmin=1)
xcoors.shape
xcoors.reshape(-1,1)
xcoors.reshape(-1,1).shape
xcoors=xcoors.reshape(-1,1)
ycoors=ycoors.reshape(-1,1)
zcoors=zcoors.reshape(-1,1)
test=numpy.hstack((xcoors, ycoors))
test.shape
test
all=numpy.hstack((test, zcoors))
all
import readline
readline.write_history_file('reload_open_coors.py')
