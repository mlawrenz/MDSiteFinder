
import mdtraj as md
indices=[]
for i in test.topology.atoms:
    if 'H' not in i.name:
        indices.append(i.index)
reduced_coors=test.xyz[0][indices,:]
