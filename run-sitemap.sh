#!/bin/bash

x=1
while [ $x -lt 5 ]
do
python  RunMDSiteFinder.py -t ../Mcl-1/analysis/trajs/strip-all-cmd.dcd  -p test/aln-apo_pocket${x}.pqr -y test/aln-apo.pdb -o cmd_fpocket_site${x} >& out${x} &
let x=($x + 1)
done
