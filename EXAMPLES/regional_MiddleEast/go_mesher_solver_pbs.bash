#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_mesher_solver
#PBS -j oe
#PBS -o OUTPUT_FILES/job.o

###########################################################
# USER PARAMETERS

## 64 CPUs ( 8*8 ), walltime 5 hour
#PBS -l nodes=8:ppn=8,walltime=5:00:00

###########################################################

cd $PBS_O_WORKDIR

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# script to run the mesher and the solver
# read DATA/Par_file to get information about the run
# compute total number of nodes needed
NPROC_XI=`grep ^NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep ^NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep ^NCHUNKS DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))

mkdir -p OUTPUT_FILES

# backup files used for this simulation
cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/

# obtain job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid

##
## mesh generation
##
sleep 2

echo
echo `date`
echo "starting MPI mesher on $numnodes processors"
echo

mpiexec -np $numnodes $PWD/bin/xmeshfem3D

echo "  mesher done: `date`"
echo

# backup important files addressing.txt and list*.txt
cp OUTPUT_FILES/*.txt $BASEMPIDIR/


##
## forward simulation
##

# set up addressing
#cp $BASEMPIDIR/addr*.txt OUTPUT_FILES/
#cp $BASEMPIDIR/list*.txt OUTPUT_FILES/

sleep 2

echo
echo `date`
echo starting run in current directory $PWD
echo

mpiexec -np $numnodes $PWD/bin/xspecfem3D

echo "finished successfully"
echo `date`

