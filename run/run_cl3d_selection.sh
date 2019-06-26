#!/bin/bash

# Script to run the HGCal L1T cluster definition
cd /afs/cern.ch/work/j/jlangfor/HGCal/L1T_new/v3_9_4/CMSSW_10_5_0/src/L1Trigger/HGCal_L1T_egammaID/run
eval `scramv1 runtime -sh`

input_type=$1
fileNumber=$2
clusteringAlgo=$3
geometry=$4

python $CMSSW_BASE/src/L1Trigger/HGCal_L1T_egammaID/python/HGCalL1T_cluster_selection.py --input $input_type --fileNumber $fileNumber --maxEvents -1 --clusteringAlgo $clusteringAlgo --geometry $geometry

