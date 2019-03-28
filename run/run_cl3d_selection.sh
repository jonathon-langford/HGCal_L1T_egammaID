#!/bin/bash

# Script to run the HGCal analysis: resource management in motherboards
cd /afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/run
eval `scramv1 runtime -sh`

input_type=$1
fileNumber=$2
clusteringAlgo=$3

python $CMSSW_BASE/src/L1Trigger/analysis/python/HGCalL1T_cluster_selection.py --input $input_type --fileNumber $fileNumber --maxEvents -1 --clusteringAlgo $clusteringAlgo

