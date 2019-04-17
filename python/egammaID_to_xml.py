# e/g ID for HGCal L1Trigger: converting model to xml format to be read by root

#usual imports
import numpy as np
import xgboost as xg
import pickle
import pandas as pd
import ROOT as r
from root_numpy import tree2array, testdata, list_branches, fill_hist
from os import system, path

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: converting e/g ID model to xml format"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#configure options
from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--modelAlgo', dest='modelAlgo', default='default', help="Clustering algorithm with which BDT was trained" )
  parser.add_option('--inputSignalType', dest='inputSignalType', default='electron', help="Input signal type" )
  parser.add_option('--inputBackgroundType', dest='inputBackgroundType', default='neutrino', help="Input background type" )
  return parser.parse_args()

(opt,args) = get_options()

modelAlgo = opt.modelAlgo

#set up global variables
modelDir = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/models/%s/"%modelAlgo

#define variables used in model
egID_vars = ['coreshowerlength','firstlayer','maxlayer','srrmean']

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load models
egID_model = xg.Booster()
egID_model.load_model( "%segID_%s_sig_%s_bkg_%s.model"%(modelDir,modelAlgo,opt.inputSignalType,opt.inputBackgroundType) )
print "  --> Loaded model: egID_%s_sig_%s_bkg_%s.model"%(modelAlgo,opt.inputSignalType,opt.inputBackgroundType)

# Define name of xml file to save
f_xml = "%segID_%s_sig_%s_bkg_%s.xml"%(modelDir,modelAlgo,opt.inputSignalType,opt.inputBackgroundType)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Convert to xml
from mlglue.tree import tree_to_tmva, BDTxgboost, BDTsklearn
target_names = ['eg_background','eg_signal']
bdt = BDTxgboost( egID_model, egID_vars, target_names, kind='binary', max_depth=6, learning_rate=0.3 )
bdt.to_tmva( f_xml )
print "  --> Converted to xml"
