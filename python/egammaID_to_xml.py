# e/g ID for HGCal L1Trigger: converting model to xml format to be read by root

#usual imports
import numpy as np
import xgboost as xg
import pickle
import pandas as pd
import ROOT as r
from root_numpy import tree2array, testdata, list_branches, fill_hist
from os import system, path
import os

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: converting e/g ID model to xml format"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#configure options
from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--modelAlgo', dest='modelAlgo', default='default', help="Clustering algorithm with which BDT was trained" )
  parser.add_option('--signal', dest='signal', default='electron', help="Input signal type" )
  parser.add_option('--background', dest='background', default='neutrino', help="Input background type" )
  parser.add_option('--release', dest='release', default='103X', help="CMSSW release (for geometry)" )
  parser.add_option('--bdtName', dest='bdt_name', default='baseline', help="String to identify BDT" )
  return parser.parse_args()

(opt,args) = get_options()

modelAlgo = opt.modelAlgo
release = opt.release
bdt_name = opt.bdt_name

#set up global variables
modelDir = os.environ['CMSSW_BASE']+"/src/L1Trigger/analysis/output/models/%s"%release

#define variables used in model
egID_var_dict = {'baseline':['cl3d_coreshowerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_srrmean'],'full':['cl3d_coreshowerlength','cl3d_showerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_szz','cl3d_srrmean','cl3d_srrtot','cl3d_seetot','cl3d_spptot']}
egID_vars = egID_var_dict[ bdt_name ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load models
egID_model = xg.Booster()
if bdt_name == "baseline": modelStr = "%s/egID_%s_%s_vs_%s.model"%(modelDir,modelAlgo,opt.signal,opt.background)
else: modelStr = "%s/egID_%s_%s_vs_%s_%s.model"%(modelDir,modelAlgo,opt.signal,opt.background,bdt_name)
egID_model.load_model( modelStr )
print "  --> Loaded model: %s"%modelStr

# Define name of xml file to save
if bdt_name == "baseline": f_xml = "%s/egID_%s_%s_vs_%s.xml"%(modelDir,modelAlgo,opt.signal,opt.background)
else: f_xml = "%s/egID_%s_%s_vs_%s_%s.xml"%(modelDir,modelAlgo,opt.signal,opt.background,bdt_name)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Convert to xml
from mlglue.tree import tree_to_tmva, BDTxgboost, BDTsklearn
target_names = ['background','signal']
bdt = BDTxgboost( egID_model, egID_vars, target_names, kind='binary', max_depth=6, learning_rate=0.3 )
bdt.to_tmva( f_xml )
print "  --> Converted to xml"
