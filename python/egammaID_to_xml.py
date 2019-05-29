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
  parser.add_option('--modelAlgo', dest='modelAlgo', default='Histomax_vardr', help="Clustering algorithm with which BDT was trained" )
  parser.add_option('--signal', dest='signal', default='electron', help="Input signal type" )
  parser.add_option('--background', dest='background', default='neutrino', help="Input background type" )
  parser.add_option('--geometry', dest='geometry', default='v9', help="HGCal Geomtry Config" )
  parser.add_option('--bdtConfig', dest='bdtConfig', default='baseline', help="String to identify BDT" )
  parser.add_option('--etaRegion', dest='etaRegion', default='total', help="Eta region of clusters used for training" )
  return parser.parse_args()

(opt,args) = get_options()

modelAlgo = opt.modelAlgo
geometry = opt.geometry
eta_region = opt.etaRegion
bdt_name = "%s_vs_%s_%s"%(opt.signal,opt.background,opt.bdtConfig)

#set up global variables
modelDir = os.environ['CMSSW_BASE']+"/src/L1Trigger/egid_analysis/HGCal_L1T_egammaID/output/models/%s"%geometry

#define variables used in model
egID_var_dict = {'electron_vs_neutrino_baseline':['cl3d_coreshowerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_srrmean'],'electron_vs_neutrino_full':['cl3d_coreshowerlength','cl3d_showerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_szz','cl3d_srrmean','cl3d_srrtot','cl3d_seetot','cl3d_spptot']}
egID_vars = egID_var_dict[ bdt_name ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load models
egID_model = xg.Booster()
if eta_region in ['low','high']: modelStr = "%s/egID_%s_%s_%seta.model"%(modelDir,modelAlgo,bdt_name,eta_region)
else: modelStr = "%s/egID_%s_%s.model"%(modelDir,modelAlgo,bdt_name)
egID_model.load_model( modelStr )
print "  --> Loaded model: %s"%modelStr

# Define name of xml file to save
if eta_region in ['low','high']: f_xml = "%s/egID_%s_%s_%seta.xml"%(modelDir,modelAlgo,bdt_name,eta_region)
else: f_xml = "%s/egID_%s_%s.xml"%(modelDir,modelAlgo,bdt_name)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Convert to xml
from mlglue.tree import tree_to_tmva, BDTxgboost, BDTsklearn
target_names = ['background','signal']
bdt = BDTxgboost( egID_model, egID_vars, target_names, kind='binary', max_depth=6, learning_rate=0.3 )
bdt.to_tmva( f_xml )
print "  --> Converted to xml"
