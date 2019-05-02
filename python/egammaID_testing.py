# Testing the egID models for HGCal L1Trigger: creates test dataframes and loads different models

# Usual imports
import ROOT
import math
import numpy as np
import pandas as pd
import xgboost as xg
from xgboost import plot_importance
import matplotlib.pyplot as plt
import pickle
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_auc_score, roc_curve
from os import path, system
import os
import sys
from array import array

#Additional functions (if needed)
from root_numpy import tree2array, fill_hist

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: testing e/g ID for HGCal L1Trigger"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#configure options
from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--modelAlgo', dest='modelAlgo', default='Histomax_vardrth10', help="Clustering algorithm used to train BDT" )
  parser.add_option('--inputs', dest='inputs', default='electron,neutrino', help="Inputs" )
  parser.add_option('--release', dest='release', default='103X', help="CMSSW release (for HGCal geometry)" )
  parser.add_option('--bdts', dest='bdt_list', default='electron_vs_neutrino,baseline', help="Colon separated list of BDTs to save value of" )
  return parser.parse_args()

(opt,args) = get_options()

#Input BDT options
modelAlgo = opt.modelAlgo
release = opt.release
inputs = opt.inputs.split(",")
bdt_list = []
for bdt in opt.bdt_list.split(":"): bdt_list.append( {'discriminator':bdt.split(",")[0],'config':bdt.split(",")[1]} )
#Extract path for bdts (.model files)
model_paths = {}
for b in bdt_list:
  bdt_name = "%s_%s"%(b['discriminator'],b['config'])
  if b['config'] == "baseline": model_paths[ bdt_name ] = os.environ['CMSSW_BASE']+"/src/L1Trigger/analysis/output/models/%s/egID_%s_%s.model"%(release,modelAlgo,b['discriminator'])
  else: model_paths[ bdt_name ] = os.environ['CMSSW_BASE']+"/src/L1Trigger/analysis/output/models/%s/egID_%s_%s_%s.model"%(release,modelAlgo,b['discriminator'],b['config'])

#Dictionaries for type mapping
typeMap = {"electron":"SingleElectron_FlatPt-2to100","photon":"SinglePhoton_FlatPt-8to150","pion":"SinglePion_FlatPt-2to100","neutrino":"SingleNeutrino"}
treeMap = {"electron":"e_sig","photon":"g_sig","pion":"pi_bkg","neutrino":"pu_bkg"}
procMap = {"electron":"signal", "photon":"signal", "pion":"background", "neutrino":"background"}

#get paths to input trees
input_paths = {}
for proc in inputs: input_paths[ proc ] = os.environ['CMSSW_BASE'] + "/src/L1Trigger/analysis/output/trees/%s/%s/%s/%s_%s_test.root"%(release,modelAlgo,typeMap[proc],typeMap[proc],modelAlgo)

# Variables to be bdt_name = "%s_%s"%(b['discriminator'],b['config'])included in test dataframe
frame_vars = ['cl3d_pt','cl3d_eta','cl3d_phi','cl3d_clusters_n','cl3d_showerlength','cl3d_coreshowerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_seetot','cl3d_seemax','cl3d_spptot','cl3d_sppmax','cl3d_szz','cl3d_srrtot','cl3d_srrmax','cl3d_srrmean','cl3d_emaxe','cl3d_bdteg']
frameDir = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/frames/%s/test"%release

# BDT input variables: HARDCODED
in_var_names = {'electron_vs_neutrino_baseline':['cl3d_coreshowerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_srrmean'],'electron_vs_neutrino_full':['cl3d_coreshowerlength','cl3d_showerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_szz','cl3d_srrmean','cl3d_srrtot','cl3d_seetot','cl3d_spptot']}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract test dataframes from inputs
testTotal = None
testFrames = {}
for proc, fileName in input_paths.iteritems():
  testFile = ROOT.TFile(fileName)
  testTree = testFile.Get( treeMap[proc] )
  # initialise new tree with only relevant variables
  _file = ROOT.TFile("tmp.root","RECREATE")
  _tree = ROOT.TTree("tmp","tmp")
  _vars = {}
  for var in frame_vars:
    _vars[ var ] = array( 'f', [-1.] )
    _tree.Branch( '%s'%var, _vars[ var ], '%s/F'%var )
  for cl3d in testTree:
    for var in frame_vars: _vars[ var ][0] = getattr( cl3d, '%s'%var )
    _tree.Fill()
  testFrames[proc] = pd.DataFrame( tree2array(_tree) )
  del _file
  del _tree
  system('rm tmp.root')

  #add columns to dataFrame to label clusters
  testFrames[proc]['proc'] = procMap[ proc ]
  testFrames[proc]['type'] = proc
  print "  --> Extracted testing dataframe from input files"

#create one total frame
testList = []
for proc in inputs: testList.append( testFrames[proc] )
testTotal = pd.concat( testList, sort=False )
del testFrames
print "  --> Created total dataframe (%s)"%opt.inputs

#save dataframe as pkl
testTotal.to_pickle("%s/testTotal_%s.pkl"%(frameDir,modelAlgo))
print "  --> Test dataframe saved as %s/testTotal_%s.pkl"%(frameDir,modelAlgo)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load models
models = {}
for b in bdt_list:
  bdt_name = "%s_%s"%(b['discriminator'],b['config'])
  models[ bdt_name ] = xg.Booster()
  models[ bdt_name ].load_model( model_paths[bdt_name] )

print "  --> Successfully loaded all models"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Configure datasets
label_encoder = LabelEncoder()
testTotal['proc'] = label_encoder.fit_transform(testTotal['proc'].values)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check performance of model
print "  --> Checking performance of model(s)"
test_X = {}
test_Y = {}
test_predY = {}
testing = {}

test_bkgEff = {}
test_sigEff = {}

bkgEff = {}
ubkgEff = {}

pT_bins = [20.,25.,30.,40.,50.,100.]
eta_bins = [1.5,2.,2.4,2.7,3.]
pT_midpoints = [22.5,27.5,35.,45.,75.]
eta_midpoints = [1.75,2.2,2.55,2.85]
epsilon_b = {}
u_epsilon_b = {}

for b in bdt_list:

  bdt_name = "%s_%s"%(b['discriminator'],b['config'])
  epsilon_b[bdt_name] = []
  u_epsilon_b[bdt_name] = []

  for bin_idx in range(len(eta_bins)-1):
  
    #Create tmp dataframe for variable bin
    tmpTotal = testTotal[abs(testTotal['cl3d_eta'])>eta_bins[bin_idx]]
    tmpTotal = tmpTotal[abs(tmpTotal['cl3d_eta'])<eta_bins[bin_idx+1]]
    total_name = bdt_name + "_eta_%sto%s"%(eta_bins[bin_idx],eta_bins[bin_idx+1])
    test_X[ total_name ] = tmpTotal[in_var_names[bdt_name]].values
    test_Y[ total_name ] = tmpTotal['proc'].values
    totalBkg = len( tmpTotal[tmpTotal['proc']==0] )
    testing[ total_name ] = xg.DMatrix( test_X[total_name], label=test_Y[total_name], feature_names=in_var_names[bdt_name] )
    test_predY[ total_name ] = models[ bdt_name ].predict( testing[total_name] ) 

    #Extract ROC curve points
    test_bkgEff[ total_name ], test_sigEff[ total_name ], nada = roc_curve( test_Y[total_name], test_predY[total_name] )
    #Find the background efficiency at 95% signal efficiency
    _idx = abs((test_sigEff[total_name]-0.95)).argmin()
    bkgEff[ total_name ] = test_bkgEff[total_name][_idx]
    epsilon_b[bdt_name].append( bkgEff[ total_name ] )
    ubkgEff[ total_name ] = math.sqrt( bkgEff[total_name]/totalBkg )
    u_epsilon_b[bdt_name].append( ubkgEff[ total_name ] )
  
  
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "  --> Performance of BDT: %s, for eta = %s - %s GeV"%(bdt_name,eta_bins[bin_idx],eta_bins[bin_idx+1])
    print "      * Total background: %g"%totalBkg
    print "      * AUC = %5.4f"%roc_auc_score( test_Y[total_name], test_predY[total_name] )
    print "      * Background eff at 95%% Signal eff = %5.4f +- %5.4f"%(bkgEff[total_name],ubkgEff[total_name])  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plot performance as a function of pT
#convert all to numpy array
colourMap = {"electron_vs_neutrino_baseline":"red","electron_vs_neutrino_full":"blue"}
plt.figure(1)
eta_midpoints = np.asarray(eta_midpoints)
eta_bin_widths = [0.25,0.2,0.15,0.15]
eta_bin_widths = np.asarray(eta_bin_widths)
for b in bdt_list:
  bdt_name = "%s_%s"%(b['discriminator'],b['config'])
  epsilon_b[bdt_name] = np.asarray(epsilon_b[bdt_name])
  u_epsilon_b[bdt_name] = np.asarray(u_epsilon_b[bdt_name]) 
  plt.errorbar( eta_midpoints, epsilon_b[bdt_name], xerr=eta_bin_widths, yerr=u_epsilon_b[bdt_name], label="%s"%bdt_name, color=colourMap[bdt_name], fmt='' )

plt.xlabel('|Cluster eta|')
plt.ylabel('Background Efficiency')
axes = plt.gca()
axes.set_xlim([1.5,3.])
axes.set_ylim([0,0.05])
plt.legend(bbox_to_anchor=(0.5,1.), loc='upper center')
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19/bdt_performance/new/performance_vs_var/eta.pdf" )
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19/bdt_performance/new/performance_vs_var/eta.png" )












