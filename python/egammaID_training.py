# e/g ID for HGCal L1Trigger: training MVA using shower shape variables

#usual imports
import ROOT
import numpy as np
import pandas as pd
import xgboost as xg
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

#Set numpy random seed
np.random.seed(123456)

#For verbose output
debug_ = True

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: training of e/g ID for HGCal L1Trigger"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#configure options
from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--clusteringAlgo', dest='clusteringAlgo', default='Histomax_vardr', help="Clustering algorithm with which to optimise BDT" )
  parser.add_option('--signal', dest='signal', default='electron', help="Input signal type" )
  parser.add_option('--background', dest='background', default='neutrino', help="Input background type" )
  parser.add_option('--geometry', dest='geometry', default='v9', help="HGCal geometry configuration" )
  parser.add_option('--bdtConfig', dest='bdtConfig', default='baseline', help="BDT config" )
  parser.add_option('--etaRegion', dest='etaRegion', default='total', help="Train bdt spearately in different eta regions" )
  parser.add_option('--dataFrame', dest='dataFrame', default=None, help="Path to dataFrame if already produced" )
  parser.add_option('--trainParams',dest='trainParams', default=None, help='Comma-separated list of colon-separated pairs corresponding to parameters for the training')
  return parser.parse_args()

(opt,args) = get_options()

#Extract options
clusteringAlgo = opt.clusteringAlgo
geometry = opt.geometry
eta_region = opt.etaRegion
if eta_region == "low": print "  --> Training on 3d clusters with: |eta| < 2.7"
elif eta_region == "high": print "  --> Training on 3d clusters with: |eta| > 2.7"
else: print "  --> No eta requirement on 3d clusters"
if opt.trainParams: trainParams = opt.trainParams.split(',')

#Dictionaries for type mapping
typeMap = {"electron":"SingleElectron_FlatPt-2to100","photon":"SinglePhoton_FlatPt-8to150","pion":"SinglePion_FlatPt-2to100","neutrino":"SingleNeutrino"}
treeMap = {"electron":"e_sig","photon":"g_sig","pion":"pi_bkg","neutrino":"pu_bkg"}
procMap = {"electron":"signal", "photon":"signal", "pion":"background", "neutrino":"background"}

#set up global variables
frameDir = os.environ['CMSSW_BASE']+"/src/L1Trigger/HGCal_L1T_egammaID/output/frames/%s"%geometry
modelDir = os.environ['CMSSW_BASE']+"/src/L1Trigger/HGCal_L1T_egammaID/output/models/%s"%geometry

# Training a validation fractions
trainFrac = 0.9
validFrac = 0.1

#get trees from files and put them in dataFrames
procFileMap = {}
#Signal
procFileMap[ opt.signal ] = os.environ['CMSSW_BASE'] + "/src/L1Trigger/HGCal_L1T_egammaID/output/trees/%s/%s/%s/%s_%s_train.root"%(geometry,clusteringAlgo,typeMap[opt.signal],typeMap[opt.signal],clusteringAlgo)
#Background
procFileMap[ opt.background ] = os.environ['CMSSW_BASE'] + "/src/L1Trigger/HGCal_L1T_egammaID/output/trees/%s/%s/%s/%s_%s_train.root"%(geometry,clusteringAlgo,typeMap[opt.background],typeMap[opt.background],clusteringAlgo)

procs = procFileMap.keys()

#define set of variables to use
egID_vars = {'electron_vs_neutrino_baseline':['cl3d_coreshowerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_srrmean'],'electron_vs_neutrino_full':['cl3d_coreshowerlength','cl3d_showerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_szz','cl3d_srrmean','cl3d_srrtot','cl3d_seetot','cl3d_spptot']}

#Define bdt name
bdt_name = "%s_vs_%s_%s"%(opt.signal,opt.background,opt.bdtConfig)
if bdt_name in egID_vars: print "  --> BDT input variables defined"
else: 
  print "  --> No input variables defined for this BDT config"
  sys.exit(1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Extract dataFrames
trainTotal = None
trainFrames = {}
#extract the trees: turn them into arrays
for proc,fileName in procFileMap.iteritems():
  print "  --> [DEBUG] process = %s, fileName = %s"%(proc,fileName)
  trainFile = ROOT.TFile("%s"%fileName)
  trainTree = trainFile.Get( treeMap[proc] )
  #initialise new tree with only relevant variables
  _file = ROOT.TFile("tmp.root","RECREATE")
  _tree = ROOT.TTree("tmp","tmp")
  _vars = {}
  for var in egID_vars[bdt_name]: 
    _vars[ var ] = array( 'f', [-1.] )
    _tree.Branch( '%s'%var, _vars[ var ], '%s/F'%var )
  for ev in trainTree:
    #Impose eta requirements
    if eta_region == "low":
      if( abs( ev.cl3d_eta ) < 2.7 ): 
        for var in egID_vars[bdt_name]: _vars[ var ][0] = getattr( ev, '%s'%var )
        _tree.Fill()
    elif eta_region == "high":
      if( abs( ev.cl3d_eta ) > 2.7 ):
        for var in egID_vars[bdt_name]: _vars[ var ][0] = getattr( ev, '%s'%var )
        _tree.Fill()
    else: 
      for var in egID_vars[bdt_name]: _vars[ var ][0] = getattr( ev, '%s'%var )
      _tree.Fill()
  trainFrames[proc] = pd.DataFrame( tree2array( _tree ) )
  del _file
  del _tree 
  system('rm tmp.root')

  #add column to dataFrame to label clusters
  trainFrames[proc]['proc'] = procMap[ proc ]
  print "  * Extracted %s dataFrame from tree"%proc

#create one total frame
trainList = []
for proc in procs: trainList.append( trainFrames[proc] )
trainTotal = pd.concat( trainList, sort=False )
del trainFrames
print "  --> Created total dataFrame (sig+bkg)"

#event filter

#add extra info to dataframes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Re-weighting
#1) assume sig and bkg equally likely, weight both to 1 (probably incorrect)
sum_sig = len( trainTotal[ trainTotal['proc']=="signal" ].index )
print "  --> [DEBUG] N_signal = %2.1f"%sum_sig
sum_bkg = len( trainTotal[ trainTotal['proc']=="background" ].index )
print "  --> [DEBUG] N_background = %2.1f"%sum_bkg
weights = list( map( lambda a: (sum_sig+sum_bkg)/sum_sig if a == "signal" else (sum_sig+sum_bkg)/sum_bkg, trainTotal['proc'] ) )
trainTotal['weight'] = weights

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Configure datasets: shuffle to get random permutation of inputs
label_encoder = LabelEncoder()

# Shuffle to get train and validation sets
theShape = trainTotal.shape[0]
egID_shuffle = np.random.permutation( theShape )
print "  --> Shuffle: ",egID_shuffle
egID_trainLimit = int(theShape*trainFrac)
egID_validLimit = int(theShape*validFrac)

#Set up various datasets for training BDT
egID_X = trainTotal[ egID_vars[bdt_name] ].values
egID_Y = label_encoder.fit_transform(trainTotal[ 'proc' ].values)
egID_w = trainTotal[ 'weight' ].values
#del trainTotal

#Shuffle 
egID_X = egID_X[ egID_shuffle ] 
egID_Y = egID_Y[ egID_shuffle ] 
egID_w = egID_w[ egID_shuffle ] 
# Define different datasets
egID_train_X, egID_valid_X, dummy_X = np.split( egID_X, [egID_trainLimit,egID_validLimit+egID_trainLimit] )
egID_train_Y, egID_valid_Y, dummy_Y = np.split( egID_Y, [egID_trainLimit,egID_validLimit+egID_trainLimit] )
egID_train_w, egID_valid_w, dummy_w = np.split( egID_w, [egID_trainLimit,egID_validLimit+egID_trainLimit] )

if( debug_ ):
  print "  --> [DEBUG] Size of training set:", len(egID_train_X)
  print "  --> [DEBUG] Size of validation set:", len(egID_valid_X) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Building the discriminant
training_egID = xg.DMatrix( egID_train_X, label=egID_train_Y, weight=egID_train_w, feature_names=egID_vars[bdt_name] )
valid_egID = xg.DMatrix( egID_valid_X, label=egID_valid_Y, weight=egID_valid_w, feature_names=egID_vars[bdt_name] )
#training_egID = xg.DMatrix( egID_train_X, label=egID_train_Y, feature_names=egID_vars )
#valid_egID = xg.DMatrix( egID_valid_X, label=egID_valid_Y, feature_names=egID_vars )

#extract training parameters for model
trainParams = {}
trainParams['objective'] = 'binary:logistic'
trainParams['nthread'] = 1
paramExt = ''
if opt.trainParams:
  paramExt = '__'
  for paramPair in trainParams:
    param = paramPair.split(":")[0]
    value = paramPair.split(":")[1]
    trainParams[param] = value
    paramExt += '%s_%s__'%(param,value)
  paramExt = paramExt[:-2]
print "  --> Training the model: %s"%trainParams
egID_model = xg.train( trainParams, training_egID )
print "  --> Done"

# Save the model here
if eta_region in ['low','high']:
  egID_model.save_model( '%s/egID_%s_%s_%seta.model'%(modelDir,clusteringAlgo,bdt_name,eta_region) )
  print "  --> Model saved: %s/egID_%s_%s_%seta.model"%(modelDir,clusteringAlgo,bdt_name,eta_region)
else:
  egID_model.save_model( '%s/egID_%s_%s.model'%(modelDir,clusteringAlgo,bdt_name) )
  print "  --> Model saved: %s/egID_%s_%s.model"%(modelDir,clusteringAlgo,bdt_name)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Check performance of model using validation set
egID_train_predY_xcheck = egID_model.predict( training_egID )
egID_valid_predY = egID_model.predict( valid_egID )

print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print " Performance: optimised w.r.t %s clustering"%clusteringAlgo
print "   * Training set   ::: AUC = %5.4f"%( roc_auc_score( egID_train_Y, egID_train_predY_xcheck ) )
print "   * Validation set ::: AUC = %5.4f"%( roc_auc_score( egID_valid_Y, egID_valid_predY ) )
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

sys.exit(1)
