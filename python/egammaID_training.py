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
  parser.add_option('--clusteringAlgo', dest='clusteringAlgo', default='default', help="Clustering algorithm with which to optimise BDT" )
  parser.add_option('--inputSignalType', dest='inputSignalType', default='electron', help="Input signal type" )
  parser.add_option('--inputBackgroundType', dest='inputBackgroundType', default='electron', help="Input background type" )
  parser.add_option('--dataFrame', dest='dataFrame', default=None, help="Path to dataFrame if already produced" )
  parser.add_option('--trainParams',dest='trainParams', default=None, help='Comma-separated list of colon-separated pairs corresponding to parameters for the training')
  return parser.parse_args()

(opt,args) = get_options()

#Define map to extract TDirectory for different clustering algo
clusteringAlgo = opt.clusteringAlgo

#set up global variables
trainDirMap = {}
#Signal
if 'hgcal_only' in opt.inputSignalType: trainDirMap['eg_signal'] = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/trees/hgcal_only/%s/"%clusteringAlgo
else: trainDirMap['eg_signal'] = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/trees/full_eta_range/%s/"%clusteringAlgo
#Background
if 'hgcal_only' in opt.inputBackgroundType: trainDirMap['eg_background'] = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/trees/hgcal_only/%s/"%clusteringAlgo
else: trainDirMap['eg_background'] = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/trees/full_eta_range/%s/"%clusteringAlgo
frameDir = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/frames/%s/"%clusteringAlgo

trainFrac = 0.9
validFrac = 0.1

#parameters for training model
if opt.trainParams: trainParams = opt.trainParams.split(',')

#get trees from files and put them in dataFrames
procFileMap = {}
#Signal
if opt.inputSignalType == "electron_hgcal_only": procFileMap['eg_signal'] = "SingleElectronPt5_100Eta1p6_2p8/SingleElectronPt5_100Eta1p6_2p8_%s_train.root"%clusteringAlgo
elif opt.inputSignalType == "electron": procFileMap['eg_signal'] = "SingleElectron_FlatPt-2to100/SingleElectron_FlatPt-2to100_%s_train.root"%clusteringAlgo
else: 
  print "  [ERROR] Invalid signal type... Leaving"
  sys.exit(1)
#Background
if opt.inputBackgroundType == "neutrino": procFileMap['eg_background'] = "SingleNeutrino/SingleNeutrino_%s_train.root"%clusteringAlgo
else:
  print "  [ERROR] Invalid background type... Leaving"
  sys.exit(1)
procs = procFileMap.keys()

#define set of variables to use
egID_vars = ['coreshowerlength','firstlayer','maxlayer','srrmean']

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Extract dataFrames if they do not already exist
trainTotal = None
if not opt.dataFrame:
  trainFrames = {}
  #extract the trees: turn them into arrays
  for proc,fileName in procFileMap.iteritems():
    trainFile = ROOT.TFile("%s%s"%(trainDirMap[proc],fileName))
    if proc == 'eg_signal': 
      trainTree = trainFile.Get('egid_signal')
      #initialise new tree with only relevant variables
      sigFile = ROOT.TFile("tmp_sig.root","RECREATE")
      sigTree = ROOT.TTree("sig","sig")
      sig_vars = {}
      for var in egID_vars: 
        sig_vars[ var ] = array( 'f', [-1.] )
        sigTree.Branch( '%s'%var, sig_vars[ var ], '%s/F'%var )
      for ev in trainTree:
        for var in egID_vars:
          sig_vars[ var ][0] = getattr( ev, 'sig_%s'%var )
        sigTree.Fill()
      trainFrames[proc] = pd.DataFrame( tree2array( sigTree ) )
      del sigFile
      del sigTree 
      system('rm tmp_sig.root')
    
    elif proc == 'eg_background' : 
      trainTree = trainFile.Get('egid_background')
      #initialise new tree with only relevant variables
      bkgFile = ROOT.TFile("tmp_bkg.root","RECREATE")
      bkgTree = ROOT.TTree("bkg","bkg")
      bkg_vars = {}
      for var in egID_vars:
        bkg_vars[ var ] = array( 'f', [-1.] )
        bkgTree.Branch( '%s'%var, bkg_vars[ var ], '%s/F'%var )
      for ev in trainTree:
        for var in egID_vars:
          bkg_vars[ var ][0] = getattr( ev, 'bkg_%s'%var )
        bkgTree.Fill()
      trainFrames[proc] = pd.DataFrame( tree2array( bkgTree ) )
      del bkgFile
      del bkgTree
      system('rm tmp_bkg.root')

    #add column to dataFrame to label clusters
    trainFrames[proc]['proc'] = proc
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
  sum_sig = len( trainTotal[ trainTotal['proc']=="eg_signal" ] )
  sum_bkg = len( trainTotal[ trainTotal['proc']=="eg_background" ] )
  weights = list( map( lambda a: (sum_sig+sum_bkg)/sum_sig if a == "eg_signal" else (sum_sig+sum_bkg)/sum_bkg, trainTotal['proc'] ) )
  trainTotal['weight'] = weights

  #save dataframe as pkl file
  trainTotal.to_pickle("%strainTotal_sig_%s_bkg_%s.pkl"%(frameDir,opt.inputSignalType,opt.inputBackgroundType))
  print "  --> dataFrame saved as %strainTotal_sig_%s_bkg_%s.pkl"%(frameDir,opt.inputSignalType,opt.inputBackgroundType)

#read in dataFrame if already made
else:
  trainTotal = pd.read_pickle( opt.dataFrame )
  print "  --> Successfully loaded dataFrame"

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
egID_X = trainTotal[ egID_vars ].values
egID_Y = label_encoder.fit_transform(trainTotal[ 'proc' ].values)
egID_w = trainTotal[ 'weight' ].values
del trainTotal

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
training_egID = xg.DMatrix( egID_train_X, label=egID_train_Y, weight=egID_train_w, feature_names=egID_vars )
valid_egID = xg.DMatrix( egID_valid_X, label=egID_valid_Y, weight=egID_valid_w, feature_names=egID_vars )
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
modelDir = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/models/%s/"%clusteringAlgo
egID_model.save_model( '%segID_%s_sig_%s_bkg_%s.model'%(modelDir,clusteringAlgo,opt.inputSignalType,opt.inputBackgroundType) )
print "  --> Model saved: %segID_%s_sig_%s_bkg_%s.model"%(modelDir,clusteringAlgo,opt.inputSignalType,opt.inputBackgroundType)
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
