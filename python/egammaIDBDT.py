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

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: training of e/g ID for HGCal L1Trigger"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#configure options
from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--clusteringAlgo', dest='clusteringAlgo', default='default', help="Clustering algorithm with which to optimise BDT" )
  parser.add_option('--dataFrame', dest='dataFrame', default=None, help="Path to dataFrame if already produced" )
  parser.add_option('--trainParams',dest='trainParams', default=None, help='Comma-separated list of colon-separated pairs corresponding to parameters for the training')
  return parser.parse_args()

(opt,args) = get_options()

#Define map to extract TDirectory for different clustering algo
clusteringAlgo = opt.clusteringAlgo

#set up global variables
trainDir = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/trees/%s/"%clusteringAlgo
frameDir = trainDir.replace('trees','frames')
trainFrac = 0.7
validFrac = 0.1
#parameters for training model
if opt.trainParams: trainParams = opt.trainParams.split(',')

#get trees from files and put them in dataFrames
procFileMap = {'eg_signal':'SingleElectronPt5_100Eta1p6_2p8/SingleElectronPt5_100Eta1p6_2p8_%s.root'%clusteringAlgo,'eg_background_PU':'SinglePionPt25Eta1p6_2p8/SinglePionPt25Eta1p6_2p8_%s.root'%clusteringAlgo}
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
    trainFile = ROOT.TFile("%s%s"%(trainDir,fileName))
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
    
    elif proc == 'eg_background_PU' : 
      trainTree = trainFile.Get('egid_background')
      #initialise new tree with only relevant variables
      bkgPUFile = ROOT.TFile("tmp_bkg.root","RECREATE")
      bkgPUTree = ROOT.TTree("bkg_PU","bkg_PU")
      bkg_vars = {}
      for var in egID_vars:
        bkg_vars[ var ] = array( 'f', [-1.] )
        bkgPUTree.Branch( '%s'%var, bkg_vars[ var ], '%s/F'%var )
      for ev in trainTree:
        for var in egID_vars:
          bkg_vars[ var ][0] = getattr( ev, 'bkg_%s'%var )
        bkgPUTree.Fill()
      trainFrames[proc] = pd.DataFrame( tree2array( bkgPUTree ) )
      del bkgPUFile
      del bkgPUTree
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
  sum_bkg_PU = len( trainTotal[ trainTotal['proc']=="eg_background_PU" ] )
  weights = list( map( lambda a: 10000./sum_sig if a == "eg_signal" else 10000./sum_bkg_PU, trainTotal['proc'] ) )
  trainTotal['weight'] = weights

  #save dataframe as pkl file
  trainTotal.to_pickle("%strainTotal.pkl"%frameDir)
  print "  --> dataFrame saved as %strainTotal.pkl"%frameDir

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
egID_train_X, egID_valid_X, egID_test_X = np.split( egID_X, [egID_trainLimit,egID_validLimit] )
egID_train_Y, egID_valid_Y, egID_test_Y = np.split( egID_Y, [egID_trainLimit,egID_validLimit] )
egID_train_w, egID_valid_w, egID_test_w = np.split( egID_w, [egID_trainLimit,egID_validLimit] )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Building the discriminant
#training_egID = xg.DMatrix( egID_train_X, label=egID_train_Y, weight=egID_train_w, feature_names=egID_vars )
#testing_egID = xg.DMatrix( egID_test_X, label=egID_test_Y, weight=egID_test_w, feature_names=egID_vars )
training_egID = xg.DMatrix( egID_train_X, label=egID_train_Y, feature_names=egID_vars )
testing_egID = xg.DMatrix( egID_test_X, label=egID_test_Y, feature_names=egID_vars )

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
modelDir = trainDir.replace('trees','models')
egID_model.save_model( '%segID_%s.model'%(modelDir,clusteringAlgo) )
print "  --> Model saved: %segID_%s.model"%(modelDir,clusteringAlgo)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Check performance of model

egID_train_predY_xcheck = egID_model.predict( training_egID )
egID_test_predY = egID_model.predict( testing_egID )
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print " Performance: optimised w.r.t %s clustering"%clusteringAlgo
print "   * Training set ::: AUC = %5.4f"%( roc_auc_score( egID_train_Y, egID_train_predY_xcheck ) )
print "   * Test set     ::: AUC = %5.4f"%( roc_auc_score( egID_test_Y, egID_test_predY ) )
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

sys.exit(1)
