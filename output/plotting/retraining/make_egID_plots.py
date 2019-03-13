# e/g ID for HGCal Trigger: plotting tools

#usual imports
import ROOT
import numpy as np
import pandas as pd
import xgboost as xg
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
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
print "HGCal L1T Analysis: plotting tools for e/g ID"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#configure options
from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--trainingAlgo', dest='trainingAlgo', default='default', help="Clustering algorithm with which BDT was trained" )
  parser.add_option('--comparisonAlgo', dest='comparisonAlgo', default='Histomax_vardrth10', help="Clustering algorithm with which to compare" )
  parser.add_option('--comparisonModel', dest='comparisonModel', default='default', help="egamma ID BDT with which to compare" )
  return parser.parse_args()

(opt,args) = get_options()

trainingAlgo = opt.trainingAlgo
trainingModel = trainingAlgo
comparisonAlgo = opt.comparisonAlgo
comparisonModel = opt.comparisonModel

print "  --> Primary: algorithm = %s, model = %s"%(trainingAlgo,trainingModel)
print "  --> Comparison: algorithm = %s, model = %s"%(comparisonAlgo,comparisonModel)

#set up global variables
egID_vars = ['coreshowerlength','firstlayer','maxlayer','srrmean']
frameDir = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/frames/"
modelDir = frameDir.replace( 'frames', 'models' )
plotDir = "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/march19/SingleElectronPt5_100/retraining/"
trainFrac = 0.7
validFrac = 0.1

#Extract dataframes and model
trainTotal = pd.read_pickle( "%s%s/trainTotal.pkl"%(frameDir,trainingAlgo) )
comparisonTotal = pd.read_pickle( "%s%s/trainTotal.pkl"%(frameDir,comparisonAlgo) )
trainModel = xg.Booster()
trainModel.load_model( "%s%s/egID_%s.model"%(modelDir,trainingAlgo,trainingAlgo) )
if( trainingAlgo != comparisonModel ):
  compModel = xg.Booster()
  compModel.load_model( "%s%s/egID_%s.model"%(modelDir,comparisonModel,comparisonModel) )
else:
  compModel = trainModel

print "  --> Successfully loaded frames and model(s)"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Configure datasets

#Data set used to train model
label_encoder = LabelEncoder()

#Extract test set using same random seed (defined above)
theShape = trainTotal.shape[0]
egID_shuffle = np.random.permutation( theShape )
print "  --> Shuffle: ",egID_shuffle
egID_trainLimit = int(theShape*trainFrac)
egID_validLimit = int(theShape*validFrac)

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

#Comparison dataset: can use full set for testing
#if using different trained model to algorithm can use full set for testing
if( comparisonAlgo != comparisonModel ):
  print "  --> Using full comparison dataset"
  egID_comparison_test_X = comparisonTotal[ egID_vars ].values
  egID_comparison_test_Y = label_encoder.fit_transform(comparisonTotal[ 'proc' ].values)
  egID_comparison_test_w = comparisonTotal[ 'weight' ].values
  del comparisonTotal
else:
  #reset random seed
  np.random.seed(123456)
  theShape = comparisonTotal.shape[0]
  egID_comparison_shuffle = np.random.permutation( theShape )
  print "  --> Using comparison test set. Shuffle (comparison): ",egID_comparison_shuffle
  egID_comparison_trainLimit = int(theShape*trainFrac)
  egID_comparison_validLimit = int(theShape*validFrac)

  egID_comparison_X = comparisonTotal[ egID_vars ].values
  egID_comparison_Y = label_encoder.fit_transform(comparisonTotal[ 'proc' ].values)
  egID_comparison_w = comparisonTotal[ 'weight' ].values

  #Define different datasets
  egID_comparison_train_X, egID_comparison_valid_X, egID_comparison_test_X = np.split( egID_comparison_X, [egID_comparison_trainLimit,egID_comparison_validLimit] )
  egID_comparison_train_Y, egID_comparison_valid_Y, egID_comparison_test_Y = np.split( egID_comparison_Y, [egID_comparison_trainLimit,egID_comparison_validLimit] )
  egID_comparison_train_w, egID_comparison_valid_w, egID_comparison_test_w = np.split( egID_comparison_w, [egID_comparison_trainLimit,egID_comparison_validLimit] )
  del comparisonTotal
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Check performance of model

print "  --> Checking performance of model(s)"
testing_egID = xg.DMatrix( egID_test_X, label=egID_test_Y, feature_names=egID_vars )
comparison_egID = xg.DMatrix( egID_comparison_test_X, label=egID_comparison_test_Y, feature_names=egID_vars )
egID_test_predY = trainModel.predict( testing_egID )
egID_comparison_test_predY = compModel.predict( comparison_egID )

print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
if( trainingModel == comparisonModel ): print " Performance: optimised w.r.t. %s clustering"%trainingAlgo
else: print " Performance:"
print "   * %s (test set)     ::: AUC = %5.4f"%(trainingAlgo,roc_auc_score( egID_test_Y, egID_test_predY ))
if( comparisonAlgo != comparisonModel ): print "   * %s (full)     ::: AUC = %5.4f"%(comparisonAlgo,roc_auc_score( egID_comparison_test_Y, egID_comparison_test_predY ))
else: print "   * %s (test set)     ::: AUC = %5.4f"%(comparisonAlgo,roc_auc_score( egID_comparison_test_Y, egID_comparison_test_predY ))
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Plotting comparison: ROC curve

train_bkgEff, train_sigEff, nada = roc_curve( egID_test_Y, egID_test_predY )
comparison_bkgEff, comparison_sigEff, nada = roc_curve( egID_comparison_test_Y, egID_comparison_test_predY )
plt.figure(1)
if( trainingAlgo == comparisonAlgo ):
  plt.plot( 1-train_bkgEff, train_sigEff, label='Optimised w.r.t %s: AUC = %5.4f'%(trainingModel,roc_auc_score( egID_test_Y, egID_test_predY )), color='darkgreen' )
  plt.plot( 1-comparison_bkgEff, comparison_sigEff, label='Optimised w.r.t. %s: AUC = %5.4f'%(comparisonModel,roc_auc_score( egID_comparison_test_Y, egID_comparison_test_predY )), color='red' )
else:
  plt.plot( 1-train_bkgEff, train_sigEff, label='%s: AUC = %5.4f'%(trainingAlgo,roc_auc_score( egID_test_Y, egID_test_predY )), color='black' )
  plt.plot( 1-comparison_bkgEff, comparison_sigEff, label='%s: AUC = %5.4f'%(comparisonAlgo,roc_auc_score( egID_comparison_test_Y, egID_comparison_test_predY )), color='darkgreen' )
plt.xlabel('1 - Background efficiency')
plt.ylabel('Signal efficiency')
axes = plt.gca()
axes.set_xlim([0.6,1.05])
axes.set_ylim([0.6,1.05])
plt.legend(bbox_to_anchor=(0.05, 0.2, 0.5, 0.6), loc='lower left')
if( trainingModel == comparisonModel ): plt.title('ID optimised w.r.t %s clustering'%trainingAlgo)
elif( trainingAlgo == comparisonAlgo ): plt.title('Performance for %s clustering'%trainingAlgo )
#plt.show()
outputString = "%sROC_alg-%s_model-%s_vs_alg-%s_model-%s"%(plotDir,trainingAlgo,trainingModel,comparisonAlgo,comparisonModel)
plt.savefig("%s.png"%outputString)
plt.savefig("%s.pdf"%outputString)





