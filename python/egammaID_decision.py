# e/g ID for HGCal L1Trigger: extracting score from model using xml file and comparing with current

# IMPORTANT: tpg models (from Jean-Baptiste) have different weights for separate eta regions. Need to initialise with
#            two xml files and evaluate according to eta region

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Function for initialising BDT: returns BDT and dictionary of input variable arrays
def initialise_egID_BDT( i_xml, input_variable_names ):
  
  #book mva reader with input variables
  in_var = {}
  for var in input_variable_names: in_var[var] = array( 'f', [0.] )
  #initialse TMVA reader and add variables
  bdt_ = ROOT.TMVA.Reader()
  for var in input_variable_names: bdt_.AddVariable( var, in_var[var] )
  #Book mva using xml file
  bdt_.BookMVA( "BDT", i_xml )

  #Return initialised BDT and dictionary of input variables
  return bdt_, in_var


# Function for evaluating the BDT score for a 3d cluster: returns score
# Takes as input: bdt, dictionary of input variables and 3D cluster
#                 also the input variable names and whether signal/background (since name of variable is different in tree)
def evaluate_egID_BDT( _bdt, _bdt_input_variables, _cl3d, input_variable_names, _identifier, training="self" ):

  # Loop over input variables and extract value from input trees
  for var in input_variable_names: 
    #if tpg trained: input variables have cl3d_ prefix. Remove to extract variable from tree
    if training == "tpg": _bdt_input_variables[ var ][0] = getattr( _cl3d, "%s_%s"%(_identifier,var.replace("cl3d_","")) )
    else: _bdt_input_variables[ var ][0] = getattr( _cl3d, "%s_%s"%(_identifier,var) )

  #return BDT score
  return _bdt.EvaluateMVA("BDT")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Main script
print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: re trained e/g ID output score + comparing to current"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#configure options
from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--modelAlgs', dest='modelAlgs', default='default,Histomax_vardrth10,tpg_Histomax_vardrth10', help="Clustering algorithms with which BDTs were trained" )
  parser.add_option('--inputAlgo', dest='inputAlgo', default='Histomax_vardrth10', help="Clustering algorithm with which to check BDT performance" )
  parser.add_option('--inputSignalType', dest='inputSignalType', default='electron', help="Input signal type" )
  parser.add_option('--inputBackgroundType', dest='inputBackgroundType', default='neutrino', help="Input background type" )
  #parser.add_option('--signalSample', dest='signalSample', default='SingleElectronPt5_100Eta1p6_2p8', help="Input signal" )
  #parser.add_option('--backgroundSample', dest='backgroundSample', default='SinglePionPt25Eta1p6_2p8', help="Input background" )
  return parser.parse_args()

(opt,args) = get_options()

inputAlgo = opt.inputAlgo
#signalSample = opt.signalSample
#backgroundSample = opt.backgroundSample

# Extract model xml file
modelAlgs = opt.modelAlgs.split(",")
model_xmls = {}
for modelAlgo in modelAlgs:
  #for tpg models: low and high eta regions
  if "tpg" in modelAlgo:
    for eta in ["loweta","higheta"]: model_xmls[ "%s_%s"%(modelAlgo,eta) ] = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/models/tpg/%s/egamma_id_histomax_352_%s_v0.xml"%(modelAlgo,eta)
  #for self trained models
  else: model_xmls[ modelAlgo ] = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/models/%s/egID_%s_sig_%s_bkg_%s.xml"%(modelAlgo,modelAlgo,opt.inputSignalType,opt.inputBackgroundType)

#Define inputs
testDirMap = {} #first extract the test file directory
#Signal
if "hgcal_only" in opt.inputSignalType: testDirMap['signal'] = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/trees/hgcal_only/%s/"%inputAlgo
else: testDirMap['signal'] = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/trees/full_eta_range/%s/"%inputAlgo
#Background
if 'hgcal_only' in opt.inputBackgroundType: testDirMap['background'] = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/trees/hgcal_only/%s/"%inputAlgo
else: testDirMap['background'] = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/trees/full_eta_range/%s/"%inputAlgo

input_map = {}
#Signal
if opt.inputSignalType == "electron_hgcal_only": input_map['signal'] = "%sSingleElectronPt5_100Eta1p6_2p8/SingleElectronPt5_100Eta1p6_2p8_%s_test.root"%(testDirMap['signal'],inputAlgo)
elif opt.inputSignalType == "electron": input_map['signal'] = "%sSingleElectron_FlatPt-2to100/SingleElectron_FlatPt-2to100_%s_test.root"%(testDirMap['signal'],inputAlgo)
else:
  print "  --> [ERROR] Invalid signal type... Leaving"
  sys.exit(1)
#Background
if opt.inputBackgroundType == "neutrino": input_map['background'] = "%sSingleNeutrino/SingleNeutrino_%s_test.root"%(testDirMap['background'],inputAlgo)
else: 
  print "  --> [ERROR] Invalid background type... Leaving"
  sys.exit(1)
#input_map = {"signal":"/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/trees/%s/%s/%s_%s.root"%(inputAlgo,signalSample,signalSample,inputAlgo),"background":"/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/trees/%s/%s/%s_%s.root"%(inputAlgo,backgroundSample,backgroundSample,inputAlgo)}
tree_map = {"signal":"egid_signal","background":"egid_background"}

#Define output
output_map = {}
output_dir = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/egID_trees/"
#Signal
if opt.inputSignalType == "electron_hgcal_only": output_map['signal'] = "%s%s/SingleElectronPt5_100Eta1p6_2p8_%s_test.root"%(output_dir,inputAlgo,inputAlgo)
elif opt.inputSignalType == "electron": output_map['signal'] = "%s%s/SingleElectron_FlatPt-2to100_%s_test.root"%(output_dir,inputAlgo,inputAlgo)
#Background
if opt.inputBackgroundType == "neutrino": output_map['background'] = "%s%s/SingleNeutrino_%s_test.root"%(output_dir,inputAlgo,inputAlgo)
#output_map = {"signal":"/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/egID_trees/%s/%s_%s.root"%(inputAlgo,signalSample,inputAlgo),"background":"/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/egID_trees/%s/%s_%s.root"%(inputAlgo,backgroundSample,inputAlgo) }
#f_out = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/egID_trees/%s"

#BDT input variables
in_var_names = ['coreshowerlength','firstlayer','maxlayer','srrmean']
tpg_in_var_names = ['cl3d_firstlayer','cl3d_coreshowerlength','cl3d_maxlayer','cl3d_srrmean']

# output variables
out_var_names = ['pt','eta','phi','clusters_n','showerlength','coreshowerlength','firstlayer','maxlayer','seetot','seemax','spptot','sppmax','szz','srrtot','srrmax','srrmean','emaxe']
for modelAlgo in modelAlgs: out_var_names.append( 'bdtscore_%s'%modelAlgo )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Extract trees from ROOT files
f_in = {}
t_in = {}
for proc in input_map:
  f_in[ proc ] = ROOT.TFile.Open( input_map[proc] )
  t_in[ proc ] = f_in[proc].Get( tree_map[proc] )
print "  --> Input trees read successfully"
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Configure output ntuple
out_var = {}
for var in out_var_names: out_var[var] = array( 'f', [0.] )

f_out = {}
t_out = {}
for proc in output_map:
  f_out[proc] = ROOT.TFile.Open( output_map[proc], "RECREATE" )
  t_out[proc] = ROOT.TTree( tree_map[proc], tree_map[proc] )

for tree in t_out.itervalues():
  for var_name, var in out_var.iteritems(): tree.Branch( "%s"%var_name, var, "%s/F"%var_name )

print "  --> Output configured"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Initialise BDTs
bdts = {}
bdt_in_var = {} #dictionary of dictionary of input variables
for modelAlgo in modelAlgs:
  #for tpg model algorithms: initialise two BDTs depending on eta region
  if "tpg" in modelAlgo:
     for eta in ["loweta","higheta"]: bdts[ "%s_%s"%(modelAlgo,eta) ], bdt_in_var[ "%s_%s"%(modelAlgo,eta) ] = initialise_egID_BDT( model_xmls[ "%s_%s"%(modelAlgo,eta) ], tpg_in_var_names )
  #for self trained models:
  else: bdts[ modelAlgo ], bdt_in_var[ modelAlgo ] = initialise_egID_BDT( model_xmls[ modelAlgo ], in_var_names )

print "  --> Initialised BDTs: %s"%modelAlgs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Evaluating the BDT

# Loop over input trees: signal and background
for proc, tree in t_in.iteritems():

  #Define identifier for signal and background: required as name of vars is different in input trees
  if proc == "signal": identifier = "sig"
  elif proc == "background": identifier = "bkg"

  # Loop over clusters and fill output trees
  for cl3d in tree:
 
    #evaluate the bdts
    for modelAlgo in modelAlgs:
      #for tpg trained: evaluate score depending on eta region
      if "tpg" in modelAlgo:
        #low eta BDT
        if getattr( cl3d, "%s_eta"%identifier ) < 2.7: out_var[ 'bdtscore_%s'%modelAlgo ][0] = evaluate_egID_BDT( bdts[ "%s_loweta"%modelAlgo ], bdt_in_var[ "%s_loweta"%modelAlgo ], cl3d, tpg_in_var_names, identifier, training="tpg" )
        #high eta BDT
        else: out_var[ 'bdtscore_%s'%modelAlgo ][0] = evaluate_egID_BDT( bdts[ "%s_higheta"%modelAlgo ], bdt_in_var[ "%s_higheta"%modelAlgo ], cl3d, tpg_in_var_names, identifier, training="tpg" )
      
      #for self trained model
      else: out_var[ 'bdtscore_%s'%modelAlgo ][0] = evaluate_egID_BDT( bdts[ modelAlgo ], bdt_in_var[ modelAlgo ], cl3d, in_var_names, identifier )


    #extract other variables from trees
    for var in out_var_names: 
      if "bdtscore" in var: continue
      out_var[ var ][0] = getattr( cl3d, "%s_%s"%(identifier,var) )

    #write cluster with BDT scores to output tree
    t_out[proc].Fill()

print "  --> Evaluated BDT scores. Written out to trees."

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Close all files
for f in f_in.itervalues(): f.Close()
for f in f_out.itervalues():
  f.Write()
  f.Close()
