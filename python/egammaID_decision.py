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
import os
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
def evaluate_egID_BDT( _bdt, _bdt_input_variables, _cl3d, input_variable_names ):

  # Loop over input variables and extract values from input tree
  for var in input_variable_names: _bdt_input_variables[ var ][0] = getattr( _cl3d, "%s"%var )
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
  parser.add_option('--modelAlgo', dest='modelAlgo', default='Histomax_vardrth10', help="Clustering algorithm with which BDTs was trained" )
  parser.add_option('--signal', dest='signal', default='electron', help="Input signal type" )
  parser.add_option('--background', dest='background', default='neutrino', help="Input background type" )
  parser.add_option('--release', dest='release', default='103X', help="CMSSW release (for HGCal geometry)" )
  parser.add_option('--discriminators', dest='discriminators', default='electron_vs_neutrino', help='BDT trained to discriminate' )
  parser.add_option('--bdtNames', dest='bdt_names', default='baseline', help="Comma separated list of BDT names to s" )
  return parser.parse_args()

(opt,args) = get_options()

# Extract model xml file
modelAlgo = opt.modelAlgo
release = opt.release
model_xml = os.environ['CMSSW_BASE']+"/src/L1Trigger/analysis/output/models/%s/egID_%s_%s_vs_%s.xml"%(release,modelAlgo,opt.signal,opt.background)

#Dictionaries for type mapping
typeMap = {"electron":"SingleElectron_FlatPt-2to100","photon":"SinglePhoton_FlatPt-8to150","pion":"SinglePion_FlatPt-2to100","neutrino":"SingleNeutrino"}
treeMap = {"electron":"e_sig","photon":"g_sig","pion":"pi_bkg","neutrino":"pu_bkg"}

# map inputs to file
inputMap = {}
#Signal
inputMap[ opt.signal ] = os.environ['CMSSW_BASE'] + "/src/L1Trigger/analysis/output/trees/%s/%s/%s/%s_%s_test.root"%(release,modelAlgo,typeMap[opt.signal],typeMap[opt.signal],modelAlgo)
#Background
inputMap[ opt.background ] = os.environ['CMSSW_BASE'] + "/src/L1Trigger/analysis/output/trees/%s/%s/%s/%s_%s_test.root"%(release,modelAlgo,typeMap[opt.background],typeMap[opt.background],modelAlgo)

#map outputs to file
outputMap = {}
output_dir = os.environ['CMSSW_BASE'] + "/src/L1Trigger/analysis/output/egID_trees/%s/%s"%(release,modelAlgo)
outputMap[ opt.signal ] = "%s/%s_%s.root"%(output_dir,typeMap[opt.signal],modelAlgo)
outputMap[ opt.background ] = "%s/%s_%s.root"%(output_dir,typeMap[opt.background],modelAlgo)

#BDT input variables
in_var_names = ['cl3d_coreshowerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_srrmean']

# output variables
out_var_names = ['pt','eta','phi','clusters_n','showerlength','coreshowerlength','firstlayer','maxlayer','seetot','seemax','spptot','sppmax','szz','srrtot','srrmax','srrmean','emaxe','bdt_tpg']
#Add new bdt score
out_var_names.append( 'bdt_%s_%s_%s_vs_%s'%(modelAlgo,release,opt.signal,opt.background) )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Extract trees from ROOT files
f_in = {}
t_in = {}
for proc in inputMap:
  f_in[ proc ] = ROOT.TFile.Open( inputMap[proc] )
  t_in[ proc ] = f_in[proc].Get( treeMap[proc] )
print "  --> Input trees read successfully"
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Configure output ntuple
out_var = {}
for var in out_var_names: out_var[var] = array( 'f', [0.] )

f_out = {}
t_out = {}
for proc in outputMap:
  f_out[proc] = ROOT.TFile.Open( outputMap[proc], "RECREATE" )
  t_out[proc] = ROOT.TTree( treeMap[proc], treeMap[proc] )

for tree in t_out.itervalues():
  for var_name, var in out_var.iteritems(): tree.Branch( "cl3d_%s"%var_name, var, "cl3d_%s/F"%var_name )

print "  --> Output configured"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Initialise self trained BDT
bdt, bdt_in_var = initialise_egID_BDT( model_xml, in_var_names )
print "  --> Initialised BDTs: %s_%s_vs_%s"%(modelAlgo,opt.signal,opt.background)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Evaluating the BDT

# Loop over input trees: signal and background
for proc, tree in t_in.iteritems():

  # Loop over clusters and fill output trees
  for cl3d in tree:
 
    #evaluate the bdt
    out_var[ 'bdt_%s_%s_%s_vs_%s'%(modelAlgo,release,opt.signal,opt.background) ][0] = evaluate_egID_BDT( bdt, bdt_in_var, cl3d, in_var_names )

    #extract other variables from trees
    for var in out_var_names: 
      if var == 'bdt_%s_%s_%s_vs_%s'%(modelAlgo,release,opt.signal,opt.background): continue
      elif var == 'bdt_tpg': out_var[ var ][0] = getattr( cl3d, "cl3d_bdteg" )
      else: out_var[ var ][0] = getattr( cl3d, "cl3d_%s"%var )

    #write cluster with BDT scores to output tree
    t_out[proc].Fill()

print "  --> Evaluated BDT scores. Written out to trees."

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Close all files
for f in f_in.itervalues(): f.Close()
for f in f_out.itervalues():
  f.Write()
  f.Close()
