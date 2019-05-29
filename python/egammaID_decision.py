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
  parser.add_option('--modelAlgo', dest='modelAlgo', default='Histomax_vardr', help="Clustering algorithm with which BDT was trained" )
  parser.add_option('--signal', dest='signal', default='electron', help="Input signal type" )
  parser.add_option('--background', dest='background', default='neutrino', help="Input background type" )
  parser.add_option('--geometry', dest='geometry', default='v9', help="HGCal geometry config" )
  parser.add_option('--bdts', dest='bdt_list', default='electron_vs_neutrino,baseline,1', help="Comma separated list of BDTs to save value of. Last int defines if BDT is trained in different eta regions" )
  parser.add_option('--dataset', dest='dataset', default='test', help="Input dataset type" )
  return parser.parse_args()

(opt,args) = get_options()

# Input bdt options
modelAlgo = opt.modelAlgo
geometry = opt.geometry
bdt_list = []
for bdt in opt.bdt_list.split(":"): bdt_list.append( {'discriminator':bdt.split(",")[0],'config':bdt.split(",")[1],'split_by_eta':int(bdt.split(",")[2])} )
#Extract xmls for bdts
model_xmls = {}
for b in bdt_list: 
  bdt_name = "%s_%s"%(b['discriminator'],b['config'])
  if b['split_by_eta']:
    for eta_region in ['low','high']: model_xmls[ "%s_%seta"%(bdt_name,eta_region) ] = os.environ['CMSSW_BASE']+"/src/L1Trigger/egid_analysis/HGCal_L1T_egammaID/output/models/%s/egID_%s_%s_%seta.xml"%(geometry,modelAlgo,bdt_name,eta_region)
  else: model_xmls[ bdt_name ] = os.environ['CMSSW_BASE']+"/src/L1Trigger/egid_analysis/HGCal_L1T_egammaID/output/models/%s/egID_%s_%s.xml"%(geometry,modelAlgo,bdt_name)

#Dictionaries for type mapping
typeMap = {"electron":"SingleElectron_FlatPt-2to100","photon":"SinglePhoton_FlatPt-8to150","pion":"SinglePion_FlatPt-2to100","neutrino":"SingleNeutrino"}
treeMap = {"electron":"e_sig","photon":"g_sig","pion":"pi_bkg","neutrino":"pu_bkg"}

# map inputs to file
inputMap = {}
#Signal
inputMap[ opt.signal ] = os.environ['CMSSW_BASE'] + "/src/L1Trigger/egid_analysis/HGCal_L1T_egammaID/output/trees/%s/%s/%s/%s_%s_%s.root"%(geometry,modelAlgo,typeMap[opt.signal],typeMap[opt.signal],modelAlgo,opt.dataset)
#Background
inputMap[ opt.background ] = os.environ['CMSSW_BASE'] + "/src/L1Trigger/egid_analysis/HGCal_L1T_egammaID/output/trees/%s/%s/%s/%s_%s_%s.root"%(geometry,modelAlgo,typeMap[opt.background],typeMap[opt.background],modelAlgo,opt.dataset)

#map outputs to file
outputMap = {}
output_dir = os.environ['CMSSW_BASE'] + "/src/L1Trigger/egid_analysis/HGCal_L1T_egammaID/output/results/%s/%s"%(geometry,modelAlgo)
outputMap[ opt.signal ] = "%s/%s_%s_%s.root"%(output_dir,typeMap[opt.signal],modelAlgo,opt.dataset)
outputMap[ opt.background ] = "%s/%s_%s_%s.root"%(output_dir,typeMap[opt.background],modelAlgo,opt.dataset)

#BDT input variables: HARDCODED 
in_var_names = {'electron_vs_neutrino_baseline':['cl3d_coreshowerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_srrmean'],'electron_vs_neutrino_full':['cl3d_coreshowerlength','cl3d_showerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_szz','cl3d_srrmean','cl3d_srrtot','cl3d_seetot','cl3d_spptot'],'electron_vs_pion_full':['cl3d_coreshowerlength','cl3d_showerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_szz','cl3d_srrmean','cl3d_srrtot','cl3d_seetot','cl3d_spptot']}

# output variables
out_var_names = ['pt','eta','phi','clusters_n','showerlength','coreshowerlength','firstlayer','maxlayer','seetot','seemax','spptot','sppmax','szz','srrtot','srrmax','srrmean','emaxe','bdt_tpg']
#Add new bdt score
for b in bdt_list: 
  if b['split_by_eta']: out_var_names.append( 'bdt_%s_%s_%s_%s'%(modelAlgo,geometry,b['discriminator'],b['config']) )
  else: out_var_names.append( 'bdt_%s_%s_%s_%s_eta_inclusive'%(modelAlgo,geometry,b['discriminator'],b['config']) )

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
# Initialise BDTs
bdts = {}
bdt_in_var = {} #dictionary of dictionary of input variables
for b in bdt_list:
  bdt_name = "%s_%s"%(b['discriminator'],b['config'])
  if b['split_by_eta']:
    for eta_region in ['low','high']:
      bdts[ "%s_%seta"%(bdt_name,eta_region) ], bdt_in_var[ "%s_%seta"%(bdt_name,eta_region) ] = initialise_egID_BDT( model_xmls["%s_%seta"%(bdt_name,eta_region)], in_var_names[bdt_name] )
      print "  --> Initialised BDT: %s_%seta for Clustering: %s, Release: %s"%(bdt_name,eta_region,modelAlgo,geometry)
  else: 
    bdts[ bdt_name ], bdt_in_var[ bdt_name ] = initialise_egID_BDT( model_xmls[bdt_name], in_var_names[bdt_name] )
    print "  --> Initialised BDT: %s (inclusive in eta) for Clustering: %s, Release: %s"%(bdt_name,modelAlgo,geometry)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Evaluating the BDT

# Loop over input trees: signal and background
for proc, tree in t_in.iteritems():

  # Loop over clusters and fill output trees
  for cl3d in tree:
 
    #evaluate the bdts
    for b in bdt_list: 
      bdt_name = "%s_%s"%(b['discriminator'],b['config'])
 
      if b['split_by_eta']:
        #Low eta region: use low eta bdt
        if abs(cl3d.cl3d_eta) < 2.7:
          out_var[ 'bdt_%s_%s_%s_%s'%(modelAlgo,geometry,b['discriminator'],b['config']) ][0] = evaluate_egID_BDT( bdts["%s_loweta"%bdt_name], bdt_in_var["%s_loweta"%bdt_name], cl3d, in_var_names[bdt_name] )
        #High eta region: use high eta bdt
        else:
          out_var[ 'bdt_%s_%s_%s_%s'%(modelAlgo,geometry,b['discriminator'],b['config']) ][0] = evaluate_egID_BDT( bdts["%s_higheta"%bdt_name], bdt_in_var["%s_higheta"%bdt_name], cl3d, in_var_names[bdt_name] )
      #if trained using clusters with all eta
      else:
        out_var[ 'bdt_%s_%s_%s_%s_eta_inclusive'%(modelAlgo,geometry,b['discriminator'],b['config']) ][0] = evaluate_egID_BDT( bdts[bdt_name], bdt_in_var[bdt_name], cl3d, in_var_names[bdt_name] )

    #extract other variables from trees
    for var in out_var_names: 
      if "bdt_%s_%s"%(modelAlgo,geometry) in var: continue
      elif "bdt_tpg" in var: out_var[ var ][0] = getattr( cl3d, "cl3d_bdteg" )
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
