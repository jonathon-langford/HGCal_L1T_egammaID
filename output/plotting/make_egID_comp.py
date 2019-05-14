# Save pandas dataframe for signal and background combination
# Output np.array giving full set of points on ROC curve

#Usual imports
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

#Function to get path from input
def get_path( _i, _proc ): return os.environ['CMSSW_BASE'] + "/src/L1Trigger/egid_analysis/HGCal_L1T_egammaID/output/trees/%s/%s/%s/%s_%s_%s.root"%(_i['geometry'],_i['cl3d_algo'],typeMap[_i[_proc]],typeMap[_i[_proc]],_i['cl3d_algo'],_i['dataset'])

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: comparing different eg ID"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#configure options
from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--inputMap', dest='input_map', default='electron,neutrino,Histomax_vardr,v9,full,electron_vs_neutrino,tpg,blue', help="Input map: signal,background,clusteringAlgo,geometry,input dataset,bdt discriminator,bdt config,plot colour" )
  return parser.parse_args()

(opt,args) = get_options()

input_list = []
for _input in opt.input_map.split(":"):
  _i = _input.split(",")
  if len(_i)!=8:
    print "  --> [ERROR] Invalid input map"
    sys.exit(1)
  input_list.append({})
  input_list[-1]['signal'] = _i[0]
  input_list[-1]['background'] = _i[1]
  input_list[-1]['cl3d_algo'] = _i[2]
  input_list[-1]['geometry'] = _i[3]
  input_list[-1]['dataset'] = _i[4]
  input_list[-1]['discriminator'] = _i[5]
  input_list[-1]['config'] = _i[6]
  input_list[-1]['colour'] = _i[7]
  input_list[-1]['name'] = "%s_%s_%s_%s_%s_%s_%s"%(_i[0],_i[1],_i[2],_i[3],_i[4],_i[5],_i[6])

# Define dictionary for type mapping
typeMap = {"electron":"SingleElectron_FlatPt-2to100","photon":"SinglePhoton_FlatPt-8to150","pion":"SinglePion_FlatPt-2to100","neutrino":"SingleNeutrino"}
treeMap = {"electron":"e_sig","photon":"g_sig","pion":"pi_bkg","neutrino":"pu_bkg"}

#Variable to go in frame
#frame_vars = ['cl3d_pt','cl3d_eta','cl3d_phi','cl3d_clusters_n','cl3d_showerlength','cl3d_coreshowerlength','cl3d_firstlayer','cl3d_maxlayer','cl3d_seetot','cl3d_seemax','cl3d_spptot','cl3d_sppmax','cl3d_szz','cl3d_srrtot','cl3d_srrmax','cl3d_srrmean','cl3d_emaxe','cl3d_bdteg']
frame_vars = ['cl3d_pt','cl3d_eta','cl3d_bdteg']

#working points
working_points = [0.975,0.95,0.9]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define dictionaries to store values of sig eff, bkg eff and output bdt score
eff_signal = {}
eff_background = {}
bdt_points = {}
working_point_idx = {}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loop over inputs
for i in input_list:

  #Initiate list to store efficiencies
  eff_signal[i['name']] = [1.]
  eff_background[i['name']] = [1.]
  bdt_points[i['name']] = [-9999.]
  working_point_idx[i['name']] = []
  
  #Dictionary to store frames for sig and bkg
  frames = {}

  for proc in ['signal','background']:
    #extract signal and background files
    iFile = ROOT.TFile( get_path(i,proc) )
    iTree = iFile.Get( treeMap[i[proc]] )
    #initialise new tree with frame variables
    _file = ROOT.TFile("tmp.root","RECREATE")
    _tree = ROOT.TTree("tmp","tmp")
    _vars = {}
    for var in frame_vars:
      _vars[var] = array('f',[-1.])
      _tree.Branch( '%s'%var, _vars[var], '%s/F'%var )
    #Loop over clusters in input tree and fill new tree
    for cl3d in iTree:
      for var in frame_vars: _vars[ var ][0] = getattr( cl3d, '%s'%var )
      _tree.Fill()
    #Convert tree to dataframe
    frames[proc] = pd.DataFrame( tree2array(_tree) )
    del _file
    del _tree
    system('rm tmp.root')

    #add columns to dataframe to label clusters  
    frames[proc]['proc'] = proc
    frames[proc]['type'] = i[proc]

  print "  --> Extracted dataframes from input files"
  # Extract total number of signal and background
  N_sig_total,N_bkg_total = float(len(frames['signal'])),float(len(frames['background']))
  # Make one combined dataframe and sort according to the bdt score
  frames_list = []
  for proc in ['signal','background']: frames_list.append( frames[proc] )
  frameTotal = pd.concat( frames_list, sort=False )

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Sort the frame according to the bdt score
  if i['config'] == "tpg": bdt_var = "cl3d_bdteg"
  else: 
    print "  --> [ERROR] Other BDTs are not yet supported"
    sys.exit(1)
  
  frameTotal = frameTotal.sort_values(bdt_var)

  #Iterate over rows in dataframe and calc eff_sig and eff_bkg for a given bdt_point
  N_sig_running, N_bkg_running = 0., 0.
  for index, row in frameTotal.iterrows():
    #add one to running counters depending on proc
    if row['proc'] == "signal": N_sig_running += 1.
    elif row['proc'] == "background": N_bkg_running += 1.
    eff_s, eff_b = 1.-(N_sig_running/N_sig_total), 1.-(N_bkg_running/N_bkg_total)
    #only add one entry for each bdt output value, i.e. if same bdt value as previous then remove last entry
    if row[bdt_var] == bdt_points[i['name']][-1]: 
      bdt_points[i['name']] = bdt_points[i['name']][:-1]
      eff_signal[i['name']] = eff_signal[i['name']][:-1]
      eff_background[i['name']] = eff_background[i['name']][:-1]
    #add entry
    bdt_points[i['name']].append( row[bdt_var] )
    eff_signal[i['name']].append( eff_s )
    eff_background[i['name']].append( eff_b )
    
  #Convert the lists into numpy arrays
  bdt_points[i['name']] = np.asarray(bdt_points[i['name']])
  eff_signal[i['name']] = np.asarray(eff_signal[i['name']])
  eff_background[i['name']] = np.asarray(eff_background[i['name']])

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Extract the indices of the working points
  for wp in working_points: working_point_idx[i['name']].append( abs((eff_signal[i['name']]-wp)).argmin() )
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Print the information
  print "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  print "  --> Input: * signal        = %s"%i['signal']
  print "             * background    = %s"%i['background']
  print "             * cl3d_algo     = %s"%i['cl3d_algo']
  print "             * geometry      = %s"%i['geometry']
  print "             * dataset       = %s"%i['dataset']
  print ""
  print "  --> BDT:   * discriminator = %s"%i['discriminator']
  print "             * config        = %s"%i['config']
  print ""
  print "  --> Working points:"
  for wp_itr in range(len(working_points)):
    wp = working_points[wp_itr]
    print "             * At epsilon_s = %4.3f ::: BDT cut = %5.4f, epsilon_b = %5.4f"%(wp,bdt_points[i['name']][working_point_idx[i['name']][wp_itr]],eff_background[i['name']][working_point_idx[i['name']][wp_itr]])
  print "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot the ROC curves for the different inputs
plt.figure(1)
for i in input_list:
  plt.plot( eff_signal[i['name']], 1-eff_background[i['name']], label="Geometry:%s"%i['geometry'], color=i['colour'])
plt.xlabel('Signal Eff.')
plt.ylabel('1 - Background Eff.')
axes = plt.gca()
axes.set_xlim([0.5,1.1])
axes.set_ylim([0.5,1.1])
plt.legend(bbox_to_anchor=(0.05,0.1), loc='lower left')
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/ROC.png" )
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/ROC.pdf" )
    
