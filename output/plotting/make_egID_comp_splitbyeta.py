# Save pandas dataframe for signal and background combination
#Output np.array giving full set of points on ROC curve

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
def get_path( _i, _proc, where='trees' ): 
  if where == 'trees': return os.environ['CMSSW_BASE'] + "/src/L1Trigger/egid_analysis/HGCal_L1T_egammaID/output/trees/%s/%s/%s/%s_%s_%s.root"%(_i['geometry'],_i['cl3d_algo'],typeMap[_i[_proc]],typeMap[_i[_proc]],_i['cl3d_algo'],_i['dataset'])
  elif where == 'results': return os.environ['CMSSW_BASE'] + "/src/L1Trigger/egid_analysis/HGCal_L1T_egammaID/output/results/%s/%s/%s_%s_%s.root"%(_i['geometry'],_i['cl3d_algo'],typeMap[_i[_proc]],_i['cl3d_algo'],_i['dataset'])
  else: 
    print "  --> Invalid location. Leaving..."
    sys.exit(1)

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: comparing different eg ID"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#configure options
from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--inputMap', dest='input_map', default='electron,neutrino,Histomax_vardr,v9,full,electron_vs_neutrino,tpg,blue', help="Input map: signal,background,clusteringAlgo,geometry,input dataset,bdt discriminator,bdt config,plot colour" )
  parser.add_option('--no-output', dest='no_output', default=0, type="int", help='Supress output plots' )
  parser.add_option('--location', dest='location', default='trees', help='Location of input ntuples: trees/results' )
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
if opt.location == "trees": frame_vars = ['cl3d_pt','cl3d_eta','cl3d_bdteg']
elif opt.location == "results": frame_vars = ['cl3d_pt','cl3d_eta','cl3d_bdt_tpg','cl3d_bdt_Histomax_vardr_v9_electron_vs_neutrino_baseline','cl3d_bdt_Histomax_vardr_v9_electron_vs_neutrino_baseline_eta_inclusive','cl3d_bdt_Histomax_vardr_v9_electron_vs_neutrino_full','cl3d_bdt_Histomax_vardr_v9_electron_vs_neutrino_full_eta_inclusive']
else:
  print "  --> Invalid location. Leaving..."
  sys.exit(1)

#working points
working_points = [0.995,0.975,0.95,0.9]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define dictionaries to store values of sig eff, bkg eff and output bdt score
eff_signal = {}
eff_background = {}
bdt_points = {}
working_point_idx = {}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loop over inputs
for i in input_list:

  #Dictionary to store frames for sig and bkg
  frames = {}

  for proc in ['signal','background']:
    #extract signal and background files
    iFile = ROOT.TFile( get_path(i,proc,where=opt.location) )
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
  #N_sig_total,N_bkg_total = float(len(frames['signal'])),float(len(frames['background']))
  # Make one combined dataframe and sort according to the bdt score
  frames_list = []
  for proc in ['signal','background']: frames_list.append( frames[proc] )
  frameTotal = pd.concat( frames_list, sort=False )
 

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Sort the frame according to the bdt score
  if i['config'] == "tpg": 
    if opt.location == "trees": bdt_var = "cl3d_bdteg"
    elif opt.location == "results": bdt_var = "cl3d_bdt_tpg"
  else: bdt_var = "cl3d_bdt_%s_%s_%s_%s"%(i['cl3d_algo'],i['geometry'],i['discriminator'],i['config'])
  
  frameTotal = frameTotal.sort_values(bdt_var)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Impose eta requirements: only consider clusters with eta in 1.5-3 (others have bdt score  = -999.)
  frameTotal = frameTotal[ abs(frameTotal['cl3d_eta'])>=1.5 ]
  frameTotal = frameTotal[ abs(frameTotal['cl3d_eta'])<3. ]

  #Create high and low eta datasets and store in dictionary
  frames_splitByEta = {'low_eta':frameTotal[ abs(frameTotal['cl3d_eta'])<2.7 ], 'high_eta':frameTotal[ abs(frameTotal['cl3d_eta'])>=2.7 ]}

  for eta_region, frame in frames_splitByEta.iteritems():
    #create key name:
    key_name = i['name'] + "_%s"%eta_region

    #Initiate list to store efficiencies
    eff_signal[key_name] = [1.]
    eff_background[key_name] = [1.]
    bdt_points[key_name] = [-9999.]
    working_point_idx[key_name] = []

    #calculate total number of signal and background
    N_sig_total, N_bkg_total = float(len(frame[frame['proc']=='signal'])), float(len(frame[frame['proc']=='background']))

    #Iterate over rows in dataframe and calc eff_sig and eff_bkg for a given bdt_point
    N_sig_running, N_bkg_running = 0., 0.
    for index, row in frame.iterrows():
      #add one to running counters depending on proc
      if row['proc'] == "signal": N_sig_running += 1.
      elif row['proc'] == "background": N_bkg_running += 1.
      eff_s, eff_b = 1.-(N_sig_running/N_sig_total), 1.-(N_bkg_running/N_bkg_total)
      #only add one entry for each bdt output value, i.e. if same bdt value as previous then remove last entry
      if row[bdt_var] == bdt_points[key_name][-1]: 
        bdt_points[key_name] = bdt_points[key_name][:-1]
        eff_signal[key_name] = eff_signal[key_name][:-1]
        eff_background[key_name] = eff_background[key_name][:-1]
      #add entry
      bdt_points[key_name].append( row[bdt_var] )
      eff_signal[key_name].append( eff_s )
      eff_background[key_name].append( eff_b )
    
    #Convert the lists into numpy arrays
    bdt_points[key_name] = np.asarray(bdt_points[key_name])
    eff_signal[key_name] = np.asarray(eff_signal[key_name])
    eff_background[key_name] = np.asarray(eff_background[key_name])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Extract the indices of the working points
    for wp in working_points: working_point_idx[key_name].append( abs((eff_signal[key_name]-wp)).argmin() )
 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Print the information
    print "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "  --> Input: * signal        = %s"%i['signal']
    print "             * background    = %s"%i['background']
    print "             * cl3d_algo     = %s"%i['cl3d_algo']
    print "             * geometry      = %s"%i['geometry']
    print "             * dataset       = %s"%i['dataset']
    if eta_region == "low_eta": print "             * region        = 1.5 < eta < 2.7"
    elif eta_region == "high_eta": print "             * region        = 2.7 < eta < 3.0"
    print ""
    print "  --> BDT:   * discriminator = %s"%i['discriminator']
    print "             * config        = %s"%i['config']
    print ""
    print "  --> Working points:"
    for wp_itr in range(len(working_points)):
      wp = working_points[wp_itr]
      print "             * At epsilon_s = %4.3f ::: BDT cut = %8.7f, epsilon_b = %5.4f"%(wp,bdt_points[key_name][working_point_idx[key_name][wp_itr]],eff_background[key_name][working_point_idx[key_name][wp_itr]])
    print "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

if opt.no_output: sys.exit(1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SANDBOX FOR PLOTS
#plt.figure(1)
#for i in input_list:
#  key_name = i['name'] + "_low_eta"
#  if "eta_inclusive" in i['config']: _label = "%s (training inclusive in $\eta$)"%(i['config'].split("_")[0])
#  else: _label = i['config']
#  if i['config'] in ['baseline','full']: plt.plot( eff_signal[key_name], 1-eff_background[key_name], label=_label, color=i['colour'], linestyle='--' )
#  else: plt.plot( eff_signal[key_name], 1-eff_background[key_name], label=_label, color=i['colour'] )
#plt.xlabel('Signal Eff. ($\epsilon_s$)')
#plt.ylabel('1 - Background Eff. ($1-\epsilon_b$)')
#plt.title('1.5$ < |\eta| < $2.7, v9 geometry')
#axes = plt.gca()
#axes.set_xlim([0.5,1.1])
#axes.set_ylim([0.5,1.1])
#plt.legend(bbox_to_anchor=(0.05,0.1), loc='lower left')
#plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_new/ROC_loweta_all.png" )
#plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_new/ROC_loweta_all.pdf" )
#
#plt.figure(2)
#for i in input_list:
#  key_name = i['name'] + "_high_eta"
#  if "eta_inclusive" in i['config']: _label = "%s (training inclusive in $\eta$)"%(i['config'].split("_")[0])
#  else: _label = i['config']
#  if i['config'] in ['baseline','full']: plt.plot( eff_signal[key_name], 1-eff_background[key_name], label=_label, color=i['colour'], linestyle='--' )
#  else: plt.plot( eff_signal[key_name], 1-eff_background[key_name], label=_label, color=i['colour'] )
#plt.xlabel('Signal Eff.')
#plt.ylabel('1 - Background Eff.')
#plt.title('2.7$ < |\eta| < $3.0, v9 geometry')
#axes = plt.gca()
#axes.set_xlim([0.5,1.1])
#axes.set_ylim([0.5,1.1])
#plt.legend(bbox_to_anchor=(0.05,0.1), loc='lower left')
#plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_new/ROC_higheta_all.png" )
#plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_new/ROC_higheta_all.pdf" )

#sys.exit(1)

#plt.figure(3)
#for i in input_list:
#  key_name = i['name'] + "_low_eta"
#  plt.plot( bdt_points[key_name], eff_signal[key_name], label="Geometry:%s"%i['geometry'], color=i['colour'] )
#plt.xlabel('$e\gamma$-ID Output Score')
#plt.ylabel('Signal Eff. ($\epsilon_s$)')
#plt.title('1.5$ < |\eta| < $2.7')
#axes = plt.gca()
#axes.set_xlim([-0.75,0.3])
#axes.set_ylim([0.8,1.05])
#plt.legend(bbox_to_anchor=(0.05,0.1), loc='lower left')
#plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_tpg/epsilon_s_loweta.png" )
#plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_tpg/epsilon_s_loweta.pdf" )

#plt.figure(4)
#for i in input_list:
#  key_name = i['name'] + "_high_eta"
#  plt.plot( bdt_points[key_name], eff_signal[key_name], label="Geometry:%s"%i['geometry'], color=i['colour'] )
#plt.xlabel('$e\gamma$-ID Output Score')
#plt.ylabel('Signal Eff. ($\epsilon_s$)')
#plt.title('2.7$ < |\eta| < $3.0')
#axes = plt.gca()
#axes.set_xlim([-0.75,0.3])
#axes.set_ylim([0.8,1.05])
#plt.legend(bbox_to_anchor=(0.05,0.1), loc='lower left')
#plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_tpg_v9_only/epsilon_s_higheta.png" )
#plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_tpg_v9_only/epsilon_s_higheta.pdf" )

#plt.figure(5)
#for i in input_list:
#  key_name = i['name'] + "_low_eta"
#  plt.plot( bdt_points[key_name], eff_background[key_name], label="Geometry:%s"%i['geometry'], color=i['colour'] )
#plt.xlabel('$e\gamma$-ID Output Score')
#plt.ylabel('Background Eff. ($\epsilon_b$)')
#plt.title('1.5$ < |\eta| < $2.7')
#axes = plt.gca()
#axes.set_xlim([-0.75,0.3])
#axes.set_ylim([0,1.05])
#plt.legend(bbox_to_anchor=(0.05,0.1), loc='lower left')
#plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_tpg_v9_only/epsilon_b_loweta.png" )
#plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_tpg_v9_only/epsilon_b_loweta.pdf" )

#plt.figure(6)
#for i in input_list:
#  key_name = i['name'] + "_high_eta"
#  plt.plot( bdt_points[key_name], eff_background[key_name], label="Geometry:%s"%i['geometry'], color=i['colour'] )
#plt.xlabel('$e\gamma$-ID Output Score')
#plt.ylabel('Background Eff. ($\epsilon_b$)')
#plt.title('2.7$ < |\eta| < $3.0')
#axes = plt.gca()
#axes.set_xlim([-0.75,0.3])
#axes.set_ylim([0,1.05])
#plt.legend(bbox_to_anchor=(0.05,0.1), loc='lower left')
#plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_tpg_v9_only/epsilon_b_higheta.png" )
#plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_tpg_v9_only/epsilon_b_higheta.pdf" )

plt.figure(1)
for i in input_list:
  key_name = i['name'] + "_low_eta"
  plt.plot( eff_signal[key_name], 1-eff_background[key_name], label="Geometry:%s"%i['geometry'], color=i['colour'] )
plt.xlabel('Signal Eff. ($\epsilon_s$)')
plt.ylabel('1 - Background Eff. ($1-\epsilon_b$)')
plt.title('1.5$ < |\eta| < $2.7')
axes = plt.gca()
axes.set_xlim([0.5,1.05])
axes.set_ylim([0.5,1.05])
plt.legend(bbox_to_anchor=(0.05,0.1), loc='lower left')
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_tpg/ROC_loweta.png" )
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_tpg/ROC_loweta.pdf" )

plt.figure(2)
for i in input_list:
  key_name = i['name'] + "_high_eta"
  plt.plot( eff_signal[key_name], 1-eff_background[key_name], label="Geometry:%s"%i['geometry'], color=i['colour'] )
plt.xlabel('Signal Eff. ($\epsilon_s$)')
plt.ylabel('1 - Background Eff. ($1-\epsilon_b$)')
plt.title('2.7$ < |\eta| < $3.0')
axes = plt.gca()
axes.set_xlim([0.5,1.05])
axes.set_ylim([0.5,1.05])
plt.legend(bbox_to_anchor=(0.05,0.1), loc='lower left')
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_tpg/ROC_higheta.png" )
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_3/egid_tpg/ROC_higheta.pdf" )

sys.exit(1)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot the ROC curves for the different inputs
plt.figure(1)
for i in input_list:
  for eta_region in ['low_eta','high_eta']:
    key_name = i['name'] + "_%s"%eta_region
    if eta_region == "low_eta": plt.plot( eff_signal[key_name], 1-eff_background[key_name], label="Geometry:%s, 1.5$< |\eta| <$2.7"%i['geometry'], color=i['colour'])
    else: plt.plot( eff_signal[key_name], 1-eff_background[key_name], label="Geometry:%s, 2.7$< |\eta| <$3"%i['geometry'], color=i['colour'], linestyle='--')
plt.xlabel('Signal Eff.')
plt.ylabel('1 - Background Eff.')
axes = plt.gca()
axes.set_xlim([0.5,1.1])
axes.set_ylim([0.5,1.1])
plt.legend(bbox_to_anchor=(0.05,0.1), loc='lower left')
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/ROC.png" )
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/ROC.pdf" )
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plot signal efficiency + background rej as a function of BDT score
plt.figure(2)
for i in input_list:
  for eta_region in ['low_eta','high_eta']:
    key_name = i['name'] + "_%s"%eta_region
    if eta_region == "low_eta": plt.plot( bdt_points[key_name], eff_signal[key_name], label="Geometry:%s, 1.5$< |\eta| <$2.7"%i['geometry'], color=i['colour'])
    else: plt.plot( bdt_points[key_name], eff_signal[key_name], label="Geometry:%s, 2.7$< |\eta| <$3"%i['geometry'], color=i['colour'], linestyle='--')
plt.xlabel('BDT Output Score')
plt.ylabel('Signal Eff.')
axes = plt.gca()
axes.set_xlim([-0.75,0.3])
axes.set_ylim([0.8,1.05])
plt.legend(bbox_to_anchor=(0.05,0.1), loc='lower left')
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/epsilon_s.png" )
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/epsilon_s.pdf" )

plt.figure(3)
for i in input_list:
  for eta_region in ['low_eta','high_eta']:
    key_name = i['name'] + "_%s"%eta_region
    if eta_region == "low_eta": plt.plot( bdt_points[key_name], eff_background[key_name], label="Geometry:%s, 1.5$< |\eta| <$2.7"%i['geometry'], color=i['colour'])
    else: plt.plot( bdt_points[key_name], eff_background[key_name], label="Geometry:%s, 2.7$< |\eta| <$3"%i['geometry'], color=i['colour'], linestyle='--')
plt.xlabel('BDT Output Score')
plt.ylabel('Background Eff.')
axes = plt.gca()
axes.set_xlim([-0.75,0.3])
axes.set_ylim([0,1.05])
plt.legend(bbox_to_anchor=(0.9,0.9), loc='upper right')
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/epsilon_b.png" )
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/epsilon_b.pdf" )

plt.figure(4)
for i in input_list:
  if i['geometry'] == "v8":
    key_name = i['name'] + "_low_eta"
    plt.plot( bdt_points[key_name], eff_signal[key_name], label="Signal (%s)"%i['geometry'], color=i['colour'])
    plt.plot( bdt_points[key_name], eff_background[key_name], label="Background (%s)"%i['geometry'], color="black")
plt.xlabel('BDT Output Score')
plt.ylabel('Efficiency')
plt.title('Geometry: v8, 1.5$< |\eta| <$2.7')
axes = plt.gca()
axes.set_xlim([-0.75,0.3])
axes.set_ylim([0.,1.05])
plt.legend(bbox_to_anchor=(.05,.1), loc='lower left')
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/s-vs-b_v8_loweta.png" )
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/s-vs-b_v8_loweta.pdf" )

plt.figure(5)
for i in input_list:
  if i['geometry'] == "v9":
    key_name = i['name'] + "_low_eta"
    plt.plot( bdt_points[key_name], eff_signal[key_name], label="Signal (%s)"%i['geometry'], color=i['colour'])
    plt.plot( bdt_points[key_name], eff_background[key_name], label="Background (%s)"%i['geometry'], color="black")
plt.xlabel('BDT Output Score')
plt.ylabel('Efficiency')
plt.title('Geometry: v9, 1.5$< |\eta| <$2.7')
axes = plt.gca()
axes.set_xlim([-0.75,0.3])
axes.set_ylim([0.,1.05])
plt.legend(bbox_to_anchor=(.05,.1), loc='lower left')
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/s-vs-b_v9_loweta.png" )
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/s-vs-b_v9_loweta.pdf" )

plt.figure(6)
for i in input_list:
  if i['geometry'] == "v8":
    key_name = i['name'] + "_high_eta"
    plt.plot( bdt_points[key_name], eff_signal[key_name], label="Signal (%s)"%i['geometry'], color=i['colour'], linestyle='--')
    plt.plot( bdt_points[key_name], eff_background[key_name], label="Background (%s)"%i['geometry'], color="black", linestyle='--')
plt.xlabel('BDT Output Score')
plt.ylabel('Efficiency')
plt.title('Geometry: v8, 2.7$< |\eta| <$3.0')
axes = plt.gca()
axes.set_xlim([-0.75,0.3])
axes.set_ylim([0.,1.05])
plt.legend(bbox_to_anchor=(.05,.1), loc='lower left')
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/s-vs-b_v8_higheta.png" )
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/s-vs-b_v8_higheta.pdf" )

plt.figure(7)
for i in input_list:
  if i['geometry'] == "v9":
    key_name = i['name'] + "_high_eta"
    plt.plot( bdt_points[key_name], eff_signal[key_name], label="Signal (%s)"%i['geometry'], color=i['colour'], linestyle='--')
    plt.plot( bdt_points[key_name], eff_background[key_name], label="Background (%s)"%i['geometry'], color="black", linestyle='--')
plt.xlabel('BDT Output Score')
plt.ylabel('Efficiency')
plt.title('Geometry: v9, 2.7$< |\eta| <$3.0')
axes = plt.gca()
axes.set_xlim([-0.75,0.3])
axes.set_ylim([0.,1.05])
plt.legend(bbox_to_anchor=(.05,.1), loc='lower left')
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/s-vs-b_v9_higheta.png" )
plt.savefig( "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/may19_2/bdt_performance/geometry_comparison/s-vs-b_v9_higheta.pdf" )

