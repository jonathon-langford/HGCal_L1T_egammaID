import ROOT
import sys
import os
import math
from array import array
from optparse import OptionParser

#Get options from option parser
def get_options():
  parser = OptionParser()
  parser = OptionParser( usage="usage: HGCalL1T_cluster_selection.py <options>" )
  parser.add_option("--input", dest="input_type", default='electron', help="Input ntuple type")
  parser.add_option("--inputPath", dest="input_path", default='/eos/home-j/jlangfor/hgcal/l1/egid/may19', help="Path to input ntuples")
  parser.add_option("--pile_up", dest="pile_up", default='200', help="PU of input ntuples")
  parser.add_option("--geometry", dest="geometry", default='v9', help="HGCal geometry configuration")
  parser.add_option("--fileNumber", dest="file_number", default='1', help="Input ntuple number")
  parser.add_option("--maxEvents", dest="maxEvents", default=10, help="Maximum number of events to process")
  parser.add_option("--clusteringAlgo", dest="clusteringAlgo", default='Histomax_vardr', help="Clustering algorithm used in ntuple production")
  return parser.parse_args()

(opt,args) = get_options()

#Define map to extract TDirectory for different clustering algo
clusteringAlgoDirectory_map = {'TDR':'Floatingpoint8ThresholdRef2dRef3dGenclustersntuple','Histomax_vardr':'Floatingpoint8ThresholdDummyHistomaxvardrClustersntuple','STC_Histomax_vardr':'Floatingpoint8SupertriggercellDummyHistomaxvardrClustersntuple'}

# Extract options
clusteringAlgo = opt.clusteringAlgo
maxEvents = int( opt.maxEvents )

#For electron and pion inputs
if( opt.input_type == "electron" ): input_type = "SingleElectron_FlatPt-2to100"
elif( opt.input_type == "photon" ): input_type = "SinglePhoton_FlatPt-8to150"
elif( opt.input_type == "pion" ): input_type = "SinglePion_FlatPt-2to100"
elif( opt.input_type == "neutrino" ): input_type = "SingleNeutrino"
else:
  print "[ERROR] Input type (%s) not supported. Exiting..."%opt.input_type
  sys.exit(1)
fNumber = int(opt.file_number)

#Construct input file name
pile_up = opt.pile_up
input_path = opt.input_path
fInput = '%s/%s/%s_%sPU/ntuple_%g.root'%(input_path,opt.geometry,input_type,pile_up,fNumber)

#PID for gen matching: depends on input type
if( opt.input_type == "electron" ): gen_pdgid = [11]
elif( opt.input_type == "photon" ): gen_pdgid = [22]
elif( opt.input_type == "pion" ): gen_pdgid = [211]
elif( opt.input_type == "neutrino" ): gen_pdgid = []

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: cl3d tree generator"
print "\n  Input:"
print "    * Input file: %s %s"%(input_type,fInput)
print "    * Max events: %g"%maxEvents
print "    * Clustering Algorithm: %s"%clusteringAlgo
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#open ttree to read from
fin = ROOT.TFile.Open( fInput )
gen_tree = fin.Get("%s/HGCalTriggerNtuple"%clusteringAlgoDirectory_map[ 'TDR' ] )
cl3d_tree = fin.Get("%s/HGCalTriggerNtuple"%clusteringAlgoDirectory_map[ clusteringAlgo ] )

#########################################################################################
# Class definitions

# class for 3d cluster variable: only initiate if cl3d passes selection
class Cluster3D:

  #Constructor method: takes event and cluster number as input
  def __init__(self, _event, _ncl3d):
    #initialise TLorentzVector
    _p4 = ROOT.TLorentzVector()
    _p4.SetPtEtaPhiE( _event.cl3d_pt[_ncl3d], _event.cl3d_eta[_ncl3d], _event.cl3d_phi[_ncl3d], _event.cl3d_energy[_ncl3d] )
    self.P4 = _p4
    self.clusters_n       = _event.cl3d_clusters_n[_ncl3d]
    self.showerlength     = _event.cl3d_showerlength[_ncl3d]
    self.coreshowerlength = _event.cl3d_coreshowerlength[_ncl3d]
    self.firstlayer       = _event.cl3d_firstlayer[_ncl3d]
    self.maxlayer         = _event.cl3d_maxlayer[_ncl3d]
    self.seetot           = _event.cl3d_seetot[_ncl3d]
    self.seemax           = _event.cl3d_seemax[_ncl3d]
    self.spptot           = _event.cl3d_spptot[_ncl3d]
    self.sppmax           = _event.cl3d_sppmax[_ncl3d]
    self.szz              = _event.cl3d_szz[_ncl3d]
    self.srrtot           = _event.cl3d_srrtot[_ncl3d]
    self.srrmax           = _event.cl3d_srrmax[_ncl3d]
    self.srrmean          = _event.cl3d_srrmean[_ncl3d]
    self.emaxe            = _event.cl3d_emaxe[_ncl3d]
    self.bdteg            = _event.cl3d_bdteg[_ncl3d]
    self.quality          = _event.cl3d_quality[_ncl3d]
    
  #functions to act on 3D cluster



#########################################################################################
# Function definitions

#function to fill tree
def fill_cl3d( _cl3d, _out_var ):
  _out_var['pt'][0] = _cl3d.P4.Pt()
  _out_var['eta'][0] = _cl3d.P4.Eta()
  _out_var['phi'][0] = _cl3d.P4.Phi()
  _out_var['clusters_n'][0] = _cl3d.clusters_n
  _out_var['showerlength'][0] = _cl3d.showerlength
  _out_var['coreshowerlength'][0] = _cl3d.coreshowerlength
  _out_var['firstlayer'][0] = _cl3d.firstlayer
  _out_var['maxlayer'][0] = _cl3d.maxlayer
  _out_var['seetot'][0] = _cl3d.seetot
  _out_var['seemax'][0] = _cl3d.seemax
  _out_var['spptot'][0] = _cl3d.spptot
  _out_var['sppmax'][0] = _cl3d.sppmax
  _out_var['szz'][0] = _cl3d.szz
  _out_var['srrtot'][0] = _cl3d.srrtot
  _out_var['srrmax'][0] = _cl3d.srrmax
  _out_var['srrmean'][0] = _cl3d.srrmean
  _out_var['emaxe'][0] = _cl3d.emaxe
  _out_var['bdteg'][0] = _cl3d.bdteg
  _out_var['quality'][0] = _cl3d.quality
  tree.Fill()
  


#########################################################################################
# Configure output
print "Configuring output ntuple..."
#output ROOT file
#fout_id = os.environ['CMSSW_BASE'] + '/src/L1Trigger/egid_analysis/HGCal_L1T_egammaID/output/trees/%s/%s/%s/%s_%s_%g.root'%(opt.geometry,clusteringAlgo,input_type,input_type,clusteringAlgo,fNumber)
fout_id = os.environ['CMSSW_BASE'] + '/src/L1Trigger/egid_analysis/HGCal_L1T_egammaID/output/trees/new_egid/%s/%s/%s/%s_%s_%g.root'%(opt.geometry,clusteringAlgo,input_type,input_type,clusteringAlgo,fNumber)
fout = ROOT.TFile( fout_id, "RECREATE" )

#Initialise ttree and define variables
procToTree = {"electron":"e_sig","photon":"g_sig","neutrino":"pu_bkg","pion":"pi_bkg"}
tree = ROOT.TTree( procToTree[opt.input_type], procToTree[opt.input_type] )

#output variables
out_var_names = ['pt','eta','phi','clusters_n','showerlength','coreshowerlength','firstlayer','maxlayer','seetot','seemax','spptot','sppmax','szz','srrtot','srrmax','srrmean','emaxe','bdteg','quality']
out_var = {}
for var in out_var_names: out_var[var] = array( 'f', [0.] )
# Create branches in output tree
for var_name, var in out_var.iteritems(): tree.Branch( "cl3d_%s"%var_name, var, "cl3d_%s/F"%var_name )

#########################################################################################
# Loop over events

for ev_idx in range(cl3d_tree.GetEntries()):

  if ev_idx == maxEvents: break
  if maxEvents == -1:
    if ev_idx % 100 == 0: print "Processing event: %g/%g"%(ev_idx+1,cl3d_tree.GetEntries())
  else:
    if ev_idx % 100 == 0: print "Processing event: %g/%g"%(ev_idx+1,maxEvents)

  #Extract event info from both gen and cluster tree
  gen_tree.GetEntry( ev_idx )
  cl3d_tree.GetEntry( ev_idx )

  #Extract number of gen particles + cl3d in event
  N_gen = gen_tree.gen_n
  N_cl3d = cl3d_tree.cl3d_n

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # GEN-MATCHED CLUSTERS
  if( 'electron' in opt.input_type )|( 'photon' in opt.input_type )|( 'pion' in opt.input_type ):

    #Loop over gen-e/gamma in event
    for gen_idx in range( N_gen ): 
      if abs( gen_tree.gen_pdgid[gen_idx] ) in gen_pdgid:
        #define TLorentzVector for gen particle
        gen_p4 = ROOT.TLorentzVector()
        gen_p4.SetPtEtaPhiE( gen_tree.gen_pt[gen_idx], gen_tree.gen_eta[gen_idx], gen_tree.gen_phi[gen_idx], gen_tree.gen_energy[gen_idx] )    
        # require gen e/g/pi pT > 20 GeV
        if gen_p4.Pt() < 20.: continue

        # loop overi 3d clusters: save index of max pt cluster if in 
        cl3d_genMatched_maxpt_idx = -1
        cl3d_genMatched_maxpt = -999
        for cl3d_idx in range( N_cl3d ):
          #requre that cluster pt > 10 GeV
          if cl3d_tree.cl3d_pt[cl3d_idx] < 10.: continue
          #define TLorentxVector for cl3d
          cl3d_p4 = ROOT.TLorentzVector()
          cl3d_p4.SetPtEtaPhiE( cl3d_tree.cl3d_pt[cl3d_idx], cl3d_tree.cl3d_eta[cl3d_idx], cl3d_tree.cl3d_phi[cl3d_idx], cl3d_tree.cl3d_energy[cl3d_idx] )
          #Require cluster to be dR < 0.2 within gen particle
          if cl3d_p4.DeltaR( gen_p4 ) < 0.2: 
            #If pT of cluster is > present max then set 
            if cl3d_p4.Pt() > cl3d_genMatched_maxpt:
               cl3d_genMatched_maxpt = cl3d_p4.Pt()
               cl3d_genMatched_maxpt_idx = cl3d_idx

        # if cl3d idx has been set then add fill cluster to tree
        if cl3d_genMatched_maxpt_idx >= 0:
          cl3d = Cluster3D( cl3d_tree, cl3d_genMatched_maxpt_idx )
          fill_cl3d( cl3d, out_var )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #BACKGROUND CLUSTERS: PU
  else:

    #Loop over 3d clusters: if pT > 20 GeV then fill as background
    for cl3d_idx in range(0, N_cl3d ):

      if cl3d_tree.cl3d_pt[cl3d_idx] > 20.: 
        cl3d = Cluster3D( cl3d_tree, cl3d_idx )
        fill_cl3d( cl3d, out_var )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    

#end of events loop
#########################################################################################

fout.Write()
fout.Close()

#raw_input("Press any key to continue...")

