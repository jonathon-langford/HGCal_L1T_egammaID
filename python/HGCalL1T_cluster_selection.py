import ROOT
import sys
import os
import math
from array import array
from optparse import OptionParser

#Get options from option parser
def get_options():
  parser = OptionParser()
  parser = OptionParser( usage="usage: HGCalL1T_cluster_selection.py -i <input ntuple type> -f <input ntuple number> -n <max events> -cl <clustering algorithm>" )
  parser.add_option("-i", "--input", dest="input_type", default='electron', help="Input ntuple type")
  parser.add_option("-f", "--fileNumber", dest="file_number", default='1', help="Input ntuple number")
  parser.add_option("-n", "--maxEvents", dest="maxEvents", default=10, help="Maximum number of events to process")
  parser.add_option("-c", "--clusteringAlgo", dest="clusteringAlgo", default='default', help="Clustering algorithm used in ntuple production")
  return parser.parse_args()

(opt,args) = get_options()

#Define map to extract TDirectory for different clustering algo
clusteringAlgoDirectory_map = {'default':'Floatingpoint8ThresholdRef2dRef3dGenclustersntuple','Histomax':'Floatingpoint8ThresholdDummyHistomaxClustersntuple','Histomax_vardrth0':'Floatingpoint8ThresholdDummyHistomaxvardrth0Clustersntuple','Histomax_vardrth10':'Floatingpoint8ThresholdDummyHistomaxvardrth10Clustersntuple','Histomax_vardrth20':'Floatingpoint8ThresholdDummyHistomaxvardrth20Clustersntuple'}
clusteringAlgo = opt.clusteringAlgo

maxEvents = int( opt.maxEvents )

#For electron and pion inputs
if( opt.input_type == "electron" ): input_type = "SingleElectronPt5_100Eta1p6_2p8"
elif( opt.input_type == "pion" ): input_type = "SinglePionPt25Eta1p6_2p8"
else:
  print "[ERROR] Input type (%s) not supported. Exiting..."%opt.input_type
  sys.exit(1)
fNumber = int(opt.file_number)

#Construct input file name
#fInput = os.environ['CMSSW_BASE'] + "/src/L1Trigger/analysis/ntuples/%s/ntuple_%g.root"%(input_type,fNumber)
fInput = '/eos/home-j/jlangfor/hgcal/l1/egid/93X/ntuples/%s_200PU/ntuple_%g.root'%(input_type,fNumber)

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: signal and background tree generator"
print "\n  Input:"
print "    * Input file: %s %s"%(input_type,fInput)
print "    * Max events: %g"%maxEvents
print "    * Clustering Algorithm: %s"%clusteringAlgo
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#open ttree to read from
fin = ROOT.TFile.Open( fInput )
gen_tree = fin.Get("%s/HGCalTriggerNtuple"%clusteringAlgoDirectory_map[ 'default' ] )
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

# function to return number of 3d clusters for a given event
def N_cl3d( _event ): return _event.cl3d_n

# function for selection of 3d cluster
def selection( _event, _ncl3d ):
  if( _event.cl3d_pt[_ncl3d] < 10. ): return False
  else: return True

#function to return true if 3D cluster is gen matched
def isGenMatched( _genEvent, _cl3d, _type="eg", min_dR=0.1 ):
  if( _type == "eg" ): gen_pids = [11,22] #pids for egamma
  else: 
    print "[ERROR] type %s is not supported. Returning false."
    return False

  #Loop over gen particles: using collection specified in options
  for gen_idx in range( _genEvent.gen_n ):
    #if pdgid in gen_pids
    if abs(_genEvent.gen_pdgid[gen_idx]) in gen_pids:
      #define TLorentzVector for gen particle
      gen_p4 = ROOT.TLorentzVector()
      gen_p4.SetPtEtaPhiE( _genEvent.gen_pt[gen_idx], _genEvent.gen_eta[gen_idx], _genEvent.gen_phi[gen_idx], _genEvent.gen_energy[gen_idx] )
      if( _cl3d.P4.DeltaR( gen_p4 ) < min_dR )&( gen_p4.Pt() > 10. ): return True

  #if no gen matched then return false
  return False


#function to fill signal/bkg
def fillSignal( _cl3d ):
  sig_pt[0] = _cl3d.P4.Pt()
  sig_eta[0] = _cl3d.P4.Eta()
  sig_phi[0] = _cl3d.P4.Phi()
  sig_clusters_n[0] = _cl3d.clusters_n
  sig_showerlength[0] = _cl3d.showerlength
  sig_coreshowerlength[0] = _cl3d.coreshowerlength
  sig_firstlayer[0] = _cl3d.firstlayer
  sig_maxlayer[0] = _cl3d.maxlayer
  sig_seetot[0] = _cl3d.seetot
  sig_seemax[0] = _cl3d.seemax
  sig_spptot[0] = _cl3d.spptot
  sig_sppmax[0] = _cl3d.sppmax
  sig_szz[0] = _cl3d.szz
  sig_srrtot[0] = _cl3d.srrtot
  sig_srrmax[0] = _cl3d.srrmax
  sig_srrmean[0] = _cl3d.srrmean
  sig_emaxe[0] = _cl3d.emaxe
  sig_bdteg[0] = _cl3d.bdteg
  sig_quality[0] = _cl3d.quality
  tree_sig.Fill()

def fillBackground( _cl3d ):
  bkg_pt[0] = _cl3d.P4.Pt()
  bkg_eta[0] = _cl3d.P4.Eta()
  bkg_phi[0] = _cl3d.P4.Phi()
  bkg_clusters_n[0] = _cl3d.clusters_n
  bkg_showerlength[0] = _cl3d.showerlength
  bkg_coreshowerlength[0] = _cl3d.coreshowerlength
  bkg_firstlayer[0] = _cl3d.firstlayer
  bkg_maxlayer[0] = _cl3d.maxlayer
  bkg_seetot[0] = _cl3d.seetot
  bkg_seemax[0] = _cl3d.seemax
  bkg_spptot[0] = _cl3d.spptot
  bkg_sppmax[0] = _cl3d.sppmax
  bkg_szz[0] = _cl3d.szz
  bkg_srrtot[0] = _cl3d.srrtot
  bkg_srrmax[0] = _cl3d.srrmax
  bkg_srrmean[0] = _cl3d.srrmean
  bkg_emaxe[0] = _cl3d.emaxe
  bkg_bdteg[0] = _cl3d.bdteg
  bkg_quality[0] = _cl3d.quality
  tree_bkg.Fill()


#########################################################################################
# Configure output
print "Configuring output ntuple..."
#output ROOT file
fout_id = os.environ['CMSSW_BASE'] + '/src/L1Trigger/analysis/output/trees/%s/%s/%s_%s_%g.root'%(clusteringAlgo,input_type,input_type,clusteringAlgo,fNumber)
fout = ROOT.TFile( fout_id, "RECREATE" )

#Initialise ttree and define variables
tree_sig = ROOT.TTree("egid_signal","egid_signal")
tree_bkg = ROOT.TTree("egid_background","egid_background")

sig_pt = array( 'f', [-1.] )
sig_eta = array( 'f', [-1.] )
sig_phi = array( 'f', [-1.] )
sig_clusters_n = array( 'f', [-1.] )
sig_showerlength = array( 'f', [-1.] )
sig_coreshowerlength = array( 'f', [-1.] )
sig_firstlayer = array( 'f', [-1.] )
sig_maxlayer = array( 'f', [-1.] )
sig_seetot = array( 'f', [-1.] )
sig_seemax = array( 'f', [-1.] )
sig_spptot = array( 'f', [-1.] )
sig_sppmax = array( 'f', [-1.] )
sig_szz = array( 'f', [-1.] )
sig_srrtot = array( 'f', [-1.] )
sig_srrmax = array( 'f', [-1.] )
sig_srrmean = array( 'f', [-1.] )
sig_emaxe = array( 'f', [-1.] )
sig_bdteg = array( 'f', [-1.] )
sig_quality = array( 'f', [-1.] )
bkg_pt = array( 'f', [-1.] )
bkg_eta = array( 'f', [-1.] )
bkg_phi = array( 'f', [-1.] )
bkg_clusters_n = array( 'f', [-1.] )
bkg_showerlength = array( 'f', [-1.] )
bkg_coreshowerlength = array( 'f', [-1.] )
bkg_firstlayer = array( 'f', [-1.] )
bkg_maxlayer = array( 'f', [-1.] )
bkg_seetot = array( 'f', [-1.] )
bkg_seemax = array( 'f', [-1.] )
bkg_spptot = array( 'f', [-1.] )
bkg_sppmax = array( 'f', [-1.] )
bkg_szz = array( 'f', [-1.] )
bkg_srrtot = array( 'f', [-1.] )
bkg_srrmax = array( 'f', [-1.] )
bkg_srrmean = array( 'f', [-1.] )
bkg_emaxe = array( 'f', [-1.] )
bkg_bdteg = array( 'f', [-1.] )
bkg_quality = array( 'f', [-1.] )

tree_sig.Branch("sig_pt", sig_pt, 'sig_pt/F')
tree_sig.Branch("sig_eta", sig_eta, 'sig_eta/F')
tree_sig.Branch("sig_phi", sig_phi, 'sig_phi/F')
tree_sig.Branch("sig_clusters_n", sig_clusters_n, 'sig_clusters_n/F')
tree_sig.Branch("sig_showerlength", sig_showerlength, 'sig_showerlength/F')
tree_sig.Branch("sig_coreshowerlength", sig_coreshowerlength, 'sig_coreshowerlength/F')
tree_sig.Branch("sig_firstlayer", sig_firstlayer, 'sig_firstlayer/F')
tree_sig.Branch("sig_maxlayer", sig_maxlayer, 'sig_maxlayer/F')
tree_sig.Branch("sig_seetot", sig_seetot, 'sig_seetot/F')
tree_sig.Branch("sig_seemax", sig_seemax, 'sig_seemax/F')
tree_sig.Branch("sig_spptot", sig_spptot, 'sig_spptot/F')
tree_sig.Branch("sig_sppmax", sig_sppmax, 'sig_sppmax/F')
tree_sig.Branch("sig_szz", sig_szz, 'sig_szz/F')
tree_sig.Branch("sig_srrtot", sig_srrtot, 'sig_srrtot/F')
tree_sig.Branch("sig_srrmax", sig_srrmax, 'sig_srrmax/F')
tree_sig.Branch("sig_srrmean", sig_srrmean, 'sig_srrmean/F')
tree_sig.Branch("sig_emaxe", sig_emaxe, 'sig_emaxe/F')
tree_sig.Branch("sig_bdteg", sig_bdteg, 'sig_bdteg/F')
tree_sig.Branch("sig_quality", sig_quality, 'sig_quality/F')
tree_bkg.Branch("bkg_pt", bkg_pt, 'bkg_pt/F')
tree_bkg.Branch("bkg_eta", bkg_eta, 'bkg_eta/F')
tree_bkg.Branch("bkg_phi", bkg_phi, 'bkg_phi/F')
tree_bkg.Branch("bkg_clusters_n", bkg_clusters_n, 'bkg_clusters_n/F')
tree_bkg.Branch("bkg_showerlength", bkg_showerlength, 'bkg_showerlength/F')
tree_bkg.Branch("bkg_coreshowerlength", bkg_coreshowerlength, 'bkg_coreshowerlength/F')
tree_bkg.Branch("bkg_firstlayer", bkg_firstlayer, 'bkg_firstlayer/F')
tree_bkg.Branch("bkg_maxlayer", bkg_maxlayer, 'bkg_maxlayer/F')
tree_bkg.Branch("bkg_seetot", bkg_seetot, 'bkg_seetot/F')
tree_bkg.Branch("bkg_seemax", bkg_seemax, 'bkg_seemax/F')
tree_bkg.Branch("bkg_spptot", bkg_spptot, 'bkg_spptot/F')
tree_bkg.Branch("bkg_sppmax", bkg_sppmax, 'bkg_sppmax/F')
tree_bkg.Branch("bkg_szz", bkg_szz, 'bkg_szz/F')
tree_bkg.Branch("bkg_srrtot", bkg_srrtot, 'bkg_srrtot/F')
tree_bkg.Branch("bkg_srrmax", bkg_srrmax, 'bkg_srrmax/F')
tree_bkg.Branch("bkg_srrmean", bkg_srrmean, 'bkg_srrmean/F')
tree_bkg.Branch("bkg_emaxe", bkg_emaxe, 'bkg_emaxe/F')
tree_bkg.Branch("bkg_bdteg", bkg_bdteg, 'bkg_bdteg/F')
tree_bkg.Branch("bkg_quality", bkg_quality, 'bkg_quality/F')


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

  ##Loop over 3D clusters in event
  #for cl3d_idx in range(0, N_cl3d(cl3d_event) ):
    
  #  #cluster selection
  #  if selection( cl3d_event, cl3d_idx ):

      #instantiate 3d cluster class
  #    cl3d = Cluster3D( cl3d_event, cl3d_idx ) 
  #    if( isGenMatched(gen_event,cl3d,min_dR=0.01) ): fillSignal( cl3d )
  #    else: fillBackground( cl3d )

  #Loop over 3D clusters in event
  for cl3d_idx in range(0, N_cl3d(cl3d_tree) ):
    
    #cluster selection
    if selection( cl3d_tree, cl3d_idx ):

      #instantiate 3d cluster class
      cl3d = Cluster3D( cl3d_tree, cl3d_idx ) 
      if( isGenMatched(gen_tree,cl3d,min_dR=0.1) ): fillSignal( cl3d )
      else: fillBackground( cl3d )



#end of events loop


#########################################################################################

fout.Write()
fout.Close()

#raw_input("Press any key to continue...")

