import ROOT
import sys
import os
import math
from array import array
from optparse import OptionParser

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: efficiency plot"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#Define map to extract TDirectory for gen and histomax
inputDirectoryMap = {'gen':'Floatingpoint8ThresholdRef2dRef3dGenclustersntuple','cl3d':'Floatingpoint8ThresholdDummyHistomaxvardrClustersntuple'}

#pile up
pile_up = "200"

# Different eg-id to compare
bdt_configs = ['full']

# Define input files
f_input = {}
cl3dTrees = {}
genTrees = {}
for b in bdt_configs:
  #if b == "full": f_input[b] = ROOT.TFile("/eos/home-j/jlangfor/hgcal/l1/egid/may19/new_egid/v9/SingleElectron_FlatPt-2to100_%sPU/unknown_tpg/SingleElectron_FlatPt-2to100_%sPU_%s.root"%(pile_up,pile_up,b))
  if b == "full": f_input[b] = ROOT.TFile("/eos/home-j/jlangfor/hgcal/l1/egid/may19/new_egid/v9/SingleElectron_FlatPt-2to100_%sPU/SingleElectron_FlatPt-2to100_%sPU_%s.root"%(pile_up,pile_up,b))
  elif b == "tpg": f_input[b] = ROOT.TFile("/eos/home-j/jlangfor/hgcal/l1/egid/may19/v9/SingleElectron_FlatPt-2to100_%sPU/SingleElectron_FlatPt-2to100_%sPU_%s.root"%(pile_up,pile_up,b))
  cl3dTrees[b] = f_input[b].Get("%s/HGCalTriggerNtuple"%inputDirectoryMap['cl3d'])
  genTrees[b] = f_input[b].Get("%s/HGCalTriggerNtuple"%inputDirectoryMap['gen'])

#PID for gen matching electron: depends on input type
gen_pdgid = 11

wp = {"full_loweta":-0.1937952,"full_higheta":0.7078400,"tpg_loweta":0.03496629,"tpg_higheta":0.13347613}

#########################################################################################
#Output file
f_out = ROOT.TFile("output_370.root","RECREATE")

h_2d_eta_vs_ptresponse = ROOT.TH2F("h_2d_eta_vs_ptresponse","",80,1.5,3.2,80,0.8,1.5)
h_2d_pass_eta_vs_ptresponse = ROOT.TH2F("h_2d_pass_eta_vs_ptresponse","",80,1.5,3.2,80,0.8,1.5)
h_2d_fail_eta_vs_ptresponse = ROOT.TH2F("h_2d_fail_eta_vs_ptresponse","",80,1.5,3.2,80,0.8,1.5)

debug_ = False

N_events = 20000
#N_events = genTrees[b].GetEntries()

#########################################################################################
# Loop over bdt configs
for b in bdt_configs:
  print "Processing for BDT config: %s"%b
  # Loop over events
  for ev_idx in range(N_events):

    if ev_idx % 1000 == 0: print "  --> Processing Event: %g/%g"%(ev_idx,genTrees[b].GetEntries())
    if(debug_): print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    if(debug_): print "Processing event: %g"%ev_idx

    #Extract event info from gen and all cl3d trees
    genTrees[b].GetEntry( ev_idx )
    cl3dTrees[b].GetEntry( ev_idx )

    #Extract number of gen particles + cl3d in event
    N_gen = genTrees[b].gen_n
    N_cl3d = cl3dTrees[b].cl3d_n

    #Loop over gen particles in event
    for gen_idx in range( N_gen ):
      #Apply selection and add to total histogram
      # if electron...
      if abs( genTrees[b].gen_pdgid[gen_idx] ) == gen_pdgid:
        if genTrees[b].gen_pt[gen_idx] > 30.:
          if( abs( genTrees[b].gen_eta[gen_idx] ) >= 1.55 )&( abs( genTrees[b].gen_eta[gen_idx] ) < 2.9 ):
          
            #Define 4 vector for gen particle
            gen_p4 = ROOT.TLorentzVector()
            gen_p4.SetPtEtaPhiE( genTrees[b].gen_pt[gen_idx], genTrees[b].gen_eta[gen_idx], genTrees[b].gen_phi[gen_idx], genTrees[b].gen_energy[gen_idx] )

            if(debug_): print "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            if(debug_): print "  GEN ELECTRON: pT = %5.4f GeV, eta = %5.4f, phi = %5.4f"%(gen_p4.Pt(),gen_p4.Eta(),gen_p4.Phi())
            if(debug_): print ""

            #Loop over clusters in event: pT > 20 GeV (choose highest pt) and dR < 0.2
            cl3d_genmatched_maxpt_idx = -1
            cl3d_genmatched_maxpt = -999.
            if(debug_): print "  CL3D: total number = %g"%N_cl3d
            for cl3d_idx in range(N_cl3d):
              #Apply selection
              if cl3dTrees[b].cl3d_pt[cl3d_idx] < 20.: continue
              #Define 4 vector for cl3d
              cl3d_p4 = ROOT.TLorentzVector()
              cl3d_p4.SetPtEtaPhiE( cl3dTrees[b].cl3d_pt[cl3d_idx], cl3dTrees[b].cl3d_eta[cl3d_idx], cl3dTrees[b].cl3d_phi[cl3d_idx], cl3dTrees[b].cl3d_energy[cl3d_idx] )
              #Require cluster to be in dR < 0.2
              if gen_p4.DeltaR( cl3d_p4 ) < 0.2:
                if(debug_): print "       --> cl3d %g: pT = %5.4f GeV, eta = %5.4f, phi = %5.4f, bdt = %5.4f, pass = %g"%(cl3d_idx,cl3d_p4.Pt(),cl3d_p4.Eta(),cl3d_p4.Phi(),cl3dTrees[b].cl3d_bdteg[cl3d_idx],cl3dTrees[b].cl3d_quality[cl3d_idx])
                #If cluster pT > present max then set
                if cl3d_p4.Pt() > cl3d_genmatched_maxpt:
                  cl3d_genmatched_maxpt = cl3d_p4.Pt()
                  cl3d_genmatched_maxpt_idx = cl3d_idx
            
            if(debug_): print "       --> Max pT cluster: %g, pass = %g"%(cl3d_genmatched_maxpt_idx,cl3dTrees[b].cl3d_quality[cl3d_genmatched_maxpt_idx])  

            if cl3d_genmatched_maxpt_idx != -1: 
              h_2d_eta_vs_ptresponse.Fill(abs(cl3dTrees[b].cl3d_eta[cl3d_genmatched_maxpt_idx]),cl3dTrees[b].cl3d_pt[cl3d_genmatched_maxpt_idx]/genTrees[b].gen_pt[gen_idx])
              if( cl3dTrees[b].cl3d_quality[cl3d_genmatched_maxpt_idx] > 0 ): h_2d_pass_eta_vs_ptresponse.Fill(abs(cl3dTrees[b].cl3d_eta[cl3d_genmatched_maxpt_idx]),cl3dTrees[b].cl3d_pt[cl3d_genmatched_maxpt_idx]/genTrees[b].gen_pt[gen_idx])
              else: h_2d_fail_eta_vs_ptresponse.Fill(abs(cl3dTrees[b].cl3d_eta[cl3d_genmatched_maxpt_idx]),cl3dTrees[b].cl3d_pt[cl3d_genmatched_maxpt_idx]/genTrees[b].gen_pt[gen_idx])
              

#end of events loop
#########################################################################################

f_out.Write()
f_out.Close()
