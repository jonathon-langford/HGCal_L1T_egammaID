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
pile_up = "0"

# Different eg-id to compare
bdt_configs = ['tpg','full']
#bdt_configs = ['tpg']

# Define input files
f_input = {}
cl3dTrees = {}
genTrees = {}
for b in bdt_configs:
  if b == "full": f_input[b] = ROOT.TFile("/eos/home-j/jlangfor/hgcal/l1/egid/v3_9_4/v9/SingleElectron_FlatPt-2to100_%sPU_retrained/SingleElectron_FlatPt-2to100_%sPU_%s.root"%(pile_up,pile_up,b))
  if b == "tpg": f_input[b] = ROOT.TFile("/eos/home-j/jlangfor/hgcal/l1/egid/v3_9_4/v9/SingleElectron_FlatPt-2to100_%sPU/SingleElectron_FlatPt-2to100_%sPU_%s.root"%(pile_up,pile_up,b))
  cl3dTrees[b] = f_input[b].Get("%s/HGCalTriggerNtuple"%inputDirectoryMap['cl3d'])
  genTrees[b] = f_input[b].Get("%s/HGCalTriggerNtuple"%inputDirectoryMap['gen'])

#PID for gen matching electron: depends on input type
gen_pdgid = 11

wp = {"full_loweta":0.5276631,"full_higheta":0.8825340,"tpg_loweta":0.03496629,"tpg_higheta":0.13347613}

#########################################################################################
# Configure output
print "Configuring output ntuple..."
#output ROOT file
f_out = ROOT.TFile( "efficiency_%sPU_v3_9_4.root"%pile_up, "RECREATE" )

# define histograms: total, cl3d_matched and passing egid for each config
h = {}

for b in bdt_configs:
  h['%s_total_genpt'%b] = ROOT.TH1F("h_%s_total_genpt"%b,"",20,0,100)
  h['%s_total_geneta'%b] = ROOT.TH1F("h_%s_total_geneta"%b,"",64,-3.2,3.2)
  h['%s_cl3d_matched_genpt'%b] = ROOT.TH1F("h_%s_cl3d_matched_genpt"%b,"",20,0,100)
  h['%s_cl3d_matched_geneta'%b] = ROOT.TH1F("h_%s_cl3d_matched_geneta"%b,"",64,-3.2,3.2)
  h['%s_egid_genpt'%b] = ROOT.TH1F("h_%s_egid_genpt"%b,"",20,0,100)
  h['%s_egid_geneta'%b] = ROOT.TH1F("h_%s_egid_geneta"%b,"",64,-3.2,3.2)
  h['%s_egid_quality_genpt'%b] = ROOT.TH1F("h_%s_egid_quality_genpt"%b,"",20,0,100)
  h['%s_egid_quality_geneta'%b] = ROOT.TH1F("h_%s_egid_quality_geneta"%b,"",64,-3.2,3.2)
  h['%s_cl3d_bdteg_loweta'%b] = ROOT.TH1F("h_%s_cl3d_bdteg_loweta"%b,"", 100,-1,1)
  h['%s_cl3d_bdteg_higheta'%b] = ROOT.TH1F("h_%s_cl3d_bdteg_higheta"%b,"",100,-1,1)


#########################################################################################
# Loop over bdt configs
for b in bdt_configs:
  print "Processing for BDT config: %s"%b
  # Loop over events
  for ev_idx in range(genTrees[b].GetEntries()):

    if ev_idx % 10000 == 0: print "  --> Processing Event: %g/%g"%(ev_idx,genTrees[b].GetEntries())

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
          
            h['%s_total_genpt'%b].Fill( genTrees[b].gen_pt[gen_idx] )
            h['%s_total_geneta'%b].Fill( genTrees[b].gen_eta[gen_idx] )
            #Define 4 vector for gen particle
            gen_p4 = ROOT.TLorentzVector()
            gen_p4.SetPtEtaPhiE( genTrees[b].gen_pt[gen_idx], genTrees[b].gen_eta[gen_idx], genTrees[b].gen_phi[gen_idx], genTrees[b].gen_energy[gen_idx] )

            #Loop over clusters in event: pT > 20 GeV (choose highest pt) and dR < 0.2
            cl3d_genmatched_maxpt_idx = -1
            cl3d_genmatched_maxpt = -999.
            for cl3d_idx in range(N_cl3d):
              #Apply selection
              if cl3dTrees[b].cl3d_pt[cl3d_idx] < 20.: continue
              #Define 4 vector for cl3d
              cl3d_p4 = ROOT.TLorentzVector()
              cl3d_p4.SetPtEtaPhiE( cl3dTrees[b].cl3d_pt[cl3d_idx], cl3dTrees[b].cl3d_eta[cl3d_idx], cl3dTrees[b].cl3d_phi[cl3d_idx], cl3dTrees[b].cl3d_energy[cl3d_idx] )
              #Require cluster to be in dR < 0.2
              if gen_p4.DeltaR( cl3d_p4 ) < 0.2:
                #If cluster pT > present max then set
                if cl3d_p4.Pt() > cl3d_genmatched_maxpt:
                  cl3d_genmatched_maxpt = cl3d_p4.Pt()
                  cl3d_genmatched_maxpt_idx = cl3d_idx

            #If genmatched idx set then fill hist
            if cl3d_genmatched_maxpt_idx >= 0:
              h['%s_cl3d_matched_genpt'%b].Fill( genTrees[b].gen_pt[gen_idx] )
              h['%s_cl3d_matched_geneta'%b].Fill( genTrees[b].gen_eta[gen_idx] )

              if abs(cl3dTrees[b].cl3d_eta[cl3d_genmatched_maxpt_idx]) < 2.7: h['%s_cl3d_bdteg_loweta'%b].Fill( cl3dTrees[b].cl3d_bdteg[cl3d_genmatched_maxpt_idx] )
              else: h['%s_cl3d_bdteg_higheta'%b].Fill( cl3dTrees[b].cl3d_bdteg[cl3d_genmatched_maxpt_idx] )

              # Using quality
              if cl3dTrees[b].cl3d_quality[cl3d_genmatched_maxpt_idx] > 0:
                h['%s_egid_quality_genpt'%b].Fill( genTrees[b].gen_pt[gen_idx] )
                h['%s_egid_quality_geneta'%b].Fill( genTrees[b].gen_eta[gen_idx] )
 
              #Using bdt score
              if abs(cl3dTrees[b].cl3d_eta[cl3d_genmatched_maxpt_idx]) < 2.7:
                if cl3dTrees[b].cl3d_bdteg[cl3d_genmatched_maxpt_idx] > wp["%s_loweta"%b]: 
                  h['%s_egid_genpt'%b].Fill( genTrees[b].gen_pt[gen_idx] )
                  h['%s_egid_geneta'%b].Fill( genTrees[b].gen_eta[gen_idx] )
              else:
                if cl3dTrees[b].cl3d_bdteg[cl3d_genmatched_maxpt_idx] > wp["%s_higheta"%b]:
                  h['%s_egid_genpt'%b].Fill( genTrees[b].gen_pt[gen_idx] )
                  h['%s_egid_geneta'%b].Fill( genTrees[b].gen_eta[gen_idx] )
            
   

#end of events loop
#########################################################################################

f_out.Write()
f_out.Close()

#raw_input("Press any key to continue...")

