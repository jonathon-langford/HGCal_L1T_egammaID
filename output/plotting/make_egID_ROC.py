import ROOT
import sys
import os
from array import array
from optparse import OptionParser

#configure options
from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--algoMap', dest='algoMap', default='default:default,default:Histomax_vardrth10', help="Comma separated list of [clustering:model] algorithms to plot" )
  parser.add_option('--inputSignalType', dest='inputSignalType', default='electron', help="Input signal type" )
  parser.add_option('--inputBackgroundType', dest='inputBackgroundType', default='neutrino', help="Input background type" )
  parser.add_option('--output', dest='output_dir', default='', help="Output directory" )
  parser.add_option('--batch', dest='batch', default=0, help="Set to 1 to supress output to screen" )
  return parser.parse_args()

(opt,args) = get_options()
batch = int( opt.batch )

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: plotting tools for e/g ID"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#Configure inputs
algoPairs = opt.algoMap.split(",")
output_dir = opt.output_dir

if opt.inputSignalType == "electron_hgcal_only": signalName = "SingleElectronPt5_100Eta1p6_2p8"
else: 
  print "  --> [ERROR] Invalid signal type. No eg ID tree in folder. Leaving..."
  sys.exit(1)
if opt.inputBackgroundType == "neutrino": backgroundName = "SingleNeutrino"
else:
  print "  --> [ERROR] Invalid background type. No eg ID tree in folder. Leaving..."
  sys.exit(1)

input_map = {}
for algoPair in algoPairs:
  clusteringAlgo = algoPair.split(":")[0]
  input_map[ "signal_%s"%clusteringAlgo ] = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/egID_trees/%s/%s_%s_test.root"%(clusteringAlgo,signalName,clusteringAlgo)
  input_map[ "background_%s"%clusteringAlgo ] = "/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/egID_trees/%s/%s_%s_test.root"%(clusteringAlgo,backgroundName,clusteringAlgo)

tree_map = {"signal":"egid_signal","background":"egid_background"}

#markerStyles = {'default':4, 'Histomax_vardrth10':34, 'tpg_Histomax_vardrth10':22}
lineColors = {'default:default':1, 'default:Histomax_vardrth10':ROOT.kOrange, 'default:tpg_Histomax_vardrth10':ROOT.kOrange+7, 'Histomax_vardrth10:Histomax_vardrth10':8,'Histomax_vardrth10:tpg_Histomax_vardrth10':9,'Histomax_vardrth10:default':2}
#lineColors = {'default':2, 'Histomax_vardrth10':8, 'tpg_Histomax_vardrth10':9}

#Configure output
ROOT.gStyle.SetOptStat(0)
if batch: ROOT.gROOT.SetBatch(ROOT.kTRUE)
canv = ROOT.TCanvas("c","c")

#Define dictionaries to store histograms for bdt score output and ROC curves
h_sig = {}
h_bkg = {}
g_roc = {}

#Loop over algo pairs and create ROC curves
for algoPair in algoPairs:

  #Extract algo pair from list
  clusteringAlgo = algoPair.split(":")[0]
  modelAlgo = algoPair.split(":")[1]

  #Open signal and background files and extract trees
  f_sig = ROOT.TFile.Open( input_map['signal_%s'%clusteringAlgo] )
  t_sig = f_sig.Get( tree_map['signal'] )
  f_bkg = ROOT.TFile.Open( input_map['background_%s'%clusteringAlgo] )
  t_bkg = f_bkg.Get( tree_map['background'] )

  #define histograms
  h_sig[ "%s_%s"%(clusteringAlgo,modelAlgo) ] = ROOT.TH1F("h_sig_%s_bdtscore_%s"%(clusteringAlgo,modelAlgo), "", 240, -1.1, 1.1 )
  h_bkg[ "%s_%s"%(clusteringAlgo,modelAlgo) ] = ROOT.TH1F("h_bkg_%s_bdtscore_%s"%(clusteringAlgo,modelAlgo), "", 240, -1.1, 1.1 )

  for ev in t_sig: h_sig[ "%s_%s"%(clusteringAlgo,modelAlgo) ].Fill( getattr( ev, "bdtscore_%s"%modelAlgo ) )
  for ev in t_bkg: h_bkg[ "%s_%s"%(clusteringAlgo,modelAlgo) ].Fill( getattr( ev, "bdtscore_%s"%modelAlgo ) )

  # Define graph to plot ROC
  g_roc[ "%s_%s"%(clusteringAlgo,modelAlgo) ] = ROOT.TGraph()

  # Loop over bins in histogram and plot points in graph
  for i in range( h_sig[ "%s_%s"%(clusteringAlgo,modelAlgo) ].GetNbinsX() ):
    eff_sig = h_sig[ "%s_%s"%(clusteringAlgo,modelAlgo) ].Integral(i+1,h_sig[ "%s_%s"%(clusteringAlgo,modelAlgo) ].GetNbinsX()+1)/h_sig[ "%s_%s"%(clusteringAlgo,modelAlgo) ].Integral()
    eff_bkg = h_bkg[ "%s_%s"%(clusteringAlgo,modelAlgo) ].Integral(i+1,h_bkg[ "%s_%s"%(clusteringAlgo,modelAlgo) ].GetNbinsX()+1)/h_bkg[ "%s_%s"%(clusteringAlgo,modelAlgo) ].Integral()
    g_roc[ "%s_%s"%(clusteringAlgo,modelAlgo) ].SetPoint(i,1-eff_bkg,eff_sig)

  #Configure ROC
  g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].GetXaxis().SetRangeUser(0.75,1.)
  g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].GetYaxis().SetRangeUser(0.6,1.1)
  g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].GetXaxis().SetTitle('1 - #epsilon_{b}')
  g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].GetXaxis().SetTitleSize(0.05)
  g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].GetXaxis().SetTitleOffset(0.9)
  g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].GetYaxis().SetTitle('#epsilon_{s}')
  g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].GetYaxis().SetTitleSize(0.05)
  g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].GetYaxis().SetTitleOffset(0.9)
  g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].SetLineWidth(2)
  #g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].SetMarkerStyle( markerStyles[modelAlgo] )
  g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].SetLineColor( lineColors[algoPair] )
  g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].SetMarkerColor( lineColors[algoPair] )  

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOTTING
for algoPair_idx in range( len(algoPairs) ):
  clusteringAlgo = algoPairs[algoPair_idx].split(":")[0]
  modelAlgo = algoPairs[algoPair_idx].split(":")[1]
  if algoPair_idx == 0: g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].Draw("AL")
  else: g_roc["%s_%s"%(clusteringAlgo,modelAlgo)].Draw("L SAME")

leg1 = ROOT.TLegend(0.15,0.15,0.6,0.45)
leg1.SetFillColor(0)
leg1.SetLineColor(0)
leg1.AddEntry(0, "(clustering algo., training algo.)", "")
for algoPair in algoPairs:
  clusteringAlgo = algoPair.split(":")[0]
  modelAlgo = algoPair.split(":")[1]
  leg1.AddEntry(g_roc[ "%s_%s"%(clusteringAlgo,modelAlgo) ], "(%s,%s)"%(clusteringAlgo,modelAlgo), "L")
leg1.Draw("Same")

lat = ROOT.TLatex()
lat.SetTextFont(42)
lat.SetLineWidth(2)
lat.SetTextAlign(11)
lat.SetNDC()
lat.SetTextSize(0.05)
lat.DrawLatex(0.1,0.92,"#bf{CMS Phase-2} #scale[0.75]{#it{Internal}}")
lat.DrawLatex(0.8,0.92,"14 TeV")
canv.Update()

if output_dir != '':
  output_file = "%s/ROC_cl3d-train"%output_dir
  for algoPair in algoPairs:
    clusteringAlgo = algoPair.split(":")[0]
    modelAlgo = algoPair.split(":")[1]
    output_file += "_%s-%s"%(clusteringAlgo,modelAlgo)  
  canv.Print( '%s.png'%output_file )
  canv.Print( '%s.pdf'%output_file )

if not batch: raw_input("Press any key to continue...")
