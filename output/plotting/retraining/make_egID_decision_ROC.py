import ROOT
import sys
import os
from array import array
from optparse import OptionParser

#configure options
from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--modelAlgs', dest='modelAlgs', default='default,Histomax_vardrth10,tpg_Histomax_vardrth10', help="Clustering algorithms with which BDTs were trained" )
  parser.add_option('--inputAlgo', dest='inputAlgo', default='Histomax_vardrth10', help="Clustering algorithm with which to check BDT performance" )
  parser.add_option('--signalSample', dest='signalSample', default='SingleElectronPt5_100Eta1p6_2p8', help="Input signal" )
  parser.add_option('--backgroundSample', dest='backgroundSample', default='SinglePionPt25Eta1p6_2p8', help="Input background" )
  return parser.parse_args()

(opt,args) = get_options()

modelAlgs = opt.modelAlgs.split(",")
inputAlgo = opt.inputAlgo
signalSample = opt.signalSample
backgroundSample = opt.backgroundSample

#Define inputs
input_map = {"signal":"/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/egID_trees/%s/%s_%s.root"%(inputAlgo,signalSample,inputAlgo),"background":"/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/egID_trees/%s/%s_%s.root"%(inputAlgo,backgroundSample,inputAlgo)}
tree_map = {"signal":"egid_signal","background":"egid_background"}

markerStyles = {'default':4, 'Histomax_vardrth10':34, 'tpg_Histomax_vardrth10':22}
lineColors = {'default':2, 'Histomax_vardrth10':8, 'tpg_Histomax_vardrth10':9}

#Configure output
ROOT.gStyle.SetOptStat(0)
canv = ROOT.TCanvas("c","c")

# extract from input files
f_sig = ROOT.TFile.Open( input_map['signal'] )
t_sig = f_sig.Get( tree_map['signal'] )
f_bkg = ROOT.TFile.Open( input_map['background'] )
t_bkg = f_bkg.Get( tree_map['background'] )
#Define dictionaries for histograms for bdts trained with different algs
h_sig = {}
h_bkg = {}
g_roc = {}
for modelAlgo in modelAlgs:

  #i
  h_sig[modelAlgo] = ROOT.TH1F("h_sig_bdtscore_%s"%(modelAlgo), "", 240, -1.1, 1.1 )
  h_bkg[modelAlgo] = ROOT.TH1F("h_bkg_bdtscore_%s"%(modelAlgo), "", 240, -1.1, 1.1 )

  for ev in t_sig: h_sig[modelAlgo].Fill( getattr( ev, "bdtscore_%s"%modelAlgo ) )
  for ev in t_bkg: h_bkg[modelAlgo].Fill( getattr( ev, "bdtscore_%s"%modelAlgo ) )

  #Define graph to plot ROC
  g_roc[modelAlgo] = ROOT.TGraph()

  #Loop over bins in histogram and plot points in graph according to signal and bkg efficiency
  for i in range( h_sig[modelAlgo].GetNbinsX() ):
    eff_sig = h_sig[modelAlgo].Integral(i+1,h_sig[modelAlgo].GetNbinsX()+1)/h_sig[modelAlgo].Integral()
    eff_bkg = h_bkg[modelAlgo].Integral(i+1,h_bkg[modelAlgo].GetNbinsX()+1)/h_bkg[modelAlgo].Integral()
    g_roc[modelAlgo].SetPoint(i, 1-eff_bkg, eff_sig)

  #Configure ROC graph
  g_roc[modelAlgo].GetXaxis().SetRangeUser(0.6,1.1)
  g_roc[modelAlgo].GetYaxis().SetRangeUser(0.6,1.1)
  g_roc[modelAlgo].GetXaxis().SetTitle('1 - #epsilon_{b}')
  g_roc[modelAlgo].GetXaxis().SetTitleSize(0.05)
  g_roc[modelAlgo].GetXaxis().SetTitleOffset(0.9)
  g_roc[modelAlgo].GetYaxis().SetTitle('#epsilon_{s}')
  g_roc[modelAlgo].GetYaxis().SetTitleSize(0.05)
  g_roc[modelAlgo].GetYaxis().SetTitleOffset(0.9)
  g_roc[modelAlgo].SetLineWidth(2)
  g_roc[modelAlgo].SetMarkerStyle( markerStyles[modelAlgo] )
  g_roc[modelAlgo].SetLineColor( lineColors[modelAlgo] )
  g_roc[modelAlgo].SetMarkerColor( lineColors[modelAlgo] )

g_roc['default'].Draw("ALP")
for modelAlgo in modelAlgs:
  if modelAlgo == "default": continue
  g_roc[ modelAlgo ].Draw('LP same')

leg1 = ROOT.TLegend(0.15,0.15,0.45,0.45)
leg1.SetFillColor(0)
leg1.SetLineColor(0)
for modelAlgo in modelAlgs:
  leg1.AddEntry(g_roc[ modelAlgo ], "%s"%modelAlgo,"LP")
leg1.Draw("Same")

lat = ROOT.TLatex()
lat.SetTextFont(42)
lat.SetLineWidth(2)
lat.SetTextAlign(11)
lat.SetNDC()
lat.SetTextSize(0.05)
lat.DrawLatex(0.1,0.92,"#bf{CMS Phase-2} #scale[0.75]{#it{Preliminary}}")
lat.DrawLatex(0.8,0.92,"14 TeV")
canv.Update()

#outputString = "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/march19/%s/ROC/ROC_%s"%(inputSignalType,var)
#canv.Print( "%s.png"%outputString )
#canv.Print( "%s.pdf"%outputString )

raw_input("Press any key to continue...")


  
