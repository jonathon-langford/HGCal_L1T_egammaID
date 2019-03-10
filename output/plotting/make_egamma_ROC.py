import ROOT
import sys
import os
from array import array
from optparse import OptionParser

#Get options from option parser
def get_options():
  parser = OptionParser()
  parser = OptionParser( usage="usage: python make_egamma_ROC.py <variable>" )
  parser.add_option("-v", "--variable", dest="variable", default="", help="Variable to plot discrimination")
  parser.add_option("-c", "--clusteringAlgs", dest="clusteringAlgs", default="default,Histomax", help="Variable to plot")
  parser.add_option("-i", "--inputSignalType", dest="inputSignalType", default="SingleElectronPt5_100", help="Input signal type")
  parser.add_option("-d", "--defaultOnly", dest="defaultOnly", default=0, help="Default clustering only plot")
  return parser.parse_args()

(opt,args) = get_options()

variables_plotting_options = {'pt':[90,5,50], 'eta':[50,-3.14,3.14], 'phi':[50,-3.14,3.14], 'clusters_n':[30,0,30], 'showerlength':[50,0,50], 'coreshowerlength':[30,0,30], 'firstlayer':[50,0,50], 'maxlayer':[50,0,50], 'seetot':[100,0,0.10], 'seemax':[100,0,0.1], 'spptot':[100,0,0.1], 'sppmax':[100,0,0.1], 'szz':[50,0,50], 'srrtot':[100,0,0.01], 'srrmax':[100,0,0.01], 'srrmean':[100,0,0.01], 'emaxe':[60,0,1.2], 'bdteg':[100,-1,1], 'quality':[4,-1,3]}

clusteringAlgs = opt.clusteringAlgs.split(",")

inputSignalType = opt.inputSignalType

markerStyles = {'default':24,'Histomax':34,'Histomax_vardrth0':32,'Histomax_vardrth10':41,'Histomax_vardrth20':29}
lineColors = {'default':1, 'Histomax':2, 'Histomax_vardrth0':5, 'Histomax_vardrth10':8, 'Histomax_vardrth20':9}

var = opt.variable
if var not in variables_plotting_options: 
  print "[ERROR] Variable (%s) not supported"%var
  sys.exit(1)

#Extact plotting options
binning = variables_plotting_options[var]
defaultOnly = int(opt.defaultOnly)

#Configure output
ROOT.gStyle.SetOptStat(0)
canv = ROOT.TCanvas("c","c")

# extract from input files
f_sig = {}
t_sig = {}
h_sig = {}
f_bkg = {}
t_bkg = {}
h_bkg = {}
g_roc = {}

for clusteringAlgo in clusteringAlgs:
  f_sig['%s'%clusteringAlgo] = ROOT.TFile.Open( os.environ['CMSSW_BASE'] + "/src/L1Trigger/analysis/output/trees/%s/%sEta1p6_2p8/%sEta1p6_2p8_%s.root"%(clusteringAlgo,inputSignalType,inputSignalType,clusteringAlgo) )
  f_bkg['%s'%clusteringAlgo] = ROOT.TFile.Open( os.environ['CMSSW_BASE'] + "/src/L1Trigger/analysis/output/trees/%s/SinglePionPt25Eta1p6_2p8/SinglePionPt25Eta1p6_2p8_%s.root"%(clusteringAlgo,clusteringAlgo) )
  t_sig['%s'%clusteringAlgo] = f_sig['%s'%clusteringAlgo].Get("egid_signal")
  t_bkg['%s'%clusteringAlgo] = f_bkg['%s'%clusteringAlgo].Get("egid_background")

  h_sig['%s'%clusteringAlgo] = ROOT.TH1F("h_sig_%s_%s"%(clusteringAlgo,var), "", binning[0], binning[1], binning[2] )
  h_bkg['%s'%clusteringAlgo] = ROOT.TH1F("h_bkg_%s_%s"%(clusteringAlgo,var), "", binning[0], binning[1], binning[2] )

  for ev in t_sig['%s'%clusteringAlgo]: h_sig['%s'%clusteringAlgo].Fill( getattr( ev, "sig_%s"%var ) )
  for ev in t_bkg['%s'%clusteringAlgo]: h_bkg['%s'%clusteringAlgo].Fill( getattr( ev, "bkg_%s"%var ) )

  #Define graph to plot ROC
  g_roc['%s'%clusteringAlgo] = ROOT.TGraph()

  #Loop over bins in histogram and plot points in graph according to signal and bkg efficiency
  for i in range( h_sig['%s'%clusteringAlgo].GetNbinsX() ):
    eff_sig = h_sig['%s'%clusteringAlgo].Integral(i+1,h_sig['%s'%clusteringAlgo].GetNbinsX()+1)/h_sig['%s'%clusteringAlgo].Integral()
    eff_bkg = h_bkg['%s'%clusteringAlgo].Integral(i+1,h_bkg['%s'%clusteringAlgo].GetNbinsX()+1)/h_bkg['%s'%clusteringAlgo].Integral()
    g_roc['%s'%clusteringAlgo].SetPoint(i, 1-eff_bkg, eff_sig)

  #Configure ROC graph
  if( var == "bdteg" ): 
    g_roc['%s'%clusteringAlgo].GetXaxis().SetRangeUser(0.5,1.1)
    g_roc['%s'%clusteringAlgo].GetYaxis().SetRangeUser(0.5,1.1)
  else: 
    g_roc['%s'%clusteringAlgo].GetXaxis().SetRangeUser(0,1.1)
    g_roc['%s'%clusteringAlgo].GetYaxis().SetRangeUser(0,1.1)
  g_roc['%s'%clusteringAlgo].GetXaxis().SetTitle('1 - #epsilon_{b}')
  g_roc['%s'%clusteringAlgo].GetXaxis().SetTitleSize(0.05)
  g_roc['%s'%clusteringAlgo].GetXaxis().SetTitleOffset(0.9)
  g_roc['%s'%clusteringAlgo].GetYaxis().SetTitle('#epsilon_{s}')
  g_roc['%s'%clusteringAlgo].GetYaxis().SetTitleSize(0.05)
  g_roc['%s'%clusteringAlgo].GetYaxis().SetTitleOffset(0.9)
  g_roc['%s'%clusteringAlgo].SetLineWidth(2)
  g_roc['%s'%clusteringAlgo].SetMarkerStyle( markerStyles['%s'%clusteringAlgo] )
  g_roc['%s'%clusteringAlgo].SetLineColor( lineColors['%s'%clusteringAlgo] )
  g_roc['%s'%clusteringAlgo].SetMarkerColor( lineColors['%s'%clusteringAlgo] )

g_roc['default'].Draw('ACP')
if not defaultOnly:
  for clusteringAlgo in clusteringAlgs:
    if clusteringAlgo == 'default': continue
    g_roc['%s'%clusteringAlgo].Draw('CP SAME')

leg1 = ROOT.TLegend(0.15,0.15,0.6,0.6)
leg1.SetFillColor(0)
leg1.SetLineColor(0)
leg1.AddEntry(g_roc['default'],"default: AUC = %5.4f"%g_roc['default'].Integral(),"LP")
if not defaultOnly: 
  for clusteringAlgo in clusteringAlgs:
    if clusteringAlgo == 'default': continue
    leg1.AddEntry(g_roc[ clusteringAlgo ], "%s: AUC = %5.4f"%(clusteringAlgo,g_roc[clusteringAlgo].Integral()) ,"LP")
leg1.Draw("Same")

lat = ROOT.TLatex()
lat.SetTextFont(42)
lat.SetLineWidth(2)
lat.SetTextAlign(11)
lat.SetNDC()
lat.SetTextSize(0.05)
lat.DrawLatex(0.1,0.92,"#bf{CMS Phase-2} #scale[0.75]{#it{Preliminary}}")
lat.DrawLatex(0.8,0.92,"14 TeV")
if( var == "bdteg" ): lat.DrawLatex(0.4,0.82,"e/#gamma BDT (optimised w.r.t. default)")
else: lat.DrawLatex(0.7,0.82, var)

canv.Update()

outputString = "/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/march19/%s/ROC/ROC_%s"%(inputSignalType,var)
canv.Print( "%s.png"%outputString )
canv.Print( "%s.pdf"%outputString )

raw_input("Press any key to continue...")


  
