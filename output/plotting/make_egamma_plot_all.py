import ROOT
import sys
import os
from array import array
from optparse import OptionParser

#Supress output
ROOT.gROOT.SetBatch(ROOT.kTRUE)

#Get options from option parser
def get_options():
  parser = OptionParser()
  parser = OptionParser( usage="usage: python make_egamma_plot.py <options>" )
  parser.add_option("-c", "--clusteringAlgs", dest="clusteringAlgs", default="default,Histomax_vardrth10", help="Variable to plot")
  parser.add_option("-i", "--inputSignalType", dest="inputSignalType", default="SingleElectronPt5_100", help="Input signal type")
  parser.add_option("-d", "--defaultOnly", dest="defaultOnly", default=0, help="Default clustering only plot")
  parser.add_option("-s", "--signalOnly", dest="signalOnly", default=0, help="Only plot signal")
  parser.add_option("-n", "--normalized", dest="normalized", default=1, help="Normalise plot")
  return parser.parse_args()

(opt,args) = get_options()

#Variable plotting options: [bins, mimnimum, maximum]
variables_plotting_options = {'pt':[90,5,50,0], 'eta':[50,-3.14,3.14,0], 'phi':[50,-3.14,3.14,0], 'clusters_n':[30,0,30,0], 'showerlength':[50,0,50,0], 'coreshowerlength':[30,0,30,0], 'firstlayer':[50,0,50,0], 'maxlayer':[50,0,50,0], 'seetot':[100,0,0.10,1], 'seemax':[100,0,0.1,1], 'spptot':[100,0,0.1,1], 'sppmax':[100,0,0.1,1], 'szz':[50,0,50,1], 'srrtot':[100,0,0.01,1], 'srrmax':[100,0,0.01,1], 'srrmean':[100,0,0.01,1], 'emaxe':[60,0,1.2,0], 'bdteg':[50,-1,0.5,0], 'quality':[6,-1,5,0]}

clusteringAlgs = opt.clusteringAlgs.split(",")

inputSignalType = opt.inputSignalType

markerStyles = {'Histomax':34,'Histomax_vardrth0':32,'Histomax_vardrth10':41,'Histomax_vardrth20':29}

#general plotting options
defaultOnly = int(opt.defaultOnly)
signalOnly = int(opt.signalOnly)
normalized = int(opt.normalized)

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: making cl3d variable plots"
print "\n  Input:"
print "    * Clustering Algorithms: %s"%clusteringAlgs
print "    * Default Only: %g"%defaultOnly
print "    * Signal Only: %g"%signalOnly
print "    * Normalized: %g"%normalized
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#Loop over variables
for var in variables_plotting_options:

  print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  print "Preparing plot: Variable = %s"%var

  #Extract plotting options
  binning = variables_plotting_options[var][:-1]
  setLogY=int(variables_plotting_options[var][-1])

  #Configure output
  ROOT.gStyle.SetOptStat(0)
  canv = ROOT.TCanvas("c","c")
  if( setLogY ): canv.SetLogy()

  # extract from input files
  f_sig = {}
  t_sig = {}
  h_sig = {}
  f_bkg = {}
  t_bkg = {}
  h_bkg = {}

  for clusteringAlgo in clusteringAlgs:
    f_sig['%s'%clusteringAlgo] = ROOT.TFile.Open( os.environ['CMSSW_BASE'] + "/src/L1Trigger/analysis/output/trees/%s/%sEta1p6_2p8/%sEta1p6_2p8_%s.root"%(clusteringAlgo,inputSignalType,inputSignalType,clusteringAlgo) )
    f_bkg['%s'%clusteringAlgo] = ROOT.TFile.Open( os.environ['CMSSW_BASE'] + "/src/L1Trigger/analysis/output/trees/%s/SinglePionPt25Eta1p6_2p8/SinglePionPt25Eta1p6_2p8_%s.root"%(clusteringAlgo,clusteringAlgo) )
    t_sig['%s'%clusteringAlgo] = f_sig['%s'%clusteringAlgo].Get("egid_signal")
    t_bkg['%s'%clusteringAlgo] = f_bkg['%s'%clusteringAlgo].Get("egid_background")

    h_sig['%s'%clusteringAlgo] = ROOT.TH1F("h_sig_%s_%s"%(clusteringAlgo,var), "", binning[0], binning[1], binning[2] )
    h_bkg['%s'%clusteringAlgo] = ROOT.TH1F("h_bkg_%s_%s"%(clusteringAlgo,var), "", binning[0], binning[1], binning[2] )

    for ev in t_sig['%s'%clusteringAlgo]: h_sig['%s'%clusteringAlgo].Fill( getattr( ev, "sig_%s"%var ) )
    for ev in t_bkg['%s'%clusteringAlgo]: h_bkg['%s'%clusteringAlgo].Fill( getattr( ev, "bkg_%s"%var ) )

    #normalise histograms
    if normalized:
      h_sig['%s'%clusteringAlgo].Scale(1./h_sig['%s'%clusteringAlgo].GetEntries())
      h_bkg['%s'%clusteringAlgo].Scale(1./h_bkg['%s'%clusteringAlgo].GetEntries())

  #Draw reference histograms
  h_sig['default'].SetLineColor(8)
  h_sig['default'].SetLineWidth(3)
  if normalized: h_sig['default'].GetYaxis().SetTitle("1/N dN/d(%s)"%var)
  else: h_sig['default'].GetYaxis().SetTitle("N")
  h_sig['default'].GetYaxis().SetTitleSize(0.05)
  h_sig['default'].GetYaxis().SetTitleOffset(0.8)
  h_sig['default'].GetXaxis().SetTitle("%s"%var)
  h_sig['default'].GetXaxis().SetTitleSize(0.05)
  h_sig['default'].GetXaxis().SetTitleOffset(0.9)
  h_bkg['default'].SetLineColor(9)
  h_bkg['default'].SetLineWidth(3)
  h_sig['default'].Draw("hist")
  if not signalOnly: h_bkg['default'].Draw("Same hist")

  #If not only plotting reference histograms, loop over other clustering algorithms
  if not defaultOnly:
    for clusteringAlgo in clusteringAlgs:
      if clusteringAlgo == 'default': continue
      h_sig['%s'%clusteringAlgo].SetMarkerStyle( markerStyles[clusteringAlgo] )
      h_sig['%s'%clusteringAlgo].SetMarkerSize( 1.4 )
      h_sig['%s'%clusteringAlgo].SetMarkerColor(8)
      h_sig['%s'%clusteringAlgo].SetLineColor(0)
      h_sig['%s'%clusteringAlgo].Draw("SAME P")
      h_bkg['%s'%clusteringAlgo].SetMarkerStyle( markerStyles[clusteringAlgo] )
      h_bkg['%s'%clusteringAlgo].SetMarkerSize( 1.4 )
      h_bkg['%s'%clusteringAlgo].SetMarkerColor(9)
      h_bkg['%s'%clusteringAlgo].SetLineColor(0)
      if not signalOnly: h_bkg['%s'%clusteringAlgo].Draw("SAME P")

  #Loop over all histograms and get maximum so including all in plot frame
  maximum_value = 0
  for hist in h_sig.itervalues():
    if hist.GetMaximum() > maximum_value: maximum_value = hist.GetMaximum()
  if not signalOnly:
    for hist in h_bkg.itervalues():
      if hist.GetMaximum() > maximum_value: maximum_value = hist.GetMaximum()
  h_sig['default'].SetMaximum( 1.1*maximum_value )
  if( setLogY ): h_sig['default'].SetMinimum( 1e-3 )

  #For legend
  if var == "bdteg": leg1 = ROOT.TLegend(0.160172,0.6,0.47,0.88)
  elif var == "eta": leg1 = ROOT.TLegend(0.360172,0.6,0.67,0.88)
  else: leg1 = ROOT.TLegend(0.560172,0.6,0.87,0.88)
  #dummy histograms and graphs for legend
  h_default = ROOT.TH1F("h_default","",10,0,10)
  h_default.SetLineWidth(3)
  h_default.SetLineColor(1)
  graphs = {}
  for clusteringAlgo in clusteringAlgs:
    if clusteringAlgo == "default": continue
    gr = ROOT.TGraph()
    graphs[ "%s"%clusteringAlgo ] = gr
  leg1.SetFillColor(0)
  leg1.SetLineColor(0)
  leg1.AddEntry("h_sig_default_%s"%var,"Signal: gen-matched e/#gamma","F")
  if not signalOnly: leg1.AddEntry("h_bkg_default_%s"%var,"Background","F")
  #Add dummy histograms and graphs to legend
  leg1.AddEntry(h_default,"default","L")
  if not defaultOnly: 
    for clusteringAlgo in clusteringAlgs:
      if clusteringAlgo == 'default': continue
      graphs[ clusteringAlgo ].SetMarkerColor(1)
      graphs[ clusteringAlgo ].SetMarkerSize(1.2)
      graphs[ clusteringAlgo ].SetMarkerStyle( markerStyles[ clusteringAlgo ] )
      leg1.AddEntry( graphs[clusteringAlgo],"%s"%clusteringAlgo,"P")
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
  if opt.clusteringAlgs == 'default,Histomax': outputString = "/eos/user/j/jlangfor/www/CMS/HGCal/L1/egID/march19/%s/cl3d_variables_HistomaxOnly/cl3d_%s"%(inputSignalType,var)
  else: outputString = "/eos/user/j/jlangfor/www/CMS/HGCal/L1/egID/march19/%s/cl3d_variables/cl3d_%s"%(inputSignalType,var)
  if( defaultOnly ): outputString += "_defaultOnly"
  if( signalOnly ): outputString += "_signalOnly"

  if not normalized: outputString += '_unnormalized'

  canv.Print( '%s.png'%outputString )
  canv.Print( '%s.pdf'%outputString )
  canv.Close()
  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#End of loop over variables
#raw_input("Press any key to continue...")


  
