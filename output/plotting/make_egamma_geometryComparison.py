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
  parser.add_option("-s", "--signalOnly", dest="signalOnly", default=0, help="Only plot signal")
  return parser.parse_args()

(opt,args) = get_options()

#Variable plotting options: [bins, mimnimum, maximum]
variables_plotting_options = {'pt':[90,5,50,0], 'eta':[50,-3.14,3.14,0], 'phi':[50,-3.14,3.14,0], 'clusters_n':[30,0,30,0], 'showerlength':[50,0,50,0], 'coreshowerlength':[30,0,30,0], 'firstlayer':[50,0,50,0], 'maxlayer':[50,0,50,0], 'seetot':[100,0,0.10,1], 'seemax':[100,0,0.1,1], 'spptot':[100,0,0.1,1], 'sppmax':[100,0,0.1,1], 'szz':[50,0,50,1], 'srrtot':[100,0,0.01,1], 'srrmax':[100,0,0.01,1], 'srrmean':[100,0,0.01,1], 'emaxe':[60,0,1.2,0], 'bdteg':[50,-1,0.5,0], 'quality':[6,-1,5,0]}

geometries = ['93X','103X']

inputStrings = {'93X':"/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/trees/Histomax_vardrth10/SingleElectronPt5_100Eta1p6_2p8/SingleElectronPt5_100Eta1p6_2p8_Histomax_vardrth10.root",'103X':"/afs/cern.ch/work/j/jlangfor/HGCal/L1/CMSSW_10_4_0/src/L1Trigger/analysis/output/trees/104X/Histomax_vardrth10/SingleElectronPt2_100Eta1p6_2p8/SingleElectronPt2_100Eta1p6_2p8_Histomax_vardrth10.root"}

geometryColors = {'93X':1,'103X':2}
geometryMarkerStyle = {'93X':41,'103X':34}

#general plotting options
signalOnly = int(opt.signalOnly)

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: making cl3d variable comparison for different geometries"
print "\n  Input:"
print "    * Clustering Algorithm: Histomax_vardrth10"
print "    * Signal Only: %g"%signalOnly
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
  #f_bkg = {}
  #t_bkg = {}
  #h_bkg = {}

  for geometry in geometries:
    f_sig['%s'%geometry] = ROOT.TFile.Open( "%s"%inputStrings[ geometry ] )
    #f_bkg['%s'%clusteringAlgo] = ROOT.TFile.Open( os.environ['CMSSW_BASE'] + "/src/L1Trigger/analysis/output/trees/%s/SinglePionPt25Eta1p6_2p8/SinglePionPt25Eta1p6_2p8_%s.root"%(clusteringAlgo,clusteringAlgo) )
    t_sig['%s'%geometry] = f_sig['%s'%geometry].Get("egid_signal")
    #t_bkg['%s'%clusteringAlgo] = f_bkg['%s'%clusteringAlgo].Get("egid_background")

    h_sig['%s'%geometry] = ROOT.TH1F("h_sig_%s_%s"%(geometry,var), "", binning[0], binning[1], binning[2] )
    #h_bkg['%s'%clusteringAlgo] = ROOT.TH1F("h_bkg_%s_%s"%(clusteringAlgo,var), "", binning[0], binning[1], binning[2] )

    for ev in t_sig['%s'%geometry]: h_sig['%s'%geometry].Fill( getattr( ev, "sig_%s"%var ) )
    #for ev in t_bkg['%s'%geometry]: h_bkg['%s'%geometry].Fill( getattr( ev, "bkg_%s"%var ) )

    #normalise histograms
    h_sig['%s'%geometry].Scale(1./h_sig['%s'%geometry].GetEntries())
    #h_bkg['%s'%geometry].Scale(1./h_bkg['%s'%geometry].GetEntries())

  #Draw reference histograms
  h_sig['93X'].SetMarkerStyle( geometryMarkerStyle['93X'] )
  h_sig['93X'].SetMarkerColor( geometryColors['93X'] )
  h_sig['93X'].SetMarkerSize( 1.5 )
  h_sig['93X'].SetLineColor(0)
  h_sig['93X'].GetYaxis().SetTitle("1/N dN/d(%s)"%var)
  h_sig['93X'].GetYaxis().SetTitleSize(0.05)
  h_sig['93X'].GetYaxis().SetTitleOffset(0.8)
  h_sig['93X'].GetXaxis().SetTitle("%s"%var)
  h_sig['93X'].GetXaxis().SetTitleSize(0.05)
  h_sig['93X'].GetXaxis().SetTitleOffset(0.9)
  #h_bkg['93X'].SetLineColor(9)
  #h_bkg['93X'].SetLineWidth(3)
  h_sig['93X'].Draw("P")
  #if not signalOnly: h_bkg['93X'].Draw("Same hist")

  #If not only plotting reference histograms, loop over other clustering algorithms
  for geometry in geometries:
    if geometry == '93X': continue
    h_sig['%s'%geometry].SetMarkerColor( geometryColors[geometry] )
    h_sig['%s'%geometry].SetMarkerStyle( geometryMarkerStyle[geometry] )
    h_sig['%s'%geometry].SetMarkerSize( 1.5 )
    h_sig['%s'%geometry].SetLineColor( 0 )
    h_sig['%s'%geometry].Draw("same P")
    #h_bkg['%s'%geometry].SetMarkerStyle( markerStyles[geometry] )
    #h_bkg['%s'%geometry].SetMarkerSize( 1.4 )
    #h_bkg['%s'%geometry].SetMarkerColor(9)
    #h_bkg['%s'%geometry].SetLineColor(0)
    #if not signalOnly: h_bkg['%s'%geometry].Draw("SAME P")

  #Loop over all histograms and get maximum so including all in plot frame
  maximum_value = 0
  for hist in h_sig.itervalues():
    if hist.GetMaximum() > maximum_value: maximum_value = hist.GetMaximum()
  #if not signalOnly:
  #  for hist in h_bkg.itervalues():
  #    if hist.GetMaximum() > maximum_value: maximum_value = hist.GetMaximum()
  h_sig['93X'].SetMaximum( 1.2*maximum_value )
  if( setLogY ): h_sig['93X'].SetMinimum( 1e-3 )

  #For legend
  if var == "bdteg": leg1 = ROOT.TLegend(0.160172,0.6,0.47,0.88)
  elif var == "eta": leg1 = ROOT.TLegend(0.360172,0.6,0.67,0.88)
  else: leg1 = ROOT.TLegend(0.560172,0.6,0.87,0.88)
  leg1.SetFillColor(0)
  leg1.SetLineColor(0)
  for geometry in geometries:
    leg1.AddEntry("h_sig_%s_%s"%(geometry,var),"%s: gen-matched e/#gamma"%geometry,"P")
  #if not signalOnly: leg1.AddEntry("h_bkg_default_%s"%var,"Background","F")
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
  outputString = "/eos/user/j/jlangfor/www/CMS/HGCal/L1/egID/march19/93X_vs_103X/cl3d_%s_93X_vs_103X"%(var)
 
  canv.Print( '%s.png'%outputString )
  canv.Print( '%s.pdf'%outputString )
  canv.Close()
  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#End of loop over variables
#raw_input("Press any key to continue...")


  
