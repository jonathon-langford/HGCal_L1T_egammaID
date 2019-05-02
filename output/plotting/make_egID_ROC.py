import ROOT
import sys
import os
from array import array
from optparse import OptionParser

#configure options
from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--inputMap', dest='input_map', default='electron,neutrino,default,103X,self,electron_vs_neutrino,baseline,1', help="Colon separated maps of ROC curves to plot, of the type: signal,background,model algo.,CMSSW release,training (self/tpg),discriminator,config,colour" )
  parser.add_option("--legend", dest="legend", default='', help="Legend entries")
  parser.add_option("--zoom", dest="zoom", default="0,0,0", help="Zooming in on plot")
  parser.add_option('--output', dest='output', default='', help="Output directory" )
  parser.add_option('--batch', dest='batch', default=0, help="Set to 1 to supress output to screen" )
  return parser.parse_args()

(opt,args) = get_options()
batch = int( opt.batch )

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "HGCal L1T Analysis: plotting tools for e/g ID"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#Configure inputs: extract map and store as list of dictionaries
input_list = []
for _input in opt.input_map.split(":"):
  inputInfo = _input.split(",")
  if len( inputInfo ) != 8:
    print "  --> [ERROR] Invalid input map"
    sys.exit(1)
  input_list.append({})
  input_list[-1]['signal'] = inputInfo[0]
  input_list[-1]['background'] = inputInfo[1]
  input_list[-1]['model_algo'] = inputInfo[2]
  input_list[-1]['release'] = inputInfo[3]
  input_list[-1]['training'] = inputInfo[4]
  input_list[-1]['discriminator'] = inputInfo[5]
  input_list[-1]['config'] = inputInfo[6]
  input_list[-1]['colour'] = inputInfo[7]

# Define dictionary for type mapping
typeMap = {"electron":"SingleElectron_FlatPt-2to100","photon":"SinglePhoton_FlatPt-8to150","pion":"SinglePion_FlatPt-2to100","neutrino":"SingleNeutrino"}
treeMap = {"electron":"e_sig","photon":"g_sig","pion":"pi_bkg","neutrino":"pu_bkg"}
  
#Configure output
output = opt.output
ROOT.gStyle.SetOptStat(0)
if batch: ROOT.gROOT.SetBatch(ROOT.kTRUE)
canv = ROOT.TCanvas("c","c")
zoom = opt.zoom.split(",")

#Define dictionary to store ROC curve
g_roc = {}

#Loop over inputs and create ROC curves
for i in input_list:

  #Open signal and background files and extract trees
  f_sig = ROOT.TFile.Open( os.environ['CMSSW_BASE'] + "/src/L1Trigger/analysis/output/egID_trees/%s/%s/%s_%s.root"%(i['release'],i['model_algo'],typeMap[i['signal']],i['model_algo']))
  f_bkg = ROOT.TFile.Open( os.environ['CMSSW_BASE'] + "/src/L1Trigger/analysis/output/egID_trees/%s/%s/%s_%s.root"%(i['release'],i['model_algo'],typeMap[i['background']],i['model_algo']))
  t_sig = f_sig.Get( treeMap[i['signal']] )
  t_bkg = f_bkg.Get( treeMap[i['background']] )

  #Define histograms
  h_sig = ROOT.TH1F("h_sig_bdt_%s_%s_%s_%s_%s"%(i['model_algo'],i['release'],i['training'],i['discriminator'],i['config']), "", 240, -1.1, 1.1 )
  h_bkg = ROOT.TH1F("h_bkg_bdt_%s_%s_%s_%s_%s"%(i['model_algo'],i['release'],i['training'],i['discriminator'],i['config']), "", 240, -1.1, 1.1 )

  #Fill histograms from trees
  for ev in t_sig:
    if i['training'] == "tpg": h_sig.Fill( getattr( ev, "cl3d_bdt_tpg" ) )
    elif i['training'] == "self": h_sig.Fill( getattr( ev, "cl3d_bdt_%s_%s_%s_%s"%(i['model_algo'],i['release'],i['discriminator'],i['config']) ) )
  for ev in t_bkg:
    if i['training'] == "tpg": h_bkg.Fill( getattr( ev, "cl3d_bdt_tpg" ) )
    elif i['training'] == "self": h_bkg.Fill( getattr( ev, "cl3d_bdt_%s_%s_%s_%s"%(i['model_algo'],i['release'],i['discriminator'],i['config']) ) )

  #Define graph to plot ROC
  rocStr = '%s_%s_%s_%s_%s_%s_%s'%(i['signal'],i['background'],i['model_algo'],i['release'],i['training'],i['discriminator'],i['config'])
  g_roc[ rocStr ] = ROOT.TGraph()

  # Loop over bins in histogram and plot points in graph
  for j in range( h_sig.GetNbinsX() ):
    eff_sig = h_sig.Integral(j+1,h_sig.GetNbinsX()+1)/h_sig.Integral()
    eff_bkg = h_bkg.Integral(j+1,h_bkg.GetNbinsX()+1)/h_bkg.Integral()
    g_roc[ rocStr ].SetPoint(j,1-eff_bkg,eff_sig)

  #Configure ROC
  if int(zoom[0]) == 1:
    g_roc[ rocStr ].GetXaxis().SetRangeUser(float(zoom[1]),1.)
    g_roc[ rocStr ].GetYaxis().SetRangeUser(float(zoom[2]),1.05)
  else:
    g_roc[ rocStr ].GetXaxis().SetRangeUser(0.75,1.)
    g_roc[ rocStr ].GetYaxis().SetRangeUser(0.6,1.05)
  g_roc[ rocStr ].GetXaxis().SetTitle('1 - #epsilon_{b}')
  g_roc[ rocStr ].GetXaxis().SetTitleSize(0.05)
  g_roc[ rocStr ].GetXaxis().SetTitleOffset(0.9)
  g_roc[ rocStr ].GetYaxis().SetTitle('#epsilon_{s}')
  g_roc[ rocStr ].GetYaxis().SetTitleSize(0.05)
  g_roc[ rocStr ].GetYaxis().SetTitleOffset(0.9)
  g_roc[ rocStr ].SetLineWidth(2)
  g_roc[ rocStr ].SetLineColor( int(i['colour']) )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loop over input list again and plot
for _idx in range( len( input_list ) ):
  i = input_list[_idx]
  rocStr = '%s_%s_%s_%s_%s_%s_%s'%(i['signal'],i['background'],i['model_algo'],i['release'],i['training'],i['discriminator'],i['config'])
  if _idx == 0: g_roc[ rocStr ].Draw("AL")
  else: g_roc[ rocStr ].Draw("L SAME")

#For legend
if len(opt.legend) > 1:
  entry_list = []
  #Entry of type: text,colour,option
  for _entry in opt.legend.split("+"):
    entryInfo = _entry.split(":")
    #print "  --> [DEBUG] entryInfo =", entryInfo, ", len(entryInfo) =", len(entryInfo)
    if len(entryInfo) != 3:
      print "  --> [ERROR] Invalid legend entry. Exiting..."
      sys.exit(1)
    entry_list.append({})
    entry_list[-1]['text'] = entryInfo[0]
    entry_list[-1]['colour'] = entryInfo[1]
    entry_list[-1]['option'] = entryInfo[2]

  graph_list = []
  #Create dummy graphs to place in legend
  for entry in entry_list:
    gr = ROOT.TGraph()
    if entry['colour'] != "empty":
      gr.SetFillColor( int(entry['colour']) )
      gr.SetLineColor( int(entry['colour']) )
      gr.SetLineWidth( 2 )
    graph_list.append( gr )

  #Create legend and add entries
  leg = ROOT.TLegend(0.15,0.15,0.6,0.45)
  leg.SetFillColor(0)
  leg.SetLineColor(0)
  for _idx in range( len( entry_list ) ):
    entry = entry_list[_idx]
    if entry['colour'] == "empty": leg.AddEntry( 0, "%s"%entry['text'], "%s"%entry['option'] )
    else: leg.AddEntry( graph_list[_idx], "%s"%entry['text'], "%s"%entry['option'] )
  leg.Draw("Same")

lat = ROOT.TLatex()
lat.SetTextFont(42)
lat.SetLineWidth(2)
lat.SetTextAlign(11)
lat.SetNDC()
lat.SetTextSize(0.05)
lat.DrawLatex(0.1,0.92,"#bf{CMS Phase-2} #scale[0.75]{#it{Internal}}")
lat.DrawLatex(0.8,0.92,"14 TeV")
canv.Update()

if output != '':
  canv.Print( '%s.png'%output )
  canv.Print( '%s.pdf'%output )

if not batch: raw_input("Press any key to continue...")
