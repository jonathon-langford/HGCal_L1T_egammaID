import ROOT
import sys
import os
from array import array
from optparse import OptionParser

#Get options from option parser
def get_options():
  parser = OptionParser()
  parser = OptionParser( usage="usage: python make_cl3dVar_plot.py <options>" )
  parser.add_option("--inputMap", dest="input_map", default="", help="List of inputs to plot, of the form: input type,clustering algo.,geometry,tpg software version,colour,marker style,plotting option:..." )
  parser.add_option("--variable", dest="variable", default="", help="Variable to plot")
  parser.add_option("--outputDir", dest="output_dir", default='', help="Output directory")
  parser.add_option("--normalized", dest="normalized", default=1, help="Normalise plot")
  parser.add_option("--legend", dest="legend", default='', help="Legend entries")
  parser.add_option("--batch", dest="batch", default=0, help="Suppress output of plots to screen")
  return parser.parse_args()

(opt,args) = get_options()

#Variable plotting options: [bins, mimnimum, maximum]
var = opt.variable
variables_plotting_options = {'pt':[150,0,150,0], 'eta':[50,-3.14,3.14,0], 'phi':[50,-3.14,3.14,0], 'clusters_n':[60,0,60,0], 'showerlength':[60,0,60,0], 'coreshowerlength':[30,0,30,0], 'firstlayer':[20,0,20,0], 'maxlayer':[50,0,50,0], 'seetot':[100,0,0.20,0], 'seemax':[50,0,0.1,0], 'spptot':[100,0,0.2,1], 'sppmax':[100,0,0.2,1], 'szz':[100,0,100,1], 'srrtot':[150,0,0.03,1], 'srrmax':[150,0,0.03,1], 'srrmean':[50,0,0.01,1], 'emaxe':[60,0,1.2,0], 'bdteg':[50,-1,1.,1], 'quality':[6,-1,5,0]}
if( var not in variables_plotting_options )&( var != 'all' ): 
  print "  --> [ERROR] Variable (%s) not supported"%var
  sys.exit(1)

# Output allowed variables
allowed_variables = "{"
for variable in variables_plotting_options: allowed_variables += "%s,"%variable
allowed_variables = allowed_variables[:-1]+"}"
#print allowed_variables

# Dictionary to store histograms
fileDict = {}
treeDict = {}
histDict = {}

# Extract input map and store in dictionary
input_list = []
for _input in opt.input_map.split(":"):
  inputInfo = _input.split(",")
  if len(inputInfo) != 7:
    print "  --> [ERROR] Invalid input. Exiting"
    sys.exit(1)
  input_list.append({})
  input_list[-1]['type'] = inputInfo[0]
  input_list[-1]['cl3d_algo'] = inputInfo[1]
  input_list[-1]['geometry'] = inputInfo[2]
  input_list[-1]['version'] = inputInfo[3]
  input_list[-1]['colour'] = inputInfo[4]
  input_list[-1]['marker'] = inputInfo[5]
  input_list[-1]['option'] = inputInfo[6]

# Define dictionary for type mapping
typeMap = {"electron":"SingleElectron_FlatPt-2to100","photon":"SinglePhoton_FlatPt-8to150","pion":"SinglePion_FlatPt-2to100","neutrino":"SingleNeutrino"}
treeMap = {"electron":"e_sig","photon":"g_sig","pion":"pi_bkg","neutrino":"pu_bkg"}

#Dictionary for input directory mapping
inputDirMap = {"v3_7_0":"/afs/cern.ch/work/j/jlangfor/HGCal/L1T_new/CMSSW_10_5_0/src/L1Trigger/egid_analysis/HGCal_L1T_egammaID/output/trees","v3_9_4":os.environ['CMSSW_BASE'] + "/src/L1Trigger/HGCal_L1T_egammaID/output/trees"}

#Store maximum value of all histograms for keeping all on same plot
maximum_value = 0

# Output options
output_dir = opt.output_dir
batch = int(opt.batch)
if batch: ROOT.gROOT.SetBatch(ROOT.kTRUE)

#Extract plotting options
binning = variables_plotting_options[var][:-1]
setLogY= variables_plotting_options[var][-1]
normalized = int(opt.normalized)

#Configure output
ROOT.gStyle.SetOptStat(0)
canv = ROOT.TCanvas("c","c")
if( setLogY ): canv.SetLogy()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loop over inputs in map and create histogram from var in tree
for i in input_list:

  hist_key = '%s_%s_%s_%s'%(i['type'],i['cl3d_algo'],i['geometry'],i['version'])
  fileDict[hist_key] = ROOT.TFile.Open( inputDirMap[i['version']]+"/%s/%s/%s/%s_%s_full.root"%(i['geometry'],i['cl3d_algo'],typeMap[i['type']],typeMap[i['type']],i['cl3d_algo']) )
  treeDict[hist_key] = fileDict[hist_key].Get( treeMap[i['type']] )

  histDict[hist_key] = ROOT.TH1F("h_%s_%s_%s"%(i['type'],i['cl3d_algo'],i['geometry']), "", binning[0], binning[1], binning[2] )

  #for ev in treeDict[hist_key]: 
  #  if( abs(getattr( ev, "cl3d_eta" ))>=2.7 )&( abs(getattr( ev, "cl3d_eta" ))<3.0 ): histDict[hist_key].Fill( getattr( ev, "cl3d_%s"%var ) )

  for ev in treeDict[hist_key]: histDict[hist_key].Fill( getattr( ev, "cl3d_%s"%var ) )

  #normalise histograms
  if normalized: histDict[hist_key].Scale(1./histDict[hist_key].GetEntries())

  #plotting options
  histDict[hist_key].SetLineWidth(2)
  histDict[hist_key].SetLineColor(int(i['colour']))
  histDict[hist_key].SetMarkerColor(int(i['colour']))
  histDict[hist_key].SetMarkerStyle(int(i['marker']))
  histDict[hist_key].SetMarkerSize(1.2)
  if normalized: histDict[hist_key].GetYaxis().SetTitle("1/N dN/d(%s)"%var)
  else: histDict[hist_key].GetYaxis().SetTitle("N")
  histDict[hist_key].GetYaxis().SetTitleSize(0.05)
  histDict[hist_key].GetYaxis().SetTitleOffset(0.8)
  histDict[hist_key].GetXaxis().SetTitle("%s"%var)
  histDict[hist_key].GetXaxis().SetTitleSize(0.05)
  histDict[hist_key].GetXaxis().SetTitleOffset(0.9)

  #if maximum > current: save new maximum
  if histDict[hist_key].GetMaximum() > maximum_value: maximum_value = histDict[hist_key].GetMaximum()
  
# Loop over input list again and plot
for _idx in range( len( input_list ) ):
  i = input_list[_idx]
  hist_key = '%s_%s_%s_%s'%(i['type'],i['cl3d_algo'],i['geometry'],i['version'])
  if( _idx == 0 ):
    if( setLogY ): 
      histDict[hist_key].SetMaximum( 1.3*maximum_value ) 
      histDict[hist_key].SetMinimum( 1e-4 )
    else:
      histDict[hist_key].SetMaximum( 1.1*maximum_value ) 
    histDict[hist_key].Draw("%s"%i['option'])
  else: 
    histDict[hist_key].Draw("SAME %s"%i['option'])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Text for canvas and legend
lat = ROOT.TLatex()
lat.SetTextFont(42)
lat.SetLineWidth(2)
lat.SetTextAlign(11)
lat.SetNDC()
lat.SetTextSize(0.05)
lat.DrawLatex(0.1,0.92,"#bf{CMS Phase-2} #scale[0.75]{#it{Internal}}")
lat.DrawLatex(0.8,0.92,"14 TeV")
#lat.DrawLatex(0.65,0.55,"2.7 < |#eta| < 3.0")

#For legend
entry_list = []
#Entry of type: text,colour,option
if opt.legend != "":
  for _entry in opt.legend.split("+"):
    entryInfo = _entry.split(",")
    if len(entryInfo) != 4:
      print "  --> [ERROR] Invalid legend entry. Exiting..."
      sys.exit(1)
    entry_list.append({})
    entry_list[-1]['text'] = entryInfo[0]
    entry_list[-1]['colour'] = entryInfo[1]
    entry_list[-1]['marker'] = entryInfo[2]
    entry_list[-1]['option'] = entryInfo[3]

  graph_list = []
  #Create dummy graphs to place in legend
  for entry in entry_list:
    gr = ROOT.TGraph()
    gr.SetFillColor( int(entry['colour']) )
    gr.SetLineColor( int(entry['colour']) )
    gr.SetLineWidth( 2 )
    gr.SetMarkerColor( int(entry['colour']) )
    gr.SetMarkerStyle( int(entry['marker']) )
    gr.SetMarkerSize( 1.2 )
    graph_list.append( gr )

  #Create legend and add entries
  if var == "eta": leg = ROOT.TLegend(0.38,0.65,0.62,0.88)
  elif var == "srrmax": leg = ROOT.TLegend(0.65,0.15,0.88,0.38)
  else: leg = ROOT.TLegend(0.65,0.65,0.88,0.88)
  #leg = ROOT.TLegend(0.38,0.38,0.62,0.62)
  leg.SetFillColor(0)
  leg.SetLineColor(0)
  for _idx in range( len( entry_list ) ):
    entry = entry_list[_idx]
    leg.AddEntry( graph_list[_idx], "%s"%entry['text'], "%s"%entry['option'] )
  leg.Draw("Same")

canv.Update()

if output_dir != '':
  output_file = "%s/cl3d_%s"%(output_dir,var)
  if not normalized: output_file += '_unnormalized'
  canv.Print( '%s.png'%output_file )
  canv.Print( '%s.pdf'%output_file )

if not batch: raw_input("Press Enter to continue...")

#Close files properly
for f in fileDict.itervalues(): f.Close()
