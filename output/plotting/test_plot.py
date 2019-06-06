import ROOT

version = "370"

f = ROOT.TFile("output_%s.root"%version)
h2d = (f.Get("h_2d_fail_eta_vs_ptresponse").RebinX(2)).RebinY(2)

canv = ROOT.TCanvas("canv","canv")
ROOT.gStyle.SetOptStat(0)

h2d.GetXaxis().SetRangeUser(1.5,3.0)
h2d.GetXaxis().SetTitle("|#eta|")
h2d.GetXaxis().SetTitleSize(0.04)
h2d.GetXaxis().SetLabelSize(0.03)
h2d.GetXaxis().SetTitleOffset(0.95)
h2d.GetYaxis().SetRangeUser(0.8,1.3)
h2d.GetYaxis().SetTitle("p_{T}^{cluster}/p_{T}^{GEN}")
h2d.GetYaxis().SetTitleSize(0.04)
h2d.GetYaxis().SetLabelSize(0.03)
h2d.GetYaxis().SetTitleOffset(1.0)

h2d.SetTitle("Processed with v%s.%s.%s"%(version[0],version[1],version[2]))

h2d.Draw("COLZ")

