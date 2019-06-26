import ROOT

var = "genpt"
pile_up = "0"
eff_min = 0.75

f_in = ROOT.TFile("efficiency_%sPU_v3_9_4.root"%pile_up)

#Define histograms

#Divide histograms
h_tpg_egid_genpt = f_in.Get("h_tpg_egid_genpt")
h_tpg_total_genpt = f_in.Get("h_tpg_total_genpt")
h_tpg_egid_geneta = f_in.Get("h_tpg_egid_geneta")
h_tpg_total_geneta = f_in.Get("h_tpg_total_geneta")

h_full_egid_genpt = f_in.Get("h_full_egid_genpt")
h_full_total_genpt = f_in.Get("h_full_total_genpt")
h_full_egid_geneta = f_in.Get("h_full_egid_geneta")
h_full_total_geneta = f_in.Get("h_full_total_geneta")

h_tpg_eff_genpt = h_tpg_egid_genpt.Clone()
h_tpg_eff_geneta = h_tpg_egid_geneta.Clone()
h_full_eff_genpt = h_full_egid_genpt.Clone()
h_full_eff_geneta = h_full_egid_geneta.Clone()

h_tpg_eff_genpt.Sumw2()
h_tpg_eff_geneta.Sumw2() 
h_full_eff_genpt.Sumw2() 
h_full_eff_geneta.Sumw2() 

h_tpg_eff_genpt.Divide(h_tpg_total_genpt)
h_tpg_eff_geneta.Divide(h_tpg_total_geneta)
h_full_eff_genpt.Divide(h_full_total_genpt)
h_full_eff_geneta.Divide(h_full_total_geneta)

#Make canvases
canv = ROOT.TCanvas("c","c")
ROOT.gStyle.SetOptStat(0)
if var == "genpt":
  #configure plot
  h_tpg_eff_genpt.GetYaxis().SetRangeUser(eff_min,1.1)
  h_tpg_eff_genpt.GetYaxis().SetTitle("(L1>thr. & matched to GEN)/GEN")
  h_tpg_eff_genpt.GetYaxis().SetLabelSize(0.03)
  h_tpg_eff_genpt.GetYaxis().SetTitleSize(0.04)
  h_tpg_eff_genpt.GetYaxis().SetTitleOffset(1.0)
  h_tpg_eff_genpt.GetXaxis().SetRangeUser(25,100)
  h_tpg_eff_genpt.GetXaxis().SetTitle("p_{T}^{GEN}  [GeV]")
  h_tpg_eff_genpt.GetXaxis().SetLabelSize(0.03)
  h_tpg_eff_genpt.GetXaxis().SetTitleSize(0.04)
  h_tpg_eff_genpt.GetXaxis().SetTitleOffset(0.95)

  h_tpg_eff_genpt.SetMarkerColor(1)
  h_tpg_eff_genpt.SetMarkerSize(1.3)
  h_tpg_eff_genpt.SetMarkerStyle(34)
  h_tpg_eff_genpt.SetLineColor(1)
  h_tpg_eff_genpt.SetLineWidth(1)

  h_full_eff_genpt.SetMarkerColor(2)
  h_full_eff_genpt.SetMarkerSize(1.3)
  h_full_eff_genpt.SetMarkerStyle(34)
  h_full_eff_genpt.SetLineColor(2)
  h_full_eff_genpt.SetLineWidth(2)

  h_tpg_eff_genpt.Draw("P Hist")
  h_full_eff_genpt.Draw("SAME P Hist")

  #Draw a line at 100% efficiency
  line1 = ROOT.TLine(25,1,100,1)
  line1.SetLineWidth(2)
  line1.SetLineStyle(2)
  line1.Draw("Same")

  #Latex
  lat = ROOT.TLatex()
  lat.SetTextFont(42)
  lat.SetLineWidth(2)
  lat.SetTextAlign(11)
  lat.SetNDC()
  lat.SetTextSize(0.04)
  lat.DrawLatex(0.12,0.85,"PU%s, Histomax_vardr, p_{T}^{L1} > 20GeV, p_{T}^{GEN} > 30GeV"%pile_up)

  #Legend
  leg1 = ROOT.TLegend(0.6,0.2,0.89,0.35)
  leg1.SetFillColor(0)
  leg1.SetLineColor(0)
  leg1.AddEntry(h_tpg_eff_genpt,"tpg (v8 3.5.2)","P")
  leg1.AddEntry(h_full_eff_genpt,"full (v9 3.9.4)","P")
  leg1.Draw("Same")

  canv.SaveAs("/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/june19_2/efficiency/efficiency_vs_genpt_%sPU.png"%pile_up)
  canv.SaveAs("/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/june19_2/efficiency/efficiency_vs_genpt_%sPU.pdf"%pile_up)

elif var == "geneta":
  #configure plot
  h_tpg_eff_geneta.GetYaxis().SetRangeUser(eff_min,1.1)
  h_tpg_eff_geneta.GetYaxis().SetTitle("(L1>thr. & matched to GEN)/GEN")
  h_tpg_eff_geneta.GetYaxis().SetLabelSize(0.03)
  h_tpg_eff_geneta.GetYaxis().SetTitleSize(0.04)
  h_tpg_eff_geneta.GetYaxis().SetTitleOffset(1.0)
  h_tpg_eff_geneta.GetXaxis().SetRangeUser(-3.2,3.2)
  h_tpg_eff_geneta.GetXaxis().SetTitle("#eta^{GEN}")
  h_tpg_eff_geneta.GetXaxis().SetLabelSize(0.03)
  h_tpg_eff_geneta.GetXaxis().SetTitleSize(0.04)
  h_tpg_eff_geneta.GetXaxis().SetTitleOffset(0.95)

  h_tpg_eff_geneta.SetMarkerColor(1)
  h_tpg_eff_geneta.SetMarkerSize(1.3)
  h_tpg_eff_geneta.SetMarkerStyle(34)
  h_tpg_eff_geneta.SetLineColor(1)
  h_tpg_eff_geneta.SetLineWidth(1)

  h_full_eff_geneta.SetMarkerColor(2)
  h_full_eff_geneta.SetMarkerSize(1.3)
  h_full_eff_geneta.SetMarkerStyle(34)
  h_full_eff_geneta.SetLineColor(2)
  h_full_eff_geneta.SetLineWidth(2)

  h_tpg_eff_geneta.Draw("P Hist")
  h_full_eff_geneta.Draw("SAME P Hist")

  #Draw a line at 100% efficiency
  line1 = ROOT.TLine(-3.2,1,3.2,1)
  line1.SetLineWidth(2)
  line1.SetLineStyle(2)
  line1.Draw("Same")

  #Also lines for different eta regions
  line2 = ROOT.TLine(-3,eff_min,-3,1.03)
  line3 = ROOT.TLine(-2.7,eff_min,-2.7,1.03)
  line4 = ROOT.TLine(-1.5,eff_min,-1.5,1.03)
  line5 = ROOT.TLine(1.5,eff_min,1.5,1.03)
  line6 = ROOT.TLine(2.7,eff_min,2.7,1.03)
  line7 = ROOT.TLine(3,eff_min,3,1.03) 
  line2.SetLineWidth(2)
  line2.SetLineStyle(2)
  line2.Draw("Same")
  line3.SetLineWidth(2)
  line3.SetLineStyle(2)
  line3.Draw("Same")
  line4.SetLineWidth(2)
  line4.SetLineStyle(2)
  line4.Draw("Same")
  line5.SetLineWidth(2)
  line5.SetLineStyle(2)
  line5.Draw("Same")
  line6.SetLineWidth(2)
  line6.SetLineStyle(2)
  line6.Draw("Same")
  line7.SetLineWidth(2)
  line7.SetLineStyle(2)
  line7.Draw("Same")

  #Latex
  lat = ROOT.TLatex()
  lat.SetTextFont(42)
  lat.SetLineWidth(2)
  lat.SetTextAlign(11)
  lat.SetNDC()
  lat.SetTextSize(0.04)
  lat.DrawLatex(0.12,0.85,"PU%s, Histomax_vardr, p_{T}^{L1} > 20GeV, p_{T}^{GEN} > 30GeV"%pile_up)

  #Legend
  leg1 = ROOT.TLegend(0.33,0.2,0.62,0.35)
  leg1.SetFillColor(0)
  leg1.SetLineColor(0)
  leg1.AddEntry(h_tpg_eff_geneta,"tpg (v8 3.5.2)","P")
  leg1.AddEntry(h_full_eff_geneta,"full (v9 3.9.4)","P")
  leg1.Draw("Same")

  canv.SaveAs("/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/june19_2/efficiency/efficiency_vs_geneta_%sPU.png"%pile_up)
  canv.SaveAs("/eos/home-j/jlangfor/www/CMS/HGCal/L1/egID/june19_2/efficiency/efficiency_vs_geneta_%sPU.pdf"%pile_up)


