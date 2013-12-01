from drawer_init import *

drawerInit = DrawerInit()
wait = True
if not wait:  gROOT.SetBatch(1)

#tfile = TFile.Open("../bin/compactified.root")
tfile = TFile.Open("../bin/compactified.0.root")
#tfile = TFile.Open("../bin/compactified.L1ETM40.0.root")
tree = tfile.Events


# Rates
h_rates_obs = TH1F("h_rates_obs", "; ; count", 30, 0, 30)
h_rates_exp = TH1F("h_rates_exp", "; ; count", 30, 0, 30)

htemp1 = TH1F("htemp1", ";", 10, 0, 10)
htemp2 = TH1F("htemp2", ";", 10, 0, 10)

ii = 1
for i, trig in enumerate(triggers):
    tree.Draw("triggerFlags[%i] >> htemp1" % i, runsel, "goff")
    tree.Draw("triggerFlags[%i] && triggerFlags[%i] >> htemp2" % (i, i_reftrig), runsel, "goff")

    print trig
    print "    %7i   %7i +/- %7i" % (htemp1.GetBinContent(2), htemp2.GetBinContent(2) * ps_reftrig, sqrt(htemp2.GetBinContent(2)) * ps_reftrig)
    print

    if not trig.startswith("HLT_L1ETM") and not trig.startswith("HLT_MET"):
        h_rates_obs.SetBinContent(ii, htemp1.GetBinContent(2))
        #h_rates_obs.SetBinError(ii, 0)
        h_rates_obs.GetXaxis().SetBinLabel(ii, trig)
        h_rates_exp.SetBinContent(ii, htemp2.GetBinContent(2) * ps_reftrig)
        h_rates_exp.SetBinError(ii, sqrt(htemp2.GetBinContent(2)) * ps_reftrig)
        h_rates_exp.GetXaxis().SetBinLabel(ii, trig)
        ii += 1

c1 = TCanvas("c1", "c1", 1100, 700)
#c1.SetLeftMargin(0.1)
c1.SetRightMargin(0.2)
c1.SetBottomMargin(0.4)
h_rates_exp.SetMaximum(h_rates_exp.GetMaximum() * 1.3)
h_rates_exp.SetMarkerSize(0)
h_rates_exp.SetLineColor(kRed)
h_rates_exp.SetLineWidth(2)
h_rates_exp.GetXaxis().SetLabelSize(0.04)
h_rates_exp.GetYaxis().SetLabelSize(0.04)
h_rates_exp.GetYaxis().SetTitleOffset(0.9)
h_rates_exp.GetYaxis().SetTitleSize(0.05)
h_rates_exp.Draw("hist")
h_rates_exp.GetXaxis().SetRangeUser(0,10)
h_rates_exp.Draw("same e1")
h_rates_obs.Draw("same p")

leg = TLegend(0.82, 0.74, 0.96, 0.96, "MET Run 207454")
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.SetShadowColor(0)
leg.AddEntry(h_rates_obs, "observed", "p")
leg.AddEntry(h_rates_exp, "predicted")
leg.Draw()

imgname = "rates_MET_R207454"
gPad.Print(imgdir + imgname + ".png")
gPad.Print(imgdir + imgname + ".pdf")
