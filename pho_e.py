"""you need to import the following modules or class...
TFile, TCanvas, TH1D, gStyle, TF1, gROOT, TLegend, gPad, TGraphErrors, TMath"""

import math
from ROOT import TFile, TCanvas, TH1D, gStyle, TF1, gROOT, TLegend, gPad
from ROOT import TGraphErrors, TMath

from mystyle import style
from usefun import *


class pho_e():
    cent_list1 = [70, 60, 50, 40, 30, 20, 10, 5, 0]
    cent_list2 = [80, 70, 60, 50, 40, 30, 20, 10, 5]
    hist_sige = "hnSigE"
    hist_bemcmatch = "hbemcmatch"
    hist_tofmatch = "htofmatch"
    hist_nphi = "hnphi"
    hist_neta = "hneta"
    hist_phiDist = "hphiDist"
    hist_zDist = "hzDist"
    hist_nHitsFit = "hnHitsFit"
    hist_nHitsDedx = "hnHitsDedx"
    hOrigin_xy = "hOrigin_xy"
    hOrigin_zy = "hOrigin_yz"
    hOrigin_zx = "hOrigin_xz"
    hpairdca = "hpairdca"

    def __init__(self, filename, canvasname):
        self.__filename__ = filename + ".root"
        self.__file__ = TFile.Open(self.__filename__)
        self.__canvas_name__ = canvasname + ".root"
        self.__canvas__ = TFile(self.__canvas_name__, "recreate")

    def __doc__(self):
        pass

    def draw_origin(self):
        file = self.__file__
        file.cd()
        h_xy = file.Get(self.hOrigin_xy)
        h_zy = file.Get(self.hOrigin_zy)
        h_zx = file.Get(self.hOrigin_zx)
        # mystyle = style("ATLAS")
        # gROOT.SetStyle("ATLAS")
        c1 = TCanvas("c1", "", 1000, 800)
        h_xy.GetXaxis().SetRangeUser(-100, 100)
        h_xy.GetXaxis().SetTitle("x/cm")
        h_xy.GetYaxis().SetTitle("y/cm")
        h_xy.GetYaxis().SetRangeUser(-100, 100)
        h_xy.Draw("colz")
        c1.SetLogz()
        c1.SaveAs("hOrigin_xy.png")
        c2 = TCanvas("c2", "", 1000, 800)
        h_zx.GetXaxis().SetTitle("z/cm")
        h_zx.GetYaxis().SetTitle("x/cm")
        c2.SetLogz()
        h_zx.Draw('colz')
        c2.SaveAs("hOrigin_zx.png")
        c3 = TCanvas("c3", "", 1000, 800)
        c3.SetLogz()
        h_zy.GetYaxis().SetTitle("z/cm")
        h_zy.GetXaxis().SetTitle("y/cm")
        h_zy.Draw("colz")
        c3.SaveAs("hOrigin_zy.png")

    def draw_nse(self, cent_low, cent_high, p_low, p_high, nse_low, nse_high):
        file = TFile.Open(self.__filename__)
        hnsige = file.Get(self.hist_sige)
        hnsige.GetAxis(0).SetRange(2, 2)
        h_sigma = hnsige.Projection(1, 2, 3, "e")
        h_sigma.SetName("h_sige_unlike")
        hnsige.GetAxis(0).SetRange(1, 1)
        h_sigma_like1 = hnsige.Projection(1, 2, 3, "e")
        h_sigma_like1.SetName("h_sige_like1")
        hnsige.GetAxis(0).SetRange(3, 3)
        h_sigma_like2 = hnsige.Projection(1, 2, 3, "e")
        h_sigma_like2.SetName("h_sige_like2")
        h_unlike = getPz(h_sigma, cent_low, cent_high, p_low, p_high)
        h1 = getPz(h_sigma_like1, cent_low, cent_high, p_low, p_high)
        h1.Sumw2()
        h2 = getPz(h_sigma_like2, cent_low, cent_high, p_low, p_high)
        h2.Sumw2()
        h_like = getavg(h1, h2)
        h = TH1D()
        h_unlike.Copy(h)
        h.Add(h_like, -1)
        canvas_name = "hnSigE_{}_{}_{:.1f}_{:.1f}".format(
            cent_low, cent_high, p_low, p_high)
        mystyle = style()
        # gROOT.SetStyle('ATLAS')
        c = TCanvas("c", "", 1000, 800)
        c.SetLogy()
        gStyle.SetOptStat(0)
        # gStyle.SetOptFit()
        y_max = h.GetMaximum()
        haxis = histo(-4, 6, 0.1, 3 * y_max, "n#sigma_{e}", "Counts")
        haxis.SetStats(0)
        haxis.Draw()
        h.SetMarkerStyle(24)
        h.Draw("esame")
        h_unlike.Draw("esame")
        h_unlike.SetMarkerStyle(25)
        h_unlike.SetMarkerColor(632)
        h_like.SetMarkerStyle(20)
        # h_like.SetMarkerColor()
        h_like.Draw("esame")
        mygaus = TF1("gaus", "gaus", -3, 3)
        h.Fit(mygaus, "r", "", -3, 3)
        chh1 = "#chi^{{2}} / ndf = {0:.2f} / {1}".format(
            mygaus.GetChisquare(), mygaus.GetNDF())
        chh5 = "N = {:.2E} #pm {:.2E}".format(
            mygaus.GetParameter(0), mygaus.GetParError(0))
        chh6 = "#mu = {:.3f} #pm {:.3f}".format(
            mygaus.GetParameter(1), mygaus.GetParError(1))
        chh7 = "#sigma = {:.3f} #pm {:.3f}".format(
            mygaus.GetParameter(2), mygaus.GetParError(2))
        cent_range = "{}-{}%".format(
            self.cent_list1[cent_high], self.cent_list2[cent_low])
        p_range = "%.1f<p<%.1f GeV/c" % (p_low, p_high)
        m = mygaus.Integral(nse_low, nse_high)
        N = mygaus.Integral(-10, 10)
        eff = "efficiency : {:.3f} #pm {:.3f}".format(
            m / N, math.sqrt(m * (1 - m / N)) / N)
        pv = pave(0.67, 0.65, 0.96, 0.9, cent_range + ' , ' +
                  p_range, chh1, chh5, chh6, chh7, eff)
        pv.Draw("same")
        l1 = drawLine(nse_low, 0, nse_low, mygaus.GetParameter(0), 2, 1, 632)
        l2 = drawLine(nse_high, 0, nse_high, mygaus.GetParameter(0), 2, 1, 632)
        leg = TLegend(0.75, 0.45, 0.9, 0.65)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.035)
        leg.SetTextFont(132)
        leg.AddEntry(h_unlike, "unlike-sign", "p")
        leg.AddEntry(h_like, "like-sign", "p")
        leg.AddEntry(h, "unlike - like", "p")
        leg.AddEntry(mygaus, "gaussian fit", "l")
        leg.Draw()
        c.Modified()
        canvas_file = self.__canvas__
        canvas_file.cd()
        c.Write(canvas_name)
        c.SaveAs(canvas_name + ".png")
        return mygaus

    def draw_nse_mean(self):
        cent_bin1 = [0, 4, 7]
        cent_bin2 = [3, 6, 8]
        p_bin1 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                  1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        p_bin2 = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                  1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        gr_mean = TGraphErrors()
        gr_width = TGraphErrors()
        gr_eff = TGraphErrors()
        for i in range(0, len(p_bin1)):
            f_gaus = self.draw_nse(cent_bin1[1], cent_bin2[
                                   1], p_bin1[i], p_bin2[i], -1, 2)
            mean = f_gaus.GetParameter(1)
            mean_err = f_gaus.GetParError(1)
            width = f_gaus.GetParameter(2)
            width_err = f_gaus.GetParError(2)
            gr_mean.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, mean)
            gr_mean.SetPointError(i, (p_bin2[i] - p_bin1[i]) * 0.5, mean_err)
            gr_width.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, width)
            gr_width.SetPointError(i, (p_bin2[i] - p_bin1[i]) * 0.5, width_err)
            p = 0.5 * (p_bin1[i] + p_bin2[i])
            if p < 1:
                nse_low = 0.5 * p - 1.5
            else:
                nse_low = -1
            m = f_gaus.Integral(nse_low, 2.0)
            N = f_gaus.Integral(-10.0, 10.0)
            gr_eff.SetPoint(i, p, m / N)
            gr_eff.SetPointError(
                i, (p_bin2[i] - p_bin1[i]) * 0.5, eff_err(m, N))
        gr_mean.SetMarkerStyle(20)

        gr_mean.SetMarkerSize(1.5)
        gr_mean.SetMarkerColor(632)
        gr_width.SetMarkerStyle(20)
        gr_width.SetMarkerSize(1.5)
        gr_width.SetMarkerColor(600)
        gr_eff.SetMarkerStyle(20)
        gr_eff.SetMarkerSize(1.5)
        gr_eff.SetMarkerColor(600)
        hx = histo(0, 10, -0.2, 1.2, "p_{T}(GeV/c)", "Mean/Width")
        canvas_file = self.__canvas__
        canvas_file.cd()
        c1 = TCanvas("c1", "c1", 1000, 800)
        hx.Draw()
        gr_width.Draw("psame")
        gr_mean.Draw("psame")
        c1.Write("nse")
        c1.SaveAs("nse.png")
        c2 = TCanvas("c2", "c2", 1000, 800)
        hx.GetYaxis().SetTitle("n#sigma_{e} cut efficiency")
        hx.Draw()
        gr_eff.Draw("psame")
        c2.Write("nse_eff")
        c2.SaveAs("nse_eff.png")

    def draw_hHitsFit(self, cent_low, cent_high, p_low, p_high):
        file = self.__file__
        hnhitsfit = file.Get(self.hist_nHitsFit)
        hnhitsfit.GetAxis(0).SetRange(2, 2)
        h_hHitsFit = hnhitsfit.Projection(1, 2, 3, "e")
        h_hHitsFit.SetName("h_nHitsFit_unlike")
        hnhitsfit.GetAxis(0).SetRange(1, 1)
        h_hHitsFit_like1 = hnhitsfit.Projection(1, 2, 3, "e")
        h_hHitsFit_like1.SetName("h_nHitsFit_like1")
        hnhitsfit.GetAxis(0).SetRange(3, 3)
        h_hHitsFit_like2 = hnhitsfit.Projection(1, 2, 3, "e")
        h_hHitsFit_like2.SetName("h_nHitsFit_like2")
        h_unlike = getPz(h_hHitsFit, cent_low, cent_high, p_low, p_high)
        h_unlike.SetName("unlike-sign")
        h1 = getPz(h_hHitsFit_like1, cent_low, cent_high, p_low, p_high)
        h1.Sumw2()
        h2 = getPz(h_hHitsFit_like2, cent_low, cent_high, p_low, p_high)
        h2.Sumw2()
        h_like = getavg(h1, h2)
        h_like.SetName("like-sign")
        h = TH1D()

        h_unlike.Copy(h)
        h.SetName("unlike - like ")
        h.Add(h_like, -1)
        ymax = h_unlike.GetMaximum()
        haxis = histo(10, 60, 0, ymax * 1.2, "nHitsFit", "Counts")
        c = TCanvas("c", "c", 1000, 800)
        haxis.SetStats(0)
        haxis.GetXaxis().SetNdivisions(520)
        haxis.Draw()
        h_cut = TH1D()
        h.Copy(h_cut)
        N = h.Integral()
        xmin = h_cut.GetXaxis().FindBin(25.0 + 1E-4)
        xmax = h_cut.GetXaxis().FindBin(1000)
        h_cut.GetXaxis().SetRange(xmin, xmax)
        m = h_cut.Integral(xmin, xmax)
        print m, N
        h_cut.SetFillColorAlpha(600, 0.5)
        h_cut.SetFillStyle(3002)
        h_cut.SetLineColor(600)
        h_cut.SetMarkerStyle(1)
        h_cut.SetMarkerColor(0)
        h_cut.Draw("histsameF2bar")

        # h_unlike.SetFillColorAlpha(1,1)
        # h_unlike.SetFillStyle(3004)
        h_unlike.SetLineColor(1)
        h_unlike.SetMarkerStyle(1)
        h_unlike.SetMarkerColor(0)
        h_unlike.Draw("histsame")

        # h_like.SetFillColorAlpha(632,0)
        # h_like.SetFillStyle(3005)
        h_like.SetLineColor(632)
        h_like.SetMarkerStyle(1)
        h_like.SetMarkerColor(0)
        h_like.Draw("histsame")
        h.SetLineColor(600)
        h.SetMarkerStyle(1)
        h.SetMarkerColor(600)
        h.Draw("histsame")
        leg = TLegend(0.65, 0.45, 0.9, 0.65)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.035)
        leg.SetTextFont(132)
        leg.AddEntry(h_unlike, "unlike-sign", "l")
        leg.AddEntry(h_like, "like-sign", "l")
        leg.AddEntry(h, "unlike - like", "l")
        leg.AddEntry(h_cut, "accepted tracks", "f")
        leg.Draw()
        data_set = "Run10 Au+Au@200GeV MinBias"
        cent_range = "{}-{}%".format(
            self.cent_list1[cent_high], self.cent_list2[cent_low])
        p_range = "%.1f<p<%.1f GeV/c" % (p_low, p_high)
        eff = "efficiency : {:.3f} #pm {:.3f}".format(m / N, eff_err(m, N))
        drawLatex(0.55, 0.77, data_set, 0.035, 132, 1)
        drawLatex(0.55, 0.72, cent_range + " , " + p_range, 0.035, 132, 1)
        drawLatex(0.55, 0.67, eff, 0.035, 132, 1)
        canvas_file = self.__canvas__
        canvas_file.cd()
        canvas_name = "hHitsFit_{}_{}_{:.1f}_{:.1f}".format(
            cent_low, cent_high, p_low, p_high)
        c.Write(canvas_name)
        c.SaveAs(canvas_name+".png")
        mn = []
        mn.append(m)
        mn.append(N)
        return mn

    def draw_HitsFit_cut_eff(self, cent_low, cent_high):
        cent_bin1 = [0, 4, 7]
        cent_bin2 = [3, 6, 8]
        p_bin1 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                  1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        p_bin2 = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                  1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        gr_eff = TGraphErrors()
        for i in range(0, len(p_bin1)):
            mn = self.draw_hHitsFit(cent_low, cent_high, p_bin1[i], p_bin2[i])
            m = mn[0]
            N = mn[1]
            gr_eff.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, m / N)
            gr_eff.SetPointError(
                i, (p_bin2[i] - p_bin1[i]) * 0.5, eff_err(m, N))
        gr_eff.SetMarkerStyle(20)
        gr_eff.SetMarkerSize(1.5)
        gr_eff.SetMarkerColor(600)
        hx = histo(0, 8, -0.2, 1.2, "p_{T}(GeV/c)", "hHitsFit cut efficiency")
        canvas_file = self.__canvas__
        canvas_file.cd()
        c2 = TCanvas("c2", "c2", 1000, 800)
        hx.Draw()
        gr_eff.Draw("psame")

        c2.Write("nHitsFit_cut_eff")
        c2.SaveAs("nHitsFit_cut.png")

    def draw_HitsDedx_cut_eff(self, cent_low, cent_high):
        cent_bin1 = [0, 4, 7]
        cent_bin2 = [3, 6, 8]
        p_bin1 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                  1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        p_bin2 = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                  1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        gr_eff = TGraphErrors()
        for i in range(0, len(p_bin1)):
            mn = self.draw_hHitsDedx(cent_low, cent_high, p_bin1[i], p_bin2[i])
            m = mn[0]
            N = mn[1]
            gr_eff.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, m / N)
            gr_eff.SetPointError(
                i, (p_bin2[i] - p_bin1[i]) * 0.5, eff_err(m, N))
        gr_eff.SetMarkerStyle(20)
        gr_eff.SetMarkerSize(1.5)
        gr_eff.SetMarkerColor(600)
        hx = histo(0, 8, -0.2, 1.2, "p_{T}(GeV/c)", "hHitsDEdx cut efficiency")
        canvas_file = self.__canvas__
        canvas_file.cd()
        c2 = TCanvas("c2", "c2", 1000, 800)
        hx.Draw()
        gr_eff.Draw("psame")
        c2.Write("nHitsDEdx_cut_eff")
        c2.SaveAs("nHitsDEdx_cut.png")

    def draw_hHitsDedx(self, cent_low, cent_high, p_low, p_high):
        file = self.__file__
        hnhitsdedx = file.Get(self.hist_nHitsDedx)
        hnhitsdedx.GetAxis(0).SetRange(2, 2)
        h_hHitsDedx = hnhitsdedx.Projection(1, 2, 3, "e")
        h_hHitsDedx.SetName("h_nHitsDedx_unlike")
        hnhitsdedx.GetAxis(0).SetRange(1, 1)
        h_hHitsDedx_like1 = hnhitsdedx.Projection(1, 2, 3, "e")
        h_hHitsDedx_like1.SetName("h_nHitsDedx_like1")
        hnhitsdedx.GetAxis(0).SetRange(3, 3)
        h_hHitsDedx_like2 = hnhitsdedx.Projection(1, 2, 3, "e")
        h_hHitsDedx_like2.SetName("h_nHitsDedx_like2")
        h_unlike = getPz(h_hHitsDedx, cent_low, cent_high, p_low, p_high)
        h_unlike.SetName("unlike-sign")
        h1 = getPz(h_hHitsDedx_like1, cent_low, cent_high, p_low, p_high)
        h1.Sumw2()
        h2 = getPz(h_hHitsDedx_like2, cent_low, cent_high, p_low, p_high)
        h2.Sumw2()
        h_like = getavg(h1, h2)
        h_like.SetName("like-sign")
        h = TH1D()

        h_unlike.Copy(h)
        h.SetName("unlike - like ")
        h.Add(h_like, -1)
        ymax = h_unlike.GetMaximum()
        haxis = histo(0, 50, 0, ymax * 1.2, "nHitsDedx", "Counts")
        c = TCanvas("c", "c", 1000, 800)
        haxis.SetStats(0)
        haxis.GetXaxis().SetNdivisions(520)
        haxis.Draw()
        h_cut = TH1D()
        h.Copy(h_cut)
        xmin = h_cut.GetXaxis().FindBin(15.0 + 1E-4)
        xmax = h_cut.GetXaxis().FindBin(1000)
        m = h_cut.Integral(xmin, xmax)
        N = h.Integral()
        h_cut.GetXaxis().SetRange(xmin, xmax)
        h_cut.SetFillColorAlpha(600, 0.5)
        h_cut.SetFillStyle(3002)
        h_cut.SetLineColor(600)
        h_cut.SetMarkerStyle(1)
        h_cut.SetMarkerColor(0)
        h_cut.Draw("histsameF2bar")

        # h_unlike.SetFillColorAlpha(1,1)
        # h_unlike.SetFillStyle(3004)
        h_unlike.SetLineColor(1)
        h_unlike.SetMarkerStyle(1)
        h_unlike.SetMarkerColor(0)
        h_unlike.Draw("histsame")

        # h_like.SetFillColorAlpha(632,0)
        # h_like.SetFillStyle(3005)
        h_like.SetLineColor(632)
        h_like.SetMarkerStyle(1)
        h_like.SetMarkerColor(0)
        h_like.Draw("histsame")
        h.SetLineColor(600)
        h.SetMarkerStyle(1)
        h.SetMarkerColor(600)
        h.Draw("histsame")
        leg = TLegend(0.65, 0.45, 0.9, 0.65)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.035)
        leg.SetTextFont(132)
        leg.AddEntry(h_unlike, "unlike-sign", "l")
        leg.AddEntry(h_like, "like-sign", "l")
        leg.AddEntry(h, "unlike - like", "l")
        leg.AddEntry(h_cut, "accepted tracks", "f")
        leg.Draw()
        data_set = "Run10 Au+Au@200GeV MinBias"
        cent_range = "{}-{}%".format(
            self.cent_list1[cent_high], self.cent_list2[cent_low])
        p_range = "%.1f<p<%.1f GeV/c" % (p_low, p_high)
        eff = "efficiency : {:.3f} #pm {:.3f}".format(m / N, eff_err(m, N))
        drawLatex(0.55, 0.77, data_set, 0.035, 132, 1)
        drawLatex(0.55, 0.72, cent_range + " , " + p_range, 0.035, 132, 1)
        drawLatex(0.55, 0.67, eff, 0.035, 132, 1)
        canvas_file = self.__canvas__
        canvas_file.cd()
        canvas_name = "hHitsDEdx_{}_{}_{:.1f}_{:.1f}".format(
            cent_low, cent_high, p_low, p_high)
        c.Write(canvas_name)
        c.SaveAs(canvas_name+".png")
        c.Write()
        c.SaveAs(canvas_name+".png")
        mn = []
        mn.append(m)
        mn.append(N)
        return mn

    def draw_mass(self, cent_low, cent_high, pt_low, pt_high):
        file = self.__file__
        hNum = file.Get("hNum")
        hDenM = file.Get("hDenM")
        hDenP = file.Get("hDenP")
        h_unlike = getPz(hNum, cent_low, cent_high, pt_low, pt_high)
        h2 = getPz(hDenM, cent_low, cent_high, pt_low, pt_high)
        h3 = getPz(hDenP, cent_low, cent_high, pt_low, pt_high)
        h_like = getavg(h2, h3)
        h = TH1D()
        h_unlike.Copy(h)
        h.Add(h_like, -1)
        c1 = TCanvas("c1", "", 1000, 800)
        c1.SetLogz()
        h.Draw()
        h_unlike.SetLineColor(632)
        h_like.SetLineColor(600)
        h_unlike.Draw("same")
        h_like.Draw("same")
        c1.SaveAs("pair_mass.png")
        file = self.__canvas__
        file.cd()
        c1.Write("mass")

    def get_subtract(self, histname, cent_low, cent_high):
        h_substract = TH1D()
        self.__file__.cd()
        h_init = self.__file__.Get(histname)
        cent_low_bin = h_init.GetAxis(1).FindBin(cent_low + 1E-4)
        cent_high_bin = h_init.GetAxis(1).FindBin(cent_high - 1E-4)
        h_init.GetAxis(1).SetRange(cent_low_bin, cent_high_bin)
        h_init.GetAxis(0).SetRange(2, 2)
        h_unlike = h_init.Projection(2, "e")
        h_unlike.SetName("h_unlike")
        h_init.GetAxis(0).SetRange(1, 1)
        h_like1 = h_init.Projection(2, "e")
        h_like1.SetName("h_like1")
        h_init.GetAxis(0).SetRange(3, 3)
        h_like2 = h_init.Projection(2, "e")
        h_like2.SetName("h_like2")
        h_like = getavg(h_like1, h_like2)
        h_like.SetName("h_like")
        h_unlike.Copy(h_substract)
        h_substract.SetName("unlike - like ")
        h_substract.Add(h_like, -1)
        h = []
        h.append(h_unlike)
        h.append(h_like)
        h.append(h_substract)
        # return h_substract
        return h

    def draw_tof_match(self, cent_low, cent_high):
        h_tpc_e = self.get_subtract(self.hist_sige, cent_low, cent_high)
        h_tpc_e.SetName("h_tpc_e")
        h_tofmatch = self.get_subtract(self.hist_tofmatch, cent_low, cent_high)
        h_tofmatch.SetName("h_tofmatch")
        ymax = h_tpc_e.GetMaximum()
        haxis = histo(0, 10, 1, ymax * 3.5, "p_{T}(GeV/c)", "Counts")
        haxis.SetNdivisions(510)

        h_tpc_e.SetLineColor(600)
        h_tpc_e.SetFillColor(600)
        h_tpc_e.SetFillStyle(1001)
        h_tpc_e.SetMarkerStyle(1)
        h_tpc_e.SetMarkerColor(600)
        h_tofmatch.SetLineColor(632)
        h_tofmatch.SetFillColorAlpha(632, 1)
        h_tofmatch.SetFillStyle(1001)
        h_tofmatch.SetMarkerStyle(1)
        h_tofmatch.SetMarkerColor(632)
        leg = TLegend(0.55, 0.5, 0.8, 0.65)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.035)
        leg.SetTextFont(132)
        leg.AddEntry(h_tpc_e, "TPC electron", "f")
        leg.AddEntry(h_tofmatch, "TOF match electon", "f")
        data_set = "Run10 Au+Au@200GeV MinBias"
        cent_range = "{}-{}% centrality".format(
            self.cent_list1[cent_high], self.cent_list2[cent_low])
        c = TCanvas("c", "c", 1000, 800)
        haxis.SetStats(0)
        haxis.Draw()
        gPad.SetLogy()
        h_tpc_e.Draw("histsamelfbar2")
        h_tofmatch.Draw("histsamelfbar2")
        leg.Draw("same")
        drawLatex(0.55, 0.72, data_set, 0.035, 132, 1)
        drawLatex(0.55, 0.67, cent_range, 0.035, 132, 1)
        c.SaveAs("tpc_e.png")
        self.__canvas__.cd()
        c.Write("tpc_e")
        c1 = TCanvas("c1", "c1", 1000, 800)
        gr = TGraphErrors()
        p_bin1 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                  1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
        p_bin2 = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2,
                  1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0]
        for i in range(0, len(p_bin1)):
            xlow = h_tpc_e.GetXaxis().FindBin(p_bin1[i] + 1E-3)
            xup = h_tpc_e.GetXaxis().FindBin(p_bin2[i] + 1E-3)
            N = h_tpc_e.Integral(xlow, xup)
            xlow = h_tofmatch.GetXaxis().FindBin(p_bin1[i] + 1E-3)
            xup = h_tofmatch.GetXaxis().FindBin(p_bin2[i] + 1E-3)
            m = h_tofmatch.Integral(xlow, xup)
            gr.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, m / N)
            gr.SetPointError(i, (p_bin2[i] - p_bin1[i]) * 0.5, eff_err(m, N))
        hx = histo(0, 10, 0, 1.5, "p_{T}(GeV/c)", "TOF match efficiency")
        hx.Draw()
        gr.Draw("psame")
        gr.SetMarkerStyle(20)
        gr.SetMarkerSize(1.5)
        gr.SetMarkerColor(632)
        canvas_file = self.__canvas__
        canvas_file.cd()
        c1.Write("tof_match_eff")
        c1.SaveAs("tof_match_eff.png")

    def draw_tof_cut(self, cent_low, cent_high, p_low, p_high):
        file = TFile.Open(self.__filename__)
        h_tofmatch = file.Get(self.hist_tofmatch)
        h_tofmatch.GetAxis(0).SetRange(2, 2)
        h_ovb = h_tofmatch.Projection(1, 2, 3, "e")
        h_ovb.SetName("h_sige_unlike")
        h_tofmatch.GetAxis(0).SetRange(1, 1)
        h_ovb_like1 = h_tofmatch.Projection(1, 2, 3, "e")
        h_ovb_like1.SetName("h_sige_like1")
        h_tofmatch.GetAxis(0).SetRange(3, 3)
        h_ovb_like2 = h_tofmatch.Projection(1, 2, 3, "e")
        h_ovb_like2.SetName("h_sige_like2")
        h_unlike = getPz(h_ovb, cent_low, cent_high, p_low, p_high)
        h1 = getPz(h_ovb_like1, cent_low, cent_high, p_low, p_high)
        h1.Sumw2()
        h2 = getPz(h_ovb_like2, cent_low, cent_high, p_low, p_high)
        h2.Sumw2()
        h_like = getavg(h1, h2)
        h = TH1D()
        h_unlike.Copy(h)
        h.Add(h_like, -1)
        canvas_name = "h_ovb_{}_{}_{:.1f}_{:.1f}".format(
            cent_low, cent_high, p_low, p_high)
        c = TCanvas("c", "", 1000, 800)
        # c.SetLogy()
        gStyle.SetOptStat(0)
        # gStyle.SetOptFit()
        y_max = h.GetMaximum()
        haxis = histo(0.5, 1.5, 0.1, 1.5 * y_max, "1/#beta", "Counts")
        haxis.SetStats(0)
        haxis.Draw()
        h.SetMarkerStyle(24)
        h.Draw("esame")
        h_unlike.Draw("esame")
        h_unlike.SetMarkerStyle(25)
        h_unlike.SetMarkerColor(632)
        h_like.SetMarkerStyle(20)
        # h_like.SetMarkerColor()
        h_like.Draw("esame")
        mygaus = TF1("gaus", "gaus", 0, 2)
        # mygaus.SetNpx(100000)
        h.Fit(mygaus, "r", "", 0.97, 1.03)
        chh1 = "#chi^{{2}} / ndf = {0:.2f} / {1}".format(
            mygaus.GetChisquare(), mygaus.GetNDF())
        chh5 = "N = {:.2E} #pm {:.2E}".format(
            mygaus.GetParameter(0), mygaus.GetParError(0))
        chh6 = "#mu = {:.3f} #pm {:.3f}".format(
            mygaus.GetParameter(1), mygaus.GetParError(1))
        chh7 = "#sigma = {:.3f} #pm {:.3f}".format(
            mygaus.GetParameter(2), mygaus.GetParError(2))
        cent_range = "{}-{}%".format(
            self.cent_list1[cent_high], self.cent_list2[cent_low])
        p_range = "%.1f<p_{T}<%.1f GeV/c" % (p_low, p_high)
        m = mygaus.Integral(-0.03, 0.03)
        N = mygaus.Integral(-10, 10)
        eff = "efficiency : {:.3f} #pm {:.3f}".format(
            m / N, math.sqrt(m * (1 - m / N)) / N)
        pv = pave(0.6, 0.65, 0.9, 0.9, cent_range + ' , ' +
                  p_range, chh1, chh5, chh6, chh7, eff)
        pv.Draw("same")
        # l1 = drawLine(0.97, 0, 0.97, mygaus.GetParameter(0), 2, 1, 632)
        # l2 = drawLine(1.03, 0, 1.03, mygaus.GetParameter(0), 2, 1, 632)
        leg = TLegend(0.7, 0.45, 0.9, 0.65)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.035)
        leg.SetTextFont(132)
        leg.AddEntry(h_unlike, "unlike-sign", "p")
        leg.AddEntry(h_like, "like-sign", "p")
        leg.AddEntry(h, "unlike - like", "p")
        leg.AddEntry(mygaus, "gaussian fit", "l")
        leg.Draw()
        c.Modified()
        canvas_file = self.__canvas__
        canvas_file.cd()
        c.Write(canvas_name)
        c.SaveAs(canvas_name + ".png")
        return mygaus

    def draw_tof_cut_eff(self, cent_low, cent_high):
        cent_bin1 = [0, 4, 7]
        cent_bin2 = [3, 6, 8]
        p_bin1 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                  1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        p_bin2 = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                  1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        gr_mean = TGraphErrors()
        gr_width = TGraphErrors()
        gr_eff = TGraphErrors()
        for i in range(0, len(p_bin1)):
            f_gaus = self.draw_tof_cut(cent_bin1[1], cent_bin2[
                                       1], p_bin1[i], p_bin2[i])
            mean = f_gaus.GetParameter(1)
            mean_err = f_gaus.GetParError(1)
            width = f_gaus.GetParameter(2)
            width_err = f_gaus.GetParError(2)
            gr_mean.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, mean)
            gr_mean.SetPointError(i, (p_bin2[i] - p_bin1[i]) * 0.5, mean_err)
            gr_width.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, width)
            gr_width.SetPointError(i, (p_bin2[i] - p_bin1[i]) * 0.5, width_err)
            m = f_gaus.Integral(0.97, 1.03)
            N = f_gaus.GetParameter(
                0) * math.sqrt(2 * math.pi) * f_gaus.GetParameter(2)
            print m, N
            gr_eff.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, m / N)
            gr_eff.SetPointError(
                i, (p_bin2[i] - p_bin1[i]) * 0.5, eff_err(m, N))
        gr_mean.SetMarkerStyle(20)
        gr_mean.SetMarkerSize(1.5)
        gr_mean.SetMarkerColor(632)
        gr_width.SetMarkerStyle(20)
        gr_width.SetMarkerSize(1.5)
        gr_width.SetMarkerColor(600)
        gr_eff.SetMarkerStyle(20)
        gr_eff.SetMarkerSize(1.5)
        gr_eff.SetMarkerColor(600)
        hx = histo(0, 10, -0.2, 1.2, "p_{T}(GeV/c)", "Mean/Width")
        canvas_file = self.__canvas__
        canvas_file.cd()
        c1 = TCanvas("c1", "c1", 1000, 800)
        hx.GetYaxis().SetTitle("1/#beta cut efficency")
        hx.Draw()
        gr_width.Draw("psame")
        gr_mean.Draw("psame")
        c1.Write("tof")
        c1.SaveAs("tof.png")
        c2 = TCanvas("c2", "c2", 1000, 800)
        hx.Draw()
        gr_eff.Draw("psame")
        c2.Write("tof_cut_eff")
        c2.SaveAs("tof_cut.png")

    def draw_bemc_match(self, cent_low, cent_high):
        h_tpc_e = self.get_subtract(self.hist_sige, cent_low, cent_high)
        h_tpc_e.SetName("h_tpc_e")
        h_bemcmatch = self.get_subtract(
            self.hist_bemcmatch, cent_low, cent_high)
        h_bemcmatch.SetName("h_bemcmatch")
        ymax = h_tpc_e.GetMaximum()
        haxis = histo(0, 10, 1, ymax * 3.5, "p_{T}(GeV/c)", "Counts")
        haxis.SetNdivisions(510)

        h_tpc_e.SetLineColor(600)
        h_tpc_e.SetFillColor(600)
        h_tpc_e.SetFillStyle(1001)
        h_tpc_e.SetMarkerStyle(1)
        h_tpc_e.SetMarkerColor(600)
        h_bemcmatch.SetLineColor(632)
        h_bemcmatch.SetFillColorAlpha(632, 1)
        h_bemcmatch.SetFillStyle(1001)
        h_bemcmatch.SetMarkerStyle(1)
        h_bemcmatch.SetMarkerColor(632)
        leg = TLegend(0.55, 0.5, 0.8, 0.65)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.035)
        leg.SetTextFont(132)
        leg.AddEntry(h_tpc_e, "TPC electron", "f")
        leg.AddEntry(h_bemcmatch, "BEMC match electon", "f")
        data_set = "Run10 Au+Au@200GeV MinBias"
        cent_range = "{}-{}% centrality".format(
            self.cent_list1[cent_high], self.cent_list2[cent_low])
        c = TCanvas("c", "c", 1000, 800)
        haxis.SetStats(0)
        haxis.Draw()
        gPad.SetLogy()
        h_tpc_e.Draw("histsamelfbar2")
        h_bemcmatch.Draw("histsamelfbar2")
        leg.Draw("same")
        drawLatex(0.55, 0.72, data_set, 0.035, 132, 1)
        drawLatex(0.55, 0.67, cent_range, 0.035, 132, 1)
        c.SaveAs("tpc_e.png")
        self.__canvas__.cd()
        c.Write("tpc_e")
        c1 = TCanvas("c1", "c1", 1000, 800)
        gr = TGraphErrors()
        p_bin1 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                  1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
        p_bin2 = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2,
                  1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0]
        for i in range(0, len(p_bin1)):
            xlow = h_tpc_e.GetXaxis().FindBin(p_bin1[i] + 1E-3)
            xup = h_tpc_e.GetXaxis().FindBin(p_bin2[i] + 1E-3)
            N = h_tpc_e.Integral(xlow, xup)
            xlow = h_bemcmatch.GetXaxis().FindBin(p_bin1[i] + 1E-3)
            xup = h_bemcmatch.GetXaxis().FindBin(p_bin2[i] + 1E-3)
            m = h_bemcmatch.Integral(xlow, xup)
            gr.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, m / N)
            gr.SetPointError(i, (p_bin2[i] - p_bin1[i]) * 0.5, eff_err(m, N))
        hx = histo(0, 10, 0, 1.5, "p_{T}(GeV/c)", "TOF match efficiency")
        hx.Draw()
        gr.Draw("psame")
        gr.SetMarkerStyle(20)
        gr.SetMarkerSize(1.5)
        gr.SetMarkerColor(632)
        canvas_file = self.__canvas__
        canvas_file.cd()
        c1.Write("bemc_match_eff")
        c1.SaveAs("bemc_match_eff.png")

    def draw_bemc_cut(self, cent_low, cent_high, p_low, p_high):
        file = TFile.Open(self.__filename__)
        h_bemcmatch = file.Get(self.hist_bemcmatch)
        h_bemcmatch.GetAxis(0).SetRange(2, 2)
        h_pve = h_bemcmatch.Projection(1, 2, 3, "e")
        h_pve.SetName("h_pve_unlike")
        h_bemcmatch.GetAxis(0).SetRange(1, 1)
        h_pve_like1 = h_bemcmatch.Projection(1, 2, 3, "e")
        h_pve_like1.SetName("h_pve_like1")
        h_bemcmatch.GetAxis(0).SetRange(3, 3)
        h_pve_like2 = h_bemcmatch.Projection(1, 2, 3, "e")
        h_pve_like2.SetName("h_pve_like2")
        h_unlike = getPz(h_pve, cent_low, cent_high, p_low, p_high)
        h1 = getPz(h_pve_like1, cent_low, cent_high, p_low, p_high)
        h1.Sumw2()
        h2 = getPz(h_pve_like2, cent_low, cent_high, p_low, p_high)
        h2.Sumw2()
        h_like = getavg(h1, h2)
        h = TH1D()
        h_unlike.Copy(h)
        h.Add(h_like, -1)
        h.SetName("unlike-like")
        canvas_name = "h_pve_{}_{}_{:.1f}_{:.1f}".format(
            cent_low, cent_high, p_low, p_high)
        c = TCanvas("c", "", 1000, 800)
        gStyle.SetOptStat(0)
        y_max = h.GetMaximum()
        haxis = histo(0, 3, 0.1, 1.5 * y_max, "P/E0", "Counts")
        haxis.SetStats(0)
        haxis.Draw()
        h.SetMarkerStyle(24)
        h.Draw("esame")
        h_unlike.Draw("esame")
        h_unlike.SetMarkerStyle(25)
        h_unlike.SetMarkerColor(632)
        h_like.SetMarkerStyle(20)
        h_like.Draw("esame")
        # h.Draw("histsameFbar")
        cent_range = "{}-{}%".format(
            self.cent_list1[cent_high], self.cent_list2[cent_low])
        p_range = "%.1f<p_{T}<%.1f GeV/c" % (p_low, p_high)
        N = h.Integral()
        bin1 = h.GetXaxis().FindBin(0 + 1E-4)
        bin2 = h.GetXaxis().FindBin(2 - 1E-4)
        m = h.Integral(bin1, bin2)
        eff = "efficiency : {:.3f} #pm {:.3f}".format(m / N, eff_err(m, N))
        pv = pave(0.6, 0.65, 0.9, 0.9, cent_range + ' , ' + p_range, eff)
        pv.Draw("same")
        leg = TLegend(0.7, 0.45, 0.9, 0.65)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.035)
        leg.SetTextFont(132)
        leg.AddEntry(h_unlike, "unlike-sign", "p")
        leg.AddEntry(h_like, "like-sign", "p")
        leg.AddEntry(h, "unlike - like", "p")
        # leg.AddEntry(h_shadow, "accepted electrons", "f")
        leg.Draw()
        c.Modified()
        canvas_file = self.__canvas__
        canvas_file.cd()
        c.Write(canvas_name)
        c.SaveAs(canvas_name + ".png")
        mn = []
        mn.append(m)
        mn.append(N)
        return mn

    def draw_bemc_cut_eff(self, cent_low, cent_high):
        cent_bin1 = [0, 4, 7]
        cent_bin2 = [3, 6, 8]
        p_bin1 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                  1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        p_bin2 = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                  1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        gr_eff = TGraphErrors()
        for i in range(0, len(p_bin1)):
            mn = self.draw_bemc_cut(cent_low, cent_high, p_bin1[i], p_bin2[i])
            m = mn[0]
            N = mn[1]
            gr_eff.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, m / N)
            gr_eff.SetPointError(
                i, (p_bin2[i] - p_bin1[i]) * 0.5, eff_err(m, N))
        gr_eff.SetMarkerStyle(20)
        gr_eff.SetMarkerSize(1.5)
        gr_eff.SetMarkerColor(600)
        hx = histo(0, 8, -0.2, 1.2, "p_{T}(GeV/c)", "BEMC cut efficiency")
        canvas_file = self.__canvas__
        canvas_file.cd()
        c2 = TCanvas("c2", "c2", 1000, 800)
        hx.Draw()
        gr_eff.Draw("psame")
        c2.Write("bemc_cut_eff")
        c2.SaveAs("bemc_cut.png")

    def draw_neta_cut(self, cent_low, cent_high, p_low, p_high):
        file = TFile.Open(self.__filename__)
        h_neta = file.Get(self.hist_neta)
        h_neta.GetAxis(0).SetRange(2, 2)
        hneta = h_neta.Projection(1, 2, 3, "e")
        hneta.SetName("hneta_unlike")
        h_neta.GetAxis(0).SetRange(1, 1)
        hneta_like1 = h_neta.Projection(1, 2, 3, "e")
        hneta_like1.SetName("hneta_like1")
        h_neta.GetAxis(0).SetRange(3, 3)
        hneta_like2 = h_neta.Projection(1, 2, 3, "e")
        hneta_like2.SetName("hneta_like2")
        h_unlike = getPz(hneta, cent_low, cent_high, p_low, p_high)
        h1 = getPz(hneta_like1, cent_low, cent_high, p_low, p_high)
        h1.Sumw2()
        h2 = getPz(hneta_like2, cent_low, cent_high, p_low, p_high)
        h2.Sumw2()
        h_like = getavg(h1, h2)
        h = TH1D()
        h_unlike.Copy(h)
        h.Add(h_like, -1)
        h.SetName("unlike-like")
        canvas_name = "hneta_{}_{}_{:.1f}_{:.1f}".format(
            cent_low, cent_high, p_low, p_high)
        c = TCanvas("c", "", 1000, 800)
        gStyle.SetOptStat(0)
        y_max = h.GetMaximum()
        haxis = histo(0, 10, 0.1, 1.5 * y_max, "n#eta", "Counts")
        haxis.SetStats(0)
        haxis.Draw()
        h.SetMarkerStyle(24)
        h.Draw("esame")
        h_unlike.Draw("esame")
        h_unlike.SetMarkerStyle(25)
        h_unlike.SetMarkerColor(632)
        h_like.SetMarkerStyle(20)
        h_like.Draw("esame")
        # h.Draw("histsameFbar")
        cent_range = "{}-{}%".format(
            self.cent_list1[cent_high], self.cent_list2[cent_low])
        p_range = "%.1f<p_{T}<%.1f GeV/c" % (p_low, p_high)
        N = h.Integral()
        bin1 = h.GetXaxis().FindBin(1 + 1E-4)
        bin2 = h.GetXaxis().FindBin(1000)
        m = h.Integral(bin1, bin2)
        eff = "efficiency : {:.3f} #pm {:.3f}".format(m / N, eff_err(m, N))
        pv = pave(0.6, 0.65, 0.9, 0.9, cent_range + ' , ' + p_range, eff)
        pv.Draw("same")
        leg = TLegend(0.7, 0.45, 0.9, 0.65)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.035)
        leg.SetTextFont(132)
        leg.AddEntry(h_unlike, "unlike-sign", "p")
        leg.AddEntry(h_like, "like-sign", "p")
        leg.AddEntry(h, "unlike - like", "p")
        # leg.AddEntry(h_shadow, "accepted electrons", "f")
        leg.Draw()
        c.Modified()
        canvas_file = self.__canvas__
        canvas_file.cd()
        c.Write(canvas_name)
        c.SaveAs(canvas_name + ".png")
        mn = []
        mn.append(m)
        mn.append(N)
        return mn

    def draw_neta_cut_eff(self, cent_low, cent_high):
        cent_bin1 = [0, 4, 7]
        cent_bin2 = [3, 6, 8]
        p_bin1 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                  1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        p_bin2 = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                  1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        gr_eff = TGraphErrors()
        for i in range(0, len(p_bin1)):
            mn = self.draw_neta_cut(cent_low, cent_high, p_bin1[i], p_bin2[i])
            m = mn[0]
            N = mn[1]
            gr_eff.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, m / N)
            gr_eff.SetPointError(
                i, (p_bin2[i] - p_bin1[i]) * 0.5, eff_err(m, N))
        gr_eff.SetMarkerStyle(20)
        gr_eff.SetMarkerSize(1.5)
        gr_eff.SetMarkerColor(600)
        hx = histo(0, 8, -0.2, 1.2, "p_{T}(GeV/c)", "n#eta cut efficiency")
        canvas_file = self.__canvas__
        canvas_file.cd()
        c2 = TCanvas("c2", "c2", 1000, 800)
        hx.Draw()
        gr_eff.Draw("psame")
        c2.Write("neta_cut_eff")
        c2.SaveAs("neta_cut.png")

    def draw_nphi_cut(self, cent_low, cent_high, p_low, p_high):
        file = TFile.Open(self.__filename__)
        h_nphi = file.Get(self.hist_nphi)
        h_nphi.GetAxis(0).SetRange(2, 2)
        hnphi = h_nphi.Projection(1, 2, 3, "e")
        hnphi.SetName("hnphi_unlike")
        h_nphi.GetAxis(0).SetRange(1, 1)
        hnphi_like1 = h_nphi.Projection(1, 2, 3, "e")
        hnphi_like1.SetName("hnphi_like1")
        h_nphi.GetAxis(0).SetRange(3, 3)
        hnphi_like2 = h_nphi.Projection(1, 2, 3, "e")
        hnphi_like2.SetName("hnphi_like2")
        h_unlike = getPz(hnphi, cent_low, cent_high, p_low, p_high)
        h1 = getPz(hnphi_like1, cent_low, cent_high, p_low, p_high)
        h1.Sumw2()
        h2 = getPz(hnphi_like2, cent_low, cent_high, p_low, p_high)
        h2.Sumw2()
        h_like = getavg(h1, h2)
        h = TH1D()
        h_unlike.Copy(h)
        h.Add(h_like, -1)
        h.SetName("unlike-like")
        canvas_name = "hnphi_{}_{}_{:.1f}_{:.1f}".format(
            cent_low, cent_high, p_low, p_high)
        c = TCanvas("c", "", 1000, 800)
        gStyle.SetOptStat(0)
        y_max = h.GetMaximum()
        haxis = histo(0, 10, 0.1, 1.5 * y_max, "n#phi", "Counts")
        haxis.SetStats(0)
        haxis.Draw()
        h.SetMarkerStyle(24)
        h.Draw("esame")
        h_unlike.Draw("esame")
        h_unlike.SetMarkerStyle(25)
        h_unlike.SetMarkerColor(632)
        h_like.SetMarkerStyle(20)
        h_like.Draw("esame")
        # h.Draw("histsameFbar")
        cent_range = "{}-{}%".format(
            self.cent_list1[cent_high], self.cent_list2[cent_low])
        p_range = "%.1f<p_{T}<%.1f GeV/c" % (p_low, p_high)
        N = h.Integral()
        bin1 = h.GetXaxis().FindBin(1 + 1E-4)
        bin2 = h.GetXaxis().FindBin(1000)
        m = h.Integral(bin1, bin2)
        eff = "efficiency : {:.3f} #pm {:.3f}".format(m / N, eff_err(m, N))
        pv = pave(0.6, 0.65, 0.9, 0.9, cent_range + ' , ' + p_range, eff)
        pv.Draw("same")
        leg = TLegend(0.7, 0.45, 0.9, 0.65)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.035)
        leg.SetTextFont(132)
        leg.AddEntry(h_unlike, "unlike-sign", "p")
        leg.AddEntry(h_like, "like-sign", "p")
        leg.AddEntry(h, "unlike - like", "p")
        # leg.AddEntry(h_shadow, "accepted electrons", "f")
        leg.Draw()
        c.Modified()
        canvas_file = self.__canvas__
        canvas_file.cd()
        c.Write(canvas_name)
        c.SaveAs(canvas_name + ".png")
        mn = []
        mn.append(m)
        mn.append(N)
        return mn

    def draw_nphi_cut_eff(self, cent_low, cent_high):
        cent_bin1 = [0, 4, 7]
        cent_bin2 = [3, 6, 8]
        p_bin1 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                  1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        p_bin2 = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                  1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        gr_eff = TGraphErrors()
        for i in range(0, len(p_bin1)):
            mn = self.draw_nphi_cut(cent_low, cent_high, p_bin1[i], p_bin2[i])
            m = mn[0]
            N = mn[1]
            gr_eff.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, m / N)
            gr_eff.SetPointError(
                i, (p_bin2[i] - p_bin1[i]) * 0.5, eff_err(m, N))
        gr_eff.SetMarkerStyle(20)
        gr_eff.SetMarkerSize(1.5)
        gr_eff.SetMarkerColor(600)
        hx = histo(0, 8, -0.2, 1.2, "p_{T}(GeV/c)", "n#phi cut efficiency")
        canvas_file = self.__canvas__
        canvas_file.cd()
        c2 = TCanvas("c2", "c2", 1000, 800)
        hx.Draw()
        gr_eff.Draw("psame")
        c2.Write("nphi_cut_eff")
        c2.SaveAs("nphi_cut.png")

    def draw_distphi_cut(self, cent_low, cent_high, p_low, p_high):
        file = TFile.Open(self.__filename__)
        h_phiDist = file.Get(self.hist_phiDist)
        h_phiDist.GetAxis(0).SetRange(2, 2)
        hphiDist = h_phiDist.Projection(1, 2, 3, "e")
        hphiDist.SetName("hphiDist_unlike")
        h_phiDist.GetAxis(0).SetRange(1, 1)
        hphiDist_like1 = h_phiDist.Projection(1, 2, 3, "e")
        hphiDist_like1.SetName("hphiDist_like1")
        h_phiDist.GetAxis(0).SetRange(3, 3)
        hphiDist_like2 = h_phiDist.Projection(1, 2, 3, "e")
        hphiDist_like2.SetName("hphiDist_like2")
        h_unlike = getPz(hphiDist, cent_low, cent_high, p_low, p_high)
        h1 = getPz(hphiDist_like1, cent_low, cent_high, p_low, p_high)
        h1.Sumw2()
        h2 = getPz(hphiDist_like2, cent_low, cent_high, p_low, p_high)
        h2.Sumw2()
        h_like = getavg(h1, h2)
        h = TH1D()
        h_unlike.Copy(h)
        h.Add(h_like, -1)
        h.SetName("unlike-like")
        canvas_name = "hphiDist_{}_{}_{:.1f}_{:.1f}".format(
            cent_low, cent_high, p_low, p_high)
        c = TCanvas("c", "", 1000, 800)
        gStyle.SetOptStat(0)
        y_max = h.GetMaximum()
        haxis = histo(-0.1, 0.1, 0.1, 1.5 * y_max, "#phi dist", "Counts")
        haxis.SetStats(0)
        haxis.Draw()
        h.SetMarkerStyle(24)
        h.Draw("esame")
        h_unlike.Draw("esame")
        h_unlike.SetMarkerStyle(25)
        h_unlike.SetMarkerColor(632)
        h_like.SetMarkerStyle(20)
        h_like.Draw("esame")
        # h.Draw("histsameFbar")
        cent_range = "{}-{}%".format(
            self.cent_list1[cent_high], self.cent_list2[cent_low])
        p_range = "%.1f<p_{T}<%.1f GeV/c" % (p_low, p_high)
        N = h.Integral()
        bin1 = h.GetXaxis().FindBin(-0.01 + 1E-4)
        bin2 = h.GetXaxis().FindBin(0.01 + 1E-4)
        m = h.Integral(bin1, bin2)
        eff = "efficiency : {:.3f} #pm {:.3f}".format(m / N, eff_err(m, N))
        pv = pave(0.6, 0.65, 0.9, 0.85, cent_range + ' , ' + p_range, eff)
        pv.Draw("same")
        leg = TLegend(0.7, 0.45, 0.9, 0.65)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.035)
        leg.SetTextFont(132)
        leg.AddEntry(h_unlike, "unlike-sign", "p")
        leg.AddEntry(h_like, "like-sign", "p")
        leg.AddEntry(h, "unlike - like", "p")
        # leg.AddEntry(h_shadow, "accepted electrons", "f")
        leg.Draw()
        c.Modified()
        canvas_file = self.__canvas__
        canvas_file.cd()
        c.Write(canvas_name)
        c.SaveAs(canvas_name + ".png")
        mn = []
        mn.append(m)
        mn.append(N)
        return mn

    def draw_distphi_cut_eff(self, cent_low, cent_high):
        cent_bin1 = [0, 4, 7]
        cent_bin2 = [3, 6, 8]
        p_bin1 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                  1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        p_bin2 = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                  1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        gr_eff = TGraphErrors()
        for i in range(0, len(p_bin1)):
            mn = self.draw_distphi_cut(
                cent_low, cent_high, p_bin1[i], p_bin2[i])
            m = mn[0]
            N = mn[1]
            gr_eff.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, m / N)
            gr_eff.SetPointError(
                i, (p_bin2[i] - p_bin1[i]) * 0.5, eff_err(m, N))
        gr_eff.SetMarkerStyle(20)
        gr_eff.SetMarkerSize(1.5)
        gr_eff.SetMarkerColor(600)
        hx = histo(0, 8, -0.2, 1.2, "p_{T}(GeV/c)", "#phi dist cut efficiency")
        canvas_file = self.__canvas__
        canvas_file.cd()
        c2 = TCanvas("c2", "c2", 1000, 800)
        hx.Draw()
        gr_eff.Draw("psame")
        c2.Write("phiDist_cut_eff")
        c2.SaveAs("phiDist_cut.png")

    def draw_distz_cut(self, cent_low, cent_high, p_low, p_high):
        file = TFile.Open(self.__filename__)
        h_zDist = file.Get(self.hist_zDist)
        h_zDist.GetAxis(0).SetRange(2, 2)
        hzDist = h_zDist.Projection(1, 2, 3, "e")
        hzDist.SetName("hzDist_unlike")
        h_zDist.GetAxis(0).SetRange(1, 1)
        hzDist_like1 = h_zDist.Projection(1, 2, 3, "e")
        hzDist_like1.SetName("hzDist_like1")
        h_zDist.GetAxis(0).SetRange(3, 3)
        hzDist_like2 = h_zDist.Projection(1, 2, 3, "e")
        hzDist_like2.SetName("hzDist_like2")
        h_unlike = getPz(hzDist, cent_low, cent_high, p_low, p_high)
        h1 = getPz(hzDist_like1, cent_low, cent_high, p_low, p_high)
        h1.Sumw2()
        h2 = getPz(hzDist_like2, cent_low, cent_high, p_low, p_high)
        h2.Sumw2()
        h_like = getavg(h1, h2)
        h = TH1D()
        h_unlike.Copy(h)
        h.Add(h_like, -1)
        h.SetName("unlike-like")
        canvas_name = "hzDist_{}_{}_{:.1f}_{:.1f}".format(
            cent_low, cent_high, p_low, p_high)
        c = TCanvas("c", "", 1000, 800)
        gStyle.SetOptStat(0)
        y_max = h.GetMaximum()
        haxis = histo(-5, 5, 0.1, 1.5 * y_max, "z dist(cm)", "Counts")
        haxis.SetStats(0)
        haxis.Draw()
        h.SetMarkerStyle(24)
        h.Draw("esame")
        h_unlike.Draw("esame")
        h_unlike.SetMarkerStyle(25)
        h_unlike.SetMarkerColor(632)
        h_like.SetMarkerStyle(20)
        h_like.Draw("esame")
        # h.Draw("histsameFbar")
        cent_range = "{}-{}%".format(
            self.cent_list1[cent_high], self.cent_list2[cent_low])
        p_range = "%.1f<p_{T}<%.1f GeV/c" % (p_low, p_high)
        N = h.Integral()
        bin1 = h.GetXaxis().FindBin(-2 + 1E-4)
        bin2 = h.GetXaxis().FindBin(2 + 1E-4)
        m = h.Integral(bin1, bin2)
        eff = "efficiency : {:.3f} #pm {:.3f}".format(m / N, eff_err(m, N))
        pv = pave(0.6, 0.65, 0.9, 0.85, cent_range + ' , ' + p_range, eff)
        pv.Draw("same")
        leg = TLegend(0.7, 0.45, 0.9, 0.65)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.035)
        leg.SetTextFont(132)
        leg.AddEntry(h_unlike, "unlike-sign", "p")
        leg.AddEntry(h_like, "like-sign", "p")
        leg.AddEntry(h, "unlike - like", "p")
        # leg.AddEntry(h_shadow, "accepted electrons", "f")
        leg.Draw()
        c.Modified()
        canvas_file = self.__canvas__
        canvas_file.cd()
        c.Write(canvas_name)
        c.SaveAs(canvas_name + ".png")
        mn = []
        mn.append(m)
        mn.append(N)
        return mn

    def draw_distz_cut_eff(self, cent_low, cent_high):
        cent_bin1 = [0, 4, 7]
        cent_bin2 = [3, 6, 8]
        p_bin1 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                  1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        p_bin2 = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                  1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        gr_eff = TGraphErrors()
        for i in range(0, len(p_bin1)):
            mn = self.draw_distz_cut(cent_low, cent_high, p_bin1[i], p_bin2[i])
            m = mn[0]
            N = mn[1]
            gr_eff.SetPoint(i, (p_bin1[i] + p_bin2[i]) * 0.5, m / N)
            gr_eff.SetPointError(
                i, (p_bin2[i] - p_bin1[i]) * 0.5, eff_err(m, N))
        gr_eff.SetMarkerStyle(20)
        gr_eff.SetMarkerSize(1.5)
        gr_eff.SetMarkerColor(600)
        hx = histo(0, 8, -0.2, 1.2, "p_{T}(GeV/c)", "zDist cut efficiency")
        canvas_file = self.__canvas__
        canvas_file.cd()
        c2 = TCanvas("c2", "c2", 1000, 800)
        hx.Draw()
        gr_eff.Draw("psame")
        c2.Write("zDist_cut_eff")
        c2.SaveAs("zDist_cut.png")


if __name__ == '__main__':
    a = pho_e("pho_e", "canvas")
    bin1 = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5]
    bin2 = [0.4, 0.6, 0.8, 1.0, 1.5, 2.5]
    cent_bin1 = [0, 4, 7]
    cent_bin2 = [3, 6, 8]
    a.draw_mass(7, 8, 0, 5)
    # h[2].Draw("same")
    # a.draw_hHitsDedx(0, 4, 0, 5)
    # a.draw_hHitsFit(0, 4, 0, 5)
    # a.draw_mass(0, 4, 0, 5)
    # a.draw_origin()
