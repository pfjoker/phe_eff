"""
version : 1.0
description : This file includes functions which are often used in plotting

"""
from ROOT import TH1D, TLatex, TPaveStats, TLine, TCanvas, gPad, gStyle
import math


def drawLatex(x=0.5, y=0.50, text="text", size=0.035, font=132, color=1):
    tl = TLatex(x, y, text)
    tl.SetTextSize(size)
    tl.SetTextFont(font)
    tl.SetTextColor(color)
    tl.DrawLatexNDC(x, y, text)


def histo(xmin=0.0, xmax=10.0, ymin=1.0, ymax=10.0, xtitle="X", ytitle="Y"):
    haxis = TH1D("haxis", "", 100, xmin, xmax)
    haxis.SetName(xtitle + ytitle)
    haxis.SetStats(0)
    haxis.SetLineColor(10)
    y = haxis.GetYaxis()
    y.SetRangeUser(ymin, ymax)
    y.SetTitle(ytitle)
    y.CenterTitle()
    y.SetTitleOffset(1.)
    y.SetTitleSize(0.04)
    # y.SetTitleFont()
    x = haxis.GetXaxis()
    x.SetTitle(xtitle)
    x.SetTitleOffset(.80)
    x.SetTitleSize(0.05)
    x.CenterTitle()
    return haxis


def frange(start, stop, step):
    i = start
    while i < stop:
        yield i
        i += step


def getPz(h, xmin, xmax, ymin, ymax, precision=1E-3):
    bin1 = h.GetXaxis().FindBin(xmin + precision)
    bin2 = h.GetXaxis().FindBin(xmax - precision)
    h.GetXaxis().SetRange(bin1, bin2)
    bin1 = h.GetYaxis().FindBin(ymin + precision)
    bin2 = h.GetYaxis().FindBin(ymax - precision)
    h.GetYaxis().SetRange(bin1, bin2)
    h2 = h.Project3D("z")
    return h2


def getavg(h1, h2, option=1):
    '''

    :param h1: 
    :param h2: 
    :param option:1 arithmetical average;2 geo average;

    :return: 
    '''
    h3 = TH1D()
    nbins1 = h1.GetNbinsX()
    nbins2 = h2.GetNbinsX()
    xmin1 = h1.GetXaxis().GetXmin()
    xmin2 = h2.GetXaxis().GetXmin()
    xmax1 = h1.GetXaxis().GetXmax()
    xmax2 = h2.GetXaxis().GetXmax()
    if (nbins1 == nbins2) and (xmin1 == xmin2) and (xmax1 == xmax2):
        h3.SetBins(nbins1, xmin1, xmax1)
        if (option == 2):
            for i in range(0, nbins1):
                a = h1.GetBinContent(i + 1)
                b = h2.GetBinContent(i + 1)
                c = 2 * math.sqrt(a * b)
                c_error = math.sqrt(a + b)
                h3.SetBinContent(i + 1, c)
                h3.SetBinError(i + 1, c_error)
        if (option == 1):
            h3.Add(h1, h2)

    return h3


def pave(x1, y1, x2, y2, *args, **text_style):
    pv = TPaveStats(x1, y1, x2, y2, "brNDC")
    pv.SetName("stats")
    pv.SetOptStat(1101)
    pv.SetBorderSize(0)
    pv.SetFillColor(10)
    pv.SetTextAlign(11)
    pv.SetLineColor(0)
    for i in args:
        pv.AddText(i)
    t_style = {'size': 0.035, 'font': 132, 'color': 1}
    if (text_style):
        t_style.update(text_style)
    pv.SetTextSize(t_style['size'])
    pv.SetTextFont(t_style['font'])
    pv.SetTextColor(t_style['color'])
    return pv


def fit_pave(x1, y1, x2, y2, fun, **text_style):
    chh1 = "#chi^{{2}} / ndf = {0:.2f} / {1}".format(
        fun.GetChisquare(), fun.GetNDF())
    chh5 = "N = {:.2E} #pm {:.2E}".format(
        fun.GetParameter(0), fun.GetParError(0))
    chh6 = "#mu = {:.3f} #pm {:.3f}".format(
        fun.GetParameter(1), fun.GetParError(1))
    chh7 = "#sigma = {:.3f} #pm {:.3f}".format(
        fun.GetParameter(2), fun.GetParError(2))
    return pave(x1, y1, x2, y2, chh1, chh5, chh6, chh7, **text_style)


def drawLine(xlow, ylow, xup, yup, lineWidth, lineStyle, lineColor):
    L1 = TLine(xlow, ylow, xup, yup)
    L1.SetLineWidth(lineWidth)
    L1.SetLineColor(lineColor)
    L1.SetLineStyle(lineStyle)
    L1.Draw("same")
    return L1


def setpad(left, right, top, bottom):
    gPad.SetFillColor(10)
    gPad.SetBorderMode(0)
    gPad.SetBorderSize(0)
    gPad.SetFrameFillColor(10)
    gPad.SetFrameBorderMode(0)
    gPad.SetFrameBorderSize(0)
    gPad.SetLeftMargin(left)
    gPad.SetRightMargin(right)
    gPad.SetTopMargin(top)
    gPad.SetBottomMargin(bottom)
    gPad.SetGridx(0)
    gPad.SetGridy(0)
    gStyle.SetOptStat(0)


def eff_err(m, N):
    v_e = (m + 1.) * (m + 2.) / ((N + 2.) * (N + 3.)) - math.pow((m + 1.) / (N + 2.), 2)
    return math.sqrt(v_e)


def shadow_hist_range(hist, xlow, xhigh, fill_style, fill_color):
    h = TH1D()
    h.SetName(hist.GetName() + "_shadow")
    x = h.GetXaxis()
    bin1 = x.FindBin(xlow)
    bin2 = x.FindBin(xhigh)
    x.SetRange(bin1, bin2)
    h.SetFillColor(fill_color)
    h.SetFillStyle(fill_style)
    return h


if __name__ == '__main__':
    c = TCanvas()
    haxis = histo()
    haxis.Draw()
    l = drawLine(1, 2, 5, 5, 1, 1, 632)
    a = 2
    c.SaveAs("xxx.png")
