# encoding:utf-8
# @Filename=plot_eff
# @Project=phe_eff
# @Date=2017-05-18.22:20
# @Author=sa
# code is far away from bugs with the god animal protecting
"""
I love animals. They taste delicious.
              ┏┓      ┏┓
            ┏┛┻━━━┛┻┓
            ┃      ☃      ┃
            ┃  ┳┛  ┗┳  ┃
            ┃      ┻      ┃
            ┗━┓      ┏━┛
                ┃      ┗━━━┓
                ┃  神兽保佑    ┣┓
                ┃　永无BUG！   ┏┛
                ┗┓┓┏━┳┓┏┛
                  ┃┫┫  ┃┫┫
                  ┗┻┛  ┗┻┛
"""
from ROOT import TFile, TCanvas, TH1D, gStyle, TF1, gROOT, TLegend, gPad, TGraphErrors, TMath

from mystyle import style
from usefun import *
from pho_e import pho_e
import math
a = pho_e("pho_e", "efficiency")
a.__doc__()
