import sys
from math import tan
from array import array
import ROOT
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--thickness', '-t', help="sensor thickness (in um)", type=float, default=200.)
parser.add_argument('--stripWidth', '-w', help="strip width (in um)", type=float, default=90.)
parser.add_argument('--refTanLorentzAngle',  help="tan Lorentz angle / B[T]", type=float, default=0.07)
parser.add_argument('--bfield', '-b', help='B field [T]', type=float, default=3.8)
parser.add_argument('--maxTanLambda',  help="max. track dx/dz", type=float, default=2.)
parser.add_argument('--nbTanLambda', help="number of bins for tan(#lambda)", type=int, default=200)
parser.add_argument('--nbPosX', help="number of bins for posigion", type=int, default=180)
args = parser.parse_args()

#thickness = 200     # sensor thickness in um
#width = 90          # strip width in um
#tanAlpha = 3.8*0.07 # tan(alpha) at 3.8T
tanAlpha = args.bfield*args.refTanLorentzAngle

ROOT.gStyle.SetOptStat(0)

hw = ROOT.TH1F("hw","charge width;tan(#lambda);width [#mum]",args.nbTanLambda,-args.maxTanLambda,args.maxTanLambda)
hc = ROOT.TH2F("hc","cluster size;#x_0 [#mum];tan(#lambda)",args.nbPosX,0.,args.stripWidth, \
                args.nbTanLambda,-args.maxTanLambda,args.maxTanLambda)
xcaxis = hc.GetXaxis()
ycaxis = hc.GetYaxis()
zcaxis = hc.GetZaxis()

dl = args.thickness*tanAlpha
for ibx in range(1,hc.GetNbinsX()+1):
    posx = xcaxis.GetBinCenter(ibx)
    for iby in range(1,hc.GetNbinsY()+1):
        tanTrack = ycaxis.GetBinCenter(iby)
        dt = args.thickness*tanTrack
        ibin = hc.GetBin(ibx,iby)
        xs = posx + dt/2.
        xbs = posx - dt/2. + dl
        x1 = min(xs,xbs)
        x2 = max(xs,xbs)
        hw.SetBinContent(iby,x2-x1)
        is1 = int(x1/args.stripWidth)
        is2 = int(x2/args.stripWidth)
        hc.SetBinContent(ibin,is2-is1+1)
        if (is2-is1+1)==6:
            print(posx,tanTrack)

cnvW = ROOT.TCanvas("cWidth","cWidth",600,600)
hw.Draw("cont")
cnvW.Update()

icmin = 1
icmax = int(hc.GetMaximum()+0.1)
hc.SetMinimum(icmin)
hc.SetMaximum(icmax)
print((icmax-icmin+1),0)
zcaxis.SetNdivisions(10100+icmax-icmin,0)
#zcaxis.CenterLabels(1)
ROOT.gStyle.SetPalette(ROOT.kBird)

cnvC = ROOT.TCanvas("cSize","cSize",600,600)
hc.Draw("zcol")
cnvC.Update()
