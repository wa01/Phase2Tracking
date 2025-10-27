import sys,os
import ROOT
import argparse

def fitHistogram(mType,h):
    f1name = "f1"+str(mType)
    f1 = ROOT.TF1(f1name,"gaus(0)")
    f1.SetParameter(0,h.GetMaximum())
    f1.SetParLimits(0,0.,2*h.GetMaximum())
    f1.SetParameter(1,0.)
    f1.SetParameter(2,h.GetRMS()/10.)
    h.Fit(f1name)
    f2name = "f2"+str(mType)
    f2 = ROOT.TF1(f2name,"gaus(0)+gaus(3)")
    f2.SetParameter(0,f1.GetParameter(0))
    f2.SetParameter(1,f1.GetParameter(1))
    f2.SetParameter(2,f1.GetParameter(2))
    f2.SetParameter(3,f1.GetParameter(0)/100.)
    f2.SetParameter(4,f1.GetParameter(1))
    f2.SetParameter(5,5*f1.GetParameter(2))
    f2.SetParLimits(5,f1.GetParameter(2),10*f1.GetParameter(2))
    h.Fit(f2name)
    ROOT.gPad.SetLogy(1)
    ROOT.gPad.Update()


def cutString(*cuts):
    #print("&&".join([ c for c in cuts if c.strip()!="" ]))
    return "&&".join([ c for c in cuts if c.strip()!="" ])

def drawCuts(canvas,cuts,effcuts=None):
    individualCuts = cuts.split("&&")
    indEffCuts = None if effcuts==None else effcuts.split("&&")
    canvas.cd(4)
    hpave = 0.05+(len(individualCuts)+1)*0.04
    if indEffCuts!=None:
        hpave += 0.05+(len(indEffCuts)+2)*0.04
        
    pave = ROOT.TPaveText(0.05,1.0-hpave,0.95,1.0)
    pave.SetBorderSize(0)
    pave.SetFillStyle(0)
    t = pave.AddText("Basic selection")
    t.SetTextFont(42)
    t.SetTextSize(0.05)
    t.SetTextAlign(13)
    t = pave.AddText("")
    t.SetTextFont(42)
    t.SetTextSize(0.04)
    t.SetTextAlign(13)
    for ic,c in enumerate(individualCuts):
        l = "  " + c
        if ic<(len(individualCuts)-1):
            l += " &&"
        t = pave.AddText(l)
        t.SetTextFont(42)
        t.SetTextSize(0.04)
        t.SetTextAlign(13)
    if indEffCuts!=None:
        t = pave.AddText("")
        t.SetTextFont(42)
        t.SetTextSize(0.04)
        t.SetTextAlign(13)
        t = pave.AddText("Efficiency selection")
        t.SetTextFont(42)
        t.SetTextSize(0.05)
        t.SetTextAlign(13)
        t = pave.AddText("")
        t.SetTextFont(42)
        t.SetTextSize(0.04)
        t.SetTextAlign(13)
        for ic,c in enumerate(indEffCuts):
            l = "  " + c
            if ic<(len(individualCuts)-1):
                l += " &&"
            t = pave.AddText(l)
            t.SetTextFont(42)
            t.SetTextSize(0.04)
            t.SetTextAlign(13)
        
    pave.Draw()
    ROOT.gPad.Update()
    return pave

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--effVar', help='variable for extra efficiency plot (format <name>,<nbins>,<min>,<max>)', \
                        type=str, default=None)
parser.add_argument('--dxMax', help='max. local dx for efficiency plots', type=float, default=0.0075)
parser.add_argument('--cuts', '-c', help="basic cut string", type=str, default="")
parser.add_argument('--output', '-o', help='output directory for graphic output', type=str, default=None)
parser.add_argument('--sampleName', help='sample label for output', type=str, default=None)
parser.add_argument('file', help='input file', type=str, nargs=1, default=None)
args = parser.parse_args()
if args.output!=None:
    assert os.path.isdir(args.output)

effVarName = None
effVarAxis = None
if args.effVar!=None:
    fields1 = args.effVar.split(";")
    assert len(fields1)==2
    effVarName = fields1[0]
    fields2 = fields1[1].split(",")
    assert len(fields2)==3
    effVarAxis = ( int(fields2[0]), float(fields2[1]), float(fields2[2]) )

extraCuts = "abs(particleType)==13"
extraCuts = "tof<12.5"
extraCuts = args.cuts

ROOT.gROOT.ProcessLine(".L setTDRStyle.C")
ROOT.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetTitleFont(42)
ROOT.gStyle.SetTitleFontSize(0.04)
ROOT.gStyle.SetTitleX(0.10)
ROOT.gStyle.SetTitleY(1.00)
ROOT.gStyle.SetTitleAlign(13)
ROOT.gStyle.SetTitleBorderSize(0)
#ROOT.gStyle.SetTitleFillColor(0)
ROOT.gStyle.SetOptTitle(1)
tf = ROOT.TFile(args.file[0])
simHitTree = simHitTree = tf.Get("analysis").Get("SimHitTree")

canvases = [ ]
histos = { }
paves = [ ]
cRes = ROOT.TCanvas("cRes","cRes",1000,1000)
canvases.append(cRes)
cRes.Divide(2,2)
ic = 0
for mType in range(23,26):
    ic += 1
    cRes.cd(ic)
    simHitTree.Draw("(localPos.x()-rhLocalPos.x())>>hRes"+str(mType)+"(200,-0.1,0.1)", \
                        cutString(extraCuts,"hasRecHit>0","moduleType=="+str(mType)))
    histos[mType] = ROOT.gDirectory.Get("hRes"+str(mType))
    histos[mType].SetTitle("Residuals module type "+str(mType))
    histos[mType].GetXaxis().SetTitle("#Delta x [cm]")
    f = fitHistogram(mType,histos[mType])
    

dxCut = "abs(localPos.x()-rhLocalPos.x())<"+str(args.dxMax)
cEff2D = ROOT.TCanvas("cEff2D","cEff2D",1000,1000)
canvases.append(cEff2D)
cEff2D.Divide(2,2)
ic = 0
hEffs = { }
for mType in range(23,26):
    ic += 1
    cEff2D.cd(ic)
    simHitTree.Draw("localPos.y():localPos.x()>>hEff1"+str(mType)+"(60,-7.5,7.5,15,-7.5,7.5)", \
                        cutString(extraCuts,"moduleType=="+str(mType)))
    hEffs[mType] = [ ROOT.gDirectory.Get("hEff1"+str(mType)), None, None, None ]
    simHitTree.Draw("localPos.y():localPos.x()>>hEff2"+str(mType)+"(60,-7.5,7.5,15,-7.5,7.5)", \
                        cutString(extraCuts,"moduleType=="+str(mType), \
                        "hasRecHit>0",dxCut))
    hEffs[mType][1] = ROOT.gDirectory.Get("hEff2"+str(mType))
    hEffs[mType][1].Divide(hEffs[mType][0])
    hEffs[mType][1].SetTitle("Efficiency 2D module type "+str(mType))
    hEffs[mType][1].GetXaxis().SetTitle("local x [cm]")
    hEffs[mType][1].GetYaxis().SetTitle("local y [cm]")
    hEffs[mType][1].GetZaxis().SetTitle("efficiency")
    hEffs[mType][1].SetMaximum(1.)
    hEffs[mType][1].SetMinimum(0.75)
    hEffs[mType][1].Draw("zcol")
    ROOT.gPad.Update()
paves.append(drawCuts(cEff2D,cutString(extraCuts),cutString("hasRecHit>0",dxCut)))

cEffX = ROOT.TCanvas("cEffX","cEffX",1000,1000)
canvases.append(cEffX)
cEffX.Divide(2,2)
ic = 0
hEffXs = { }
for mType in range(23,26):
    ic += 1
    cEffX.cd(ic)
    ROOT.gPad.SetGridx(1)
    ROOT.gPad.SetGridy(1)
    simHitTree.Draw("localPos.x()>>hEffX1"+str(mType)+"(240,-7.5,7.5)", \
                        cutString(extraCuts,"moduleType=="+str(mType)))
    hEffXs[mType] = [ ROOT.gDirectory.Get("hEffX1"+str(mType)), None, None, None ]
    simHitTree.Draw("localPos.x()>>hEffX2"+str(mType)+"(240,-7.5,7.5)", \
                        cutString(extraCuts,"moduleType=="+str(mType), \
                            "hasRecHit>0",dxCut))
    hEffXs[mType][1] = ROOT.gDirectory.Get("hEffX2"+str(mType))
    hEffXs[mType][2] = ROOT.gPad.DrawFrame(hEffXs[mType][0].GetXaxis().GetXmin(),0.5, \
                                           hEffXs[mType][0].GetXaxis().GetXmax(),1.05)
    hEffXs[mType][2].SetTitle("Efficiency 1D module type "+str(mType))
    hEffXs[mType][2].GetXaxis().SetTitle("local x [cm]")
    hEffXs[mType][2].GetYaxis().SetTitle("efficiency")
    #hEffXs[mType][1].Divide(hEffXs[mType][0])
    hEffXs[mType][3] = ROOT.TEfficiency(hEffXs[mType][1],hEffXs[mType][0])
    hEffXs[mType][3].SetMarkerSize(0.3)
    #hEffXs[mType][2].SetTitle("Efficiency 1D module type "+str(mType)+";local x [cm];efficiency")
    #hEffXs[mType][2].SetTitle("Efficiency 1D module type "+str(mType))
    #hEffXs[mType][2].GetXaxis().SetTitle("local x [cm]")
    #hEffXs[mType][2].GetYaxis().SetTitle("efficiency")
    #hEffXs[mType][2].SetMaximum(1.05)
    #hEffXs[mType][2].SetMinimum(0.5)
    hEffXs[mType][3].Draw("same Z")
    ROOT.gPad.Update()
paves.append(drawCuts(cEffX,cutString(extraCuts),cutString("hasRecHit>0",dxCut)))

if args.effVar!=None:
    cEffV = ROOT.TCanvas("cEffV","cEffV",1000,1000)
    canvases.append(cEffV)
    cEffV.Divide(2,2)
    ic = 0
    hEffVs = { }
    for mType in range(23,26):
        ic += 1
        cEffV.cd(ic)
        ROOT.gPad.SetGridx(1)
        ROOT.gPad.SetGridy(1)
        simHitTree.Draw(effVarName+">>hEffV1"+str(mType)+ \
                            "("+str(effVarAxis[0])+","+str(effVarAxis[1])+","+str(effVarAxis[2])+")", \
                            cutString(extraCuts,"moduleType=="+str(mType)))
        hEffVs[mType] = [ ROOT.gDirectory.Get("hEffV1"+str(mType)), None, None, None ]
        simHitTree.Draw(effVarName+">>hEffV2"+str(mType)+ \
                            "("+str(effVarAxis[0])+","+str(effVarAxis[1])+","+str(effVarAxis[2])+")", \
                            cutString(extraCuts,"moduleType=="+str(mType),"hasRecHit>0",dxCut))
        hEffVs[mType][1] = ROOT.gDirectory.Get("hEffV2"+str(mType))
        hEffVs[mType][2] = ROOT.gPad.DrawFrame(hEffVs[mType][0].GetXaxis().GetXmin(),0.0, \
                                               hEffVs[mType][0].GetXaxis().GetXmax(),1.05)
        hEffVs[mType][2].SetTitle("Efficiency module type "+str(mType))
        hEffVs[mType][2].GetXaxis().SetTitle(effVarName)
        hEffVs[mType][2].GetYaxis().SetTitle("efficiency")
        hEffVs[mType][3] = ROOT.TEfficiency(hEffVs[mType][1],hEffVs[mType][0])
        hEffVs[mType][3].SetMarkerSize(0.3)
        #hEffVs[mType][1].Divide(hEffVs[mType][0])
        #hEffVs[mType][1].SetTitle("Efficiency module type "+str(mType))
        #hEffVs[mType][1].GetXaxis().SetTitle(effVarName)
        #hEffVs[mType][1].GetYaxis().SetTitle("efficiency")
        #hEffVs[mType][1].SetMaximum(1.05)
        #hEffVs[mType][1].SetMinimum(0.)
        #hEffVs[mType][1].GetXaxis().SetTitle(effVarName)
        hEffVs[mType][3].Draw("same Z")
        ROOT.gPad.Update()
    paves.append(drawCuts(cEffV,cutString(extraCuts),cutString("hasRecHit>0",dxCut)))

if args.output!=None:
    for c in canvases:
        basename = os.path.join(args.output,c.GetName())
        if args.sampleName!=None:
            basename += "_" + args.sampleName
        c.SaveAs(basename+".pdf")
        c.SaveAs(basename+".png")
