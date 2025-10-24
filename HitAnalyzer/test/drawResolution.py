import sys
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
    #print(cuts)
    #print("&&".join([ c for c in cuts if c.strip()!="" ]))
    return "&&".join([ c for c in cuts if c.strip()!="" ])

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--effVar', help='variable for extra efficiency plot (format <name>,<nbins>,<min>,<max>)', \
                        type=str, default=None)
parser.add_argument('--dxMax', help='max. local dx for efficiency plots', type=float, default=0.0075)
parser.add_argument('--cuts', '-c', help="basic cut string", type=str, default="")
parser.add_argument('file', help='input file', type=str, nargs=1, default=None)
args = parser.parse_args()

effVarName = None
effVarAxis = None
if args.effVar!=None:
    fields = args.effVar.split(",")
    assert len(fields)==4
    effVarName = fields[0]
    effVarAxis = ( int(fields[1]), float(fields[2]), float(fields[3]) )

extraCuts = "abs(particleType)==13"
extraCuts = "tof<12.5"
extraCuts = args.cuts
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
tf = ROOT.TFile(args.file[0])
simHitTree = simHitTree = tf.Get("analysis").Get("SimHitTree")

histos = { }
cRes = ROOT.TCanvas("c","c",1000,1000)
cRes.Divide(2,2)
ic = 0
for mType in range(23,26):
    ic += 1
    cRes.cd(ic)
    simHitTree.Draw("(localPos.x()-rhLocalPos.x())>>hRes"+str(mType)+"(200,-0.1,0.1)", \
                        cutString(extraCuts,"hasRecHit>0","moduleType=="+str(mType)))
    histos[mType] = ROOT.gDirectory.Get("hRes"+str(mType))
    f = fitHistogram(mType,histos[mType])
    

dxCut = "abs(localPos.x()-rhLocalPos.x())<"+str(args.dxMax)
cEff = ROOT.TCanvas("cEff2D","cEff2D",1000,1000)
cEff.Divide(2,2)
ic = 0
hEffs = { }
for mType in range(23,26):
    ic += 1
    cEff.cd(ic)
    simHitTree.Draw("localPos.y():localPos.x()>>hEff1"+str(mType)+"(60,-7.5,7.5,15,-7.5,7.5)", \
                        cutString(extraCuts,"moduleType=="+str(mType)))
    hEffs[mType] = [ ROOT.gDirectory.Get("hEff1"+str(mType)), None ]
    simHitTree.Draw("localPos.y():localPos.x()>>hEff2"+str(mType)+"(60,-7.5,7.5,15,-7.5,7.5)", \
                        cutString(extraCuts,"moduleType=="+str(mType), \
                        "hasRecHit>0",dxCut))
    hEffs[mType][1] = ROOT.gDirectory.Get("hEff2"+str(mType))
    hEffs[mType][1].Divide(hEffs[mType][0])
    hEffs[mType][1].SetMaximum(1.)
    hEffs[mType][1].SetMinimum(0.75)
    hEffs[mType][1].Draw("zcol")
    ROOT.gPad.Update()

cEffX = ROOT.TCanvas("cEffX","cEffX",1000,1000)
cEffX.Divide(2,2)
ic = 0
hEffXs = { }
for mType in range(23,26):
    ic += 1
    cEffX.cd(ic)
    simHitTree.Draw("localPos.x()>>hEffX1"+str(mType)+"(240,-7.5,7.5)", \
                        cutString(extraCuts,"moduleType=="+str(mType)))
    hEffXs[mType] = [ ROOT.gDirectory.Get("hEffX1"+str(mType)), None ]
    simHitTree.Draw("localPos.x()>>hEffX2"+str(mType)+"(240,-7.5,7.5)", \
                        cutString(extraCuts,"moduleType=="+str(mType), \
                            "hasRecHit>0",dxCut))
    hEffXs[mType][1] = ROOT.gDirectory.Get("hEffX2"+str(mType))
    hEffXs[mType][1].Divide(hEffXs[mType][0])
    hEffXs[mType][1].SetMaximum(1.05)
    hEffXs[mType][1].SetMinimum(0.)
    hEffXs[mType][1].Draw()
    ROOT.gPad.Update()

if args.effVar!=None:
    cEffV = ROOT.TCanvas("cEffV","cEffV",1000,1000)
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
        hEffVs[mType] = [ ROOT.gDirectory.Get("hEffV1"+str(mType)), None ]
        simHitTree.Draw(effVarName+">>hEffV2"+str(mType)+ \
                            "("+str(effVarAxis[0])+","+str(effVarAxis[1])+","+str(effVarAxis[2])+")", \
                            cutString(extraCuts,"moduleType=="+str(mType),"hasRecHit>0",dxCut))
        hEffVs[mType][1] = ROOT.gDirectory.Get("hEffV2"+str(mType))
        hEffVs[mType][1].Divide(hEffVs[mType][0])
        hEffVs[mType][1].SetMaximum(1.05)
        hEffVs[mType][1].SetMinimum(0.)
        hEffVs[mType][1].GetXaxis().SetTitle(effVarName)
        hEffVs[mType][1].Draw()
        ROOT.gPad.Update()

