import sys
import ROOT

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


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
tf = ROOT.TFile(sys.argv[1])
simHitTree = simHitTree = tf.Get("analysis").Get("SimHitTree")

histos = { }
cRes = ROOT.TCanvas("c","c",1000,1000)
cRes.Divide(2,2)
ic = 0
for mType in range(23,26):
    ic += 1
    cRes.cd(ic)
    simHitTree.Draw("(localPos.x()-rhLocalPos.x())>>hRes"+str(mType)+"(200,-0.1,0.1)", \
                        "hasRecHit>0&&abs(particleType)==13&&moduleType=="+str(mType))
    histos[mType] = ROOT.gDirectory.Get("hRes"+str(mType))
    f = fitHistogram(mType,histos[mType])
    
    
cEff = ROOT.TCanvas("c","c",1000,1000)
cEff.Divide(2,2)
ic = 0
hEffs = { }
for mType in range(23,26):
    ic += 1
    cEff.cd(ic)
    simHitTree.Draw("localPos.y():localPos.x()>>hEff1"+str(mType)+"(60,-7.5,7.5,15,-7.5,7.5)", \
                        "abs(particleType)==13&&moduleType=="+str(mType))
    hEffs[mType] = [ ROOT.gDirectory.Get("hEff1"+str(mType)), None ]
    simHitTree.Draw("localPos.y():localPos.x()>>hEff2"+str(mType)+"(60,-7.5,7.5,15,-7.5,7.5)", \
                        "abs(particleType)==13&&moduleType=="+str(mType)+ \
                        "&&hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.005")
    hEffs[mType][1] = ROOT.gDirectory.Get("hEff2"+str(mType))
    hEffs[mType][1].Divide(hEffs[mType][0])
    hEffs[mType][1].SetMaximum(1.)
    hEffs[mType][1].SetMinimum(0.75)
    hEffs[mType][1].Draw("zcol")
    ROOT.gPad.Update()

cEffX = ROOT.TCanvas("c","c",1000,1000)
cEffX.Divide(2,2)
ic = 0
hEffXs = { }
for mType in range(23,26):
    ic += 1
    cEffX.cd(ic)
    simHitTree.Draw("localPos.x()>>hEffX1"+str(mType)+"(240,-7.5,7.5)", \
                        "abs(particleType)==13&&moduleType=="+str(mType))
    hEffXs[mType] = [ ROOT.gDirectory.Get("hEffX1"+str(mType)), None ]
    simHitTree.Draw("localPos.x()>>hEffX2"+str(mType)+"(240,-7.5,7.5)", \
                        "abs(particleType)==13&&moduleType=="+str(mType)+ \
                        "&&hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.005")
    hEffXs[mType][1] = ROOT.gDirectory.Get("hEffX2"+str(mType))
    hEffXs[mType][1].Divide(hEffXs[mType][0])
    hEffXs[mType][1].SetMaximum(1.)
    hEffXs[mType][1].SetMinimum(0.9)
    hEffXs[mType][1].Draw("zcol")
    ROOT.gPad.Update()
