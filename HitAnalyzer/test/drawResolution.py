import sys,os
import ROOT
import argparse

class HistogramDefinition:

    requiredFields = [ 'canvasName', 'histogramName', 'histogramTitle', 'variable', 'baseCuts', \
                        'xNbins', 'xMin', 'xMax' ]
    allFields = requiredFields + [ 'effCuts', 'xTitle', 'yTitle', 'yNbins', 'yMin', 'yMax', 'zMin', 'zMax' ]

    def __init__(self,name,inputDict):
        #
        self.name = name
        self.parameters = { x:None for x in HistogramDefinition.allFields }
        for k,v in inputDict.items():
            if k.startswith('__'):
                continue
            if k in HistogramDefinition.allFields:
                self.parameters[k] = v
            else:
                print("Warning: key",k,"is not a standard field name - ignoring the entry")
        #
        if self.parameters['canvasName']==None:
            self.parameters['canvasName'] = "c" + self.name[0].upper() + self.name[1:]
        if self.parameters['histogramName']==None:
            self.parameters['histogramName'] = "h" + self.name[0].upper() + self.name[1:]
        #
        for f in HistogramDefinition.requiredFields:
            assert ( f in self.parameters ) and self.parameters[f]!=None

    def __getitem__(self,field):
        if field in self.parameters:
            return self.parameters[field]
        return None

        

class HistogramDefinitions:

    def __init__(self):
        self.allDefinitions = { }
        self.allHistoNames = set()
        self.allCanvases = set()

    def add(self,hdef):
        assert not hdef.name in self.allDefinitions
        assert not hdef['histogramName'] in self.allHistoNames
        assert not hdef['canvasName'] in self.allCanvases
        self.allDefinitions[hdef.name] = hdef
        self.allHistoNames.add(hdef['histogramName'])
        self.allCanvases.add(hdef['canvasName'])
        

    def __getitem__(self,name):
        if name in self.allDefinitions:
            return self.allDefinitions[name]
        return None
        
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
    #print("&&".join([ c for c in cuts if ( c!=None and c.strip()!="" ) ]))
    return "&&".join([ c for c in cuts if ( c!=None and c.strip()!="" ) ])

def drawCuts(pave,title,cuts):
    t = pave.AddText(title)
    t.SetTextFont(42)
    t.SetTextSize(0.05)
    t.SetTextAlign(13)
    t = pave.AddText("")
    for ic,c in enumerate(cuts):
        l = "  " + c
        if ic<(len(cuts)-1):
            l += " &&"
        t = pave.AddText(l)

    
def drawCutPave(canvas,cuts,effcuts=None):
    indBaseCuts = cuts.split("&&")
    indEffCuts = None if effcuts==None else effcuts.split("&&")
    #print(indBaseCuts)
    #print(indEffCuts)
    canvas.cd(4)
    hpave = 0.05+(len(indBaseCuts)+1)*0.04
    if indEffCuts!=None:
        hpave += 0.05+(len(indEffCuts)+2)*0.04
        
    pave = ROOT.TPaveText(0.05,1.0-hpave,0.95,1.0)
    pave.SetBorderSize(0)
    pave.SetFillStyle(0)
    pave.SetTextFont(42)
    pave.SetTextSize(0.04)
    pave.SetTextAlign(13)
    drawCuts(pave,"Basic selection",indBaseCuts)
    if indEffCuts!=None:
        t = pave.AddText("")
        drawCuts(pave,"Efficiency selection",indEffCuts)
        
    pave.Draw()
    ROOT.gPad.Update()
    return pave

def drawHistoByDef(tree,hDef,extraCuts):
    result = { 'cnv' : None, 'histos' : { }, 'pave' : None }
    histos = result['histos']

    cnv = ROOT.TCanvas(hDef['canvasName'],hDef['canvasName'],1000,1000)
    result['cnv'] = cnv
    cnv.Divide(2,2)

    is1D = hDef['yNbins']==None
    #print(hDef['canvasName'],'is',is1D)
    
    ic = 0
    hEffVs = { }
    effCuts = hDef['effCuts']
    variable = hDef['variable']
    for mType in range(23,26):
        ic += 1
        cnv.cd(ic)
        ROOT.gPad.SetGridx(1)
        ROOT.gPad.SetGridy(1)
        if not is1D:
            ROOT.gPad.SetRightMargin(0.125)
        hName = hDef['histogramName'] + str(mType)
        hTitle = hDef['histogramTitle'] + " module type " +str(mType)
        nbx = hDef['xNbins']
        xmin = hDef['xMin']
        xmax = hDef['xMax']
        if is1D:
            tree.Draw("("+variable+")>>"+hName+"_1("+str(nbx)+","+str(xmin)+","+str(xmax)+")",
                        cutString(extraCuts,hDef['baseCuts'],"moduleType=="+str(mType)))
            if effCuts!=None:
                tree.Draw("("+variable+")>>"+hName+"_2("+str(nbx)+","+str(xmin)+","+str(xmax)+")",
                            cutString(extraCuts,hDef['baseCuts'],"moduleType=="+str(mType),effCuts))
        else:
            nby = hDef['yNbins']
            ymin = hDef['yMin']
            ymax = hDef['yMax']
            tree.Draw(variable+">>"+hName+"_1("+str(nbx)+","+str(xmin)+","+str(xmax)+","+ \
                          str(nby)+","+str(ymin)+","+str(ymax)+")",
                          cutString(extraCuts,hDef['baseCuts'],"moduleType=="+str(mType)))
            if effCuts!=None:
                tree.Draw(variable+">>"+hName+"_2("+str(nbx)+","+str(xmin)+","+str(xmax)+","+ \
                            str(nby)+","+str(ymin)+","+str(ymax)+")",
                            cutString(extraCuts,hDef['baseCuts'],"moduleType=="+str(mType),effCuts))
        histos[mType] = [ ROOT.gDirectory.Get(hName+"_1"), None, None, None ]
        if effCuts!=None:
            histos[mType][1] = ROOT.gDirectory.Get(hName+"_2")
        xtitle = hDef['xTitle'] if hDef['xTitle'] else variable
        ytitle = hDef['yTitle'] if hDef['yTitle'] else ""
        if is1D:
            ymin = hDef['yMin'] if hDef['yMin']!=None else 0.
            ymax = hDef['yMax'] if hDef['yMax']!=None else 1.05
            if effCuts!=None:
                histos[mType][2] = ROOT.gPad.DrawFrame(xmin,ymin,xmax,ymax)
                histos[mType][2].SetTitle(hTitle)
                histos[mType][2].GetXaxis().SetTitle(xtitle)
                histos[mType][2].GetYaxis().SetTitle(ytitle)
                histos[mType][3] = ROOT.TEfficiency(histos[mType][1],histos[mType][0])
                histos[mType][3].SetMarkerSize(0.3)
                histos[mType][3].Draw("same Z")
            else:
                histos[mType][0].SetTitle(hTitle)
                histos[mType][0].GetXaxis().SetTitle(xtitle)
                histos[mType][0].GetYaxis().SetTitle(ytitle)
                histos[mType][0].Draw()
        else:
            zmin = hDef['zMin'] if hDef['zMin']!=None else 0.
            zmax = hDef['zMax'] if hDef['zMax']!=None else 1.05
            if effCuts!=None:
                histos[mType][1].Divide(histos[mType][0])
                histos[mType][1].SetTitle(hTitle)
                histos[mType][1].GetXaxis().SetTitle(xtitle)
                histos[mType][1].GetYaxis().SetTitle(ytitle)
                histos[mType][1].SetMinimum(zmin)
                histos[mType][1].SetMaximum(zmax)
                histos[mType][1].Draw("ZCOL")
            else:
                histos[mType][0].SetTitle(hTitle)
                histos[mType][0].GetXaxis().SetTitle(xtitle)
                histos[mType][0].GetYaxis().SetTitle(ytitle)
                histos[mType][0].Draw("ZCOL")
        ROOT.gPad.Update()
    result['pave'] = drawCutPave(cnv,cutString(extraCuts,hDef['baseCuts']),cutString(hDef['effCuts']))
    return result

    
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--definitions', '-d', help='python module with dictionaries defining efficiency histograms', \
                        type=str, default=None)
parser.add_argument('--effVar', help='variable for extra efficiency plot (format <name>,<nbins>,<min>,<max>)', \
                        type=str, default=None)
parser.add_argument('--dxMax', help='max. local dx for efficiency plots', type=float, default=0.0075)
parser.add_argument('--cuts', '-c', help="basic cut string", type=str, default="")
parser.add_argument('--output', '-o', help='output directory for graphic output', type=str, default=None)
parser.add_argument('--sampleName', help='sample label for output', type=str, default=None)
parser.add_argument('--fitResiduals', '-f', \
                        help='comma-separated list of names of histogram sets with residuals to be fit', \
                        type=str, default=None)
parser.add_argument('file', help='input file', type=str, nargs=1, default=None)
args = parser.parse_args()
if args.output!=None:
    assert os.path.isdir(args.output)
fitResiduals = args.fitResiduals.split(",") if args.fitResiduals else [ ]
#
# load histogram definitions
#
allHDefs = HistogramDefinitions()
if args.definitions!=None:
    module = __import__(args.definitions)
    for n in dir(module):
        if n.startswith('__'):
            continue
        hDict = getattr(module,n)
        #print(n,type(hDict))
        #sys.exit()
        assert type(hDict)==dict
        hDef = HistogramDefinition(n,hDict)
        allHDefs.add(hDef)
        print("Added",hDef['canvasName'])

if args.effVar!=None:
    effVarDict = { }
    # split into string defining the variable(s) and (1 or 2) axis definition(s)
    fields1 = args.effVar.split(";")
    assert len(fields1)<=3
    effVarDict['variable'] = fields1[0]
    #effVarDict['canvasName'] = "cEffArg"
    #effVarDict['histogramName'] = "hEffArg"
    effVarDict['histogramTitle'] = "hEffArg"
    # x-axis
    fields2 = fields1[1].split(",")
    assert len(fields2)==3 
    effVarDict['xNbins'] = int(fields2[0])
    effVarDict['xMin'] = float(fields2[1])
    effVarDict['xMax'] = float(fields2[2])
    # check for info on y axis (== presence of 2nd variable)
    if len(fields1)==3:
        assert ":" in effVarDict['variable']
        fields3 = fields1[2].split(",")
        effVarDict['yNbins'] = int(fields3[0])
        effVarDict['yMin'] = float(fields3[1])
        effVarDict['yMax'] = float(fields3[2])
        effVarDict['xTitle'] = effVarDict['variable'].split(":")[1]
        effVarDict['yTitle'] = effVarDict['variable'].split(":")[2]
    else:
        effVarDict['yMin'] = 0.
        effVarDict['yMax'] = 1.05
        effVarDict['xTitle'] = effVarDict['variable']
        effVarDict['yTitle'] = 'efficiency'
    effVarDict['baseCuts'] = args.cuts
    effVarDict['effCuts'] = cutString("hasRecHit>0","abs(localPos.x()-rhLocalPos.x())<"+str(args.dxMax))
    #xxx = HistogramDefinition("effV",effVarDict)
    #print("xxx",xxx)
    #print("xxx",xxx.parameters)
    allHDefs.add(HistogramDefinition("effV",effVarDict))
    #print(allHDefs.allDefinitions.keys())
    #print(allHDefs.allCanvases)
    #print(allHDefs['hEffArg'])
        
#effVarName = None
#effVarAxis = None
#if args.effVar!=None:
#    fields1 = args.effVar.split(";")
#    assert len(fields1)==2
#    effVarName = fields1[0]
#    fields2 = fields1[1].split(",")
#    assert len(fields2)==3
#    effVarAxis = ( int(fields2[0]), float(fields2[1]), float(fields2[2]) )

#extraCuts = "abs(particleType)==13"
#extraCuts = "tof<12.5"
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


#cRes = ROOT.TCanvas("cRes","cRes",1000,1000)
#canvases.append(cRes)
#cRes.Divide(2,2)
#ic = 0
#for mType in range(23,26):
#    ic += 1
#    cRes.cd(ic)
#    simHitTree.Draw("(localPos.x()-rhLocalPos.x())>>hRes"+str(mType)+"(200,-0.1,0.1)", \
#                        cutString(extraCuts,"hasRecHit>0","moduleType=="+str(mType)))
#    histos[mType] = ROOT.gDirectory.Get("hRes"+str(mType))
#    histos[mType].SetTitle("Residuals module type "+str(mType))
#    histos[mType].GetXaxis().SetTitle("#Delta x [cm]")
#    f = fitHistogram(mType,histos[mType])

allObjects = [ ]
for hdef in allHDefs.allDefinitions.values():
    #
    # draw histograms according to definition
    #
    allObjects.append(drawHistoByDef(simHitTree,hdef,extraCuts))
    #
    # perform fit of resolution histogram
    #
    if hdef.name in fitResiduals:
        objects = allObjects[-1]
        cnv = objects['cnv']
        # fit and redraw each panel
        ic = 0
        for mType in range(23,26):
            ic += 1
            cnv.cd(ic)
            f = fitHistogram(mType,objects['histos'][mType][0])
#allObjects.append(drawHistoByDef(simHitTree,allHDefs['effX1'],extraCuts))
#allObjects.append(drawHistoByDef(simHitTree,allHDefs['eff2D1'],extraCuts))

#!# dxCut = "abs(localPos.x()-rhLocalPos.x())<"+str(args.dxMax)
#!# cEff2D = ROOT.TCanvas("cEff2D","cEff2D",1000,1000)
#!# canvases.append(cEff2D)
#!# cEff2D.Divide(2,2)
#!# ic = 0
#!# hEffs = { }
#!# for mType in range(23,26):
#!#     ic += 1
#!#     cEff2D.cd(ic)
#!#     simHitTree.Draw("localPos.y():localPos.x()>>hEff1"+str(mType)+"(60,-7.5,7.5,15,-7.5,7.5)", \
#!#                         cutString(extraCuts,"moduleType=="+str(mType)))
#!#     hEffs[mType] = [ ROOT.gDirectory.Get("hEff1"+str(mType)), None, None, None ]
#!#     simHitTree.Draw("localPos.y():localPos.x()>>hEff2"+str(mType)+"(60,-7.5,7.5,15,-7.5,7.5)", \
#!#                         cutString(extraCuts,"moduleType=="+str(mType), \
#!#                         "hasRecHit>0",dxCut))
#!#     hEffs[mType][1] = ROOT.gDirectory.Get("hEff2"+str(mType))
#!#     hEffs[mType][1].Divide(hEffs[mType][0])
#!#     hEffs[mType][1].SetTitle("Efficiency 2D module type "+str(mType))
#!#     hEffs[mType][1].GetXaxis().SetTitle("local x [cm]")
#!#     hEffs[mType][1].GetYaxis().SetTitle("local y [cm]")
#!#     hEffs[mType][1].GetZaxis().SetTitle("efficiency")
#!#     hEffs[mType][1].SetMaximum(1.)
#!#     hEffs[mType][1].SetMinimum(0.75)
#!#     hEffs[mType][1].Draw("zcol")
#!#     ROOT.gPad.Update()
#!# paves.append(drawCutPave(cEff2D,cutString(extraCuts),cutString("hasRecHit>0",dxCut)))

#!# cEffX = ROOT.TCanvas("cEffX","cEffX",1000,1000)
#!# canvases.append(cEffX)
#!# cEffX.Divide(2,2)
#!# ic = 0
#!# hEffXs = { }
#!# for mType in range(23,26):
#!#     ic += 1
#!#     cEffX.cd(ic)
#!#     ROOT.gPad.SetGridx(1)
#!#     ROOT.gPad.SetGridy(1)
#!#     simHitTree.Draw("localPos.x()>>hEffX1"+str(mType)+"(240,-7.5,7.5)", \
#!#                         cutString(extraCuts,"moduleType=="+str(mType)))
#!#     hEffXs[mType] = [ ROOT.gDirectory.Get("hEffX1"+str(mType)), None, None, None ]
#!#     simHitTree.Draw("localPos.x()>>hEffX2"+str(mType)+"(240,-7.5,7.5)", \
#!#                         cutString(extraCuts,"moduleType=="+str(mType), \
#!#                             "hasRecHit>0",dxCut))
#!#     hEffXs[mType][1] = ROOT.gDirectory.Get("hEffX2"+str(mType))
#!#     hEffXs[mType][2] = ROOT.gPad.DrawFrame(hEffXs[mType][0].GetXaxis().GetXmin(),0.5, \
#!#                                            hEffXs[mType][0].GetXaxis().GetXmax(),1.05)
#!#     hEffXs[mType][2].SetTitle("Efficiency 1D module type "+str(mType))
#!#     hEffXs[mType][2].GetXaxis().SetTitle("local x [cm]")
#!#     hEffXs[mType][2].GetYaxis().SetTitle("efficiency")
#!#     #hEffXs[mType][1].Divide(hEffXs[mType][0])
#!#     hEffXs[mType][3] = ROOT.TEfficiency(hEffXs[mType][1],hEffXs[mType][0])
#!#     hEffXs[mType][3].SetMarkerSize(0.3)
#!#     #hEffXs[mType][2].SetTitle("Efficiency 1D module type "+str(mType)+";local x [cm];efficiency")
#!#     #hEffXs[mType][2].SetTitle("Efficiency 1D module type "+str(mType))
#!#     #hEffXs[mType][2].GetXaxis().SetTitle("local x [cm]")
#!#     #hEffXs[mType][2].GetYaxis().SetTitle("efficiency")
#!#     #hEffXs[mType][2].SetMaximum(1.05)
#!#     #hEffXs[mType][2].SetMinimum(0.5)
#!#     hEffXs[mType][3].Draw("same Z")
#!#     ROOT.gPad.Update()
#!# paves.append(drawCutPave(cEffX,cutString(extraCuts),cutString("hasRecHit>0",dxCut)))

#!# if args.effVar!=None:
#!#     cEffV = ROOT.TCanvas("cEffV","cEffV",1000,1000)
#!#     canvases.append(cEffV)
#!#     cEffV.Divide(2,2)
#!#     ic = 0
#!#     hEffVs = { }
#!#     for mType in range(23,26):
#!#         ic += 1
#!#         cEffV.cd(ic)
#!#         ROOT.gPad.SetGridx(1)
#!#         ROOT.gPad.SetGridy(1)
#!#         simHitTree.Draw(effVarName+">>hEffV1"+str(mType)+ \
#!#                             "("+str(effVarAxis[0])+","+str(effVarAxis[1])+","+str(effVarAxis[2])+")", \
#!#                             cutString(extraCuts,"moduleType=="+str(mType)))
#!#         hEffVs[mType] = [ ROOT.gDirectory.Get("hEffV1"+str(mType)), None, None, None ]
#!#         simHitTree.Draw(effVarName+">>hEffV2"+str(mType)+ \
#!#                             "("+str(effVarAxis[0])+","+str(effVarAxis[1])+","+str(effVarAxis[2])+")", \
#!#                             cutString(extraCuts,"moduleType=="+str(mType),"hasRecHit>0",dxCut))
#!#         hEffVs[mType][1] = ROOT.gDirectory.Get("hEffV2"+str(mType))
#!#         hEffVs[mType][2] = ROOT.gPad.DrawFrame(hEffVs[mType][0].GetXaxis().GetXmin(),0.0, \
#!#                                                hEffVs[mType][0].GetXaxis().GetXmax(),1.05)
#!#         hEffVs[mType][2].SetTitle("Efficiency module type "+str(mType))
#!#         hEffVs[mType][2].GetXaxis().SetTitle(effVarName)
#!#         hEffVs[mType][2].GetYaxis().SetTitle("efficiency")
#!#         hEffVs[mType][3] = ROOT.TEfficiency(hEffVs[mType][1],hEffVs[mType][0])
#!#         hEffVs[mType][3].SetMarkerSize(0.3)
#!#         #hEffVs[mType][1].Divide(hEffVs[mType][0])
#!#         #hEffVs[mType][1].SetTitle("Efficiency module type "+str(mType))
#!#         #hEffVs[mType][1].GetXaxis().SetTitle(effVarName)
#!#         #hEffVs[mType][1].GetYaxis().SetTitle("efficiency")
#!#         #hEffVs[mType][1].SetMaximum(1.05)
#!#         #hEffVs[mType][1].SetMinimum(0.)
#!#         #hEffVs[mType][1].GetXaxis().SetTitle(effVarName)
#!#         hEffVs[mType][3].Draw("same Z")
#!#         ROOT.gPad.Update()
#!#     paves.append(drawCutPave(cEffV,cutString(extraCuts),cutString("hasRecHit>0",dxCut)))

if args.output!=None:
    for c in canvases:
        basename = os.path.join(args.output,c.GetName())
        if args.sampleName!=None:
            basename += "_" + args.sampleName
        c.SaveAs(basename+".pdf")
        c.SaveAs(basename+".png")
