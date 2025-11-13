import sys,os
from histogramDefinition import *
import ROOT
import argparse
from fnmatch import fnmatch
        
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

def drawString(pave,title,s):
    t = pave.AddText(title)
    t.SetTextFont(42)
    t.SetTextSize(0.05)
    t.SetTextAlign(13)
    t = pave.AddText("")
    t = pave.AddText("  "+s)
    t = pave.AddText("")

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

    
def drawCutPave(canvas,variable,cuts,effcuts=None):
    indBaseCuts = cuts.split("&&")
    indEffCuts = None if effcuts==None else effcuts.split("&&")
    #print(indBaseCuts)
    #print(indEffCuts)
    canvas.cd(4)
    hpave = 0.05+2*0.04+0.05+(len(indBaseCuts)+1)*0.04
    if indEffCuts!=None:
        hpave += 0.05+(len(indEffCuts)+2)*0.04
        
    pave = ROOT.TPaveText(0.05,1.0-hpave,0.95,1.0)
    pave.SetBorderSize(0)
    pave.SetFillStyle(0)
    pave.SetTextFont(42)
    pave.SetTextSize(0.04)
    pave.SetTextAlign(13)
    drawString(pave,"Variable(s)",variable)
    drawCuts(pave,"Basic selection",indBaseCuts)
    if indEffCuts!=None:
        t = pave.AddText("")
        drawCuts(pave,"Efficiency selection",indEffCuts)
        
    pave.Draw()
    ROOT.gPad.Update()
    return pave

def fillHistoByDef(tree,hDef,extraCuts):
    histos = { }

    savedDir = ROOT.gDirectory
    ROOT.gROOT.cd()
    #cnv = ROOT.TCanvas(hDef.getParameter('canvasName'),hDef.getParameter('canvasName'),1000,1000)
    #result['cnv'] = cnv
    #cnv.Divide(2,2)

    is1D = hDef.getParameter('yNbins')==None
    isProfile = hDef.getParameter('profile')!=None and hDef.getParameter('profile')
    #print(hDef['canvasName'],'is',is1D)
    
    ic = 0
    hEffVs = { }
    effCuts = hDef.getParameter('effCuts')
    variable = hDef.getParameter('variable')
    for mType in range(23,26):
        #print("Checking mType",mType,"for",hDef.name)
        ic += 1
        #
        # draw histogram?
        #
        if hDef.vetoMType(mType):
            continue
        ##
        #cnv.cd(ic)
        #ROOT.gPad.SetGridx(1)
        #ROOT.gPad.SetGridy(1)
        #if not is1D:
        #    ROOT.gPad.SetRightMargin(0.125)
        hName = hDef.getParameter('histogramName') + str(mType)
        hTitle = hDef.getParameter('histogramTitle') + " module type " +str(mType)
        nbx = hDef.getParameter('xNbins',mType)
        xmin = hDef.getParameter('xMin',mType)
        xmax = hDef.getParameter('xMax',mType)
        #print("Starting for ",hDef.name,hName,hTitle)
        if is1D and ( not isProfile ):
            histos[mType] = [ ROOT.TH1F(hName+"_1",hName+"_1",nbx,xmin,xmax), None, None, None ]
            tree.Project(hName+"_1",variable, \
                        cutString(extraCuts,hDef.getParameter('baseCuts'),"moduleType=="+str(mType)))
            if effCuts!=None:
                histos[mType][1] = ROOT.TH1F(hName+"_2",hName+"_2",nbx,xmin,xmax)
                tree.Project(hName+"_2",variable, \
                            cutString(extraCuts,hDef.getParameter('baseCuts'),"moduleType=="+str(mType),effCuts))
                histos[mType][3] = ROOT.TEfficiency(histos[mType][1],histos[mType][0])
            else:
                # always keep final histogram in 4th position
                histos[mType][3] = histos[mType][0]
        elif isProfile:
            ymin = hDef.getParameter('yMin',mType)
            ymax = hDef.getParameter('yMax',mType)
            histos[mType] = [ ROOT.TProfile(hName+"_1",hName+"_1",nbx,xmin,xmax,ymin,ymax,'S'), None, None, None ]
            tree.Project(hName+"_1",variable, \
                          cutString(extraCuts,hDef.getParameter('baseCuts'),"moduleType=="+str(mType)))
            # always keep final histogram in 4th position
            histos[mType][3] = histos[mType][0]
        else:
            nby = hDef.getParameter('yNbins',mType)
            ymin = hDef.getParameter('yMin',mType)
            ymax = hDef.getParameter('yMax',mType)
            histos[mType] = [ ROOT.TH2F(hName+"_1",hName+"_1",nbx,xmin,xmax,nby,ymin,ymax), None, None, None ]
            tree.Project(hName+"_1",variable, \
                          cutString(extraCuts,hDef.getParameter('baseCuts'),"moduleType=="+str(mType)))
            if effCuts!=None:
                histos[mType][1] = ROOT.TH2F(hName+"_2",hName+"_2",nbx,xmin,xmax,nby,ymin,ymax)
                #tree.Draw(variable+">>"+hName+"_2("+str(nbx)+","+str(xmin)+","+str(xmax)+","+ \
                #            str(nby)+","+str(ymin)+","+str(ymax)+")",
                #            cutString(extraCuts,hDef.getParameter('baseCuts'),"moduleType=="+str(mType),effCuts))
                tree.Project(hName+"_2",variable, \
                            cutString(extraCuts,hDef.getParameter('baseCuts'),"moduleType=="+str(mType),effCuts))
                histos[mType][1].Divide(histos[mType][0])
                # always keep final histogram in 4th position
                histos[mType][3] = histos[mType][1]
            else:
                # always keep final histogram in 4th position
                histos[mType][3] = histos[mType][0]
        #print("Ending for ",hDef.name,hName,hTitle)
                
    savedDir.cd()
    return histos


def drawHistoByDef(histos,hDef,logY=False,same=False):
    result = { 'cnv' : None, 'histos' : histos, 'pave' : None }

    savedDir = ROOT.gDirectory
    ROOT.gROOT.cd()
    cnv = ROOT.TCanvas(hDef.getParameter('canvasName'),hDef.getParameter('canvasName'),1000,1000)
    result['cnv'] = cnv
    cnv.Divide(2,2)

    is1D = hDef.getParameter('yNbins')==None
    isProfile = hDef.getParameter('profile')!=None and hDef.getParameter('profile')
    #print(hDef['canvasName'],'is',is1D)
    
    ic = 0
    hEffVs = { }
    effCuts = hDef.getParameter('effCuts')
    variable = hDef.getParameter('variable')
    for mType in range(23,26):
        #print("Checking mType",mType,"for",hDef.name)
        ic += 1
        #
        # draw histogram?
        #
        if hDef.vetoMType(mType):
            continue
        #
        cnv.cd(ic)
        ROOT.gPad.SetGridx(1)
        ROOT.gPad.SetGridy(1)
        if not is1D:
            ROOT.gPad.SetRightMargin(0.125)
        hName = hDef.getParameter('histogramName') + str(mType)
        hTitle = hDef.getParameter('histogramTitle') + " module type " +str(mType)

        xtitle = hDef.getParameter('xTitle',mType) if hDef.getParameter('xTitle',mType) else variable
        nbx = hDef.getParameter('xNbins',mType)
        xmin = hDef.getParameter('xMin',mType)
        xmax = hDef.getParameter('xMax',mType)

        ytitle = hDef.getParameter('yTitle',mType) if hDef.getParameter('yTitle',mType) else ""
        if is1D and ( not isProfile ):
            ymin = hDef.getParameter('yMin',mType) if hDef.getParameter('yMin',mType)!=None else 0.
            ymax = hDef.getParameter('yMax',mType) if hDef.getParameter('yMax',mType)!=None else 1.05
            if effCuts!=None:
                if not same:
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
                histos[mType][0].Draw("same" if same else "")
        elif isProfile:
            histos[mType][0].SetTitle(hTitle)
            histos[mType][0].GetXaxis().SetTitle(xtitle)
            histos[mType][0].GetYaxis().SetTitle(ytitle)
            histos[mType][0].SetMarkerSize(0.5)
            #histos[mType][0].SetFillColor(ROOT.TColor.GetColorBright(ROOT.kGray))
            #histos[mType][0].SetFillColor(ROOT.kGray)
            histos[mType][0].Draw("same" if same else "")
        else:
            assert not same
            zmin = hDef.getParameter('zMin',mType) if hDef.getParameter('zMin',mType)!=None else 0.
            zmax = hDef.getParameter('zMax',mType) if hDef.getParameter('zMax',mType)!=None else 1.05
            if effCuts!=None:
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
        if logY or hDef.getParameter('logY'):
            ROOT.gPad.SetLogy(1)
        ROOT.gPad.Update()
    result['pave'] = drawCutPave(cnv,hDef.getParameter('variable'), \
                                     cutString(extraCuts,hDef.getParameter('baseCuts')), \
                                     cutString(hDef.getParameter('effCuts')))

    savedDir.cd()
    return result

def addHistogram(varString,cuts,effCuts=None,name='userHist'):
    extraHDict = { }
    # split into string defining the variable(s) and (1 or 2) axis definition(s)
    fields1 = varString.split(";")
    assert len(fields1)<=3
    extraHDict['variable'] = fields1[0]
    #extraHDict['canvasName'] = "cEffArg"
    #extraHDict['histogramName'] = "hEffArg"
    extraHDict['histogramTitle'] = name
    # x-axis
    fields2 = fields1[1].split(",")
    assert len(fields2)==3 
    extraHDict['xNbins'] = int(fields2[0])
    extraHDict['xMin'] = float(fields2[1])
    extraHDict['xMax'] = float(fields2[2])
    # check for info on y axis (== presence of 2nd variable)
    if len(fields1)==3:
        assert ":" in extraHDict['variable']
        fields3 = fields1[2].split(",")
        extraHDict['yNbins'] = int(fields3[0])
        extraHDict['yMin'] = float(fields3[1])
        extraHDict['yMax'] = float(fields3[2])
        extraHDict['xTitle'] = extraHDict['variable'].split(":")[1]
        extraHDict['yTitle'] = extraHDict['variable'].split(":")[0]
    else:
        extraHDict['yMin'] = 0.
        extraHDict['yMax'] = 1.05
        extraHDict['xTitle'] = extraHDict['variable']
        extraHDict['yTitle'] = 'efficiency' if effCuts!=None else 'events/bin'
    extraHDict['baseCuts'] = cuts
    if effCuts!=None:
        extraHDict['effCuts'] = effCuts
    #xxx = HistogramDefinition("effV",extraHDict)
    #print("xxx",xxx)
    #print("xxx",xxx.parameters)
    #allHDefs.add(HistogramDefinition("effV",extraHDict))
    #print(allHDefs.allDefinitions.keys())
    #print(allHDefs.allCanvases)
    #print(allHDefs['hEffArg'])
    return HistogramDefinition(name,extraHDict)


    
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--definitions', '-d', help='python module with dictionaries defining efficiency histograms', \
                        type=str, default=None)
parser.add_argument('--histogram', \
                    help='definition of extra histogram (format <variable>;<nbx>,<xmin>,<xmax>[;<nny>,<ymin>,<ymax>)', \
                    action='append', type=str, default=[])
parser.add_argument('--dxMax', help='max. local dx for efficiency plots', type=float, default=0.0075)
parser.add_argument('--cuts', '-c', help="basic cut string", type=str, default="")
parser.add_argument('--effCuts', '-e', help="basic cut string", type=str, default=None)
parser.add_argument('--output', '-o', help='output directory for graphic output', type=str, default=None)
parser.add_argument('--sampleName', help='sample label for output', type=str, default=None)
parser.add_argument('--fitResiduals', '-f', \
                        help='comma-separated list of names of histogram sets with residuals to be fit', \
                        type=str, default=None)
parser.add_argument('--selectedHistograms', help='comma-separated names of histogram definitions to be used', \
                        type=str, default='*')
parser.add_argument('--vetoedHistograms', help='comma-separated names of histogram definitions not to be used',
                        type=str, default='')
parser.add_argument('--logY', help='use log scale', action='store_true', default=False)
parser.add_argument('--list', '-l', help='list typle contents', action='store_true', default=False)
parser.add_argument('file', help='input file', type=str, nargs=1, default=None)
args = parser.parse_args()
if args.output!=None:
    assert os.path.isdir(args.output)
    print("***",args.output)
fitResiduals = args.fitResiduals.split(",") if args.fitResiduals else [ ]
selectedHistoNames = args.selectedHistograms.split(",")
vetoedHistoNames = args.vetoedHistograms.split(",")
#
# load histogram definitions
#
allHDefs = loadHistogramDefinitions(args.definitions,selectedHistoNames,vetoedHistoNames)

for ih,h in enumerate(args.histogram):
    allHDefs.add(addHistogram(h,args.cuts,args.effCuts,name="userH"+str(ih+1)))
#!#     varEffDict = { }
#!#     # split into string defining the variable(s) and (1 or 2) axis definition(s)
#!#     fields1 = args.varEff.split(";")
#!#     assert len(fields1)<=3
#!#     varEffDict['variable'] = fields1[0]
#!#     #varEffDict['canvasName'] = "cEffArg"
#!#     #varEffDict['histogramName'] = "hEffArg"
#!#     varEffDict['histogramTitle'] = "hEffArg"
#!#     # x-axis
#!#     fields2 = fields1[1].split(",")
#!#     assert len(fields2)==3 
#!#     varEffDict['xNbins'] = int(fields2[0])
#!#     varEffDict['xMin'] = float(fields2[1])
#!#     varEffDict['xMax'] = float(fields2[2])
#!#     # check for info on y axis (== presence of 2nd variable)
#!#     if len(fields1)==3:
#!#         assert ":" in varEffDict['variable']
#!#         fields3 = fields1[2].split(",")
#!#         varEffDict['yNbins'] = int(fields3[0])
#!#         varEffDict['yMin'] = float(fields3[1])
#!#         varEffDict['yMax'] = float(fields3[2])
#!#         varEffDict['xTitle'] = varEffDict['variable'].split(":")[1]
#!#         varEffDict['yTitle'] = varEffDict['variable'].split(":")[0]
#!#     else:
#!#         varEffDict['yMin'] = 0.
#!#         varEffDict['yMax'] = 1.05
#!#         varEffDict['xTitle'] = varEffDict['variable']
#!#         varEffDict['yTitle'] = 'efficiency'
#!#     varEffDict['baseCuts'] = args.cuts
#!#     varEffDict['effCuts'] = cutString("hasRecHit>0","abs(localPos.x()-rhLocalPos.x())<"+str(args.dxMax))
#!#     #xxx = HistogramDefinition("effV",varEffDict)
#!#     #print("xxx",xxx)
#!#     #print("xxx",xxx.parameters)
#!#     allHDefs.add(HistogramDefinition("effV",varEffDict))
#!#     #print(allHDefs.allDefinitions.keys())
#!#     #print(allHDefs.allCanvases)
#!#     #print(allHDefs['hEffArg'])
        
#extraCuts = "abs(particleType)==13"
#extraCuts = "tof<12.5"
extraCuts = args.cuts

ROOT.gROOT.ProcessLine(".L setTDRStyle.C")
ROOT.gROOT.ProcessLine(".L floatMod.C+")
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
if args.list:
    simHitTree.Print()
    sys.exit()

canvases = [ ]
histos = { }
paves = [ ]

allObjects = [ ]
for cName in allHDefs.canvasNames():
    same = False
    cHistos = { }
    for hName in allHDefs.byCanvas[cName]:
        cHistos[hName] = fillHistoByDef(simHitTree,allHDefs.byCanvas[cName][hName],extraCuts)
        allObjects.append(drawHistoByDef(cHistos[hName],allHDefs.byCanvas[cName][hName],logY=args.logY,same=same))
        same = True
#    yMin = min([ x.GetMinimum() for x in cHistos
#sys.exit()

#!# allObjects = [ ]
#!# for hdef in allHDefs.allDefinitions.values():
#!#     #
#!#     # draw histograms according to definition
#!#     #
#!#     histos = fillHistoByDef(simHitTree,hdef,extraCuts)
#!#     allObjects.append(drawHistoByDef(histos,hdef))
#!#     #
#!#     # perform fit of resolution histogram
#!#     #
#!#     if hdef.name in fitResiduals:
#!#         objects = allObjects[-1]
#!#         cnv = objects['cnv']
#!#         # fit and redraw each panel
#!#         ic = 0
#!#         for mType in range(23,26):
#!#             ic += 1
#!#             cnv.cd(ic)
#!#             f = fitHistogram(mType,objects['histos'][mType][0])


if args.output!=None:
    for c in [ x['cnv'] for x in allObjects ]:
        basename = os.path.join(args.output,c.GetName())
        if args.sampleName!=None:
            basename += "_" + args.sampleName
        print(basename)
        c.SaveAs(basename+".pdf")
        c.SaveAs(basename+".png")
