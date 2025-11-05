import sys,os
import ROOT
import argparse
from fnmatch import fnmatch

class HistogramDefinition:

    reqGenFields = [ 'canvasName', 'histogramName', 'histogramTitle', 'variable', 'baseCuts' ]
    reqHistFields = [ ]
    requiredFields = reqGenFields + reqHistFields
    optGenFields =  [ 'effCuts', 'logY' ]
    optHistFields = [ 'xNbins', 'xMin', 'xMax', 'xTitle', 'yTitle', 'yNbins', 'yMin', 'yMax', \
                          'zMin', 'zMax', 'display' ]
    optionalFields = optGenFields + optHistFields
    allFields = requiredFields + optionalFields
    allHistFields = reqHistFields + optHistFields

    def __init__(self,name,inputDict):
        #
        self.name = name
        self.parameters = { x:None for x in HistogramDefinition.allFields }
        for k,v in inputDict.items():
            if k.startswith('__'):
                continue
            if k in HistogramDefinition.allFields:
                #
                # general variable
                #
                self.parameters[k] = v
            elif k.startswith("mType"):
                #
                # mType-specific histogram parameters
                #
                assert type(v)==dict and ( not k in self.parameters ) and len(k)>5 and k[5:].isdigit()
                self.parameters[k] = { x:None for x in HistogramDefinition.allHistFields }
                for kh,vh in v.items():
                    if kh in HistogramDefinition.allHistFields:
                        self.parameters[k][kh] = vh
                    else:
                        print("Warning: key",kh, \
                                  "is not a standard histogram field name - ignoring the entry in", \
                                  self.name)
            else:
                print("Warning: key",k,"is not a standard field name - ignoring the entry in",self.name)
        #
        if self.parameters['canvasName']==None:
            self.parameters['canvasName'] = "c" + self.name[0].upper() + self.name[1:]
        if self.parameters['histogramName']==None:
            self.parameters['histogramName'] = "h" + self.name[0].upper() + self.name[1:]
        #
        for f in HistogramDefinition.requiredFields:
            assert ( f in self.parameters ) and self.parameters[f]!=None

    #def __getitem__(self,field):
    #    if field in self.parameters:
    #        return self.parameters[field]
    #    return None

    def getParameter(self,name,mType=None):
        result = None
        #
        # give priority to parameter specific to a module type
        #
        mTName = "mType"+str(mType) if mType!=None else None
        if ( mTName in self.parameters ) and ( name in self.parameters[mTName] ):
          result = self.parameters[mTName][name]
          if result!=None:
            return self.parameters[mTName][name]
        #
        # not found: use general parameter
        #
        if name in self.parameters:
            return self.parameters[name]
        return None

    def vetoMType(self,mType):
        ''' Check for an mType entry with display = False
        '''
        #
        # try to get 'display' parameter
        #
        name = 'display'
        mTName = "mType"+str(mType) if mType!=None else None
        if ( mTName in self.parameters ) and ( name in self.parameters[mTName] ):
            if self.parameters[mTName][name]!=None:
                return not self.parameters[mTName][name]
        return False

class HistogramDefinitions:

    def __init__(self):
        self.allDefinitions = { }
        self.allHistoNames = set()
        self.allCanvases = set()

    def add(self,hdef):
        assert not hdef.name in self.allDefinitions
        assert not hdef.getParameter('histogramName') in self.allHistoNames
        assert not hdef.getParameter('canvasName') in self.allCanvases
        self.allDefinitions[hdef.name] = hdef
        self.allHistoNames.add(hdef.getParameter('histogramName'))
        self.allCanvases.add(hdef.getParameter('canvasName'))
        

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

def drawHistoByDef(tree,hDef,extraCuts):
    result = { 'cnv' : None, 'histos' : { }, 'pave' : None }
    histos = result['histos']

    cnv = ROOT.TCanvas(hDef.getParameter('canvasName'),hDef.getParameter('canvasName'),1000,1000)
    result['cnv'] = cnv
    cnv.Divide(2,2)

    is1D = hDef.getParameter('yNbins')==None
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
        nbx = hDef.getParameter('xNbins',mType)
        xmin = hDef.getParameter('xMin',mType)
        xmax = hDef.getParameter('xMax',mType)
        if is1D:
            tree.Draw("("+variable+")>>"+hName+"_1("+str(nbx)+","+str(xmin)+","+str(xmax)+")",
                        cutString(extraCuts,hDef.getParameter('baseCuts'),"moduleType=="+str(mType)))
            if effCuts!=None:
                tree.Draw("("+variable+")>>"+hName+"_2("+str(nbx)+","+str(xmin)+","+str(xmax)+")",
                            cutString(extraCuts,hDef.getParameter('baseCuts'),"moduleType=="+str(mType),effCuts))
        else:
            nby = hDef.getParameter('yNbins',mType)
            ymin = hDef.getParameter('yMin',mType)
            ymax = hDef.getParameter('yMax',mType)
            tree.Draw(variable+">>"+hName+"_1("+str(nbx)+","+str(xmin)+","+str(xmax)+","+ \
                          str(nby)+","+str(ymin)+","+str(ymax)+")",
                          cutString(extraCuts,hDef.getParameter('baseCuts'),"moduleType=="+str(mType)))
            if effCuts!=None:
                tree.Draw(variable+">>"+hName+"_2("+str(nbx)+","+str(xmin)+","+str(xmax)+","+ \
                            str(nby)+","+str(ymin)+","+str(ymax)+")",
                            cutString(extraCuts,hDef.getParameter('baseCuts'),"moduleType=="+str(mType),effCuts))
        histos[mType] = [ ROOT.gDirectory.Get(hName+"_1"), None, None, None ]
        if effCuts!=None:
            histos[mType][1] = ROOT.gDirectory.Get(hName+"_2")
        xtitle = hDef.getParameter('xTitle',mType) if hDef.getParameter('xTitle',mType) else variable
        ytitle = hDef.getParameter('yTitle',mType) if hDef.getParameter('yTitle',mType) else ""
        if is1D:
            ymin = hDef.getParameter('yMin',mType) if hDef.getParameter('yMin',mType)!=None else 0.
            ymax = hDef.getParameter('yMax',mType) if hDef.getParameter('yMax',mType)!=None else 1.05
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
            zmin = hDef.getParameter('zMin',mType) if hDef.getParameter('zMin',mType)!=None else 0.
            zmax = hDef.getParameter('zMax',mType) if hDef.getParameter('zMax',mType)!=None else 1.05
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
        if hDef.getParameter('logY'):
            ROOT.gPad.SetLogy(1)
        ROOT.gPad.Update()
    result['pave'] = drawCutPave(cnv,hDef.getParameter('variable'), \
                                     cutString(extraCuts,hDef.getParameter('baseCuts')), \
                                     cutString(hDef.getParameter('effCuts')))
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
parser.add_argument('--selectedHistograms', help='comma-separated names of histogram definitions to be used', \
                        type=str, default='*')
parser.add_argument('--vetoedHistograms', help='comma-separated names of histogram definitions not to be used',
                        type=str, default='')
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
        #
        # check if in list of histograms to be displayed
        #
        selFlg = False
        for p in selectedHistoNames:
            if fnmatch(n,p):
                selFlg = True
                break
        if not selFlg:
            continue
        #
        # check if in list of histograms to be vetoed
        #
        selFlg = True
        for p in vetoedHistoNames:
            if fnmatch(n,p):
                selFlg = False
                break
        if not selFlg:
            continue
        #
        # add histogram
        #
        hDef = HistogramDefinition(n,hDict)
        allHDefs.add(hDef)
        print("Added",hDef.getParameter('canvasName'))

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
        effVarDict['yTitle'] = effVarDict['variable'].split(":")[0]
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


if args.output!=None:
    for c in [ x['cnv'] for x in allObjects ]:
        basename = os.path.join(args.output,c.GetName())
        if args.sampleName!=None:
            basename += "_" + args.sampleName
        print(basename)
        c.SaveAs(basename+".pdf")
        c.SaveAs(basename+".png")
