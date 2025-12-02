import sys
from math import sqrt,log
import ROOT
import argparse

class FitHistogram:

    def interpolate(y,y1,y2,x1=0.,x2=1.):
        ''' Linear interpolation between two points (x1,y1) and (x2,y2) 
            to yield x corresponding to the target value y.
        '''
        return x1 + (y-y1)/(y2-y1)*(x2-x1)
        
    def __init__(self,histogram):
        self.hist = histogram
        self.graph = None
        self.intGraph = None

    def fwhm(self):
        ''' Return lowest / highest x corresponding to ymax/2, and ymax/2
        '''
        #
        # target y value (1/2 maximum) and bin width (assume constant bin width)
        #
        y = self.hist.GetBinContent(self.hist.GetMaximumBin())/2.
        dx = self.hist.GetBinWidth(1)
        #
        # preset result (low / high x values)
        #
        xlow = None
        xhigh = None
        #
        # loop over pairs of adjacent histogram bins
        #
        x1 = self.hist.GetBinLowEdge(1)
        y1 = self.hist.GetBinContent(1)
        for i in range(2,self.hist.GetNbinsX()):
            # check if target value between contents of neighbouring bins
            x2 = self.hist.GetBinLowEdge(i)
            y2 = self.hist.GetBinContent(i)
            first = None
            if y>=y1 and y<y2:
                first = True
            elif y<=y1 and y>y2:
                first = False

            if first!=None:
                # calculate interpolated x value
                x = FitHistogram.interpolate(y,y1,y2,x1,x2)
                # store lowest and highest x corresponding to y
                if first and xlow==None:
                    xlow = x
                elif not first:
                    xhigh = x
            # move to next bin
            x1 = x2
            y1 = y2

        #
        # require result ( assumes that first and last bins are < ymax/2 ) and
        # correct for bin width / 2 ( assumes that bin value corresponds to center of bin )
        assert xlow!=None and xhigh!=None
        return (xlow+dx,xhigh+dx,y)

    def quantile(self,q):
        
class FitCanvas:

    def __init__(self,canvas,mtype):
        self.canvas = canvas
        self.mtype = mtype
        self.pad = None
        self.fhist = None
        
        padName = mainCnv.GetName()+str(mtype-22)
        for o in mainCnv.GetListOfPrimitives():
            if o.InheritsFrom(ROOT.TPad.Class()) and o.GetName()==padName:
                assert self.pad==None
                self.pad = o
        assert self.pad!=None

        histNameBase = "h"+mainCnv.GetName()[1:]+str(mtype)
        self.fhist = None
        for o in self.pad.GetListOfPrimitives():
            if o.InheritsFrom(ROOT.TH1.Class()) and o.GetName().startswith(histNameBase):
                assert self.fhist==None
                self.fhist = FitHistogram(o)
        assert self.fhist!=None

    def histogram(self):
        print("fhist",self.fhist,self.fhist.hist)
        return self.fhist.hist
            
    def fwhm(self):
        return self.fhist.fwhm()
    
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--moduleType', '-m', help='module type', type=int, choices=[23,24,25], default=23)
parser.add_argument('--fwhm', help='determining FWHM', action='store_true', default=False)
parser.add_argument('--slices', '-s', help='fit in slices in y', action='store_true', default=False)
parser.add_argument('file', help='input file', type=str, nargs=1, default=None)
args = parser.parse_args()

tf = ROOT.TFile(args.file[0])
tf.ls()
mainCnv = None
for k in tf.GetListOfKeys():
    o = k.ReadObj()
    if o.IsA()==ROOT.TCanvas.Class():
        assert mainCnv==None
        mainCnv = k.ReadObj()
assert mainCnv!=None

canvases = [ ]
slices = [ ]

mainCnv.Draw()
fitCanvas = FitCanvas(mainCnv,args.moduleType)
if not args.slices:
    assert fitCanvas.histogram().GetDimension()==1
    #print(fitCanvas.fwhm())
    x1,x2,y = fitCanvas.fwhm()
    print("<x> = {:6.1f}um, dx = {:6.1f}um, sig = {:6.1f}um ( interval {:6.4f} - {:6.4f}cm )".format( \
            10000*(x1+x2)/2.,10000*(x2-x1),10000*(x2-x1)/2/sqrt(2*log(2)),x1,x2))

    fitCanvas.pad.cd()
    fwhmArrow = ROOT.TArrow()
    fwhmArrow.SetLineColor(2)
    fwhmArrow.SetLineWidth(2)
    fwhmArrow.DrawArrow(x1,y,x2,y,0.005,"<>")
    fitCanvas.pad.Update()

else:
    hist2D = fitCanvas.histogram()
    assert hist2D.GetDimension()==2
    yaxis = hist2D.GetYaxis()
    nby = hist2D.GetNbinsY()
    ncol = int(sqrt(nby))
    while ncol*ncol<nby:
        ncol += 1
    canvases = [ ]
    #cnv = ROOT.TCanvas("cslice","cslice",1200,1200)
    #cnv.Divide(ncol,ncol)
    #slices = [ ]
    fwhmArrow = ROOT.TArrow()
    fwhmArrow.SetLineColor(2)
    fwhmArrow.SetLineWidth(2)
    hNameBase = hist2D.GetName()
    hTitleBase = hist2D.GetTitle()
    for iby in range(nby):
        #print(fitCanvas.histogram())
        #cnv.cd(iby+1)
        hName = hNameBase+"_"+str(iby+1)
        hTitle = hTitleBase+" (slice "+str(iby+1)+")"
        #print("new hName",hName)
        hProj = hist2D.ProjectionX(hName,iby+1,iby+1)
        #slices.append(hProj)
        #print("xxx",hProj)
        slices.append(FitHistogram(hProj))
        #print(slices[-1].hist)
        hProj.SetTitle(hTitle)
        #slices[-1].hist.Draw()
        #if hProj.GetEntries()<100:
        if hProj.GetMaximum()<100:
            continue
        canvases.append(ROOT.TCanvas("c"+hName[1:],"c"+hName[1:],800,800))
        canvases[-1].SetGridx(1)
        canvases[-1].SetGridx(2)
        hProj.Draw()
        x1,x2,y = slices[-1].fwhm()
        fmt = "{:5.2f} < dx/dz < {:5.2f} :"
        fmt += " <x> = {:6.1f}um, dx = {:6.1f}um, sig = {:6.1f}um ( interval {:6.4f} - {:6.4f}cm )"
        print(fmt.format(yaxis.GetBinLowEdge(iby+1),yaxis.GetBinUpEdge(iby+1), \
                10000*(x1+x2)/2.,10000*(x2-x1),10000*(x2-x1)/2/sqrt(2*log(2)),x1,x2))
        fwhmArrow.DrawArrow(x1,y,x2,y,0.005,"<>")
        ROOT.gPad.Update()
    #cnv.Update()
    
