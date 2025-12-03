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
        self.graph_ = None
        self.cumulativeGraph_ = None
        self.cumNormGraph_ = None

    def graph(self):
        ''' Histogram converted to list of (x,y) coordinates with x = bin center and y = bin contents.
            Ignores under- / overflow.
        '''
        if self.graph_==None:
            self.graph_ = [ ]

            for i in range(self.hist.GetNbinsX()):
                xbin = self.hist.GetBinLowEdge(i+1) + self.hist.GetBinWidth(i+1)/2.
                self.graph_.append( (xbin,self.hist.GetBinContent(i+1)) )

        return self.graph_

    def cumulativeGraph(self):
        ''' Histogram converted to list of (x,y) coordinates with x = bin center and y = cumulated bin contents.
            Ignores under- / overflow.
        '''
        if self.cumulativeGraph_==None:
            self.cumulativeGraph_ = [ ]
            sum = 0.
            for i in range(self.hist.GetNbinsX()):
                xbin = self.hist.GetBinLowEdge(i+1) + self.hist.GetBinWidth(i+1)/2.
                sum += self.hist.GetBinContent(i+1)
                self.cumulativeGraph_.append( (xbin,sum) )

        return self.cumulativeGraph_

    def cumulativeNormGraph(self):
        ''' Histogram converted to list of (x,y) coordinates with x = bin center and y = cumulated bin contents.
            Normalized to total histogram contents (including under- / overflow).
        '''
        if self.cumNormGraph_==None:
            cg = self.cumulativeGraph()
            sum = self.hist.GetSumOfWeights()
            self.cumNormGraph_ = [ ( x,y/sum ) for x,y in cg ]
            
        return self.cumNormGraph_
        
        
    def intersects(self,value,cumulative=False,norm=False,direction=0):
        ''' Calculate x-coordinates for intersection(s) of a graph defined by a list of (x,y) points sorted in x
              with y==value. Uses linear interpolation.
            Arguments:
               value ....... target value
               cumulative .. if true, use cumulative graph
               norm ........ if true, use normalized cumulative graph
               direction ... three possible values: 0 = any intersection, +1/-1 = consider only segments
                             with positive / negative slope.
        '''
        #
        # Get graph and do basic check
        #
        graph = None
        if cumulative:
            graph = self.cumulativeNormGraph() if norm else self.cumulativeGraph()
        else:
            assert norm==False
            graph = self.graph()
            
        result = [ ]
        if len(graph)<2:
            return result
        #
        # loop over adjacent pairs of points
        #
        x1 = None
        y1 = None
        for x2,y2 in graph:
            #
            # start checking at 2nd point
            #
            if x1!=None:
                assert x2>x1
                #
                # value in interval?
                #
                if value>=min(y1,y2) and value<=max(y1,y2):
                    # check if dy is positive or negative, and compare with required sign of direction
                    if direction==0 or direction*(y2-y1)>0:
                        result.append(FitHistogram.interpolate(value,y2,y1,x2,x1))
            #
            # move to next point
            #
            x1 = x2
            y1 = y2

        return result

    def fwhm(self):
        ''' Return lowest / highest x corresponding to ymax/2, and ymax/2
        '''
        #
        # target y value (1/2 maximum)
        #
        y = self.hist.GetBinContent(self.hist.GetMaximumBin())/2.
        #
        # use first intersection with upward slope
        #
        xups = self.intersects(y,cumulative=False,direction=1)
        print("xups",xups)
        xlow = xups[0] if xups else None
        #
        # use last intersection with downward slope
        #
        xdowns = self.intersects(y,cumulative=False,direction=-1)
        print("xdowns",xdowns)
        xhigh = xdowns[-1] if xdowns else None
        #dx = self.hist.GetBinWidth(1)
        ###
        ## preset result (low / high x values)
        ##
        #xlow = None
        #xhigh = None
        ##
        ## loop over pairs of adjacent histogram bins
        ##
        #x1 = self.hist.GetBinLowEdge(1)
        #y1 = self.hist.GetBinContent(1)
        #for i in range(2,self.hist.GetNbinsX()):
        #    # check if target value between contents of neighbouring bins
        #    x2 = self.hist.GetBinLowEdge(i)
        #    y2 = self.hist.GetBinContent(i)
        #    first = None
        #    if y>=y1 and y<y2:
        #        first = True
        #    elif y<=y1 and y>y2:
        #        first = False
        #
        #    if first!=None:
        #        # calculate interpolated x value
        #        x = FitHistogram.interpolate(y,y1,y2,x1,x2)
        #        # store lowest and highest x corresponding to y
        #        if first and xlow==None:
        #            xlow = x
        #        elif not first:
        #            xhigh = x
        #    # move to next bin
        #    x1 = x2
        #    y1 = y2

        #
        # require result ( assumes that first and last bins are < ymax/2 ) and
        # correct for bin width / 2 ( assumes that bin value corresponds to center of bin )
        assert xlow!=None and xhigh!=None
        return (xlow,xhigh,y)
        #return (xlow+dx,xhigh+dx,y)

    def quantile(self,prob):
        ''' Return x-value to quantile q. 
        '''
        result = self.intersects(prob,cumulative=True,norm=True,direction=1)
        #print(prob,result)
        assert len(result)<2

        return result[0] if result else None
        
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
        #print("fhist",self.fhist,self.fhist.hist)
        return self.fhist.hist
            
    def fwhm(self):
        return self.fhist.fwhm()
    
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--moduleType', '-m', help='module type', type=int, choices=[23,24,25], default=23)
parser.add_argument('--fwhm', help='determining FWHM', action='store_true', default=False)
parser.add_argument('--slices', '-s', help='fit in slices in y', action='store_true', default=False)
parser.add_argument('--quantile', '-q', help='calculate quantiles corresponding to +- x sigma (can be repeated)', \
                        action='append', type=float, default=[])
parser.add_argument('--dbgQuantiles', help='modX,dxdz pairs defining individual bins for the quantile determination', \
                        action='append', type=str, default=[ ])
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

dbgQuantPoints = [ ]
for qp in args.dbgQuantiles:
    dbgQuantPoints.append( tuple([float(x.strip()) for x in qp.strip().split(",")]) )

canvases = [ ]
slices = [ ]

mainCnv.Draw()
fitCanvas = FitCanvas(mainCnv,args.moduleType)
fwhmArrow = ROOT.TArrow()
fwhmArrow.SetLineColor(2)
fwhmArrow.SetLineWidth(2)
quantLine = ROOT.TLine()
quantLine.SetLineColor(4)
quantLine.SetLineWidth(2)
if fitCanvas.histogram().GetDimension()==1 and ( not args.slices ):
    assert fitCanvas.histogram().GetDimension()==1
    #print(fitCanvas.fwhm())
    x1,x2,y = fitCanvas.fwhm()
    print("<x> = {:6.1f}um, dx = {:6.1f}um, sig = {:6.1f}um ( interval {:6.4f} - {:6.4f}cm )".format( \
            10000*(x1+x2)/2.,10000*(x2-x1),10000*(x2-x1)/2/sqrt(2*log(2)),x1,x2))

    fitCanvas.pad.cd()
    fwhmArrow.DrawArrow(x1,y,x2,y,0.005,"<>")
    
    quantiles = [ ]
    qlmax = slices[-1].hist.GetMaximum()/10.
    for sigma in args.quantile:
        qs = [ ]
        for sgn in [-1,1]:
            p = ROOT.TMath.Freq(sgn*sigma)
            qs.append(slices[-1].quantile(p))
            #print(sigma,sgn,p,qs[-1])
        quantiles.append(qs)
        if not ( None in qs ):
            quantLine.DrawLine(qs[0],0.,qs[0],qlmax)
            quantLine.DrawLine(qs[1],0.,qs[1],qlmax)

    fitCanvas.pad.Update()

elif fitCanvas.histogram().GetDimension()==2 and args.slices:
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

        quantiles = [ ]
        qlmax = slices[-1].hist.GetMaximum()/10.
        for sigma in args.quantile:
            qs = [ ]
            for sgn in [-1,1]:
                p = ROOT.TMath.Freq(sgn*sigma)
                qs.append(slices[-1].quantile(p))
                #print(sigma,sgn,p,qs[-1])
            quantiles.append(qs)
            if not ( None in qs ):
                quantLine.DrawLine(qs[0],0.,qs[0],qlmax)
                quantLine.DrawLine(qs[1],0.,qs[1],qlmax)
                
        ROOT.gPad.Update()
    
elif fitCanvas.histogram().GetDimension()==3:
    #
    # calculate summary of quantiles in x
    #
    h = fitCanvas.histogram()
    nby = h.GetNbinsY()
    ymin = h.GetYaxis().GetXmin()
    ymax = h.GetYaxis().GetXmax()
    nbz = h.GetNbinsZ()
    zmin = h.GetZaxis().GetXmin()
    zmax = h.GetZaxis().GetXmax()
    print(nby,ymin,ymax,nbz,zmin,zmax)
    #
    # define summary histograms
    #
    summaries = [ ]
    for isigma,sigma in enumerate(args.quantile):
        y1 = ROOT.TMath.Freq(-sigma)
        y2 = ROOT.TMath.Freq(sigma)
        if nby==1:
            hcentre = ROOT.TH1F(h.GetName()+"Qc"+str(isigma), \
                                h.GetTitle()+" (mean from quantiles)".format(y1,y2), \
                                nbz,zmin,zmax)
            hcentre.GetXaxis().SetTitle(h.GetZaxis().GetTitle())
            hcentre.GetYaxis().SetTitle("mean of {:4.1%}/{:4.1%} quantiles".format(y1,y2))
            hhalfwidth = ROOT.TH1F(h.GetName()+"Qhw"+str(isigma), \
                                   h.GetTitle()+" (#sigma from quantiles)".format(y1,y2), \
                                   nbz,zmin,zmax)
            hhalfwidth.GetXaxis().SetTitle(h.GetZaxis().GetTitle())
            hhalfwidth.GetYaxis().SetTitle("#sigma from {:4.1%}/{:4.1%} quantiles".format(y1,y2))
        elif nbz==1:
            hcentre = ROOT.TH1F(h.GetName()+"Qc"+str(isigma), \
                                h.GetTitle()+" (mean from quantiles)".format(y1,y2), \
                                nby,ymin,ymax)
            hcentre.GetXaxis().SetTitle(h.GetYaxis().GetTitle())
            hcentre.GetYaxis().SetTitle("mean of {:4.1%}/{:4.1%} quantiles".format(y1,y2))
            hhalfwidth = ROOT.TH2F(h.GetName()+"Qhw"+str(isigma), \
                                   h.GetTitle()+" (#sigma from quantiles)".format(y1,y2), \
                                   nby,ymin,ymax)
            hhalfwidth.GetXaxis().SetTitle(h.GetYaxis().GetTitle())
            hhalfwidth.GetYaxis().SetTitle("#sigma from {:4.1%}/{:4.1%} quantiles".format(y1,y2))
        else:
            hcentre = ROOT.TH2F(h.GetName()+"Qc"+str(isigma), \
                                h.GetTitle()+" (mean from quantiles)".format(y1,y2), \
                                nby,ymin,ymax,nbz,zmin,zmax)
            hcentre.GetXaxis().SetTitle(h.GetYaxis().GetTitle())
            hcentre.GetYaxis().SetTitle(h.GetZaxis().GetTitle())
            hcentre.GetZaxis().SetTitle("mean of {:4.1%}/{:4.1%} quantiles".format(y1,y2))
            hhalfwidth = ROOT.TH2F(h.GetName()+"Qhw"+str(isigma), \
                                   h.GetTitle()+" (#sigma from quantiles)".format(y1,y2), \
                                   nby,ymin,ymax,nbz,zmin,zmax)
            hhalfwidth.GetXaxis().SetTitle(h.GetYaxis().GetTitle())
            hhalfwidth.GetYaxis().SetTitle(h.GetZaxis().GetTitle())
            hhalfwidth.GetZaxis().SetTitle("#sigma from {:4.1%}/{:4.1%} quantiles".format(y1,y2))
        summaries.append((hcentre,hhalfwidth))

    dbgBins = [ ]
    for y,z in dbgQuantPoints:
        dbgBins.append( (h.GetYaxis().FindBin(y), h.GetZaxis().FindBin(z)) )

    allDbgObjects = { x:[ ] for x in range(len(args.quantile)) }
    for iy in range(nby):
        for iz in range(nbz):
            htmp = h.ProjectionX(iymin=iy+1,iymax=iy+1,izmin=iz+1,izmax=iz+1)
            if htmp.GetSumOfWeights()<100:
                continue
            fhtmp = FitHistogram(htmp)
            for isigma,sigma in enumerate(args.quantile):
                q1 = fhtmp.quantile(ROOT.TMath.Freq(-sigma))
                q2 = fhtmp.quantile(ROOT.TMath.Freq(sigma))
                if q1!=None and q2!=None:
                    hsum = summaries[isigma][0]
                    if nby==1:
                        ibin = hsum.GetBin(iz+1)
                    elif nbz==1:
                        ibin = hsum.GetBin(iy+1)
                    else:
                        ibin = hsum.GetBin(iy+1,iz+1)
                    hsum.SetBinContent(ibin,(q1+q2)/2.)
                    hsum = summaries[isigma][1]
                    hsum.SetBinContent(ibin,(q2-q1)/2./sigma)

                #allDbgObjects[isigma] = [ ]
                if (iy+1,iz+1) in dbgBins:
                    dbgObjects = [ ]
                    hResDbg =htmp.Clone("hResDbg_"+str(isigma)+"_"+str(len(allDbgObjects[isigma])))
                    hResDbg.SetTitle("hResDbg "+str(isigma)+" y"+str(iy+1)+" z"+str(iz+1))
                    hResDbg.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
                    savePad = ROOT.gPad
                    dbgObjects.append(ROOT.TCanvas("c"+hResDbg.GetName()[1:],"c"+hResDbg.GetName()[1:],500,500))
                    dbgObjects[0].SetGridx(1)
                    dbgObjects[0].SetGridy(1)
                    dbgObjects.append(hResDbg)
                    hResDbg.Scale(1./hResDbg.GetMaximum())
                    hResDbg.Draw("hist")
                    if q1!=None and q2!=None:
                        quantLine.DrawLine(q1,0.,q1,hResDbg.GetMaximum()/1.)
                        quantLine.DrawLine(q2,0.,q2,hResDbg.GetMaximum()/1.)
                        g = ROOT.TGraph()
                        g.SetLineColor(ROOT.kMagenta)
                        g.SetLineStyle(2)
                        #g.SetMarkerStyle(20)
                        for x,y in fhtmp.cumulativeNormGraph():
                            g.SetPoint(g.GetN(),x,y)
                        g.Draw("L")
                        dbgObjects.append(g)
                        
                    dbgObjects[0].Update()
                    allDbgObjects[isigma].append(dbgObjects)
                    savePad.cd()

    nsigma = len(args.quantile)
    c = ROOT.TCanvas("cSum","cSum",500*nsigma,1000)
    c.Divide(nsigma,2)
    for i in range(nsigma):
        c.cd(i+1)
        if summaries[i][0].GetDimension()==2:
            ROOT.gPad.SetRightMargin(0.15)
            summaries[i][0].Draw("ZCOL")
        else:
            summaries[i][0].Draw("HIST")
        ROOT.gPad.Update()
        c.cd(nsigma+i+1)
        if summaries[i][1].GetDimension()==2:
            summaries[i][1].Draw("ZCOL")
        else:
            summaries[i][1].Draw("HIST")
        ROOT.gPad.Update()
    c.Update()
    
