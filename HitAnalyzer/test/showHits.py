import sys,os
import re
import ROOT
import argparse

    
class SimHit:
    ''' Debug information for a SimHit
    '''
    def __init__(self,trackId,pabs,loss,pdgId,x1,x2,y1,y2,closest):
        self.trackId = trackId
        self.pabs = pabs
        self.energyLoss = loss
        self.pdgId = pdgId
        self.entry = ( x1, y1 )
        self.exit = ( x2, y2 )
        self.closest = closest

    def __str__(self):
        line = "SimHit: track id "
        if self.trackId!=None:
            line += "{:4d}".format(self.trackId)
        else:
            line += "{:4s}".format("None")
        line += ", pabs {:8.4f}, eloss {:8.6f}".format(self.pabs,self.energyLoss)
        line += ", pid {:6d}".format(self.pdgId)
        line += ", x/y {:8.4f}/{:8.4f}".format((self.entry[0]+self.exit[0])/2.,(self.entry[1]+self.exit[1])/2.)
        if self.closest:
            line += " (closest)"
        return line
    
class RecHit:
    ''' Debug information for a RecHit
    '''
    def __init__(self,detId,layer,mType):
        self.detId = detId
        self.layer = layer
        self.mType = mType
        self.rhPos = None
        self.rhErr = None
        self.detDims = None
        self.clusterSize = None
        self.trackIds = None
        self.simHits = [ ]

    def setPosition(self,x,y):
        self.rhPos = ( x, y )

    def setError(self,ex,ey,rho):
        self.rhErr = ( ex, ey, rho )

    def __str__(self):
        line = "RecHit: raw id {:10d}, layer {:4d}, mType {:2d}, clusterSize {:2d}".format(self.detId, \
                    self.layer,self.mType,self.clusterSize)
        line += ", trackIds "+",".join([ str(x) for x in self.trackIds]) + "\n"
        line += "   pos(x/y) {:8.4f}/{:8.4f}".format(self.rhPos[0],self.rhPos[1])
        line += ",   err(x/y) {:8.4f}/{:8.4f}".format(self.rhErr[0],self.rhErr[1])
        return line

def updateRange(range,x,y):
    ''' update a range (xmin,xmax,ymin,ymax) for a new point (x,y)
    '''
    if range[0]==None or range[0]>x:
        range[0] = x
    if range[1]==None or range[1]<x:
        range[1] = x
    if range[2]==None or range[2]>y:
        range[2] = y
    if range[3]==None or range[3]<y:
        range[3] = y

    return range
        
if __name__=="__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--allHits', '-a', help='show all hits on a module', action='store_true', default=False)
    parser.add_argument('--fullSize', '-f', help='do not zoom in on hits', action='store_true', default=False)
    parser.add_argument('--output', '-o', help='directory for storing canvases', type=str, default=None)
    parser.add_argument('--maxFrames', help='max. number of frames shown', type=int, default=None)
    parser.add_argument('file', help='input log file', type=str, nargs=1, default=None)
    args = parser.parse_args()
    if args.output!=None:
        assert args.maxFrames!=None and args.maxFrames>0
        args.output = os.path.expanduser(args.output)
        assert os.path.isdir(args.output)
    #
    # loop over input file and store information related to each RecHit
    #
    recHits = [ ]
    with open(args.file[0]) as txtFile:
        rh = None
        for l in txtFile:
            line = l[:-1].strip()
            fields = line.replace(","," ").split()
            if line.startswith("RecHit"):
                #
                # Start new RecHit with a "RecHit ..." line
                #
                if rh!=None:
                    recHits.append(rh)
                assert fields[2]=="rawid" and fields[4]=="layer" and fields[6]=="mType"
                rh = RecHit(int(fields[3]),int(fields[5]),int(fields[7]))
            elif rh!=None:
                #
                # complete info for current RecHit
                #
                if fields[0]=="localPos":
                    assert fields[2]=="/"
                    rh.setPosition(float(fields[1]),float(fields[3]))
                elif fields[0]=="localErr":
                    assert fields[2]=="/" and fields[4]=="/"
                    rh.setError(float(fields[1]),float(fields[3]),float(fields[5]))
                elif fields[0]=="width":
                    assert fields[2]=="length" and fields[4]=="thickness"
                    rh.detDims = ( float(fields[1]),float(fields[3]),float(fields[5]) )
                elif line.startswith("cluster size"):
                    rh.clusterSize = int(fields[3])
                elif line.startswith("track ids:"):
                    rh.trackIds = sorted([ int(x) for x in fields[2:] ])
                elif line.startswith("SimHit"):
                    closest = "*" in fields[0]
                    if fields[1]=="track":
                        trackId = int(fields[3])
                        fields = fields[4:]
                    else:
                        trackId = None
                        fields = fields[1:]
                    rh.simHits.append(SimHit(trackId,float(fields[2]),float(fields[4]),int(fields[6]), \
                                             float(fields[9]),float(fields[11]),float(fields[14]),float(fields[16]), \
                                             closest))
                else:
                    #
                    # no more recognised line: store RecHit
                    #
                    recHits.append(rh)
                    rh = None
    #
    # after EOF, store last open RecHit (if any)
    #
    if rh!=None:
        recHits.append(rh)
#
# ROOT objects and styles for drawing
#
cnv = ROOT.TCanvas("c","c",1000,800)
rhMarker = ROOT.TMarker()
rhMarker.SetMarkerStyle(29)
rhMarker.SetMarkerColor(2)
rhEllipse = ROOT.TEllipse()
rhEllipse.SetLineColor(2)
rhEllipse.SetFillStyle(0)
rhEllipse.SetLineWidth(2)
shMarker = ROOT.TMarker()
shMarker.SetMarkerStyle(20)
shMarker.SetMarkerColor(1)
shMarker.SetMarkerSize(0.5)
shcMarker = ROOT.TMarker()
shcMarker.SetMarkerStyle(20)
shcMarker.SetMarkerColor(4)
shcMarker.SetMarkerSize(1.)
shLine = ROOT.TLine()
shLine.SetLineColor(1)
shLine.SetLineWidth(3)
shcLine = ROOT.TLine()
shcLine.SetLineColor(4)
shcLine.SetLineWidth(3)
#
# loop over RecHits
#
detId = None
frame = None
closestSHs = [ ]
cnvIdx = 0
for rh in recHits:
    #print("DetId",rh.detId)
    if detId==None or (not args.allHits ) or detId!=rh.detId:
        #
        # redraw all SimHits with "closest" attribute in current frame
        #print(len(closestSHs),"closest SimHits")
        for sh in closestSHs:
            shcLine.DrawLine(sh.entry[0],sh.entry[1],sh.exit[0],sh.exit[1])
            shcMarker.DrawMarker(*sh.entry)
            shcMarker.DrawMarker(*sh.exit)
        if closestSHs:
            cnv.Update()
        #
        # if a frame has been displayed: expect input before continuing
        #
        if detId!=None and args.output==None:
            try:
                input("Next? ")
            except EOFError:
                # end RecHit loop on Ctrl-D
                break
        #
        # delete existing frame
        #
        if frame!=None:
            if args.output!=None:
                fnout = "RecHit"
                if args.allHits:
                    fnout += "s"
                fnout += "_{:04d}_l{:04d}_m{:02d}.png".format(cnvIdx+1,rh.layer,rh.mType)
                cnv.SaveAs(os.path.join(args.output,fnout))
            cnvIdx += 1
            if args.maxFrames!=None and cnvIdx>=args.maxFrames:
                break
            del frame
            
        #
        # print RecHit info and draw frame with module dimensions
        #
        print(rh)
        width = rh.detDims[0]
        length = rh.detDims[1]
        frame = cnv.DrawFrame(-width/2.,-length/2.,width/2.,length/2.)
        frame.SetTitle("Layer {:d} ModuleType {:d} DetId {:d}".format(rh.layer,rh.mType,rh.detId))
        frange = [ None, None, None, None ]
        closestSHs = [ ]
    elif args.allHits:
        #
        # print next RecHit in case all hits / module where requested
        #
        print(rh)
    #
    # draw RecHit with error ellipse
    #
    detId = rh.detId
    rhMarker.DrawMarker(*rh.rhPos)
    # for the time being, ignore correlation (if present)
    rhEllipse.DrawEllipse(rh.rhPos[0],rh.rhPos[1],rh.rhErr[0],rh.rhErr[1],0,360,0)
    frange = updateRange(frange,rh.rhPos[0]+rh.rhErr[0],rh.rhPos[1]+rh.rhErr[1])
    frange = updateRange(frange,rh.rhPos[0]+rh.rhErr[0],rh.rhPos[1]-rh.rhErr[1])
    frange = updateRange(frange,rh.rhPos[0]-rh.rhErr[0],rh.rhPos[1]+rh.rhErr[1])
    frange = updateRange(frange,rh.rhPos[0]-rh.rhErr[0],rh.rhPos[1]+rh.rhErr[1])
    #
    # loop over all associated SimHits (a single SimHit can show up several times
    #   in case all RecHits / module where requested)
    #
    hasClosest = False
    for sh in rh.simHits:
        #
        # print and draw hit
        #
        print(sh)
        frange = updateRange(frange,sh.entry[0],sh.entry[1])
        frange = updateRange(frange,sh.exit[0],sh.exit[1])
        if sh.closest:
            # "closest hit" style
            hasClosest = True
            shcLine.DrawLine(sh.entry[0],sh.entry[1],sh.exit[0],sh.exit[1])
            shcMarker.DrawMarker(*sh.entry)
            shcMarker.DrawMarker(*sh.exit)
            closestSHs.append(sh)
        else:
            # standard style
            shLine.DrawLine(sh.entry[0],sh.entry[1],sh.exit[0],sh.exit[1])
            shMarker.DrawMarker(*sh.entry)
            shMarker.DrawMarker(*sh.exit)

    if rh.simHits and ( not hasClosest ):
        print("No closest hit")

    #
    # readjust axis range unless full module size was requested
    #
    if not args.fullSize:
        fdx = frange[1] - frange[0]
        frame.GetXaxis().SetRangeUser(frange[0]-0.1*fdx,frange[1]+0.1*fdx)
        fdy = frange[3] - frange[2]
        frame.GetYaxis().SetRangeUser(frange[2]-0.1*fdy,frange[3]+0.1*fdy)
    frame.GetXaxis().SetTitle("x [cm]")
    frame.GetYaxis().SetTitle("y [cm]")
    #
    # update canvas
    #
    cnv.SetGridx(1)
    cnv.SetGridy(1)
    cnv.Update()
    print()
    #break

    
    
    
