import sys
from math import pi,atan2,sin,cos,sqrt
import ROOT

class HitPos4D:

    def __init__(self,t,x,y,z):
        self.t = t
        self.x = x
        self.y = y
        self.z = z

def interpolate(x,xyzs):
    ''' assume xyzs sorted in x
    '''
    for i in range(len(xyzs)-1):
        x1,y1,z1 = xyzs[i]
        x2,y2,z2 = xyzs[i+1]
        if (x-x1)*(x-x2)<=0.:
            return (y1 + (x-x1)/(x2-x1)*(y2-y1)), (z1 + (x-x1)/(x2-x1)*(z2-z1))
    # outside range
    minxyz = ( 999999, 999999, 999999 )
    maxxyz = ( -999999, -999999, -999999 )
    for i in range(len(xyzs)):
        if xyzs[i][0]<minxyz[0]:
            minxyz = xysz[i]
        if xyzs[i][0]>maxxyz[0]:
            maxxyz = xysz[i]
        
    #xs = [ v[0] for v in xysz ]
    #ys = [ v[1] for v in xysz ]
    #xs = [ v[2] for v in xysz ]
    if x1<sum([ v[0] for v in xyzs ])/len(xyzs):
        return minxyz[0],minxyz[2]
    else:
        return maxxyz[0],maxxyz[2]
        

rootObjects = [ ]
def drawHits(hitPositions,direction,inTimeWindow):
    cnvXY.cd()
    f = cnvXY.DrawFrame(-120.,-120.,120.,120.)
    cnvRZ.cd()
    f = cnvRZ.DrawFrame(0.,0.,350.,120.)

    cnvXY.cd()
    rootObjects.append(f)
    e = ROOT.TEllipse()
    e.SetLineColor(4)
    e.SetFillStyle(0)
    rootObjects.append(e)
    for r in [23.5,36.5,51.5]:
        e.DrawEllipse(0.,0.,r,r,0.,360.,0.)
    cnvXY.SetGridx(1)
    cnvXY.SetGridy(1)
    mout1 = ROOT.TMarker()
    mout1.SetMarkerStyle(20)
    mout1.SetMarkerColor(3)
    mout2 = ROOT.TMarker()
    mout2.SetMarkerStyle(24)
    mout2.SetMarkerColor(3)
    minw1 = ROOT.TMarker()
    minw1.SetMarkerStyle(20)
    minw1.SetMarkerColor(2)
    minw2 = ROOT.TMarker()
    minw2.SetMarkerStyle(24)
    minw2.SetMarkerColor(2)
    rootObjects.append(mout1)
    rootObjects.append(mout2)
    rootObjects.append(minw1)
    rootObjects.append(minw2)
    for i in range(len(hitPositions)):
        if inTimeWindow[i]:
            if direction[i]:
                marker = mout1
                #mout1.DrawMarker(hitPositions[i].x,hitPositions[i].y)
            else:
                marker = minw1
                #minw1.DrawMarker(hitPositions[i].x,hitPositions[i].y)
        else:
            if direction[i]:
                marker = mout2
                #mout2.DrawMarker(hitPositions[i].x,hitPositions[i].y)
            else:
                marker = minw2
                #minw2.DrawMarker(hitPositions[i].x,hitPositions[i].y)
        cnvXY.cd()
        marker.DrawMarker(hitPositions[i].x,hitPositions[i].y)
        cnvRZ.cd()
        marker.DrawMarker(abs(hitPositions[i].z),sqrt(hitPositions[i].x**2+hitPositions[i].y**2))

    xc = sum([ v.x for v in hitPositions ])/len(hitPositions)
    yc = sum([ v.y for v in hitPositions ])/len(hitPositions)
    print("center",xc,yc)
    gtxyzs = [ ]
    for i in range(len(hitPositions)):
        phi = atan2(hitPositions[i].y-yc,hitPositions[i].x-xc)
        #if phis and sgn==None:
        #    sgn = 1 if phi>phis[-1] else -1
        #if phis and sgn*(phi-phis[-1])<0:
        #    phi += sgn*2*pi
        gtxyzs.append((hitPositions[i].t,phi,sqrt((hitPositions[i].x-xc)**2+(hitPositions[i].y-yc)**2), \
                       hitPositions[i].z))
    gtxyzs.sort()

    phis = [ ]
    sgn = None
    gxyzs = [ ]
    for i,txyz in enumerate(gtxyzs):
        t,x,y,z = txyz
        #print(i,txy,y*cos(x)+xc,y*sin(x)+yc)
        if gxyzs and sgn==None:
            sgn = 1 if x>gxyzs[-1][0] else -1
            #print("Setting sgn to",sgn,"at point",i+1,y*cos(x)+xc,y*sin(x)+yc)
        if gxyzs:
            #print(sgn*(x-gxys[-1][0]))
            if sgn*(x-gxyzs[-1][0])<(-pi):
                #print("Correcting phi by 2pi from",x,"to",x+sgn*2*pi)
                x += sgn*2*pi
            if sgn*(x-gxyzs[-1][0])<0:
                #print("Starting new loop at point",i+1,sgn,x,gxys[-1][0])
                x += sgn*2*pi
        gxyzs.append((x,y,z))
        phis.append(x)

    #g = ROOT.TGraph()
    #rootObjects.append(g)
    #for i,gxy in enumerate(gxys):
    #    g.SetPoint(i,gxy[0],gxy[1])
    #spline = ROOT.TSpline3("myspline",g)
    #rootObjects.append(spline)
    phiMin = min(phis)
    phiMax = max(phis)
    #print(phiMin,phiMax)
    polyLineXY = ROOT.TPolyLine()
    polyLineRZ = ROOT.TPolyLine()
    phi = phiMin
    while phi<phiMax:
        #r = spline.Eval(phi)
        r,zglob = interpolate(phi,gxyzs)
        xglob = xc+cos(phi)*r
        yglob = yc+sin(phi)*r
        polyLineXY.SetNextPoint(xglob,yglob)
        polyLineRZ.SetNextPoint(abs(zglob),sqrt(xglob**2+yglob**2))
        #print("-",phi,r,xc+cos(phi)*r,yc+sin(phi)*r)
        phi += 0.1
    rootObjects.append(polyLineXY)
    cnvXY.cd()
    polyLineXY.Draw()
    cnvXY.Update()
    cnvRZ.cd()
    polyLineRZ.Draw()
    cnvRZ.Update()

    
    input("Enter")

tf = ROOT.TFile(sys.argv[1])
analysis = tf.Get("analysis")
simHitTree = analysis.Get("SimHitTree")
simTrackTree = analysis.Get("SimTrackTree")
print(simHitTree.GetEntries())
#simHitTree.Print()
print(simTrackTree.GetEntries())
assert simHitTree.GetEntries()==simTrackTree.GetEntries()
#simTrackTree.Print()
#cnv = ROOT.TCanvas()
#simTrackTree.Draw("SimTrack_type","abs(SimTrack_type)<200")
#cnv.Update()

tMin = 0.
tMax = 999999.
if len(sys.argv)==4:
    tMin = float(sys.argv[2])
    tMax = float(sys.argv[3])
    
cnvXY = ROOT.TCanvas("cxy","cxy",1000,1000)
cnvRZ = ROOT.TCanvas("crz","crz",1200,400)

for evHits,evTracks in zip(simHitTree,simTrackTree):
    #
    # look for central pion track with pt>0.5GeV
    #
    itkSel = None
    tkId = None
    for itk in range(len(evTracks.SimTrack_trackId)):
        if evTracks.SimTrack_momentum[itk].Rho()<0.5 or abs(evTracks.SimTrack_momentum[itk].Eta())>1.5 or \
          abs(evTracks.SimTrack_type[itk])!=211:
            continue
        itkSel = itk
        tkId = evTracks.SimTrack_trackId[itk]
        break
    if tkId==None:
        continue

    print("Found track:")
    print("  pid, pt, eta, phi = {:3d}, {:5.1f}, {:5.2f}, {:5.2f}".format(evTracks.SimTrack_type[itkSel], \
                evTracks.SimTrack_momentum[itkSel].Rho(), \
                evTracks.SimTrack_momentum[itkSel].Eta(), \
                evTracks.SimTrack_momentum[itkSel].Phi()))
              
    ishSel = [ ]
    for ish in range(len(evHits.pabs)):
        if evHits.trackId[ish]==tkId and evHits.particleType[ish]==evTracks.SimTrack_type[itkSel]:
            ishSel.append(ish)
    ishSel.sort(key=lambda x:evHits.tof[x])

    print("Found",len(ishSel),"SimHits for this track")
    print("   {:>5s}, {:>6s}, {:>5s}, {:>6s}, {:>5s}, {:>5s} , {:<5s} , {:<5s}".format("R","z","tof", \
            "pt","eta","phi","out","timeW"))

    nshInTime = sum([ (evHits.tof[i]>tMin and evHits.tof[i]<tMax) for i in ishSel ])
    print("Found",nshInTime,"hits in time window")
    if nshInTime==0:
        continue
    
    for i in ishSel:
        print("   {:5.1f}, {:6.1f}, {:5.2f}, {:6.2f}, {:5.2f}, {:5.2f}".format( \
                evHits.globalPos[i].Rho(),evHits.globalPos[i].z(),evHits.tof[i], \
                evHits.pabs[i]*evHits.globalDir[i].Rho(),evHits.pabs[i]*evHits.globalDir[i].Eta(), \
                evHits.pabs[i]*evHits.globalDir[i].Phi()),",", \
                (evHits.globalDir[i].x()*evHits.globalPos[i].x()+evHits.globalDir[i].y()*evHits.globalPos[i].y())>0,
                ",",(evHits.tof[i]>tMin and evHits.tof[i]<tMax))
    if len(ishSel)>0:
        hitPositions = [ ]
        directions = [ ]
        timeWindows = [ ]
        for i in ishSel:
            hitPositions.append(HitPos4D(evHits.tof[i],evHits.globalPos[i].x(), \
                                       evHits.globalPos[i].y(),evHits.globalPos[i].z()))
            directions.append((evHits.globalDir[i].x()*evHits.globalPos[i].x()+ \
                               evHits.globalDir[i].y()*evHits.globalPos[i].y())>0)
            timeWindows.append((evHits.tof[i]>tMin and evHits.tof[i]<tMax))
        drawHits(hitPositions,directions,timeWindows)
        #drawHits([ evHits.tof[i] for i in ishSel ], \
        #         [ evHits.globalPos[i].x() for i in ishSel ], [ evHits.globalPos[i].y() for i in ishSel ], \
        #         [ (evHits.globalDir[i].x()*evHits.globalPos[i].x()+evHits.globalDir[i].y()*evHits.globalPos[i].y())>0 \
        #               for i in ishSel ], \
        #         [ ( evHits.tof[i]>tMin and evHits.tof[i]<tMax ) for i in ishSel ]
        #             )
        print([ ( evHits.tof[i]>tMin and evHits.tof[i]<tMax ) for i in ishSel ])
        #continue
