#
# efficiency histogram definitions for drawResolution
#
#template : {
#    'canvasName' : '', \
#    'histogramName' : '', \
#    'histogramTitle' : '', \
#    'variable' : '', \
#    'baseCut' : '', \
#    'effCut' : '', \
#    'xTitle' : '', \
#    'xNbins' : 0, \
#    'xMin' : 0, \
#    'xMax' : 0, \
#    'yTitle' : '', \
#    'yMin' : 0, \
#    'yMax' : 0, \
#}

#
# sensors:
#   2S: 2 sensors; sensor 10 x 10cm2, pitch 90um, length 5cm, 2x2016 strips / sensor
#   PS: 1 sensor, 5 x 10 cm2, pitch 100um
#   PS-s: length 2.5cm, 2 x 960 strips / sensor
#   PS-p: length 1.5mm, 32 x 960 strips / sensor
#
resX = {
    #'canvasName' : 'cEffX1', \
    #'histogramName' : 'effX1', \
    'histogramTitle' : 'residuals (x)', \
    'variable' : 'localPos.x()-rhLocalPos.x()', \
    'baseCuts' : 'tof<12.5&&hasRecHit>0', \
    'xTitle' : '#Delta x [cm]', \
    'xNbins' : 300, \
    'xMin' : -0.075, \
    'xMax' : 0.075, \
    'logY' : True
#    #'yTitle' : 'efficiency', \
#    #'yMin' : 0.8, \
#    #'yMax' : 1.05
}

resY = {
    #'canvasName' : 'cEffX1', \
    #'histogramName' : 'effX1', \
    'histogramTitle' : 'residuals (y)', \
    'variable' : 'localPos.y()-rhLocalPos.y()', \
    'baseCuts' : 'tof<12.5&&hasRecHit>0', \
    'xTitle' : '#Delta y [cm]', \
    'xNbins' : 200, \
    'xMin' : -5., \
    'xMax' : 5., \
    'mType23' : {
        'xNbins' : 250, \
        'xMin' : -0.25, \
        'xMax' : 0.25, \
        }, \
    'mType24' : {
        'xNbins' : 250, \
        'xMin' : -2.5, \
        'xMax' : 2.5, \
        }, \
    'mType25' : {
        'xNbins' : 300, \
        'xMin' : -4.5, \
        'xMax' : 4.5, \
        }, \
    'logY' : True
#    #'yTitle' : 'efficiency', \
#    #'yMin' : 0.8, \
#    #'yMax' : 1.05
}

pullX = {
    'histogramTitle' : 'pulls (x)', \
    'variable' : '(localPos.x()-rhLocalPos.x())/rhLocalErr.x()', \
    'baseCuts' : 'tof<12.5&&hasRecHit>0', \
    'xTitle' : 'pull local x', \
    'xNbins' : 200, \
    'xMin' : -5, \
    'xMax' : 5, \
    'logY' : True
#    #'yTitle' : 'efficiency', \
#    #'yMin' : 0.8, \
#    #'yMax' : 1.05
}

resY = {
    'histogramTitle' : 'pulls (y)', \
    'variable' : '(localPos.y()-rhLocalPos.y())/rhLocalErr.y()', \
    'baseCuts' : 'tof<12.5&&hasRecHit>0', \
    'xTitle' : 'pull local y', \
    'xNbins' : 200, \
    'xMin' : -5, \
    'xMax' : 5, \
    'logY' : True
#    #'yTitle' : 'efficiency', \
#    #'yMin' : 0.8, \
#    #'yMax' : 1.05
}


stdRes2D = {
    #'canvasName' : 'cEffX1', \
    #'histogramName' : 'effX1', \
    'histogramTitle' : 'my resolution', \
    'variable' : 'abs(path.x()):localPos.x()-rhLocalPos.x()', \
    'baseCuts' : 'tof<12.5&&hasRecHit>0', \
    #'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'residual x [cm]', \
    'xNbins' : 200, \
    'xMin' : -0.1, \
    'xMax' : 0.1, \
    'yTitle' : 'dx [cm]', \
    'yNbins' : 20,
    'yMin' : 0., \
    'yMax' : 0.1
}

widthVsPath = {
    #'canvasName' : 'cEffX1', \
    #'histogramName' : 'effX1', \
    'histogramTitle' : 'width vs path', \
    'variable' : 'clusterSize:abs(path.x())', \
    'baseCuts' : 'tof<12.5&&hasRecHit>0', \
    #'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'SimHit path (x)', \
    'xNbins' : 100, \
    'xMin' : 0, \
    'xMax' : 0.2, \
    'yTitle' : 'clustersize', \
    'yNbins' : 25,
    'yMin' : 0., \
    'yMax' : 25,
    'profile' : True
}

widthVsModX90 = {
    'histogramTitle' : 'width vs x modula 90 #mu m', \
    'variable' : 'clusterSize:floatMod(10000*localPos.x(),90)', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'xTitle' : 'x modulo 90 #mu m [#mu m]', \
    'xNbins' : 100, \
    'xMin' : -5, \
    'xMax' : 95, \
    'yTitle' : 'cluster width', \
#    'yNbins' : 10,
    'yMin' : 0, \
    'yMax' : 25, \
    'mType23' : { 'display' : False },
    'mType24' : { 'display' : False },
    'profile' : True
}

widthVsMod2DX90 = {
    'histogramTitle' : 'width vs x modula 90 #mu m', \
    'variable' : 'clusterSize:floatMod(10000*localPos.x(),90)', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'xTitle' : 'x modulo 90 #mu m [#mu m]', \
    'xNbins' : 100, \
    'xMin' : -5, \
    'xMax' : 95, \
    'yTitle' : 'cluster width', \
    'yNbins' : 10, \
    'yMin' : 0, \
    'yMax' : 10, \
    'mType23' : { 'display' : False },
    'mType24' : { 'display' : False }
}

widthVsModX100 = {
    'histogramTitle' : 'width vs x modula 100 #mu m', \
    'variable' : 'clusterSize:floatMod(10000*localPos.x(),100)', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'xTitle' : 'x modulo 100 #mu m [#mu m]', \
    'xNbins' : 110, \
    'xMin' : -5, \
    'xMax' : 105, \
    'yTitle' : 'cluster width', \
#    'yNbins' : 10,
    'yMin' : 0, \
    'yMax' : 25, \
    'mType25' : { 'display' : False }, \
    'profile' : True
}

widthVsMod2DX100 = {
    'histogramTitle' : 'width vs x modula 100 #mu m', \
    'variable' : 'clusterSize:floatMod(10000*localPos.x(),100)', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'xTitle' : 'x modulo 100 #mu m [#mu m]', \
    'xNbins' : 110, \
    'xMin' : -5, \
    'xMax' : 105, \
    'yTitle' : 'cluster width', \
    'yNbins' : 10, \
    'yMin' : 0, \
    'yMax' : 10, \
    'mType25' : { 'display' : False }
}

widthVsModY015 =  {
    'histogramTitle' : 'Width vs. y % 1.5 mm', \
    'variable' : 'clusterSize:floatMod(100*localPos.y(),15)/10.', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'xTitle' : 'x modulo 1.5mm [mm]', \
    #'xNbins' : 340, \
    #'xMin' : -0.1, \
    #'xMax' : 1.6, \
    'xNbins' : 200, \
    'xMin' : -0.25, \
    'xMax' : 1.75, \
    'yTitle' : 'cluster width', \
#    'yNbins' : 10,
    'yMin' : 0, \
    'yMax' : 25, \
    'mType24' : { 'display' : False }, \
    'mType25' : { 'display' : False }, \
    'profile' : True
}

widthVsModY250 =  {
    'histogramTitle' : 'Width vs. y % 2.5 cm', \
    'variable' : 'clusterSize:floatMod(10*localPos.y(),25)/10.', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'xTitle' : 'x modulo 2.5cm [cm]', \
    'xNbins' : 300, \
    'xMin' : -0.25, \
    'xMax' : 2.75, \
    'yTitle' : 'cluster width', \
#    'yNbins' : 10,
    'yMin' : 0, \
    'yMax' : 25, \
    'mType23' : { 'display' : False }, \
    'mType25' : { 'display' : False }, \
    'profile' : True
}


widthVsModY500 =  {
    'histogramTitle' : 'Width vs. y % 5 cm', \
    'variable' : 'clusterSize:floatMod(10*localPos.y(),50)/10.', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'xTitle' : 'x modulo 5cm [cm]', \
    'xNbins' : 275, \
    'xMin' : -0.25, \
    'xMax' : 5.25, \
    'yTitle' : 'cluster width', \
#    'yNbins' : 10,
    'yMin' : 0, \
    'yMax' : 25, \
    'mType23' : { 'display' : False }, \
    'mType24' : { 'display' : False }, \
    'profile' : True
}


pathX = {
    'histogramTitle' : 'path (x)', \
    'variable' : '10000*abs(path.x())', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    #'baseCuts' : 'tof<12.5&&hasRecHit>0', \
    #'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'SimHit path (x) [#mum]', \
    'xNbins' : 200, \
    'xMin' : 0, \
    'xMax' : 1000
}

pathY = {
    'histogramTitle' : 'path (y)', \
    'variable' : '10000*abs(path.y())', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    #'baseCuts' : 'tof<12.5&&hasRecHit>0', \
    #'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'SimHit path (y [#mum])', \
    'xNbins' : 200, \
    'xMin' : 0, \
    'xMax' : 2000
}

effX = {
    #'canvasName' : 'cEffX1', \
    #'histogramName' : 'effX1', \
    'histogramTitle' : 'Efficiency 1D', \
    'variable' : 'localPos.x()', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'local x [cm]', \
    'mType25' : {
        'xNbins' : 250, \
        'xMin' : -5, \
        'xMax' : 5, \
    }, \
    'mType24' : { \
        'xNbins' : 250, \
        'xMin' : -5.0, \
        'xMax' : 5.0, \
    }, \
    'mType23' : { \
        'xNbins' : 550, \
        'xMin' : -5.5, \
        'xMax' : 5.5, \
    }, \
    'yTitle' : 'efficiency', \
    'yMin' : 0.8, \
    'yMax' : 1.05
}
    
    
effY = {
    #'canvasName' : 'cEffX1', \
    #'histogramName' : 'effX1', \
    'histogramTitle' : 'Efficiency Y 1D', \
    'variable' : 'localPos.y()', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'local y [cm]', \
    'xNbins' : 550, \
    'xMin' : -5.5, \
    'xMax' : 5.5, \
    'yTitle' : 'efficiency', \
    'yMin' : 0.8, \
    'yMax' : 1.05
}

effAlphaLoose = {
    #'canvasName' : 'cEffX1', \
    #'histogramName' : 'effX1', \
    'histogramTitle' : 'Efficiency vs. alpha(xz) tight', \
    'variable' : 'atan2(abs(path.x()),abs(path.z()))/3.1415*180', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'alpha(xz) [deg]', \
    'xNbins' : 100, \
    'xMin' : 0, \
    'xMax' : 100., \
    'yTitle' : 'efficiency', \
    'yMin' : 0.5, \
    'yMax' : 1.05
}

effAlphaTight = {
    #'canvasName' : 'cEffX1', \
    #'histogramName' : 'effX1', \
    'histogramTitle' : 'Efficiency vs. alpha(xz) loose', \
    'variable' : 'atan2(abs(path.x()),abs(path.z()))/3.1415*180', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.1', \
    'xTitle' : 'alpha(xz) [deg]', \
    'xNbins' : 100, \
    'xMin' : 0, \
    'xMax' : 100., \
    'yTitle' : 'efficiency', \
    'yMin' : 0.5, \
    'yMax' : 1.05
}

#effDx1 = {
#    #'canvasName' : 'cEffX1', \
#    #'histogramName' : 'effX1', \
#    'histogramTitle' : 'Efficiency vs. dx tight', \
#    'variable' : 'abs(path.x())', \
#    'baseCuts' : 'tof<12.5&&pabs>0.3', \
#    'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
#    'xTitle' : 'dx [cm]', \
#    'xNbins' : 100, \
#    'xMin' : 0, \
#    'xMax' : 0.1, \
#    'yTitle' : 'efficiency', \
#    'yMin' : 0.5, \
#    'yMax' : 1.05
#}

#effDx2 = {
#    #'canvasName' : 'cEffX1', \
#    #'histogramName' : 'effX1', \
#    'histogramTitle' : 'Efficiency vs. dx loose', \
#    'variable' : 'abs(path.x())', \
#    'baseCuts' : 'tof<12.5&&pabs>0.3', \
#    'effCuts' : 'hasRecHit>0', \
#    'xTitle' : 'dx [cm]', \
#    'xNbins' : 100, \
#    'xMin' : 0, \
#    'xMax' : 0.1, \
#    'yTitle' : 'efficiency', \
#    'yMin' : 0.5, \
#    'yMax' : 1.05
#}

#widthVsModX = {
#    'histogramTitle' : 'Cluster width vs. x%90um', \
#    'variable' : 'clusterSize:(10000*(localPos.x()+15))-90*int(10000*(localPos.x()+15)/90.)', \
#    'baseCuts' : 'tof<12.5&&pabs>0.3&abs(path.x()/path.z())<0.1', \
##    'effCuts' : 'hasRecHit>0', \
#    'xTitle' : 'x mod 90 #mu m [#mu m]', \
#    'xNbins' : 120, \
#    'xMin' : -10, \
#    'xMax' : 110, \
#    'yTitle' : 'cluster width', \
#    'yNbins' : 10,
#    'yMin' : -0.5, \
#    'yMax' : 9.5
#    }

effVsModX100 =  {
    'histogramTitle' : 'Efficiency vs. x % 100 #mu m', \
    'variable' : 'floatMod(10000*localPos.x(),100)', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'x modulo 100 #mu m [#mu m]', \
    'xNbins' : 220, \
    'xMin' : -5, \
    'xMax' : 105, \
#    'yTitle' : 'cluster width', \
#    'yNbins' : 10,
    'yMin' : 0, \
    'yMax' : 1.05,
    'mType25' : { 'display' : False }
}

effVsModX090 =  {
    'histogramTitle' : 'Efficiency vs. x % 90 #mu m', \
    'variable' : 'floatMod(10000*localPos.x(),90)', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'x modulo 90 #mu m [#mu m]', \
    'xNbins' : 200, \
    'xMin' : -5, \
    'xMax' : 95, \
#    'yTitle' : 'cluster width', \
#    'yNbins' : 10,
    'yMin' : 0, \
    'yMax' : 1.05,
    'mType23' : { 'display' : False },
    'mType24' : { 'display' : False }
}
    
effVsModY015 =  {
    'histogramTitle' : 'Efficiency vs. y % 1.5 mm', \
    'variable' : 'floatMod(100*localPos.y(),15)/10.', \
    #'variable' : 'floatMod(100*localPos.y(),15)', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'x modulo 1.5mm [mm]', \
    #'xNbins' : 340, \
    #'xMin' : -0.1, \
    #'xMax' : 1.6, \
    'xNbins' : 400, \
    #'xNbins' : 200, \
    'xMin' : -0.25, \
    'xMax' : 1.75, \
#    'yTitle' : 'cluster width', \
#    'yNbins' : 10,
    'yMin' : 0, \
    'yMax' : 1.05,
    'mType24' : { 'display' : False },
    'mType25' : { 'display' : False }
}

effVsModY250 =  {
    'histogramTitle' : 'Efficiency vs. y % 2.5 cm', \
    'variable' : 'floatMod(10*localPos.y(),25)/10.', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'x modulo 2.5cm [cm]', \
    'xNbins' : 600, \
    #'xNbins' : 300, \
    'xMin' : -0.25, \
    'xMax' : 2.75, \
#    'yTitle' : 'cluster width', \
#    'yNbins' : 10,
    'yMin' : 0, \
    'yMax' : 1.05,
    'mType23' : { 'display' : False },
    'mType25' : { 'display' : False }
}


effVsModY500 =  {
    'histogramTitle' : 'Efficiency vs. y % 5 cm', \
    'variable' : 'floatMod(10*localPos.y(),50)/10.', \
    'baseCuts' : 'tof<12.5&&pabs>0.3', \
    'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'x modulo 5cm [cm]', \
    'xNbins' : 550, \
    #'xNbins' : 275, \
    'xMin' : -0.25, \
    'xMax' : 5.25, \
#    'yTitle' : 'cluster width', \
#    'yNbins' : 10,
    'yMin' : 0, \
    'yMax' : 1.05,
    'mType23' : { 'display' : False },
    'mType24' : { 'display' : False }
}

