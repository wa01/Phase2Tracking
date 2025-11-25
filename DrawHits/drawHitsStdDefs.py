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
# sensor dimensions:
#
# 2S: 2 modules with 1016 strips at 90 μm x 5 cm
# PS: 2 x 960 strips at 100 μm x 2.5 cm
# PSp: 32 x 960 pixel at 100 μm x 1.5 mm

stdRes = {
    #'canvasName' : 'cEffX1', \
    #'histogramName' : 'effX1', \
    'histogramTitle' : 'my resolution', \
    'variable' : 'localPos.x()-rhLocalPos.x()', \
    'baseCuts' : 'tof<12.5&&hasRecHit>0', \
    #'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'local x [cm]', \
    'xNbins' : 200, \
    'xMin' : -0.1, \
    'xMax' : 0.1, \
    #'yTitle' : 'efficiency', \
    #'yMin' : 0.8, \
    #'yMax' : 1.05
}
    
effX1 = {
    'canvasName' : 'cEffX1', \
    'histogramName' : 'effX1', \
    'histogramTitle' : 'Efficiency 1D', \
    'variable' : 'localPos.x()', \
    'baseCuts' : 'tof<12.5', \
    'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'local x [cm]', \
    'xNbins' : 240, \
    'xMin' : -7.5, \
    'xMax' : 7.5, \
    'yTitle' : 'efficiency', \
    'yMin' : 0.8, \
    'yMax' : 1.05
}
    
eff2D1 = {
    'canvasName' : 'cEff2D1', \
    'histogramName' : 'eff2D1', \
    'histogramTitle' : 'Efficiency 2D', \
    'variable' : 'localPos.y():localPos.x()', \
    'baseCuts' : 'tof<12.5', \
    'effCuts' : 'hasRecHit>0&&abs(localPos.x()-rhLocalPos.x())<0.0075', \
    'xTitle' : 'local x [cm]', \
    'xNbins' : 75, \
    'xMin' : -7.5, \
    'xMax' : 7.5, \
    'yTitle' : 'efficiency', \
    'yNbins' : 15, \
    'yMin' : -7.5, \
    'yMax' : 7.5, \
    'zMin' : 0.75, \
    'zMax' : 1.05
}
