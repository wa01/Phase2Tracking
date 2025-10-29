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
