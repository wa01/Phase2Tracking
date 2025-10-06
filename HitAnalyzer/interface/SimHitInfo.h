#ifndef HITANALYZER_SIMHITINFO_H
#define HITANALYZER_SIMHITINFO_H

// This is largely copied from https://github.com/cms-sw/cmssw/blob/master/RecoLocalTracker/Phase2TrackerRecHits/test/RecHitsValidation.cc

// #include "CommonTools/UtilAlgos/interface/TFileService.h"
// #include "CommonTools/Utils/interface/TFileDirectory.h"

// #include "FWCore/ServiceRegistry/interface/Service.h"
// #include "FWCore/Framework/interface/ConsumesCollector.h"
// #include "FWCore/Framework/interface/ESHandle.h"
// #include "FWCore/Framework/interface/ESHandle.h"
// #include "FWCore/Framework/interface/Event.h"
// #include "FWCore/Framework/interface/EventSetup.h"
// #include "FWCore/Framework/interface/MakerMacros.h"
// #include "FWCore/Framework/interface/one/EDAnalyzer.h"
// #include "FWCore/MessageLogger/interface/MessageLogger.h"
// #include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
// #include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "FWCore/PluginManager/interface/ModuleDef.h"
// #include "FWCore/Utilities/interface/InputTag.h"

// #include "DataFormats/Common/interface/DetSetVector.h"
// #include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
// #include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
// #include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

// #include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
// #include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
// #include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
// #include "SimDataFormats/Track/interface/SimTrackContainer.h"
// #include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
// #include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "TTree.h"

class SimHitInfo {

public:
  struct SimHitData {
    std::vector<ROOT::Math::XYZPointF> localPos;
    std::vector<ROOT::Math::XYZPointF> globalPos;
    std::vector<ROOT::Math::XYZVectorF> localDir;
    std::vector<ROOT::Math::XYZVectorF> globalDir;
    std::vector<ROOT::Math::XYZVectorF> path;
    std::vector<float> theta;
    std::vector<float> phi;
    std::vector<float> pabs;
    std::vector<float> tof;
    std::vector<float> energyLoss;
    std::vector<unsigned short> processType;
    std::vector<int>   particleType;
    std::vector<unsigned short> layer;
    std::vector<unsigned short> moduleType;
    std::vector<ROOT::Math::XYZVectorF> detNormal;
    std::vector<unsigned int> trackId;
  };
  SimHitData simHitData;

  SimHitInfo() {
    tTopo = 0;
    tkGeom = 0;
  };

  ~SimHitInfo() {}

  void setBranches(TTree& tree);
  
  void fillSimHitInfo(const PSimHit& simHit);
  
  void clear();
  
  void setTopology(const TrackerTopology* topo) {
    tTopo = topo;
  }

  void setGeometry(const TrackerGeometry* geom) {
    tkGeom = geom;
  }  
  

 private:
  const TrackerTopology* tTopo;
  const TrackerGeometry* tkGeom;
};

#endif
