#ifndef HITANALYZER_RECHITINFO_H
#define HITANALYZER_RECHITINFO_H

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

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
// #include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

// #include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
// #include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
// #include "SimDataFormats/Vertex/interface/SimVertexContainer.h"


#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "TTree.h"

class RecHitInfo {

public:
  struct RecHitData
  {
    std::vector<float> Hit_cluster_global_x;
    std::vector<float> Hit_cluster_global_y;
    std::vector<float> Hit_cluster_global_z;
    std::vector<unsigned short> Hit_layer;
    std::vector<unsigned short> Hit_ModuleType;
    std::vector<unsigned short> Hit_cluster_size;
    std::vector<int> Hit_cluster_SimTrack_size;
    std::vector<float> Hit_cluster_local_x;
    std::vector<float> Hit_cluster_local_y;
    std::vector<float> Hit_cluster_local_z;
    std::vector<bool> Hit_cluster_haveSimHit;
    std::vector<float> Hit_cluster_closestSimHit_local_x;
    std::vector<float> Hit_cluster_closestSimHit_local_y;
    std::vector<float> Hit_cluster_closestSimHit_local_z;
    std::vector<unsigned int> Hit_det_rawid;
    std::vector<unsigned short> Hit_cluster_firstStrip;
    std::vector<unsigned short> Hit_cluster_firstRow;
    std::vector<unsigned short> Hit_cluster_column;
    std::vector<unsigned short> Hit_cluster_edge;
    std::vector<unsigned short> Hit_cluster_threshold;
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

  RecHitData recHitData;

  RecHitInfo() {
    tTopo = 0;
    tkGeom = 0;
  };

  ~RecHitInfo() {}

  void setBranches(TTree& tree);
  
  std::vector<unsigned int> getSimTrackId(
    edm::Handle<edm::DetSetVector<PixelDigiSimLink> >& pixelSimLinks,
    const DetId& detId, unsigned int channel);
  
  void fillSimHitInfo(const PSimHit& simHit);

  void fillRecHitInfo(const Phase2TrackerRecHit1D& recHit, unsigned int rawid, const GeomDetUnit* geomDetUnit,
		      edm::Handle<edm::DetSetVector<PixelDigiSimLink> >* pixelSimLinks,
		      std::map<unsigned int, SimTrack>& simTracks, edm::Handle<edm::PSimHitContainer>* simHitsRaw);
  
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
  //const edm::Handle<edm::DetSetVector<PixelDigiSimLink> >* pixelSimLinks_;
};

#endif
