#ifndef HITANALYZER_SIMTRACKINFO_H
#define HITANALYZER_SIMTRACKINFO_H

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "FWCore/Framework/interface/ConsumesCollector.h"
//#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/EventSetup.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/Framework/interface/one/EDAnalyzer.h"
//#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/PluginManager/interface/ModuleDef.h"
//#include "FWCore/Utilities/interface/InputTag.h"

// #include "DataFormats/Common/interface/DetSetVector.h"
// #include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
// #include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
// #include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
// #include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

// #include "DataFormats/DetId/interface/DetId.h"
// #include "DataFormats/GeometrySurface/interface/LocalError.h"
// #include "DataFormats/GeometryVector/interface/LocalPoint.h"
// #include "Geometry/CommonDetUnit/interface/GeomDet.h"
// #include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
// #include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
// #include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
// #include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
// #include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
// #include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "TTree.h"

class SimTrackInfo {

public:
  struct SimTrackData {
    std::vector<ROOT::Math::XYZPointF> SimTrack_posTk;
    std::vector<ROOT::Math::XYZVectorF> SimTrack_momTk;
    std::vector<ROOT::Math::XYZVectorF> SimTrack_momentum;    
    /* std::vector<float> SimTrack_xTk; */
    /* std::vector<float> SimTrack_yTk; */
    /* std::vector<float> SimTrack_zTk; */
    /* std::vector<float> SimTrack_ptTk; */
    /* std::vector<float> SimTrack_etaTk; */
    /* std::vector<float> SimTrack_phiTk; */
    /* std::vector<float> SimTrack_pt; */
    /* std::vector<float> SimTrack_eta; */
    /* std::vector<float> SimTrack_phi; */
    std::vector<int> SimTrack_type;
    std::vector<int> SimTrack_charge;
    std::vector<unsigned short> SimTrack_trackInfo;
    std::vector<unsigned int> SimTrack_trackId;
  };
  SimTrackData simTrackData;

  SimTrackInfo() {}

  ~SimTrackInfo() {}

  void setBranches(TTree& tree);
  
  void fillSimTrackInfo(const SimTrack& simTrack);
  
  void clear();
  
};

#endif
