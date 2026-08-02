#ifndef HITANALYZER_SIMTRACKINFO_H
#define HITANALYZER_SIMTRACKINFO_H

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
// #include "SimDataFormats/Track/interface/SimTrackFwd.h"
// #include "SimDataFormats/Vertex/interface/SimVertexFwd.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "TTree.h"

class TrackInfo {

public:
  struct TrackData {
    float Track_chi2;
    unsigned short Track_ndof;
    unsigned short Track_validHits;
    unsigned short Track_lostHits;
    unsigned short Track_missInnerHits;
    unsigned short Track_missOuterHits;
    unsigned short Track_algo;
    // std::vector<short> Track_quality;
    short Track_charge;
    float Track_dxy;
    float Track_dz;
    float Track_pt;
    float Track_eta;
    float Track_phi;
    float Track_dxyErr;
    float Track_dzErr;
    float Track_ptErr;
    float Track_etaErr;
    float Track_phiErr;
    /* ROOT::Math::XYZVectorF Track_momTk; */
    ROOT::Math::XYZVectorF Track_momentum;
    ROOT::Math::XYZPointF Track_refPos;

    std::vector<short int> Hit_clusterSize;
    std::vector<unsigned short> Hit_layer;
    /* std::vector<ROOT::Math::XYZPointF> Hit_globPos; */
    /* std::vector<ROOT::Math::XYZPointF> Hit_locPos; */
    
    /* std::vector<ROOT::Math::XYZPointF> SimTrack_posTk; */
    /* std::vector<ROOT::Math::XYZVectorF> SimTrack_momTk; */
    ROOT::Math::XYZVectorF SimTrack_momentum;    
    int SimTrack_type;
    int SimTrack_charge;
    unsigned short SimTrack_trackInfo;
    int SimTrack_trackId;
  };
  TrackData trackData;

  const TrackerGeometry* tkGeom_;
  const TrackerTopology* tkTopo_;

  TrackInfo() {
    tkGeom_ = 0;
    tkTopo_ = 0;
  }

  ~TrackInfo() {}

  void setupEvent(const TrackerGeometry* geom, const TrackerTopology* topo) {
    tkGeom_ = geom;
    tkTopo_ = topo;
  }
  
  void setBranches(TTree& tree);
  
  void fillTrackInfo(const reco::Track* track, const SimTrack* simTrack, const SimVertex* simVertex,
		     TTree& tree);
  
  void clear();
  
};

#endif
