#ifndef HITANALYZER_SIMHITINFO_H
#define HITANALYZER_SIMHITINFO_H

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

// #include "DataFormats/GeometrySurface/interface/LocalError.h"
// #include "DataFormats/GeometryVector/interface/LocalPoint.h"
// #include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "TTree.h"

class SimHitInfo {

  typedef std::map< unsigned int, std::vector<const PSimHit*> > DetSimHitsMap;
  typedef std::pair< unsigned int, std::vector<const PSimHit*> > DetSimHitsPair;
  typedef std::map< unsigned int, std::vector<const Phase2TrackerRecHit1D*> > DetRecHitsMap;
  typedef std::pair< unsigned int, std::vector<const Phase2TrackerRecHit1D*> > DetRecHitsPair;
  
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
    std::vector<bool> hasRecHit;

    std::vector<ROOT::Math::XYZPointF> rhLocalPos;
    std::vector<ROOT::Math::XYZPointF> rhGlobalPos;
    // std::vector<float> Hit_cluster_global_x;
    // std::vector<float> Hit_cluster_global_y;
    // std::vector<float> Hit_cluster_global_z;
    /* std::vector<unsigned short> Hit_layer; */
    /* std::vector<unsigned short> Hit_ModuleType; */
    std::vector<unsigned short> clusterSize;
    std::vector<int> clusterNSimTracks;
    // std::vector<float> Hit_cluster_local_x;
    // std::vector<float> Hit_cluster_local_y;
    // std::vector<float> Hit_cluster_local_z;
    /* std::vector<bool> Hit_cluster_haveSimHit; */
    /* std::vector<float> Hit_cluster_closestSimHit_local_x; */
    /* std::vector<float> Hit_cluster_closestSimHit_local_y; */
    /* std::vector<float> Hit_cluster_closestSimHit_local_z; */
    std::vector<unsigned int> detRawId;
    std::vector<unsigned short> clusterFirstStrip;
    std::vector<unsigned short> clusterFirstRow;
    std::vector<unsigned short> clusterColumn;
    std::vector<unsigned short> clusterEdge;
    std::vector<unsigned short> clusterThreshold;
  };
  SimHitData simHitData;

  SimHitInfo() {
    tTopo_ = 0;
    tkGeom_ = 0;
  };

  ~SimHitInfo() {}

  void setBranches(TTree& tree);
  
  std::vector<unsigned int> getSimTrackId(const DetId&, unsigned int);
  
  void fillSimHitsPerDet(edm::Handle<edm::PSimHitContainer> *simHitsRaw);

  void fillRecHitsPerDet(const Phase2TrackerRecHit1DCollectionNew& rechits);

  const Phase2TrackerRecHit1D* matchRecHitOnDet(const PSimHit* simHit, const DetId& detId,
  						const std::vector<const Phase2TrackerRecHit1D*>& detRecHits);

 
  void fillSimHitInfo(const PSimHit& simHit);
  
  void clear();

  void setupEvent(const TrackerTopology* topo, const TrackerGeometry* geom,
		  const edm::DetSetVector<PixelDigiSimLink>* links,
		  edm::Handle<edm::PSimHitContainer> *simHitsRaw,
		  const Phase2TrackerRecHit1DCollectionNew& rechits) {
    tTopo_ = topo;
    tkGeom_ = geom;
    pixelSimLinks = links;
    //
    // fill SimHit map
    //
    // simHitsPerDet_.clear();
    fillSimHitsPerDet(simHitsRaw);
    //
    // fill RecHit map (uses SimHit map!)
    //
    // recHitsPerDet_.clear();
    fillRecHitsPerDet(rechits);
  }
  
 private:
  const TrackerTopology* tTopo_;
  const TrackerGeometry* tkGeom_;
  const edm::DetSetVector<PixelDigiSimLink>* pixelSimLinks;

  DetSimHitsMap simHitsPerDet_;
  DetRecHitsMap recHitsPerDet_;
};

#endif

