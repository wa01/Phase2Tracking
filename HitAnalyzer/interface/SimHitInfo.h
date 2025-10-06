#ifndef HITANALYZER_SIMHITINFO_H
#define HITANALYZER_SIMHITINFO_H

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

// #include "DataFormats/GeometrySurface/interface/LocalError.h"
// #include "DataFormats/GeometryVector/interface/LocalPoint.h"
// #include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

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
