#include "UserAnalysis/HitAnalyzer/interface/SimHitInfo.h"
#include "DataFormats/DetId/interface/DetId.h"

void SimHitInfo::setBranches(TTree& tree) {
  tree.Branch("localPos",&simHitData.localPos);
  tree.Branch("globalPos",&simHitData.globalPos);
  tree.Branch("localDir",&simHitData.localDir);
  tree.Branch("globalDir",&simHitData.globalDir);
  tree.Branch("path",&simHitData.path);
  tree.Branch("theta",&simHitData.theta);
  tree.Branch("phi",&simHitData.phi);
  tree.Branch("pabs",&simHitData.pabs);
  tree.Branch("tof",&simHitData.tof);
  tree.Branch("energyLoss",&simHitData.energyLoss);
  tree.Branch("processType",&simHitData.processType);
  tree.Branch("particleType",&simHitData.particleType);
  tree.Branch("layer",&simHitData.layer);
  tree.Branch("moduleType",&simHitData.moduleType);
  tree.Branch("detNormal",&simHitData.detNormal);
  tree.Branch("trackId",&simHitData.trackId);
};
  
void SimHitInfo::fillSimHitInfo(const PSimHit& simHit) {
  //
  // Get the detector unit's id
  //
  DetId detId(simHit.detUnitId());
  unsigned int layer = (tTopo->side(detId) != 0) * 1000;  // don't split up endcap sides
  layer += tTopo->layer(detId);
  TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
  // Get the geomdet
  const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
  if (!geomDetUnit) {
    std::cout << "*** did not find geomDetUnit ***" << std::endl;
    return;
  }
  ROOT::Math::XYZPointF localPos(simHit.localPosition().x(),simHit.localPosition().y(),
				 simHit.localPosition().z());
  simHitData.localPos.push_back(localPos);
  GlobalPoint globalPosition(geomDetUnit->toGlobal(simHit.localPosition()));
  ROOT::Math::XYZPointF globalPos(globalPosition.x(),globalPosition.y(),globalPosition.z());
  simHitData.globalPos.push_back(globalPos);
  ROOT::Math::XYZVectorF localDir(simHit.localDirection().x(),simHit.localDirection().y(),
				  simHit.localDirection().z());
  simHitData.localDir.push_back(localDir);
  GlobalVector globalDirection(geomDetUnit->toGlobal(simHit.localDirection()));
  ROOT::Math::XYZVectorF globalDir(globalDirection.x(),globalDirection.y(),globalDirection.z());
  simHitData.globalDir.push_back(globalDir);
  ROOT::Math::XYZVectorF path(simHit.exitPoint().x()-simHit.entryPoint().x(),
			      simHit.exitPoint().y()-simHit.entryPoint().y(),
			      simHit.exitPoint().z()-simHit.entryPoint().z());
  simHitData.path.push_back(path);
  simHitData.theta.push_back(simHit.thetaAtEntry());
  simHitData.phi.push_back(simHit.phiAtEntry());
  simHitData.pabs.push_back(simHit.pabs());
  simHitData.tof.push_back(simHit.timeOfFlight());
  simHitData.energyLoss.push_back(simHit.energyLoss());
  simHitData.processType.push_back(simHit.processType());
  simHitData.particleType.push_back(simHit.particleType());
  simHitData.layer.push_back(layer);
  simHitData.moduleType.push_back((unsigned short)mType);
  GlobalVector detNormalGlobal(geomDetUnit->toGlobal(LocalVector(0.,0.,1.)));
  ROOT::Math::XYZVectorF detNormal(detNormalGlobal.x(),detNormalGlobal.y(),detNormalGlobal.z());
  simHitData.detNormal.push_back(detNormal);
  simHitData.trackId.push_back(simHit.trackId());
};
  
void SimHitInfo::clear() {
  simHitData.localPos.clear();
  simHitData.globalPos.clear();
  simHitData.localDir.clear();
  simHitData.globalDir.clear();
  simHitData.path.clear();
  simHitData.theta.clear();
  simHitData.phi.clear();
  simHitData.pabs.clear();
  simHitData.tof.clear();
  simHitData.energyLoss.clear();
  simHitData.processType.clear();
  simHitData.particleType.clear();
  simHitData.layer.clear();
  simHitData.moduleType.clear();
  simHitData.detNormal.clear();
  simHitData.trackId.clear();
};
