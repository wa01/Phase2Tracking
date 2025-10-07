#include "UserAnalysis/HitAnalyzer/interface/SimTrackInfo.h"
// #include "DataFormats/DetId/interface/DetId.h"
#include "SimDataFormats/Track/interface/SimTrack.h"

void SimTrackInfo::setBranches(TTree& tree) {
  tree.Branch("SimTrack_posTk", &simTrackData.SimTrack_posTk);
  tree.Branch("SimTrack_momTk", &simTrackData.SimTrack_momTk);
  tree.Branch("SimTrack_momentum", &simTrackData.SimTrack_momentum);
  // tree.Branch("SimTrack_xTk",         &simTrackData.SimTrack_xTk);
  // tree.Branch("SimTrack_yTk",         &simTrackData.SimTrack_yTk);
  // tree.Branch("SimTrack_zTk",         &simTrackData.SimTrack_zTk);
  // tree.Branch("SimTrack_ptTk",        &simTrackData.SimTrack_ptTk);
  // tree.Branch("SimTrack_etaTk",       &simTrackData.SimTrack_etaTk);
  // tree.Branch("SimTrack_phiTk",       &simTrackData.SimTrack_phiTk);
  // tree.Branch("SimTrack_pt",        &simTrackData.SimTrack_pt);
  // tree.Branch("SimTrack_eta",       &simTrackData.SimTrack_eta);
  // tree.Branch("SimTrack_phi",       &simTrackData.SimTrack_phi);
  tree.Branch("SimTrack_type",      &simTrackData.SimTrack_type);
  tree.Branch("SimTrack_charge",    &simTrackData.SimTrack_charge);
  tree.Branch("SimTrack_trackInfo", &simTrackData.SimTrack_trackInfo);
};
  
void SimTrackInfo::fillSimTrackInfo(const SimTrack& simTrack) {
  ROOT::Math::XYZPointF tkPos(simTrack.trackerSurfacePosition().x(),
			      simTrack.trackerSurfacePosition().y(),
			      simTrack.trackerSurfacePosition().z());
  simTrackData.SimTrack_posTk.push_back(tkPos);
  
  // ROOT::Math::RhoEtaPhiVectorF tkMom(simTrack.trackerSurfaceMomentum().pt(),
  // 				     simTrack.trackerSurfaceMomentum().eta(),
  // 				     simTrack.trackerSurfaceMomentum().phi());
  ROOT::Math::XYZVectorF tkMom(simTrack.trackerSurfaceMomentum().px(),
			       simTrack.trackerSurfaceMomentum().py(),
			       simTrack.trackerSurfaceMomentum().pz());
  simTrackData.SimTrack_momTk.push_back(tkMom);
  // std::cout << "Adding simtrack mom pt " << simTrack.trackerSurfaceMomentum().pt() 
  // 	    << " " << simTrack.momentum().pt() << std::endl;
  
  // ROOT::Math::RhoEtaPhiVectorF mom(simTrack.momentum().pt(),simTrack.momentum().eta(),simTrack.momentum().phi());
  ROOT::Math::XYZVectorF mom(simTrack.momentum().px(),simTrack.momentum().py(),simTrack.momentum().pz());
  simTrackData.SimTrack_momentum.push_back(mom);
  
  // simTrackData.SimTrack_xTk.push_back(simTrack.trackerSurfacePosition().x());
  // simTrackData.SimTrack_yTk.push_back(simTrack.trackerSurfacePosition().y());
  // simTrackData.SimTrack_zTk.push_back(simTrack.trackerSurfacePosition().z());
  // simTrackData.SimTrack_ptTk.push_back(simTrack.trackerSurfaceMomentum().pt());
  // simTrackData.SimTrack_etaTk.push_back(simTrack.trackerSurfaceMomentum().eta());
  // simTrackData.SimTrack_phiTk.push_back(simTrack.trackerSurfaceMomentum().phi());
  // simTrackData.SimTrack_pt.push_back(simTrack.momentum().pt());
  // simTrackData.SimTrack_eta.push_back(simTrack.momentum().eta());
  // simTrackData.SimTrack_phi.push_back(simTrack.momentum().phi());
  simTrackData.SimTrack_type.push_back(simTrack.type());
  simTrackData.SimTrack_charge.push_back(simTrack.charge());
  simTrackData.SimTrack_trackInfo.push_back(simTrack.getTrackInfo());
};


  
void SimTrackInfo::clear() {
  simTrackData.SimTrack_posTk.clear();
  simTrackData.SimTrack_momTk.clear();
  simTrackData.SimTrack_momentum.clear();
  // simTrackData.SimTrack_xTk.clear();
  // simTrackData.SimTrack_yTk.clear();
  // simTrackData.SimTrack_zTk.clear();
  // simTrackData.SimTrack_ptTk.clear();
  // simTrackData.SimTrack_etaTk.clear();
  // simTrackData.SimTrack_phiTk.clear();
  // simTrackData.SimTrack_pt.clear();
  // simTrackData.SimTrack_eta.clear();
  // simTrackData.SimTrack_phi.clear();
  simTrackData.SimTrack_type.clear();
  simTrackData.SimTrack_charge.clear();
  simTrackData.SimTrack_trackInfo.clear();
};
