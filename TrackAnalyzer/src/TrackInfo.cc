#include "Phase2Tracking/TrackAnalyzer/interface/TrackInfo.h"
// #include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

void TrackInfo::setBranches(TTree& tree) {
  //
  // Branch definitions for output tree
  //
  // reco::Track-related variables
  //
  tree.Branch("Track_chi2", &trackData.Track_chi2);
  tree.Branch("Track_ndof", &trackData.Track_ndof);
  tree.Branch("Track_validHits", &trackData.Track_validHits);
  tree.Branch("Track_lostHits", &trackData.Track_lostHits);
  tree.Branch("Track_missInnerHits", &trackData.Track_missInnerHits);
  tree.Branch("Track_missOuterHits", &trackData.Track_missOuterHits);
  tree.Branch("Track_algo", &trackData.Track_algo);
  // tree.Branch("Track_quality", &trackData.Track_quality);
  tree.Branch("Track_charge", &trackData.Track_charge);
  tree.Branch("Track_dxy", &trackData.Track_dxy);
  tree.Branch("Track_dz", &trackData.Track_dz);
  tree.Branch("Track_pt", &trackData.Track_pt);
  tree.Branch("Track_eta", &trackData.Track_eta);
  tree.Branch("Track_phi", &trackData.Track_phi);
  tree.Branch("Track_dxyErr", &trackData.Track_dxyErr);
  tree.Branch("Track_dzErr", &trackData.Track_dzErr);
  tree.Branch("Track_ptErr", &trackData.Track_ptErr);
  tree.Branch("Track_etaErr", &trackData.Track_etaErr);
  tree.Branch("Track_phiErr", &trackData.Track_phiErr);
  tree.Branch("Track_refPos", &trackData.Track_refPos);
  tree.Branch("Track_momentum", &trackData.Track_momentum);
  //
  // Vectors of hit-related variables
  //
  tree.Branch("Hit_clusterSize", &trackData.Hit_clusterSize);
  tree.Branch("Hit_layer", &trackData.Hit_layer);
  tree.Branch("Hit_mType", &trackData.Hit_mType);
  // tree.Branch("Hit_globPos", &trackData.Hit_globPos);
  tree.Branch("Hit_locPos", &trackData.Hit_locPos);
  //
  // SimTrack-related variables
  //
  // tree.Branch("SimTrack_posTk", &trackData.SimTrack_posTk);
  // tree.Branch("SimTrack_momTk", &trackData.SimTrack_momTk);
  tree.Branch("SimTrack_momentum", &trackData.SimTrack_momentum);    
  tree.Branch("SimTrack_type", &trackData.SimTrack_type);
  tree.Branch("SimTrack_charge", &trackData.SimTrack_charge);
  tree.Branch("SimTrack_trackInfo", &trackData.SimTrack_trackInfo);
  tree.Branch("SimTrack_trackId", &trackData.SimTrack_trackId);

};

void TrackInfo::fillTrackInfo(const reco::Track* track, const SimTrack* simTrack, const SimVertex* simVertex,
			      TTree& tree) {
  //
  // fill TTree "tree" with information of a reco::Track, possibly with associated SimTrack / SimVertex
  //
  //
  // define & store reference point for impacts (default = origin)
  // TO BE IMPROVED:
  //    reference point (=simVertex) only used for dxy and dz ( not for the uncertainties or momentum) !!!!
  //
  math::XYZPoint refPoint(0.,0.,0.);
  if ( simVertex!=nullptr ) {
    refPoint.SetX(simVertex->position().x());
    refPoint.SetY(simVertex->position().y());
    refPoint.SetZ(simVertex->position().z());
  }
  trackData.Track_refPos.SetXYZ(refPoint.x(),refPoint.y(),refPoint.z());
  trackData.Track_momentum.SetXYZ(track->momentum().x(),track->momentum().y(),track->momentum().z()); 
  // trackData.Track_chi2.push_back(track->chi2());
  trackData.Track_chi2 = track->chi2();
  trackData.Track_ndof = track->ndof();
  trackData.Track_validHits = track->numberOfValidHits();
  trackData.Track_lostHits = track->numberOfLostHits();
  trackData.Track_missInnerHits = track->missingInnerHits();
  trackData.Track_missOuterHits = track->missingOuterHits();
  trackData.Track_algo = track->algo();
  // trackData.Track_quality = track->quality());
  trackData.Track_charge = track->charge();
  trackData.Track_dxy = track->dxy(refPoint);
  trackData.Track_dz = track->dz(refPoint);
  trackData.Track_pt = track->pt();
  trackData.Track_eta = track->eta();
  trackData.Track_phi = track->phi();
  trackData.Track_dxyErr = track->dxyError();
  trackData.Track_dzErr = track->dzError();
  trackData.Track_ptErr = track->ptError();
  trackData.Track_etaErr = track->etaError();
  trackData.Track_phiErr = track->phiError();
  //
  // Variables / TrackingRecHit
  //
  for ( auto rh=track->recHitsBegin(); rh!=track->recHitsEnd(); ++rh ) {
    //
    // Initialization of hit-related variables
    //
    short int hitClusterSize = -1;
    unsigned short hitMtype(0);
    unsigned short hitLayer = 0;
    // std::cout << "RH valid " << (**rh).isValid() << std::endl;
    // ROOT::Math::XYZPointF hitGlobPos(0.,0.,0.);
    ROOT::Math::XYZPointF hitLocPos(0.,0.,0.);
    // ROOT::Math::XYZPointF hitGlobPos((**rh).globalPosition().x(),
    // 				       (**rh).globalPosition().y(),(**rh).globalPosition().z());
     
    //
    // Det-related information
    //
    DetId detId = (**rh).geographicalId();
    // std::cout << "detId " << detId.det() << " " << detId.subdetId() << std::endl;
    //
    // Only use hits in the Phase2 outer tracker
    if ( detId.det()==DetId::Tracker && detId.subdetId()==5 ) {
      // std::cout << "Found Tracker" << std::endl;
      // std::cout << "layer " << (tkTopo_->side(detId) != 0) * 1000 << std::endl;
      // Det-related information
      hitLayer = (tkTopo_->side(detId) != 0) * 1000 + tkTopo_->layer(detId);
      hitMtype = (unsigned short)tkGeom_->getDetectorType(detId);
      //
      // for valid hits, make sure that hits are of type Phase2TrackerRecHit1D
      //
      if ( (**rh).isValid() ) {
	//
	// verify via dynamic cast
	//
	const Phase2TrackerRecHit1D* phase2rh = dynamic_cast<const Phase2TrackerRecHit1D*>(&(**rh));
	hitLocPos.SetXYZ(phase2rh->localPosition().x(),phase2rh->localPosition().y(),phase2rh->localPosition().z());
	
	if ( phase2rh ) {
	  //
	  // access to cluster
	  //
	  Phase2TrackerRecHit1D::ClusterRef const& cluster = (*phase2rh).cluster();
	  if ( cluster.isNonnull() ) {
	    // std::cout << typeid(*cluster).name() << " " << cluster->size() << std::endl << std::endl;
	    hitClusterSize = cluster->size();
	  }
	  // else {
	  // 	std::cout << "Null cluster" << std::endl;
	  // }

	}
      }
      //
      // add to vectors with Hit-related information
      //
      trackData.Hit_clusterSize.push_back(hitClusterSize);
      trackData.Hit_layer.push_back(hitLayer);
      trackData.Hit_mType.push_back(hitMtype);
      // trackData.Hit_globPos.push_back(hitGlobPos);
      trackData.Hit_locPos.push_back(hitLocPos);
	      
    }
    // else {
    //   std::cout << "No tracker !!!" << std::endl;
    //   continue;
    // }
  }
  //
  // SimTrack information
  //
  if ( simTrack!=nullptr ) {
    trackData.SimTrack_momentum.SetXYZ(simTrack->momentum().px(),simTrack->momentum().py(),simTrack->momentum().pz());
  
    trackData.SimTrack_type = simTrack->type();
    trackData.SimTrack_charge = simTrack->charge();
    trackData.SimTrack_trackInfo = simTrack->getTrackInfo();
    trackData.SimTrack_trackId = simTrack->trackId();
  }
  //
  // Fill TTree
  //
  tree.Fill();      
};

void TrackInfo::clear() {
  //
  // Clear vectors at each event
  //
  trackData.Hit_clusterSize.clear();
  trackData.Hit_layer.clear();
  trackData.Hit_mType.clear();
  // trackData.Hit_globPos.clear();
  trackData.Hit_locPos.clear();

};
