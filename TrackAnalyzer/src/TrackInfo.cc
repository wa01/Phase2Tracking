#include "Phase2Tracking/TrackAnalyzer/interface/TrackInfo.h"
// #include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

void TrackInfo::setBranches(TTree& tree) {
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

    tree.Branch("Hit_clusterSize", &trackData.Hit_clusterSize);
    tree.Branch("Hit_layer", &trackData.Hit_layer);
    // tree.Branch("Hit_globPos", &trackData.Hit_globPos);
    // tree.Branch("Hit_locPos", &trackData.Hit_locPos);

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
    // trackData.Track_refPos = ROOT::Math::XYZPointF(refPoint.x(),refPoint.y(),refPoint.z()));
    trackData.Track_refPos.SetXYZ(refPoint.x(),refPoint.y(),refPoint.z());
    //
    // store Track momentum vector
    //
    // ROOT::Math::XYZVectorF tkMom(track->momentum().x(),track->momentum().y(),track->momentum().z());
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

    for ( auto rh=track->recHitsBegin(); rh!=track->recHitsEnd(); ++rh ) {
      DetId detId = (**rh).geographicalId();
      // std::cout << "detId " << detId.det() << " " << detId.subdetId() << std::endl;

      short int hitClusterSize = -1;
      unsigned short hitLayer = 0;
      // std::cout << "RH valid " << (**rh).isValid() << std::endl;
      // ROOT::Math::XYZPointF hitGlobPos(0.,0.,0.);
      // ROOT::Math::XYZPointF hitLocPos(0.,0.,0.);
      // ROOT::Math::XYZPointF hitGlobPos((**rh).globalPosition().x(),
      // 				       (**rh).globalPosition().y(),(**rh).globalPosition().z());
     
      if ( detId.det()==DetId::Tracker && detId.subdetId()==5 ) {
	// std::cout << "Found Tracker" << std::endl;
	// std::cout << "layer " << (tkTopo_->side(detId) != 0) * 1000 << std::endl;
	hitLayer = (tkTopo_->side(detId) != 0) * 1000;
	if ( (**rh).isValid() ) {
	  TrackerGeometry::ModuleType mType = tkGeom_->getDetectorType(detId);
	  const Phase2TrackerRecHit1D* phase2rh = dynamic_cast<const Phase2TrackerRecHit1D*>(&(**rh));
	  if ( phase2rh ) {
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
	trackData.Hit_clusterSize.push_back(hitClusterSize);
	trackData.Hit_layer.push_back(hitLayer);
	// trackData.Hit_globPos.push_back(hitGlobPos);
	// trackData.Hit_locPos.push_back(hitLocPos);
	      
      }
      else {
	std::cout << "No tracker !!!" << std::endl;
	continue;
      }
      tree.Fill();      
    }
    // ROOT::Math::XYZPointF tkPos(simTrack->trackerSurfacePosition().x(),
    // 				simTrack->trackerSurfacePosition().y(),
    // 				simTrack->trackerSurfacePosition().z());
    // trackData.SimTrack_posTk.push_back(tkPos);
    
    // ROOT::Math::XYZVectorF tkMom(simTrack->trackerSurfaceMomentum().px(),
    // 				 simTrack->trackerSurfaceMomentum().py(),
    // 				 simTrack->trackerSurfaceMomentum().pz());
    // trackData.SimTrack_momTk.push_back(tkMom);
  // std::cout << "Adding simtrack mom pt " << simTrack.trackerSurfaceMomentum().pt() 
  // 	    << " " << simTrack.momentum().pt() << std::endl;
  
    trackData.SimTrack_momentum.SetXYZ(simTrack->momentum().px(),simTrack->momentum().py(),simTrack->momentum().pz());
  
    trackData.SimTrack_type = simTrack->type();
    trackData.SimTrack_charge = simTrack->charge();
    trackData.SimTrack_trackInfo = simTrack->getTrackInfo();
    trackData.SimTrack_trackId = simTrack->trackId();
};

void TrackInfo::clear() {
    // trackData.Track_chi2.clear();
    // trackData.Track_ndof.clear();
    // trackData.Track_validHits.clear();
    // trackData.Track_lostHits.clear();
    // trackData.Track_missInnerHits.clear();
    // trackData.Track_missOuterHits.clear();
    // trackData.Track_algo.clear();
    // // trackData.Track_quality.clear();
    // trackData.Track_charge.clear();
    // trackData.Track_dxy.clear();
    // trackData.Track_dz.clear();
    // trackData.Track_pt.clear();
    // trackData.Track_eta.clear();
    // trackData.Track_phi.clear();
    // trackData.Track_dxyErr.clear();
    // trackData.Track_dzErr.clear();
    // trackData.Track_ptErr.clear();
    // trackData.Track_etaErr.clear();
    // trackData.Track_phiErr.clear();
    // trackData.Track_refPos.clear();
    // trackData.Track_momTk.clear();

  trackData.Hit_clusterSize.clear();
  trackData.Hit_clusterSize.clear();
  trackData.Hit_layer.clear();
  // trackData.Hit_globPos.clear();
  // trackData.Hit_locPos.clear();

    // // trackData.SimTrack_posTk.clear();
    // // trackData.SimTrack_momTk.clear();
    // trackData.SimTrack_momentum.clear();    
    // trackData.SimTrack_type.clear();
    // trackData.SimTrack_charge.clear();
    // trackData.SimTrack_trackInfo.clear();
    // trackData.SimTrack_trackId.clear();
    

};

