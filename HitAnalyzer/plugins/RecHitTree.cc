// This is largely copied from https://github.com/cms-sw/cmssw/blob/master/RecoLocalTracker/Phase2TrackerRecHits/test/RecHitsValidation.cc

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"

struct HitInfo
{
  std::vector<float> Hit_cluster_global_x;
  std::vector<float> Hit_cluster_global_y;
  std::vector<float> Hit_cluster_global_z;
  std::vector<int> Hit_layer;
  std::vector<int> Hit_ModuleType;
  std::vector<int> Hit_cluster_size;
  std::vector<int> Hit_cluster_SimTrack_size;
  std::vector<float> Hit_cluster_local_x;
  std::vector<float> Hit_cluster_local_y;
  std::vector<float> Hit_cluster_local_z;
  std::vector<bool> Hit_cluster_haveSimHit;
  std::vector<float> Hit_cluster_closestSimHit_local_x;
  std::vector<float> Hit_cluster_closestSimHit_local_y;
  std::vector<float> Hit_cluster_closestSimHit_local_z;
  std::vector<float> Hit_cluster_closestSimHit_p;
  std::vector<float> Hit_cluster_closestSimHit_tof;
  std::vector<float> Hit_cluster_closestSimHit_E;
  std::vector<int> Hit_SimTrack_idx;
  std::vector<int> Hit_SimHit_idx;

  std::vector<int> SimHit_barrelendcap; // 0: Barrel; 1: Endcap
  std::vector<int> SimHit_layer;
  std::vector<int> SimHit_ModuleType;
  std::vector<float> SimHit_local_x;
  std::vector<float> SimHit_local_y;
  std::vector<float> SimHit_local_z;
  std::vector<float> SimHit_global_x;
  std::vector<float> SimHit_global_y;
  std::vector<float> SimHit_global_z;
  std::vector<float> SimHit_p;
  std::vector<float> SimHit_tof;
  std::vector<float> SimHit_E;
  std::vector<int> SimHit_pdgId;
  std::vector<int> SimHit_processType;
  std::vector<int> SimHit_SimTrack_idx;
  std::vector<int> SimHit_n_matched_Hit;
  std::vector<int> SimHit_Hit_idx;

  std::vector<float> SimTrack_pt;
  std::vector<float> SimTrack_eta;
  std::vector<float> SimTrack_phi;
  std::vector<int> SimTrack_pdgId;
  std::vector<int> SimTrack_n_matched_hits;
  std::vector<int> SimTrack_n_simhits;
  //std::vector<float> SimTrack_Hit_idx;
};

class RecHitTree : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit RecHitTree(const edm::ParameterSet&);
    void analyze(const edm::Event&, const edm::EventSetup&);

  private:
    virtual void beginJob() override;
    virtual void endJob() override;
    void initEventStructure();
    std::vector<unsigned int> getSimTrackId(const edm::Handle<edm::DetSetVector<PixelDigiSimLink> >&, const DetId&, unsigned int);

    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> esTokenGeom_;
    const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> esTokenTopo_;
    const edm::EDGetTokenT<Phase2TrackerRecHit1DCollectionNew> tokenRecHits_;
    const edm::EDGetTokenT<Phase2TrackerCluster1DCollectionNew> tokenClusters_;
    const edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink> > tokenLinks_;
    const edm::EDGetTokenT<edm::PSimHitContainer> tokenSimHitsB_;
    const edm::EDGetTokenT<edm::PSimHitContainer> tokenSimHitsE_;
    const edm::EDGetTokenT<edm::SimTrackContainer> tokenSimTracks_;
    const edm::EDGetTokenT<std::vector<reco::GenParticle>> tokenGenParticles_;

    const double simtrackminpt_;

    TTree* hitTree;
    HitInfo* hitInfo;

};

RecHitTree::RecHitTree(const edm::ParameterSet& cfg)
  : esTokenGeom_(esConsumes()),
    esTokenTopo_(esConsumes()),
    tokenRecHits_(consumes<Phase2TrackerRecHit1DCollectionNew>(cfg.getParameter<edm::InputTag>("rechits"))),
    tokenClusters_(consumes<Phase2TrackerCluster1DCollectionNew>(cfg.getParameter<edm::InputTag>("clusters"))),
    tokenLinks_(consumes<edm::DetSetVector<PixelDigiSimLink> >(cfg.getParameter<edm::InputTag>("links"))),
    tokenSimHitsB_(consumes<edm::PSimHitContainer>(cfg.getParameter<edm::InputTag>("simhitsbarrel"))),
    tokenSimHitsE_(consumes<edm::PSimHitContainer>(cfg.getParameter<edm::InputTag>("simhitsendcap"))),
    tokenSimTracks_(consumes<edm::SimTrackContainer>(cfg.getParameter<edm::InputTag>("simtracks"))),
    tokenGenParticles_(consumes<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genparticles"))),
    simtrackminpt_(cfg.getParameter<double>("SimTrackMinPt"))
{
  hitInfo = new HitInfo;
}

void RecHitTree::analyze(const edm::Event& event, const edm::EventSetup& eventSetup) {
  initEventStructure();

  // Get the geometry
  const TrackerGeometry* tkGeom = &eventSetup.getData(esTokenGeom_);
  const TrackerTopology* tTopo = &eventSetup.getData(esTokenTopo_);

  // Get the RecHits
  edm::Handle<Phase2TrackerRecHit1DCollectionNew> rechits;
  event.getByToken(tokenRecHits_, rechits);

  // Get the Clusters
  edm::Handle<Phase2TrackerCluster1DCollectionNew> clusters;
  event.getByToken(tokenClusters_, clusters);

  // Get the PixelDigiSimLinks
  edm::Handle<edm::DetSetVector<PixelDigiSimLink> > pixelSimLinks;
  event.getByToken(tokenLinks_, pixelSimLinks);

  // Get the SimHits
  edm::Handle<edm::PSimHitContainer> simHitsRaw[2];
  event.getByToken(tokenSimHitsB_, simHitsRaw[0]);
  event.getByToken(tokenSimHitsE_, simHitsRaw[1]);

  // Get the SimTracks
  edm::Handle<edm::SimTrackContainer> simTracksRaw;
  event.getByToken(tokenSimTracks_, simTracksRaw);

  // Get the GenParticles
  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByToken(tokenGenParticles_, genParticles);

  // Rearrange the simTracks for ease of use <simTrackID, simTrack>
  std::map<unsigned int, SimTrack> simTracks;
  // Make a map of <trackId, SimTrack index in the tree>
  std::map<unsigned int, unsigned int> simTrackIdmap;
  for (edm::SimTrackContainer::const_iterator simTrackIt(simTracksRaw->begin()); simTrackIt != simTracksRaw->end();
       ++simTrackIt) {
    if (simTrackIt->momentum().pt() > simtrackminpt_) {
      simTracks.insert(std::pair<unsigned int, SimTrack>(simTrackIt->trackId(), *simTrackIt));
      hitInfo->SimTrack_pt.push_back(simTrackIt->momentum().pt());
      hitInfo->SimTrack_eta.push_back(simTrackIt->momentum().eta());
      hitInfo->SimTrack_phi.push_back(simTrackIt->momentum().phi());
      // Index seems to be 1-based, convert to 0 based
      //int genpartidx = simTrackIt->genpartIndex() - 1;
      //int genParticle_size = genParticles->size();
      //int pdgId = 0;
      //if ( !(genpartidx<genParticle_size) ){
      //  std::cout << "Warning! SimTrack genParticle index out of range!" << std::endl;
      //}
      //if ( (genpartidx>-1) && (genpartidx<genParticle_size) ){
      //  const reco::GenParticle& genpart = genParticles->at(genpartidx);
      //  double delta2 = pow(simTrackIt->momentum().pt()-genpart.pt(),2)+pow(simTrackIt->momentum().eta()-genpart.eta(),2)+pow(simTrackIt->momentum().phi()-genpart.phi(),2);
      //  if (delta2>0.01){
      //    std::cout << "Warning! SimTrack genParticle match seems wrong!" << std::endl;
      //  }
      //  else {
      //    pdgId = genpart.pdgId();
      //  }
      //}
      hitInfo->SimTrack_pdgId.push_back(simTrackIt->type());
      simTrackIdmap.insert(std::pair<unsigned int, unsigned int>(simTrackIt->trackId(), hitInfo->SimTrack_pt.size()-1));
      hitInfo->SimTrack_n_matched_hits.push_back(0);
      hitInfo->SimTrack_n_simhits.push_back(0);
    }
  }

  // Get the simhits
  std::vector<std::map<edm::PSimHitContainer::const_iterator, unsigned int>> simHitIdmaps;
  for (unsigned int simhitidx = 0; simhitidx < 2; ++simhitidx) {  // loop over both barrel and endcap hits
    std::map<edm::PSimHitContainer::const_iterator, unsigned int> hitidmap;
    for (edm::PSimHitContainer::const_iterator simhitIt(simHitsRaw[simhitidx]->begin());
         simhitIt != simHitsRaw[simhitidx]->end();
         ++simhitIt) {
      // Get the detector unit's id
      unsigned int rawid(simhitIt->detUnitId());
      DetId detId(rawid);
      unsigned int layer = (tTopo->side(detId) != 0) * 1000;  // don't split up endcap sides
      layer += tTopo->layer(detId);
      //if (!layer) {
      //  layer += tTopo->layer(detId);
      //} else {
      //  layer += (catECasRings_ ? tTopo->tidRing(detId) * 10 : tTopo->layer(detId));
      //}
      
      // determine the detector we are in
      TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
      unsigned int det = 0;
      if (mType == TrackerGeometry::ModuleType::Ph2PSP) {
        det = 1;
      } else if (mType == TrackerGeometry::ModuleType::Ph2PSS) {
        det = 2;
      } else if (mType == TrackerGeometry::ModuleType::Ph2SS) {
        det = 3;
      } else { // it seems the det == 0 means inner tracker
        //std::cout << "UNKNOWN DETECTOR TYPE!" << std::endl;
        continue;
      }

      // Get the geomdet
      const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
      if (!geomDetUnit)
        continue;
      Local3DPoint localPosSimHit = simhitIt->localPosition();
      Global3DPoint globalPosSimHit = geomDetUnit->surface().toGlobal(localPosSimHit);
      hitidmap.insert(std::pair<edm::PSimHitContainer::const_iterator, unsigned int>(simhitIt,hitInfo->SimHit_barrelendcap.size()));
      hitInfo->SimHit_barrelendcap.push_back(simhitidx);
      hitInfo->SimHit_layer.push_back(layer);
      hitInfo->SimHit_ModuleType.push_back(det);
      hitInfo->SimHit_local_x.push_back(localPosSimHit.x());
      hitInfo->SimHit_local_y.push_back(localPosSimHit.y());
      hitInfo->SimHit_local_z.push_back(localPosSimHit.z());
      hitInfo->SimHit_global_x.push_back(globalPosSimHit.x());
      hitInfo->SimHit_global_y.push_back(globalPosSimHit.y());
      hitInfo->SimHit_global_z.push_back(globalPosSimHit.z());
      hitInfo->SimHit_p.push_back(simhitIt->pabs());
      hitInfo->SimHit_tof.push_back(simhitIt->timeOfFlight());
      hitInfo->SimHit_E.push_back(simhitIt->energyLoss());
      hitInfo->SimHit_n_matched_Hit.push_back(0);
      hitInfo->SimHit_Hit_idx.push_back(-1);
      hitInfo->SimHit_pdgId.push_back(simhitIt->particleType());

      // Connect SimHits with SimTracks
      unsigned int tkId = simhitIt->trackId();
      if (simTracks.find(tkId)!=simTracks.end()){
      }
      std::map<unsigned int, unsigned int>::const_iterator istfind(simTrackIdmap.find(tkId));
      if (istfind != simTrackIdmap.end()) {
        hitInfo->SimHit_SimTrack_idx.push_back(simTrackIdmap[tkId]);
        hitInfo->SimTrack_n_simhits[simTrackIdmap[tkId]] += 1;
      }
      else {
        hitInfo->SimHit_SimTrack_idx.push_back(-1);
      }
    }
    simHitIdmaps.push_back(hitidmap);
  }

  for (Phase2TrackerRecHit1DCollectionNew::const_iterator DSViter = rechits->begin(); DSViter != rechits->end(); ++DSViter) {
    // Get the detector unit's id
    unsigned int rawid(DSViter->detId());
    DetId detId(rawid);
    unsigned int layer = (tTopo->side(detId) != 0) * 1000;  // don't split up endcap sides
    layer += tTopo->layer(detId);
    //if (!layer) {
    //  layer += tTopo->layer(detId);
    //} else {
    //  layer += (catECasRings_ ? tTopo->tidRing(detId) * 10 : tTopo->layer(detId));
    //}
    
    // determine the detector we are in
    TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
    unsigned int det = 0;
    if (mType == TrackerGeometry::ModuleType::Ph2PSP) {
      det = 1;
    } else if (mType == TrackerGeometry::ModuleType::Ph2PSS) {
      det = 2;
    } else if (mType == TrackerGeometry::ModuleType::Ph2SS) {
      det = 3;
    } else {
      std::cout << "UNKNOWN DETECTOR TYPE!" << std::endl;
    }

    // Get the geomdet
    const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
    if (!geomDetUnit)
      continue;

    // Loop over the rechits in the detector unit
    for (edmNew::DetSet<Phase2TrackerRecHit1D>::const_iterator rechitIt = DSViter->begin(); rechitIt != DSViter->end(); ++rechitIt) {
      // determine the position
      LocalPoint localPosClu = rechitIt->localPosition();
      Global3DPoint globalPosClu = geomDetUnit->surface().toGlobal(localPosClu);

      // Get the cluster from the rechit
      const Phase2TrackerCluster1D* clustIt = &*rechitIt->cluster();

      // Get all the simTracks that form the cluster
      std::vector<unsigned int> clusterSimTrackIds;
      for (unsigned int i(0); i < clustIt->size(); ++i) {
        unsigned int channel(Phase2TrackerDigi::pixelToChannel(clustIt->firstRow() + i, clustIt->column()));
        std::vector<unsigned int> simTrackIds_unselected(getSimTrackId(pixelSimLinks, detId, channel));
        std::vector<unsigned int> simTrackIds;
        for (auto istId : simTrackIds_unselected) {
          std::map<unsigned int, SimTrack>::const_iterator istfind(simTracks.find(istId));
          if (istfind != simTracks.end())
            simTrackIds.push_back(istId);
        }
        for (unsigned int i = 0; i < simTrackIds.size(); ++i) {
          bool add = true;
          for (unsigned int j = 0; j < clusterSimTrackIds.size(); ++j) {
            // only save simtrackids that are not present yet
            if (simTrackIds.at(i) == clusterSimTrackIds.at(j))
              add = false;
          }
          if (add)
            clusterSimTrackIds.push_back(simTrackIds.at(i));
        }
      }

      // pick the leading sim track in case there are more than 1 by taking the one with highest pT
      float simtkpt_max = -1;
      unsigned int matched_simtrack_idx = -1;
      for (unsigned int isimtk : clusterSimTrackIds) {
        float simtkpt = -1;
        std::map<unsigned int, SimTrack>::const_iterator istfind(simTracks.find(isimtk));
        if (istfind != simTracks.end()) {
          hitInfo->SimTrack_n_matched_hits[simTrackIdmap[isimtk]] += 1;
          simtkpt = simTracks[isimtk].momentum().pt();
          if (simtkpt>simtkpt_max) {
            simtkpt_max = simtkpt;
            matched_simtrack_idx = simTrackIdmap[isimtk];
          }
        }
      }
      hitInfo->Hit_SimTrack_idx.push_back(matched_simtrack_idx);

      // find the closest simhit
      // this is needed because otherwise you get cases with simhits and clusters being swapped
      // when there are more than 1 cluster with common simtrackids
      const PSimHit* simhit = 0;  // bad naming to avoid changing code below. This is the closest simhit in x
      float minx = 10000;
      int matched_simhitid = -1;
      for (unsigned int simhitidx = 0; simhitidx < 2; ++simhitidx) {  // loop over both barrel and endcap hits
        for (edm::PSimHitContainer::const_iterator simhitIt(simHitsRaw[simhitidx]->begin());
             simhitIt != simHitsRaw[simhitidx]->end();
             ++simhitIt) {
          // only check SimHits that passes selections (in Outer tracker)
          if (simHitIdmaps[simhitidx].find(simhitIt)==simHitIdmaps[simhitidx].end()) {
            continue;
          }
          // check SimHit detId is the same with the RecHit
          if (rawid == simhitIt->detUnitId()) {
            //std::cout << "=== " << rawid << " " << &*simhitIt << " " << simhitIt->trackId() << " " << simhitIt->localPosition().x() << " " << simhitIt->localPosition().y() << std::endl;
            auto it = std::lower_bound(clusterSimTrackIds.begin(), clusterSimTrackIds.end(), simhitIt->trackId());
            // check SimHit track id is included in the cluster
            if (it != clusterSimTrackIds.end() && *it == simhitIt->trackId()) {
              if (!simhit || fabs(simhitIt->localPosition().x() - localPosClu.x()) < minx) {
                minx = fabs(simhitIt->localPosition().x() - localPosClu.x());
                simhit = &*simhitIt;
                matched_simhitid = simHitIdmaps[simhitidx][simhitIt];
              }
            }
          }
        }
      }

      if (matched_simhitid>=0) {
        hitInfo->SimHit_n_matched_Hit[matched_simhitid] += 1;
        if (hitInfo->SimHit_Hit_idx[matched_simhitid]==-1)
          hitInfo->SimHit_Hit_idx[matched_simhitid] = hitInfo->Hit_cluster_global_x.size();
      }

      hitInfo->Hit_cluster_global_x.push_back(globalPosClu.x());
      hitInfo->Hit_cluster_global_y.push_back(globalPosClu.y());
      hitInfo->Hit_cluster_global_z.push_back(globalPosClu.z());
      hitInfo->Hit_layer.push_back(layer);
      hitInfo->Hit_ModuleType.push_back(det);
      hitInfo->Hit_cluster_size.push_back(clustIt->size());
      hitInfo->Hit_cluster_SimTrack_size.push_back(clusterSimTrackIds.size());
      hitInfo->Hit_cluster_local_x.push_back(localPosClu.x());
      hitInfo->Hit_cluster_local_y.push_back(localPosClu.y());
      hitInfo->Hit_cluster_local_z.push_back(localPosClu.z());
      if (!simhit){
        hitInfo->Hit_cluster_haveSimHit.push_back(false);
        hitInfo->Hit_cluster_closestSimHit_local_x.push_back(999);
        hitInfo->Hit_cluster_closestSimHit_local_y.push_back(999);
        hitInfo->Hit_cluster_closestSimHit_local_z.push_back(999);
        hitInfo->Hit_cluster_closestSimHit_p.push_back(-1);
        hitInfo->Hit_cluster_closestSimHit_tof.push_back(999);
        hitInfo->Hit_cluster_closestSimHit_E.push_back(-1);
        hitInfo->Hit_SimHit_idx.push_back(-1);
      }
      else {
        hitInfo->Hit_cluster_haveSimHit.push_back(true);
        hitInfo->Hit_cluster_closestSimHit_local_x.push_back(simhit->localPosition().x());
        hitInfo->Hit_cluster_closestSimHit_local_y.push_back(simhit->localPosition().y());
        hitInfo->Hit_cluster_closestSimHit_local_z.push_back(simhit->localPosition().z());
        hitInfo->Hit_cluster_closestSimHit_p.push_back(simhit->pabs());
        hitInfo->Hit_cluster_closestSimHit_tof.push_back(simhit->timeOfFlight());
        hitInfo->Hit_cluster_closestSimHit_E.push_back(simhit->energyLoss());
        hitInfo->Hit_SimHit_idx.push_back(matched_simhitid);
      }
    }
  }
  hitTree->Fill();
}

void RecHitTree::beginJob()
{
  edm::Service<TFileService> fs;
  hitTree = fs->make<TTree>( "HitTree", "HitTree" );

  hitTree->Branch("Hit_cluster_global_x",                   &hitInfo->Hit_cluster_global_x);
  hitTree->Branch("Hit_cluster_global_y",                   &hitInfo->Hit_cluster_global_y);
  hitTree->Branch("Hit_cluster_global_z",                   &hitInfo->Hit_cluster_global_z);
  hitTree->Branch("Hit_layer",                              &hitInfo->Hit_layer);
  hitTree->Branch("Hit_ModuleType",                         &hitInfo->Hit_ModuleType);
  hitTree->Branch("Hit_cluster_size",                       &hitInfo->Hit_cluster_size);
  hitTree->Branch("Hit_cluster_SimTrack_size",              &hitInfo->Hit_cluster_SimTrack_size);
  hitTree->Branch("Hit_cluster_local_x",                    &hitInfo->Hit_cluster_local_x);
  hitTree->Branch("Hit_cluster_local_y",                    &hitInfo->Hit_cluster_local_y);
  hitTree->Branch("Hit_cluster_local_z",                    &hitInfo->Hit_cluster_local_z);
  hitTree->Branch("Hit_cluster_haveSimHit",                 &hitInfo->Hit_cluster_haveSimHit);
  hitTree->Branch("Hit_cluster_closestSimHit_local_x",      &hitInfo->Hit_cluster_closestSimHit_local_x);
  hitTree->Branch("Hit_cluster_closestSimHit_local_y",      &hitInfo->Hit_cluster_closestSimHit_local_y);
  hitTree->Branch("Hit_cluster_closestSimHit_local_z",      &hitInfo->Hit_cluster_closestSimHit_local_z);
  hitTree->Branch("Hit_cluster_closestSimHit_p",            &hitInfo->Hit_cluster_closestSimHit_p);
  hitTree->Branch("Hit_cluster_closestSimHit_tof",          &hitInfo->Hit_cluster_closestSimHit_tof);
  hitTree->Branch("Hit_cluster_closestSimHit_E",            &hitInfo->Hit_cluster_closestSimHit_E);
  hitTree->Branch("Hit_SimTrack_idx",                       &hitInfo->Hit_SimTrack_idx);
  hitTree->Branch("Hit_SimHit_idx",                         &hitInfo->Hit_SimHit_idx);
  hitTree->Branch("SimHit_barrelendcap",                    &hitInfo->SimHit_barrelendcap);
  hitTree->Branch("SimHit_layer",                           &hitInfo->SimHit_layer);
  hitTree->Branch("SimHit_ModuleType",                      &hitInfo->SimHit_ModuleType);
  hitTree->Branch("SimHit_local_x",                         &hitInfo->SimHit_local_x);
  hitTree->Branch("SimHit_local_y",                         &hitInfo->SimHit_local_y);
  hitTree->Branch("SimHit_local_z",                         &hitInfo->SimHit_local_z);
  hitTree->Branch("SimHit_global_x",                        &hitInfo->SimHit_global_x);
  hitTree->Branch("SimHit_global_y",                        &hitInfo->SimHit_global_y);
  hitTree->Branch("SimHit_global_z",                        &hitInfo->SimHit_global_z);
  hitTree->Branch("SimHit_p",                               &hitInfo->SimHit_p);
  hitTree->Branch("SimHit_tof",                             &hitInfo->SimHit_tof);
  hitTree->Branch("SimHit_E",                               &hitInfo->SimHit_E);
  hitTree->Branch("SimHit_SimTrack_idx",                    &hitInfo->SimHit_SimTrack_idx);
  hitTree->Branch("SimHit_n_matched_Hit",                   &hitInfo->SimHit_n_matched_Hit);
  hitTree->Branch("SimHit_Hit_idx",                         &hitInfo->SimHit_Hit_idx);
  hitTree->Branch("SimHit_pdgId",                           &hitInfo->SimHit_pdgId);
  hitTree->Branch("SimTrack_pt",                            &hitInfo->SimTrack_pt);
  hitTree->Branch("SimTrack_eta",                           &hitInfo->SimTrack_eta);
  hitTree->Branch("SimTrack_phi",                           &hitInfo->SimTrack_phi);
  hitTree->Branch("SimTrack_pdgId",                         &hitInfo->SimTrack_pdgId);
  hitTree->Branch("SimTrack_n_matched_hits",                &hitInfo->SimTrack_n_matched_hits);
  hitTree->Branch("SimTrack_n_simhits",                     &hitInfo->SimTrack_n_simhits);


}

void RecHitTree::endJob()
{}

void RecHitTree::initEventStructure()
{
  hitInfo->Hit_cluster_global_x.clear();
  hitInfo->Hit_cluster_global_y.clear();
  hitInfo->Hit_cluster_global_z.clear();
  hitInfo->Hit_layer.clear();
  hitInfo->Hit_ModuleType.clear();
  hitInfo->Hit_cluster_size.clear();
  hitInfo->Hit_cluster_SimTrack_size.clear();
  hitInfo->Hit_cluster_local_x.clear();
  hitInfo->Hit_cluster_local_y.clear();
  hitInfo->Hit_cluster_local_z.clear();
  hitInfo->Hit_cluster_haveSimHit.clear();
  hitInfo->Hit_cluster_closestSimHit_local_x.clear();
  hitInfo->Hit_cluster_closestSimHit_local_y.clear();
  hitInfo->Hit_cluster_closestSimHit_local_z.clear();
  hitInfo->Hit_cluster_closestSimHit_p.clear();
  hitInfo->Hit_cluster_closestSimHit_tof.clear();
  hitInfo->Hit_cluster_closestSimHit_E.clear();
  hitInfo->Hit_SimTrack_idx.clear();
  hitInfo->Hit_SimHit_idx.clear();
  hitInfo->SimHit_barrelendcap.clear(); 
  hitInfo->SimHit_layer.clear();
  hitInfo->SimHit_ModuleType.clear();
  hitInfo->SimHit_local_x.clear();
  hitInfo->SimHit_local_y.clear();
  hitInfo->SimHit_local_z.clear();
  hitInfo->SimHit_global_x.clear();
  hitInfo->SimHit_global_y.clear();
  hitInfo->SimHit_global_z.clear();
  hitInfo->SimHit_p.clear();
  hitInfo->SimHit_tof.clear();
  hitInfo->SimHit_E.clear();
  hitInfo->SimHit_SimTrack_idx.clear();
  hitInfo->SimHit_n_matched_Hit.clear();
  hitInfo->SimHit_Hit_idx.clear();
  hitInfo->SimHit_pdgId.clear();
  hitInfo->SimTrack_pt.clear();
  hitInfo->SimTrack_eta.clear();
  hitInfo->SimTrack_phi.clear();
  hitInfo->SimTrack_pdgId.clear();
  hitInfo->SimTrack_n_matched_hits.clear();
  hitInfo->SimTrack_n_simhits.clear();
}

std::vector<unsigned int> RecHitTree::getSimTrackId(
    const edm::Handle<edm::DetSetVector<PixelDigiSimLink> >& pixelSimLinks, const DetId& detId, unsigned int channel) {
  std::vector<unsigned int> retvec;
  edm::DetSetVector<PixelDigiSimLink>::const_iterator DSViter(pixelSimLinks->find(detId));
  if (DSViter == pixelSimLinks->end())
    return retvec;
  for (edm::DetSet<PixelDigiSimLink>::const_iterator it = DSViter->data.begin(); it != DSViter->data.end(); ++it) {
    if (channel == it->channel()) {
      retvec.push_back(it->SimTrackId());
    }
  }
  return retvec;
}

DEFINE_FWK_MODULE(RecHitTree);
