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

#include "TTree.h"

struct HitInfo
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
};

struct SimHitInfo
{
  std::vector<float> local_x;
  std::vector<float> local_y;
  std::vector<float> local_z;
  std::vector<float> global_x;
  std::vector<float> global_y;
  std::vector<float> global_z;
  std::vector<float> local_dx;
  std::vector<float> local_dy;
  std::vector<float> local_dz;
  std::vector<float> global_dx;
  std::vector<float> global_dy;
  std::vector<float> global_dz;
  std::vector<float> theta;
  std::vector<float> phi;
  std::vector<float> pabs;
  std::vector<float> tof;
  std::vector<float> energyLoss;
  std::vector<unsigned short> processType;
  std::vector<int>   particleType;
  //std::vector<unsigned int> detUnitId;
  std::vector<unsigned short> layer;
  std::vector<unsigned short> moduleType;
  std::vector<unsigned int> trackId;
  std::vector<unsigned int> originalTrackId;
};

struct SimTrackInfo
{
  std::vector<float> SimTrack_xTk;
  std::vector<float> SimTrack_yTk;
  std::vector<float> SimTrack_zTk;
  std::vector<float> SimTrack_ptTk;
  std::vector<float> SimTrack_etaTk;
  std::vector<float> SimTrack_phiTk;
  std::vector<float> SimTrack_pt;
  std::vector<float> SimTrack_eta;
  std::vector<float> SimTrack_phi;
  std::vector<int> SimTrack_type;
  std::vector<int> SimTrack_charge;
  std::vector<int> SimTrack_trackInfo;
};

class RecHitTreeWA : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit RecHitTreeWA(const edm::ParameterSet&);
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

    const double simtrackminpt_;

    TTree* hitTree;
    HitInfo* hitInfo;
    TTree* simHitTree;
    SimHitInfo* simHitInfo;
    TTree* simTrackTree;
    SimTrackInfo* simTrackInfo;
};

RecHitTreeWA::RecHitTreeWA(const edm::ParameterSet& cfg)
  : esTokenGeom_(esConsumes()),
    esTokenTopo_(esConsumes()),
    tokenRecHits_(consumes<Phase2TrackerRecHit1DCollectionNew>(cfg.getParameter<edm::InputTag>("rechits"))),
    tokenClusters_(consumes<Phase2TrackerCluster1DCollectionNew>(cfg.getParameter<edm::InputTag>("clusters"))),
    tokenLinks_(consumes<edm::DetSetVector<PixelDigiSimLink> >(cfg.getParameter<edm::InputTag>("links"))),
    tokenSimHitsB_(consumes<edm::PSimHitContainer>(cfg.getParameter<edm::InputTag>("simhitsbarrel"))),
    tokenSimHitsE_(consumes<edm::PSimHitContainer>(cfg.getParameter<edm::InputTag>("simhitsendcap"))),
    tokenSimTracks_(consumes<edm::SimTrackContainer>(cfg.getParameter<edm::InputTag>("simtracks"))),
    simtrackminpt_(cfg.getParameter<double>("SimTrackMinPt"))
{
  hitInfo = new HitInfo;
  simHitInfo = new SimHitInfo;
  simTrackInfo = new SimTrackInfo;
}

void RecHitTreeWA::analyze(const edm::Event& event, const edm::EventSetup& eventSetup) {
  initEventStructure();

  // Get the geometry
  const TrackerGeometry* tkGeom = &eventSetup.getData(esTokenGeom_);
  const TrackerTopology* tTopo = &eventSetup.getData(esTokenTopo_);

  // Get the RecHits
  edm::Handle<Phase2TrackerRecHit1DCollectionNew> rechits;
  event.getByToken(tokenRecHits_, rechits);
  //std::cout << "Phase2TrackerRecHit1DCollectionNew: " << rechits->size() << std::endl;
  
  // Get the Clusters
  edm::Handle<Phase2TrackerCluster1DCollectionNew> clusters;
  event.getByToken(tokenClusters_, clusters);
  if ( clusters->size() != rechits->size() ) {
    std::cout << "*** Different sizes for clusters (" << clusters->size()
	      << ") and rechits (" << rechits->size() << ")" << std::endl;
  }

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

  // Rearrange the simTracks for ease of use <simTrackID, simTrack>
  std::map<unsigned int, SimTrack> simTracks;
  for (edm::SimTrackContainer::const_iterator simTrackIt(simTracksRaw->begin()); simTrackIt != simTracksRaw->end();
       ++simTrackIt) {
    if (simTrackIt->momentum().pt() > simtrackminpt_) {
      simTracks.insert(std::pair<unsigned int, SimTrack>(simTrackIt->trackId(), *simTrackIt));
    }
    simTrackInfo->SimTrack_xTk.push_back(simTrackIt->trackerSurfacePosition().x());
    simTrackInfo->SimTrack_yTk.push_back(simTrackIt->trackerSurfacePosition().y());
    simTrackInfo->SimTrack_zTk.push_back(simTrackIt->trackerSurfacePosition().z());
    simTrackInfo->SimTrack_ptTk.push_back(simTrackIt->trackerSurfaceMomentum().pt());
    simTrackInfo->SimTrack_etaTk.push_back(simTrackIt->trackerSurfaceMomentum().eta());
    simTrackInfo->SimTrack_phiTk.push_back(simTrackIt->trackerSurfaceMomentum().phi());
    simTrackInfo->SimTrack_pt.push_back(simTrackIt->momentum().pt());
    simTrackInfo->SimTrack_eta.push_back(simTrackIt->momentum().eta());
    simTrackInfo->SimTrack_phi.push_back(simTrackIt->momentum().phi());
    simTrackInfo->SimTrack_type.push_back(simTrackIt->type());
    simTrackInfo->SimTrack_charge.push_back(simTrackIt->charge());
    simTrackInfo->SimTrack_trackInfo.push_back(simTrackIt->getTrackInfo());
  }
  simTrackTree->Fill();

  for (unsigned int simhitidx = 0; simhitidx < 2; ++simhitidx) {  // loop over both barrel and endcap hits
    for (edm::PSimHitContainer::const_iterator simhitIt(simHitsRaw[simhitidx]->begin());
	 simhitIt != simHitsRaw[simhitidx]->end(); ++simhitIt) {
          // Get the detector unit's id
      DetId detId(simhitIt->detUnitId());
      unsigned int layer = (tTopo->side(detId) != 0) * 1000;  // don't split up endcap sides
      layer += tTopo->layer(detId);    
      hitInfo->Hit_ModuleType.push_back((unsigned short)(tkGeom->getDetectorType(detId)));
      // Get the geomdet
      const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
      if (!geomDetUnit) {
	std::cout << "*** did not find geomDetUnit ***" << std::endl;
	continue;
      }
      GlobalPoint globalPosition(geomDetUnit->toGlobal(simhitIt->localPosition()));
      GlobalVector globalDirection(geomDetUnit->toGlobal(simhitIt->localDirection()));

      simHitInfo->local_x.push_back(simhitIt->localPosition().x());
      simHitInfo->local_y.push_back(simhitIt->localPosition().y());
      simHitInfo->local_z.push_back(simhitIt->localPosition().z());
      simHitInfo->global_x.push_back(globalPosition.x());
      simHitInfo->global_y.push_back(globalPosition.y());
      simHitInfo->global_z.push_back(globalPosition.z());
      simHitInfo->local_dx.push_back(simhitIt->localDirection().x());
      simHitInfo->local_dy.push_back(simhitIt->localDirection().y());
      simHitInfo->local_dz.push_back(simhitIt->localDirection().z());
      simHitInfo->global_dx.push_back(globalDirection.x());
      simHitInfo->global_dy.push_back(globalDirection.y());
      simHitInfo->global_dz.push_back(globalDirection.z());
      simHitInfo->theta.push_back(simhitIt->thetaAtEntry());
      simHitInfo->phi.push_back(simhitIt->phiAtEntry());
      simHitInfo->pabs.push_back(simhitIt->pabs());
      simHitInfo->tof.push_back(simhitIt->timeOfFlight());
      simHitInfo->energyLoss.push_back(simhitIt->energyLoss());
      simHitInfo->processType.push_back(simhitIt->processType());
      simHitInfo->particleType.push_back(simhitIt->particleType());
      //simHitInfo->detUnitId.push_back(simhitIt->detUnitId());
      simHitInfo->layer.push_back(layer);
      simHitInfo->moduleType.push_back((unsigned short)(tkGeom->getDetectorType(detId)));
      simHitInfo->trackId.push_back(simhitIt->trackId());
      simHitInfo->originalTrackId.push_back(simhitIt->originalTrackId());
    }
  }
  simHitTree->Fill();

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
    //TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
    //unsigned int det = 0;
    //if (mType == TrackerGeometry::ModuleType::Ph2PSP) {
    //  det = 1;
    //} else if (mType == TrackerGeometry::ModuleType::Ph2PSS || mType == TrackerGeometry::ModuleType::Ph2SS) {
    //  det = 2;
    //} else {
    //std::cout << "UNKNOWN DETECTOR TYPE!" << std::endl;
    //}
    // std::cout << "Rawid " << rawid << " module type " << (int)mType << std::endl;

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
      //std::cout << "  cluster size " << clustIt->size() << " firstStrip " << clustIt->firstStrip() << " firstRow " << clustIt->firstRow() << " edge " << clustIt->edge()
      //   << " column " << clustIt->column() << " threshold " << clustIt->threshold() << std::endl;

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

      // find the closest simhit
      // this is needed because otherwise you get cases with simhits and clusters being swapped
      // when there are more than 1 cluster with common simtrackids
      const PSimHit* simhit = 0;  // bad naming to avoid changing code below. This is the closest simhit in x
      float minx = 10000;
      for (unsigned int simhitidx = 0; simhitidx < 2; ++simhitidx) {  // loop over both barrel and endcap hits
        for (edm::PSimHitContainer::const_iterator simhitIt(simHitsRaw[simhitidx]->begin());
             simhitIt != simHitsRaw[simhitidx]->end();
             ++simhitIt) {
          // check SimHit detId is the same with the RecHit
          if (rawid == simhitIt->detUnitId()) {
            //std::cout << "=== " << rawid << " " << &*simhitIt << " " << simhitIt->trackId() << " " << simhitIt->localPosition().x() << " " << simhitIt->localPosition().y() << std::endl;
            auto it = std::lower_bound(clusterSimTrackIds.begin(), clusterSimTrackIds.end(), simhitIt->trackId());
            // check SimHit track id is included in the cluster
            if (it != clusterSimTrackIds.end() && *it == simhitIt->trackId()) {
              if (!simhit || fabs(simhitIt->localPosition().x() - localPosClu.x()) < minx) {
                minx = fabs(simhitIt->localPosition().x() - localPosClu.x());
                simhit = &*simhitIt;
              }
            }
          }
        }
      }

      hitInfo->Hit_cluster_global_x.push_back(globalPosClu.x());
      hitInfo->Hit_cluster_global_y.push_back(globalPosClu.y());
      hitInfo->Hit_cluster_global_z.push_back(globalPosClu.z());
      hitInfo->Hit_layer.push_back(layer);
      hitInfo->Hit_ModuleType.push_back((unsigned short)(tkGeom->getDetectorType(detId)));
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
      }
      else {
        hitInfo->Hit_cluster_haveSimHit.push_back(true);
        hitInfo->Hit_cluster_closestSimHit_local_x.push_back(simhit->localPosition().x());
        hitInfo->Hit_cluster_closestSimHit_local_y.push_back(simhit->localPosition().y());
        hitInfo->Hit_cluster_closestSimHit_local_z.push_back(simhit->localPosition().z());
      }
      hitInfo->Hit_det_rawid.push_back(rawid);
      hitInfo->Hit_cluster_firstStrip.push_back(clustIt->firstStrip());
      hitInfo->Hit_cluster_firstRow.push_back(clustIt->firstRow());
      hitInfo->Hit_cluster_column.push_back(clustIt->column());
      hitInfo->Hit_cluster_edge.push_back(clustIt->edge());
      hitInfo->Hit_cluster_threshold.push_back(clustIt->threshold());
    }
  }
  hitTree->Fill();
}

void RecHitTreeWA::beginJob()
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
  hitTree->Branch("Hit_det_rawid",                          &hitInfo->Hit_det_rawid);
  hitTree->Branch("Hit_cluster_firstStrip",                 &hitInfo->Hit_cluster_firstStrip);
  hitTree->Branch("Hit_cluster_firstRow",                   &hitInfo->Hit_cluster_firstRow);
  hitTree->Branch("Hit_cluster_column",                     &hitInfo->Hit_cluster_column);
  hitTree->Branch("Hit_cluster_edge",                       &hitInfo->Hit_cluster_edge);
  hitTree->Branch("Hit_cluster_threshold",                  &hitInfo->Hit_cluster_threshold);

  simHitTree = fs->make<TTree>( "SimHitTree", "SimHitTree" );
  simHitTree->Branch("local_x",	&simHitInfo->local_x);
  simHitTree->Branch("local_y",	&simHitInfo->local_y);
  simHitTree->Branch("local_z",	&simHitInfo->local_z);
  simHitTree->Branch("global_x",	&simHitInfo->global_x);
  simHitTree->Branch("global_y",	&simHitInfo->global_y);
  simHitTree->Branch("global_z",	&simHitInfo->global_z);
  simHitTree->Branch("local_dx",	&simHitInfo->local_dx);
  simHitTree->Branch("local_dy",	&simHitInfo->local_dy);
  simHitTree->Branch("local_dz",	&simHitInfo->local_dz);
  simHitTree->Branch("global_dx",	&simHitInfo->global_dx);
  simHitTree->Branch("global_dy",	&simHitInfo->global_dy);
  simHitTree->Branch("global_dz",	&simHitInfo->global_dz);
  simHitTree->Branch("theta",	&simHitInfo->theta);
  simHitTree->Branch("phi",	&simHitInfo->phi);
  simHitTree->Branch("pabs",	&simHitInfo->pabs);
  simHitTree->Branch("tof",	&simHitInfo->tof);
  simHitTree->Branch("energyLoss",	&simHitInfo->energyLoss);
  simHitTree->Branch("processType",	&simHitInfo->processType);
  simHitTree->Branch("particleType",	&simHitInfo->particleType);
  //simHitTree->Branch("detUnitId",	&simHitInfo->detUnitId);
  simHitTree->Branch("layer",   &simHitInfo->layer);
  simHitTree->Branch("moduleType",      &simHitInfo->moduleType);
  simHitTree->Branch("trackId",	&simHitInfo->trackId);
  simHitTree->Branch("originalTrackId",	&simHitInfo->originalTrackId);

  
  simTrackTree = fs->make<TTree>( "SimTrackTree", "SimTrackTree" );
  simTrackTree->Branch("SimTrack_xTk",         &simTrackInfo->SimTrack_xTk);
  simTrackTree->Branch("SimTrack_yTk",         &simTrackInfo->SimTrack_yTk);
  simTrackTree->Branch("SimTrack_zTk",         &simTrackInfo->SimTrack_zTk);
  simTrackTree->Branch("SimTrack_ptTk",        &simTrackInfo->SimTrack_ptTk);
  simTrackTree->Branch("SimTrack_etaTk",       &simTrackInfo->SimTrack_etaTk);
  simTrackTree->Branch("SimTrack_phiTk",       &simTrackInfo->SimTrack_phiTk);
  simTrackTree->Branch("SimTrack_pt",        &simTrackInfo->SimTrack_pt);
  simTrackTree->Branch("SimTrack_eta",       &simTrackInfo->SimTrack_eta);
  simTrackTree->Branch("SimTrack_phi",       &simTrackInfo->SimTrack_phi);
  simTrackTree->Branch("SimTrack_type",      &simTrackInfo->SimTrack_type);
  simTrackTree->Branch("SimTrack_charge",    &simTrackInfo->SimTrack_charge);
  simTrackTree->Branch("SimTrack_trackInfo", &simTrackInfo->SimTrack_trackInfo);

}

void RecHitTreeWA::endJob()
{}

void RecHitTreeWA::initEventStructure()
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
  hitInfo->Hit_det_rawid.clear();
  hitInfo->Hit_cluster_firstStrip.clear();
  hitInfo->Hit_cluster_firstRow.clear();
  hitInfo->Hit_cluster_column.clear();
  hitInfo->Hit_cluster_edge.clear();
  hitInfo->Hit_cluster_threshold.clear();

  simHitInfo->local_x.clear();
  simHitInfo->local_y.clear();
  simHitInfo->local_z.clear();
  simHitInfo->global_x.clear();
  simHitInfo->global_y.clear();
  simHitInfo->global_z.clear();
  simHitInfo->local_dx.clear();
  simHitInfo->local_dy.clear();
  simHitInfo->local_dz.clear();
  simHitInfo->global_dx.clear();
  simHitInfo->global_dy.clear();
  simHitInfo->global_dz.clear();
  simHitInfo->theta.clear();
  simHitInfo->phi.clear();
  simHitInfo->pabs.clear();
  simHitInfo->tof.clear();
  simHitInfo->energyLoss.clear();
  simHitInfo->processType.clear();
  simHitInfo->particleType.clear();
  // simHitInfo->detUnitId.clear();
  simHitInfo->layer.clear();
  simHitInfo->moduleType.clear();
  simHitInfo->trackId.clear();
  simHitInfo->originalTrackId.clear();

  simTrackInfo->SimTrack_xTk.clear();
  simTrackInfo->SimTrack_yTk.clear();
  simTrackInfo->SimTrack_zTk.clear();
  simTrackInfo->SimTrack_ptTk.clear();
  simTrackInfo->SimTrack_etaTk.clear();
  simTrackInfo->SimTrack_phiTk.clear();
  simTrackInfo->SimTrack_pt.clear();
  simTrackInfo->SimTrack_eta.clear();
  simTrackInfo->SimTrack_phi.clear();
  simTrackInfo->SimTrack_type.clear();
  simTrackInfo->SimTrack_charge.clear();
  simTrackInfo->SimTrack_trackInfo.clear();

}

std::vector<unsigned int> RecHitTreeWA::getSimTrackId(
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

DEFINE_FWK_MODULE(RecHitTreeWA);
