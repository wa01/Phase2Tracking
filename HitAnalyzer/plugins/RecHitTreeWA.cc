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

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "UserAnalysis/HitAnalyzer/interface/SimTrackInfo.h"
#include "UserAnalysis/HitAnalyzer/interface/SimHitInfo.h"
#include "UserAnalysis/HitAnalyzer/interface/RecHitInfo.h"

#include "TTree.h"

class RecHitTreeWA : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit RecHitTreeWA(const edm::ParameterSet&);
    void analyze(const edm::Event&, const edm::EventSetup&);

  private:
    virtual void beginJob() override;
    virtual void endJob() override;
    void initEventStructure();
    std::vector<unsigned int> getSimTrackId(edm::Handle<edm::DetSetVector<PixelDigiSimLink> >&,
					    const DetId&, unsigned int);


    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> esTokenGeom_;
    const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> esTokenTopo_;
    const edm::EDGetTokenT<Phase2TrackerRecHit1DCollectionNew> tokenRecHits_;
    const edm::EDGetTokenT<Phase2TrackerCluster1DCollectionNew> tokenClusters_;
    const edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink> > tokenLinks_;
    const edm::EDGetTokenT<edm::PSimHitContainer> tokenSimHitsB_;
    const edm::EDGetTokenT<edm::PSimHitContainer> tokenSimHitsE_;
    const edm::EDGetTokenT<edm::SimTrackContainer> tokenSimTracks_;

    const double simtrackminpt_;

    SimHitInfo simHitInfo_;
    RecHitInfo recHitInfo_;
    SimTrackInfo simTrackInfo_;
  
    TTree* recHitTree;
    TTree* simTrackTree;
    TTree* simHitTree;

    bool debugHitMatch_;
};

std::vector<unsigned int> RecHitTreeWA::getSimTrackId(edm::Handle<edm::DetSetVector<PixelDigiSimLink> >& pixelSimLinks,
						      const DetId& detId, unsigned int channel) {
  std::vector<unsigned int> retvec;
  edm::DetSetVector<PixelDigiSimLink>::const_iterator DSViter(pixelSimLinks->find(detId));
  if (DSViter == pixelSimLinks->end())  return retvec;
  for (edm::DetSet<PixelDigiSimLink>::const_iterator it = DSViter->data.begin(); it != DSViter->data.end(); ++it) {
    if (channel == it->channel()) {
      retvec.push_back(it->SimTrackId());
    }
  }
  return retvec;
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
    simtrackminpt_(cfg.getParameter<double>("SimTrackMinPt")),
    debugHitMatch_(cfg.getParameter<bool>("debugHitMatch"))
{
  //recHitInfo_ = new RecHitInfo;
  //simTrackInfo = new SimTrackInfo;
  //simHitInfo_ = new SimHitInfo();
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
  // std::cout << "#simTracks " << simTracksRaw.product()->size() << std::endl;

  // Rearrange the simTracks for ease of use <simTrackID, simTrack>
  std::map<unsigned int, SimTrack> simTracks;
  for (edm::SimTrackContainer::const_iterator simTrackIt(simTracksRaw->begin()); simTrackIt != simTracksRaw->end();
       ++simTrackIt) {
    // std::cout << simTrackIt->type() << "  " << simTrackIt->trackerSurfaceMomentum().pt() << " "
    // 	      << simTrackIt->momentum().pt() << std::endl;
    if (simTrackIt->momentum().pt() > simtrackminpt_) {
      simTracks.insert(std::pair<unsigned int, SimTrack>(simTrackIt->trackId(), *simTrackIt));
    }
    simTrackInfo_.fillSimTrackInfo(*simTrackIt);
  }
  simTrackTree->Fill();

  simHitInfo_.setupEvent(tTopo,tkGeom,pixelSimLinks.product(),simHitsRaw,*rechits.product(),&simTracks);
  for (unsigned int simhitidx = 0; simhitidx < 2; ++simhitidx) {  // loop over both barrel and endcap hits
    for (edm::PSimHitContainer::const_iterator simhitIt(simHitsRaw[simhitidx]->begin());
	 simhitIt != simHitsRaw[simhitidx]->end(); ++simhitIt) {
      // Get the detector unit's id
      DetId detId(simhitIt->detUnitId());
      unsigned int layer = (tTopo->side(detId) != 0) * 1000;  // don't split up endcap sides
      layer += tTopo->layer(detId);    
      //hitInfo->Hit_ModuleType.push_back((unsigned short)(tkGeom->getDetectorType(detId)));
      TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
      // Restrict to Phase2 OT
      if ( mType!=TrackerGeometry::ModuleType::Ph2PSP && 
	   mType!=TrackerGeometry::ModuleType::Ph2PSS &&
	   mType!=TrackerGeometry::ModuleType::Ph2SS )  continue;

      // Get the geomdet
      const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
      if (!geomDetUnit) {
	std::cout << "*** did not find geomDetUnit ***" << std::endl;
	continue;
      }

      // std::cout << "Filling sim hit info" << std::endl;
      simHitInfo_.fillSimHitInfo(*simhitIt);
      
    }
  }
  simHitTree->Fill();

  for (Phase2TrackerRecHit1DCollectionNew::const_iterator DSViter = rechits->begin();
       DSViter != rechits->end(); ++DSViter) {
    // Get the detector unit's id
    unsigned int rawid(DSViter->detId());
    DetId detId(rawid);
    unsigned int layer = (tTopo->side(detId) != 0) * 1000;  // don't split up endcap sides
    layer += tTopo->layer(detId);

    // Get the geomdet
    const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
    if (!geomDetUnit)
      continue;

    recHitInfo_.setTopology(tTopo);
    recHitInfo_.setGeometry(tkGeom);

    // Loop over the rechits in the detector unit
    for (edmNew::DetSet<Phase2TrackerRecHit1D>::const_iterator rechitIt = DSViter->begin();
	 rechitIt != DSViter->end(); ++rechitIt) {
      recHitInfo_.fillRecHitInfo(*rechitIt,rawid,geomDetUnit,&pixelSimLinks,simTracks,simHitsRaw,debugHitMatch_);
    }
  }
  recHitTree->Fill();

}


void RecHitTreeWA::beginJob()
{
  edm::Service<TFileService> fs;
  recHitTree = fs->make<TTree>( "RecHitTree", "RecHitTree" );
  recHitInfo_.setBranches(*recHitTree);

  simHitTree = fs->make<TTree>( "SimHitTree", "SimHitTree" );
  simHitInfo_.setBranches(*simHitTree);
  
  simTrackTree = fs->make<TTree>( "SimTrackTree", "SimTrackTree" );
  simTrackInfo_.setBranches(*simTrackTree);

}

void RecHitTreeWA::endJob()
{}

void RecHitTreeWA::initEventStructure()
{
  recHitInfo_.clear();

  simHitInfo_.clear();

  simTrackInfo_.clear();

}

DEFINE_FWK_MODULE(RecHitTreeWA);

