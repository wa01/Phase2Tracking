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
  std::vector<float> Hit_x;
  std::vector<float> Hit_y;
  std::vector<float> Hit_z;
};

class RecHitTree : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit RecHitTree(const edm::ParameterSet&);
    void analyze(const edm::Event&, const edm::EventSetup&);

  private:
    virtual void beginJob() override;
    virtual void endJob() override;
    void initEventStructure();

    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> esTokenGeom_;
    const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> esTokenTopo_;
    const edm::EDGetTokenT<Phase2TrackerRecHit1DCollectionNew> tokenRecHits_;
    const edm::EDGetTokenT<Phase2TrackerCluster1DCollectionNew> tokenClusters_;

    TTree* hitTree;
    HitInfo* hitInfo;
};

RecHitTree::RecHitTree(const edm::ParameterSet& cfg)
  : esTokenGeom_(esConsumes()),
    esTokenTopo_(esConsumes()),
    tokenRecHits_(consumes<Phase2TrackerRecHit1DCollectionNew>(cfg.getParameter<edm::InputTag>("rechits"))),
    tokenClusters_(consumes<Phase2TrackerCluster1DCollectionNew>(cfg.getParameter<edm::InputTag>("clusters")))
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

  for (Phase2TrackerRecHit1DCollectionNew::const_iterator DSViter = rechits->begin(); DSViter != rechits->end(); ++DSViter) {
    // Get the detector unit's id
    unsigned int rawid(DSViter->detId());
    DetId detId(rawid);

    // Get the geomdet
    const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
    if (!geomDetUnit)
      continue;

    // Loop over the rechits in the detector unit
    for (edmNew::DetSet<Phase2TrackerRecHit1D>::const_iterator rechitIt = DSViter->begin(); rechitIt != DSViter->end(); ++rechitIt) {
      // determine the position
      LocalPoint localPosClu = rechitIt->localPosition();
      Global3DPoint globalPosClu = geomDetUnit->surface().toGlobal(localPosClu);

      hitInfo->Hit_x.push_back(globalPosClu.x());
      hitInfo->Hit_y.push_back(globalPosClu.y());
      hitInfo->Hit_z.push_back(globalPosClu.z());
    }
  }
  hitTree->Fill();
}

void RecHitTree::beginJob()
{
  edm::Service<TFileService> fs;
  hitTree = fs->make<TTree>( "HitTree", "HitTree" );

  hitTree->Branch("Hit_x",    &hitInfo->Hit_x);
  hitTree->Branch("Hit_y",    &hitInfo->Hit_y);
  hitTree->Branch("Hit_z",    &hitInfo->Hit_z);
}

void RecHitTree::endJob()
{}

void RecHitTree::initEventStructure()
{
  hitInfo->Hit_x.clear();
  hitInfo->Hit_y.clear();
  hitInfo->Hit_z.clear();
}

DEFINE_FWK_MODULE(RecHitTree);
