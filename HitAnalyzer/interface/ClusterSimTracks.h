#ifndef HITANALYZER_CLUSTERSIMTRACKS_H
#define HITANALYZER_CLUSTERSIMTRACKS_H

#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
/* #include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h" */
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"

class ClusterSimTracks {
  
 public:

 ClusterSimTracks(const Phase2TrackerCluster1D& cluster, const DetId& detId,
		  const edm::DetSetVector<PixelDigiSimLink>& pixelSimLinks):
  cluster_(cluster), detId_(detId), pixelSimLinks_(pixelSimLinks), simTracksValid_(false) {};

  std::set<unsigned int> simTrackIdsPerChannel(unsigned int channel) const;

  bool simTrackInCluster(unsigned int simTrackId);
  
  const std::set<unsigned int>& simTrackIds();
  
 private:
  const Phase2TrackerCluster1D& cluster_;
  const DetId& detId_;
  const edm::DetSetVector<PixelDigiSimLink>& pixelSimLinks_;
  bool simTracksValid_;
  std::set<unsigned int> simTrackIds_;
};
#endif
