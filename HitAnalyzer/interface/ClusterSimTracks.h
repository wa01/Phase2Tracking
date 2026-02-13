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

 ClusterSimTracks(const Phase2TrackerCluster1D& cluster, const DetId& detId, const edm::DetSetVector<PixelDigiSimLink>& pixelSimLinks):
  cluster_(cluster), detId_(detId), pixelSimLinks_(pixelSimLinks), channelsChecked_(0) {};

  std::set<unsigned int> getSimTrackId(unsigned int channel) const;

  bool listIsComplete() {
    //
    // Return true if simTrackIds have been collected from all cluster channels
    //
    return channelsChecked_==cluster_.size();
  };

  bool simTrackInCluster(unsigned int simTrackId);
  
  const std::set<unsigned int>& simTrackIds() {
    //
    // List of all simTrackIds associated to cluster.
    //   Only available if listIsComplete() is true !!!
    //
    if ( not listIsComplete() )
      throw std::logic_error("ClusterSimTracks: trying to retrieve incomplete list of SimTrackIds");

    return simTrackIds_;
  };
  
 private:
  const Phase2TrackerCluster1D& cluster_;
  const DetId& detId_;
  const edm::DetSetVector<PixelDigiSimLink>& pixelSimLinks_;
  unsigned short channelsChecked_;
  std::set<unsigned int> simTrackIds_;
};
#endif
