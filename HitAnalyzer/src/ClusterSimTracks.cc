#include "Phase2Tracking/HitAnalyzer/interface/ClusterSimTracks.h"
// #include "Phase2Tracking/HitAnalyzer/interface/CommonHitInfo.h"
// #include "DataFormats/DetId/interface/DetId.h"
//
// return SimTrack ids for a given DetId / channel
//
std::set<unsigned int> ClusterSimTracks::getSimTrackId(unsigned int channel) const {
  std::set<unsigned int> result;
  edm::DetSetVector<PixelDigiSimLink>::const_iterator DSViter(pixelSimLinks_.find(detId_));
  if (DSViter == pixelSimLinks_.end())  return result;
  for (edm::DetSet<PixelDigiSimLink>::const_iterator it = DSViter->data.begin(); it != DSViter->data.end(); ++it) {
    if (channel == it->channel()) {
      result.insert(it->SimTrackId());
    }
  }
  return result;
};

bool ClusterSimTracks::simTrackInCluster(unsigned int simTrackId) {
  //
  // Return true if <simTrackId> contributed to one or more channels of <cluster>
  //   Keep a (partial) cache of SimTrackIds
  //
  //
  // Start with existing list
  //
  bool matched = std::find(simTrackIds_.begin(),simTrackIds_.end(),simTrackId)!=simTrackIds_.end();
  //
  // This is the final result if the list is complete
  //
  if ( listIsComplete() )  return matched;
  //
  // Otherwise: loop over missing channels
  //
  //
  while ( ++channelsChecked_<cluster_.size() ) {
    // get channel number
    unsigned int channel(Phase2TrackerDigi::pixelToChannel(cluster_.firstRow() + channelsChecked_, cluster_.column()));
    // get SimTrackIds for this channel
    std::set<unsigned int> simTrackIds = getSimTrackId(channel);
    matched = std::find(simTrackIds.begin(),simTrackIds.end(),simTrackId)!=simTrackIds.end();
    if ( matched ) break;
    simTrackIds_.insert(simTrackIds.begin(),simTrackIds.end());
  }
  return matched; 
};

