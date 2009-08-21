#include "RES_TrackFitter.hh"

#include "RES_Event.hh"

RES_TrackFitter::RES_TrackFitter()
{
}

RES_TrackFitter::~RES_TrackFitter()
{
}

RES_Event RES_TrackFitter::Fit(RES_Event genEvent)
{
  RES_Event recEvent = genEvent;
  recEvent.SetEventType(reconstructed);
  return recEvent;
}
