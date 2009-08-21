#ifndef RES_TrackFitter_hh
#define RES_TrackFitter_hh

#include "RES_Event.hh"

class RES_TrackFitter
{

public:
  RES_TrackFitter();
  ~RES_TrackFitter();

public:
  void SetCurrentRecEvent(RES_Event event) {m_currentRecEvent = event;}
  RES_Event Fit(RES_Event genEvent);

  RES_Event m_currentRecEvent;

};

#endif /* RES_TrackFitter_hh */
