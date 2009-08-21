#ifndef RES_TrackFitter_hh
#define RES_TrackFitter_hh

class RES_Event;

class RES_TrackFitter
{

public:
  RES_TrackFitter();
  ~RES_TrackFitter();

public:
  RES_Event Fit(RES_Event genEvent);

};

#endif /* RES_TrackFitter_hh */
