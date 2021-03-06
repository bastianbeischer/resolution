// $Id: RES_EventActionReconstruction.hh,v 1.3 2009/10/14 09:24:31 beischer Exp $

#ifndef RES_EventActionReconstruction_hh
#define RES_EventActionReconstruction_hh

#include "G4UserEventAction.hh"

class RES_TrackFitter;

class RES_EventActionReconstruction : public G4UserEventAction
{

public:
  RES_EventActionReconstruction(RES_TrackFitter* trackFitter);
  ~RES_EventActionReconstruction();

  void BeginOfEventAction(const G4Event* event);
  void EndOfEventAction(const G4Event* event);

private:
  RES_TrackFitter* m_trackFitter;

};

#endif /* RES_EventActionReconstruction_hh */
