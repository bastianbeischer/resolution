// $Id: RES_EventActionGeneration.hh,v 1.4 2010/07/19 20:20:09 beischer Exp $

#ifndef RES_EventActionGeneration_hh
#define RES_EventActionGeneration_hh

#include "G4UserEventAction.hh"

class RES_Event;

class RES_EventActionGeneration : public G4UserEventAction
{

public:
  RES_EventActionGeneration();
  ~RES_EventActionGeneration();

  void BeginOfEventAction(const G4Event* event);
  void EndOfEventAction(const G4Event* event);

private:
  void SmearHits(RES_Event* event);
  void AddNoiseHits(RES_Event* event);

};

#endif /* RES_EventActionGeneration_hh */
