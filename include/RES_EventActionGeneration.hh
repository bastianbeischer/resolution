// $Id: RES_EventActionGeneration.hh,v 1.3 2009/10/14 16:51:35 beischer Exp $

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

};

#endif /* RES_EventActionGeneration_hh */
