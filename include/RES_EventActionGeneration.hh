// $Id: RES_EventActionGeneration.hh,v 1.2 2009/10/14 09:24:31 beischer Exp $

#ifndef RES_EventActionGeneration_hh
#define RES_EventActionGeneration_hh

#include "G4UserEventAction.hh"

class RES_EventActionGeneration : public G4UserEventAction
{

public:
  RES_EventActionGeneration();
  ~RES_EventActionGeneration();

  void BeginOfEventAction(const G4Event* event);
  void EndOfEventAction(const G4Event* event);

};

#endif /* RES_EventActionGeneration_hh */
