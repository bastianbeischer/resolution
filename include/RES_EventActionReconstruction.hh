#ifndef RES_EventActionReconstruction_hh
#define RES_EventActionReconstruction_hh

#include "G4UserEventAction.hh"

class RES_EventActionReconstruction : public G4UserEventAction
{

public:
  RES_EventActionReconstruction();
  ~RES_EventActionReconstruction();

  void BeginOfEventAction(const G4Event* event);
  void EndOfEventAction(const G4Event* event);

};

#endif /* RES_EventActionReconstruction_hh */
