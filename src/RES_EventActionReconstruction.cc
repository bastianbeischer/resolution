#include "RES_EventActionReconstruction.hh"

#include "G4Event.hh"

RES_EventActionReconstruction::RES_EventActionReconstruction()
{
}

RES_EventActionReconstruction::~RES_EventActionReconstruction()
{
}

void RES_EventActionReconstruction::BeginOfEventAction(const G4Event* event)
{
  if( event->GetEventID() % 100 == 0 )
    G4cout << ">>> Event " << event->GetEventID() << G4endl;
}

void RES_EventActionReconstruction::EndOfEventAction(const G4Event* /*event*/)
{
}
