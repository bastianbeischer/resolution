#include "RES_EventActionReconstruction.hh"

#include "RES_TrackFitter.hh"

#include "G4Event.hh"
#include "G4SDManager.hh"
#include "RES_FiberSD.hh"

RES_EventActionReconstruction::RES_EventActionReconstruction(RES_TrackFitter* trackFitter)
{
  m_trackFitter = trackFitter;
}

RES_EventActionReconstruction::~RES_EventActionReconstruction()
{
}

void RES_EventActionReconstruction::BeginOfEventAction(const G4Event* event)
{
  if( (event->GetEventID() > 0) && (event->GetEventID() % 100 == 0) )
    G4cout << ">>> Event " << event->GetEventID() << G4endl;
}

void RES_EventActionReconstruction::EndOfEventAction(const G4Event* event)
{
  G4HCofThisEvent* HCofTE = event->GetHCofThisEvent();
  G4String collectionName = "fiberHitsCollection";
  G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName);
  RES_FiberHitsCollection* fiberHC = (RES_FiberHitsCollection*) HCofTE->GetHC(HCID);
  
  RES_Event recEvent;

  G4int NbHits = fiberHC->entries();
  for (int i = 0; i < NbHits; i++) {
    RES_FiberHit* hit = (*fiberHC)[i];
    recEvent.AddHit(hit->GetPosition().x()/cm,hit->GetPosition().y()/cm, hit->GetPosition().z()/cm);
  }

  recEvent.SetEventType(reconstructed);

  G4PrimaryParticle* primary = event->GetPrimaryVertex()->GetPrimary();
  G4double momentum = primary->GetMomentum().mag();
  recEvent.SetMomentum(momentum/GeV);

  m_trackFitter->SetCurrentRecEvent(recEvent);
}
