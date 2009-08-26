#include "RES_EventActionReconstruction.hh"

#include "RES_TrackFitter.hh"
#include "RES_FiberSD.hh"
#include "RES_Event.hh"

#include "G4Event.hh"
#include "G4SDManager.hh"


RES_EventActionReconstruction::RES_EventActionReconstruction(RES_TrackFitter* trackFitter)
{
  m_trackFitter = trackFitter;
}

RES_EventActionReconstruction::~RES_EventActionReconstruction()
{
}

void RES_EventActionReconstruction::BeginOfEventAction(const G4Event* /*event*/)
{
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
    recEvent.AddHit(hit->GetModuleID(),hit->GetFiberID(),hit->GetPosition().x(),hit->GetPosition().y(),hit->GetPosition().z());
  }

  recEvent.SetEventType(reconstructed);

  G4PrimaryParticle* primary = event->GetPrimaryVertex()->GetPrimary();
  G4double momentum = primary->GetMomentum().mag();
  recEvent.SetMomentum(momentum);

  m_trackFitter->SetCurrentRecEvent(recEvent);
}
