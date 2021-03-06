// $Id: RES_EventActionReconstruction.cc,v 1.8 2009/12/11 12:52:24 beischer Exp $

#include "RES_EventActionReconstruction.hh"

#include "RES_TrackFitter.hh"
#include "RES_SD.hh"
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
  RES_HitsCollection* HC = (RES_HitsCollection*) HCofTE->GetHC(HCID);
  
  RES_Event recEvent; 

  G4int NbHits = HC->entries();

  for (int i = 0; i < NbHits; i++) {
    RES_Hit* hit = (*HC)[i];
    recEvent.AddHit(hit->GetModuleID(),hit->GetLayerID(),hit->GetPosition().x(),hit->GetPosition().y(),hit->GetPosition().z());
  }

  recEvent.SetEventType(reconstructed);

  G4PrimaryParticle* primary = event->GetPrimaryVertex()->GetPrimary();
  G4ThreeVector momentum = primary->GetMomentum();
  recEvent.SetMomentum(momentum.mag());
  recEvent.SetTransverseMomentum(sqrt(momentum.y()*momentum.y() + momentum.z()*momentum.z()));

  m_trackFitter->SetCurrentRecEvent(recEvent);
}
