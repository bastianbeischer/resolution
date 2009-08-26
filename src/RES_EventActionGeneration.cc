#include "RES_EventActionGeneration.hh"

#include "RES_RunManager.hh"
#include "RES_DataHandler.hh"
#include "RES_FiberSD.hh"

#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

RES_EventActionGeneration::RES_EventActionGeneration()
{
}

RES_EventActionGeneration::~RES_EventActionGeneration()
{
}

void RES_EventActionGeneration::BeginOfEventAction(const G4Event* event)
{
  if( (event->GetEventID() > 0) && (event->GetEventID() % 100 == 0) )
    G4cout << ">>> Event " << event->GetEventID() << G4endl;
}

void RES_EventActionGeneration::EndOfEventAction(const G4Event* event)
{
  RES_RunManager* runManager = (RES_RunManager*) G4RunManager::GetRunManager();

  if (runManager->GetStoreResults()) {
    RES_DataHandler* dataHandler = runManager->GetDataHandler();

    G4HCofThisEvent* HCofTE = event->GetHCofThisEvent();
    G4String collectionName = "fiberHitsCollection";
    G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName);
    RES_FiberHitsCollection* fiberHC = (RES_FiberHitsCollection*) HCofTE->GetHC(HCID);
  
    RES_Event newEvent;

    G4int NbHits = fiberHC->entries();
    for (int i = 0; i < NbHits; i++) {
      RES_FiberHit* hit = (*fiberHC)[i];
      newEvent.AddHit(hit->GetModuleID(),hit->GetFiberID(),hit->GetPosition().x(),hit->GetPosition().y(),hit->GetPosition().z());
    }

    newEvent.SetEventType(generated);

    G4PrimaryParticle* primary = event->GetPrimaryVertex()->GetPrimary();
    G4double momentum = primary->GetMomentum().mag();
    newEvent.SetMomentum(momentum);

    dataHandler->AddEvent(newEvent);
  }
}
