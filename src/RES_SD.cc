// $Id: RES_SD.cc,v 1.1 2009/12/11 12:52:25 beischer Exp $

#include "RES_SD.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"

RES_SD::RES_SD(G4String name) :
  G4VSensitiveDetector(name)
{
  G4String HCname = "fiberHitsCollection";
  collectionName.insert(HCname);
}

RES_SD::~RES_SD()
{
}

void RES_SD::Initialize(G4HCofThisEvent* HCE)
{
  fiberHitsCollection = new RES_HitsCollection(SensitiveDetectorName, collectionName[0]);
  static G4int HCID = -1;
  if (HCID < 0) HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  HCE->AddHitsCollection( HCID, fiberHitsCollection );
}

G4bool RES_SD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  RES_Hit* newHit = new RES_Hit();

  const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int ownCopyNb = touchable->GetReplicaNumber(0);
  G4int motherCopyNb = touchable->GetReplicaNumber(1);
  newHit->SetModuleID(motherCopyNb);
  newHit->SetLayerID(ownCopyNb);
  newHit->SetPosition(aStep->GetPreStepPoint()->GetPosition());
  fiberHitsCollection->insert(newHit);
  newHit->Draw();

  return true;
}

void RES_SD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel>0) { 
    G4int NbHits = fiberHitsCollection->entries();
    G4cout << "\n-------->Hits Collection: in this event there are " << NbHits 
	   << " hits in the tracker fibers: " << G4endl;
    for (G4int i=0;i<NbHits;i++) (*fiberHitsCollection)[i]->Print();
  } 
}
