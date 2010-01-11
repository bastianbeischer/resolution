// $Id: RES_SD.cc,v 1.3 2010/01/11 15:25:50 beischer Exp $

#include "RES_SD.hh"

#include "RES_DetectorConstruction.hh"
#include "RES_Module.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"

#include "CLHEP/Random/RandFlat.h"

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
  const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int ownCopyNb = touchable->GetReplicaNumber(0);
  G4int motherCopyNb = touchable->GetReplicaNumber(1);

  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  RES_Module* module = det->GetModule(motherCopyNb);
  G4double eff = ownCopyNb == 0 ? module->GetUpperEfficiency() : module->GetLowerEfficiency();
  G4double rand = CLHEP::RandFlat::shoot(0.,1.);
  
  if (rand < eff) {
    RES_Hit* newHit = new RES_Hit();
    newHit->SetModuleID(motherCopyNb);
    newHit->SetLayerID(ownCopyNb);
    newHit->SetPosition(aStep->GetPreStepPoint()->GetPosition());
    fiberHitsCollection->insert(newHit);
    newHit->Draw();
  }
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
