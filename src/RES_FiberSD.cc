#include "RES_FiberSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"

RES_FiberSD::RES_FiberSD(G4String name) :
  G4VSensitiveDetector(name)
{
  G4String HCname = "fiberHitsCollection";
  collectionName.insert(HCname);
}

RES_FiberSD::~RES_FiberSD()
{
}

void RES_FiberSD::Initialize(G4HCofThisEvent* HCE)
{
  fiberHitsCollection = new RES_FiberHitsCollection(SensitiveDetectorName, collectionName[0]);
  static G4int HCID = -1;
  if (HCID < 0) HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  HCE->AddHitsCollection( HCID, fiberHitsCollection );
}

G4bool RES_FiberSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  RES_FiberHit* newHit = new RES_FiberHit();
  newHit->SetPosition(aStep->GetPreStepPoint()->GetPosition());
  fiberHitsCollection->insert(newHit);
  newHit->Draw();

  return true;
}

void RES_FiberSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel>=0) { 
    G4int NbHits = fiberHitsCollection->entries();
    G4cout << "\n-------->Hits Collection: in this event there are " << NbHits 
	   << " hits in the tracker fibers: " << G4endl;
    for (G4int i=0;i<NbHits;i++) (*fiberHitsCollection)[i]->Print();
  } 
}
