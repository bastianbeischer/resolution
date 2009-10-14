// $Id: RES_FiberSD.hh,v 1.2 2009/10/14 09:24:31 beischer Exp $

#ifndef RES_FiberSD_hh
#define RES_FiberSD_hh

#include "G4VSensitiveDetector.hh"
#include "G4THitsCollection.hh"

#include "RES_FiberHit.hh"

class G4Step;
class G4HCofThisEvent;

typedef G4THitsCollection<RES_FiberHit> RES_FiberHitsCollection;

class RES_FiberSD : public G4VSensitiveDetector
{

public:
  RES_FiberSD(G4String name);
  ~RES_FiberSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  RES_FiberHitsCollection* fiberHitsCollection;

};

#endif /* RES_FiberSD_hh */
