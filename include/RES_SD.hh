// $Id: RES_SD.hh,v 1.1 2009/12/11 12:52:23 beischer Exp $

#ifndef RES_SD_hh
#define RES_SD_hh

#include "G4VSensitiveDetector.hh"
#include "G4THitsCollection.hh"

#include "RES_Hit.hh"

class G4Step;
class G4HCofThisEvent;

typedef G4THitsCollection<RES_Hit> RES_HitsCollection;

class RES_SD : public G4VSensitiveDetector
{

public:
  RES_SD(G4String name);
  ~RES_SD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  RES_HitsCollection* fiberHitsCollection;

};

#endif /* RES_SD_hh */
