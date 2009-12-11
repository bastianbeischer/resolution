// $Id: RES_EventActionGeneration.cc,v 1.17 2009/12/11 12:52:24 beischer Exp $

#include "RES_EventActionGeneration.hh"

#include "RES_DetectorConstruction.hh"
#include "RES_RunManager.hh"
#include "RES_DataHandler.hh"
#include "RES_AlignmentManager.hh"
#include "RES_SD.hh"

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
    RES_HitsCollection* HC = (RES_HitsCollection*) HCofTE->GetHC(HCID);
  
    RES_Event newEvent;

    G4int NbHits = HC->entries();
    for (int i = 0; i < NbHits; i++) {
      RES_Hit* hit = (*HC)[i];
      newEvent.AddHit(hit->GetModuleID(),hit->GetLayerID(),hit->GetPosition().x(),hit->GetPosition().y(),hit->GetPosition().z());
    }

    newEvent.SetEventType(generated);

    SmearHits(&newEvent);

    G4PrimaryParticle* primary = event->GetPrimaryVertex()->GetPrimary();
    G4ThreeVector momentum = primary->GetMomentum();
    newEvent.SetMomentum(momentum.mag());
    newEvent.SetTransverseMomentum(sqrt(momentum.y()*momentum.y() + momentum.z()*momentum.z()));

    dataHandler->AddEvent(newEvent);
  }
}

void RES_EventActionGeneration::SmearHits(RES_Event* event)
{
  G4int nHits = event->GetNbOfHits();
  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();

  for (int i = 0; i < nHits; i++) {
    G4int iModule = event->GetModuleID(i);
    G4int iLayer  = event->GetLayerID(i);
    G4double angle = det->GetModuleAngle(iModule);
    if (iLayer > 0) angle += det->GetModuleInternalAngle(iModule);

    RES_AlignmentManager* alignMgr = RES_AlignmentManager::GetInstance();

    angle += alignMgr->GetAngleShift(iModule);
    
    // collect hit information
    G4double x = event->GetHitPosition(i).x();
    G4double y = event->GetHitPosition(i).y();
    G4double z = event->GetHitPosition(i).z();
    G4ThreeVector hit(x,y,z);

    // transform from detector reference frame (x,y,z) to module frame (u,v,z), smear the hits, transform back to detector frame
    G4RotationMatrix forwardRotation(angle, 0., 0.);
    G4RotationMatrix backwardRotation(-angle, 0., 0.);

    // if (m_fitMethod == transverse)
    //   hit.setX(0.);

    G4double sigmaU, sigmaV, sigmaZ;
    if (iLayer == 0) {
      sigmaU = det->GetModuleUpperSigmaU(iModule);
      sigmaV = det->GetModuleUpperSigmaV(iModule);
      sigmaZ = det->GetModuleUpperSigmaZ(iModule);
    }
    else {
      sigmaU = det->GetModuleLowerSigmaU(iModule);
      sigmaV = det->GetModuleLowerSigmaV(iModule);
      sigmaZ = det->GetModuleLowerSigmaZ(iModule);
    }



    hit = forwardRotation*hit;
    //    hit.setX(CLHEP::RandFlat::shoot());
    //    hit.setX(CLHEP::RandGauss::shoot(0., sigmaU));
    hit.setX(0.);
    hit.setY(CLHEP::RandGauss::shoot(hit.y(), sigmaV));
    hit.setZ(CLHEP::RandGauss::shoot(hit.z(), sigmaZ));
    hit = backwardRotation*hit;

    hit.setX(hit.x() + alignMgr->GetXshift(iModule));
    hit.setY(hit.y() + alignMgr->GetYshift(iModule));

    event->AddSmearedHit(hit.x(), hit.y(), hit.z());
  }
}
