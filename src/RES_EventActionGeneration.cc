// $Id: RES_EventActionGeneration.cc,v 1.28 2010/07/22 15:49:05 beischer Exp $

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

  RES_RunManager* mgr = RES_RunManager::GetRunManager();
  G4int nDof = mgr->GetFixedDof();

  if (nDof >= 0) {
    RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    std::vector<G4int> used;
    
    // NASTY HARDCODED!!! for testbeam!
    used.push_back(0);
    used.push_back(1);
    used.push_back(10);
    used.push_back(11);

    G4int nModules = det->GetNumberOfModules();
    G4int iMin = 0;
    G4int iMax = 2*nModules;
    G4int nPar = 4;
    G4int numberOfLayersToRemove = 2*nModules - nPar - nDof;
    for (G4int i = 0; i < numberOfLayersToRemove; i++) {
      G4int layer = floor(CLHEP::RandFlat::shoot(iMin, iMax));
      while(std::find(used.begin(), used.end(), layer) != used.end()) {
        layer = floor(CLHEP::RandFlat::shoot(iMin, iMax));
      }
      used.push_back(layer);
      RES_Module* module = det->GetModule(layer/2);
      layer%2 ? module->SetUpperEfficiency(0.) : module->SetLowerEfficiency(0.);
    }
  }

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
    for (G4int i = 0; i < NbHits; i++) {
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

  //NASTY HARDCODED AGAIN
  RES_RunManager* mgr = RES_RunManager::GetRunManager();
  G4int nDof = mgr->GetFixedDof();
  if (nDof > 0) {
    RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    for (unsigned int i = 0; i < det->GetNumberOfModules(); i++) {
      RES_Module* module = det->GetModule(i);
      module->SetUpperEfficiency(1.);
      module->SetLowerEfficiency(1.);
    }
  }
}

void RES_EventActionGeneration::SmearHits(RES_Event* event)
{
  G4int nHits = event->GetNbOfHits();
  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();

  for (G4int i = 0; i < nHits; i++) {
    G4int iModule = event->GetModuleID(i);
    G4int iLayer  = event->GetLayerID(i);
    RES_Module* module = det->GetModule(iModule);
    G4double angle = module->GetAngle();
    if (iLayer > 0) angle += module->GetInternalAngle();

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
      sigmaU = module->GetUpperSigmaU();
      sigmaV = module->GetUpperSigmaV();
      sigmaZ = module->GetUpperSigmaZ();
    }
    else {
      sigmaU = module->GetLowerSigmaU();
      sigmaV = module->GetLowerSigmaV();
      sigmaZ = module->GetLowerSigmaZ();
    }

    // if (iModule == 0 || iModule == 5)
    //   sigmaV += 10*um;

    //    sigmaV += 10*um;

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
