// $Id: RES_RunManager.cc,v 1.20.2.1 2010/09/08 17:45:00 beischer Exp $

#include "RES_RunManager.hh"

#include "G4RunManager.hh"

#include "RES_DetectorConstruction.hh"
#include "RES_RunMessenger.hh"
#include "RES_DataHandler.hh"
#include "RES_EventActionGeneration.hh"
#include "RES_EventActionReconstruction.hh"
#include "RES_TrackFitter.hh"
#include "RES_Event.hh"
#include "RES_PrimaryGeneratorAction.hh"

RES_RunManager* RES_RunManager::m_instance = 0;

RES_RunManager* RES_RunManager::GetRunManager()
{
  if (!m_instance) m_instance = new RES_RunManager();
  return m_instance;
}

RES_RunManager::RES_RunManager() :
  G4RunManager()
{
  m_messenger = new RES_RunMessenger(this);
  m_dataHandler = new RES_DataHandler();
  m_trackFitter = RES_TrackFitter::GetInstance();
  m_eventActionGen = new RES_EventActionGeneration();
  m_eventActionRec = new RES_EventActionReconstruction(m_trackFitter);
  m_storeResults = false;
  m_fixedDof = -1;

  SetActionsForGeneration();
}

RES_RunManager::~RES_RunManager()
{
  delete m_dataHandler;
}

void RES_RunManager::SetStoreResults(G4bool value)
{
  m_storeResults = value;
  m_dataHandler->Initialize();
}

void RES_RunManager::StartGenerationRun(G4int nEvents)
{
  m_dataHandler->Initialize();
  SetActionsForGeneration();
  BeamOn(nEvents);
  if (m_storeResults) m_dataHandler->WriteFile();
  G4cout << G4endl
         << " ---------------------------------------------- " << G4endl
         << " Successfully generated " << nEvents << " events" << G4endl
         << " ---------------------------------------------- " << G4endl
         << G4endl;
}

void RES_RunManager::StartReconstructionRun()
{
  m_dataHandler->Initialize();
  SetActionsForReconstruction();

  RES_PrimaryGeneratorAction* primaryGeneratorAction = (RES_PrimaryGeneratorAction*) GetUserPrimaryGeneratorAction();

  G4bool randOriginValue = primaryGeneratorAction->GetRandomOrigin();
  G4bool randDirectionValue = primaryGeneratorAction->GetRandomDirection();
  primaryGeneratorAction->SetRandomOrigin(false);
  primaryGeneratorAction->SetRandomDirection(false);

  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  std::vector<double> backupEfficiencies;
  for (unsigned int i = 0; i < det->GetNumberOfModules(); i++) {
    RES_Module* module = det->GetModule(i);
    backupEfficiencies.push_back(module->GetUpperEfficiency());
    backupEfficiencies.push_back(module->GetLowerEfficiency());
    module->SetUpperEfficiency(1.0);
    module->SetLowerEfficiency(1.0);
  }

  int Nevents = m_dataHandler->GetNumberOfGeneratedEvents();
  for (int i = 0; i < Nevents; i++) {

    if( (i > 0) && (i % 100 == 0) )
      G4cout << ">>> Event " << i << G4endl;

    m_dataHandler->LoadGeneratedEntry(i);
    RES_Event genEvent = m_dataHandler->GetCurrentEvent();
    m_trackFitter->SetCurrentGenEvent(genEvent);
    RES_Event recEvent;
    recEvent = m_trackFitter->Fit();

    if (m_storeResults) {
      m_dataHandler->AddEvent(recEvent);
    }
  }
  if (m_storeResults) {
    m_dataHandler->WriteFile();
  }

  for (unsigned int i = 0; i < det->GetNumberOfModules(); i++) {
    RES_Module* module = det->GetModule(i);
    module->SetUpperEfficiency(backupEfficiencies.at(2*i));
    module->SetLowerEfficiency(backupEfficiencies.at(2*i+1));
  }

  primaryGeneratorAction->SetRandomOrigin(randOriginValue);  
  primaryGeneratorAction->SetRandomDirection(randDirectionValue);  

  G4int nEvents = m_dataHandler->GetNumberOfReconstructedEvents();
  G4cout << G4endl
         << " -------------------------------------------------- " << G4endl
         << " Successfully reconstructed " << nEvents << " events" << G4endl
         << " -------------------------------------------------- " << G4endl
         << G4endl;
}

void RES_RunManager::StartReconstructionRunWithoutLayer(G4int layer)
{
  if (layer == 1) {
    m_trackFitter->AddLayerToBeSkipped(0);
    m_trackFitter->AddLayerToBeSkipped(2);
    m_trackFitter->AddLayerToBeSkipped(4);
  }
  else if (layer == 2) {
    m_trackFitter->AddLayerToBeSkipped(1);
    m_trackFitter->AddLayerToBeSkipped(3);
    m_trackFitter->AddLayerToBeSkipped(5);
  }
  else if (layer == 3) {
    m_trackFitter->AddLayerToBeSkipped(6);
    m_trackFitter->AddLayerToBeSkipped(8);
  }
  else if (layer == 4) {
    m_trackFitter->AddLayerToBeSkipped(7);
    m_trackFitter->AddLayerToBeSkipped(9);
  }
  else if (layer == 5) {
    m_trackFitter->AddLayerToBeSkipped(10);
    m_trackFitter->AddLayerToBeSkipped(12);
  }
  else if (layer == 6) {
    m_trackFitter->AddLayerToBeSkipped(11);
    m_trackFitter->AddLayerToBeSkipped(13);
  }
  else if (layer == 7) {
    m_trackFitter->AddLayerToBeSkipped(14);
    m_trackFitter->AddLayerToBeSkipped(16);
    m_trackFitter->AddLayerToBeSkipped(18);
  }
  else if (layer == 8) {
    m_trackFitter->AddLayerToBeSkipped(15);
    m_trackFitter->AddLayerToBeSkipped(17);
    m_trackFitter->AddLayerToBeSkipped(19);
  }
  //  m_trackFitter->AddLayerToBeSkipped(layer);
  StartReconstructionRun();
  m_trackFitter->ClearLayersToBeSkipped();
  //  m_trackFitter->RemoveLayerToBeSkipped(layer);
}

void RES_RunManager::ScanChi2Function(G4int iPar, G4int jPar, G4String filename)
{
  m_dataHandler->Initialize();
  SetActionsForReconstruction();

  int Nevents = m_dataHandler->GetNumberOfGeneratedEvents();
  for (int i = 0; i < Nevents; i++) {

    if( (i > 0) && (i % 100 == 0) )
      G4cout << ">>> Event " << i << G4endl;

    m_dataHandler->LoadGeneratedEntry(i);
    RES_Event genEvent = m_dataHandler->GetCurrentEvent();
    m_trackFitter->SetCurrentGenEvent(genEvent);
    m_trackFitter->ScanChi2Function(iPar,jPar,filename);
  }
}

void RES_RunManager::SetActionsForGeneration()
{
  SetUserAction(m_eventActionGen);
}

void RES_RunManager::SetActionsForReconstruction()
{
  SetUserAction(m_eventActionRec);
}
