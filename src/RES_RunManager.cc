#include "RES_RunManager.hh"

#include "RES_RunMessenger.hh"
#include "RES_DataHandler.hh"
#include "RES_EventActionGeneration.hh"
#include "RES_EventActionReconstruction.hh"
#include "RES_TrackFitter.hh"
#include "RES_Event.hh"
#include "RES_PrimaryGeneratorAction.hh"

RES_RunManager::RES_RunManager() :
  G4RunManager()
{
  m_messenger = new RES_RunMessenger(this);
  m_dataHandler = new RES_DataHandler();
  m_trackFitter = RES_TrackFitter::GetInstance();
  m_eventActionGen = new RES_EventActionGeneration();
  m_eventActionRec = new RES_EventActionReconstruction(m_trackFitter);
  
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
  SetActionsForGeneration();
  BeamOn(nEvents);
  if (m_storeResults) m_dataHandler->WriteFile();
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

  int Nevents = m_dataHandler->GetNumberOfGeneratedEvents();
  for (int i = 0; i < Nevents; i++) {

    if( (i > 0) && (i % 100 == 0) )
      G4cout << ">>> Event " << i << G4endl;

    m_dataHandler->LoadGeneratedEntry(i);
    RES_Event genEvent = m_dataHandler->GetCurrentEvent();
    m_trackFitter->SetCurrentGenEvent(genEvent);
    RES_Event recEvent = m_trackFitter->Fit();
    if (m_storeResults) {
      m_dataHandler->AddEvent(recEvent);
    }
  }
  if (m_storeResults) {
    m_dataHandler->WriteFile();
  }

  primaryGeneratorAction->SetRandomOrigin(randOriginValue);  
  primaryGeneratorAction->SetRandomDirection(randDirectionValue);  
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
