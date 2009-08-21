#include "RES_RunManager.hh"

#include "RES_RunMessenger.hh"
#include "RES_DataHandler.hh"
#include "RES_EventActionGeneration.hh"
#include "RES_EventActionReconstruction.hh"

#include "RES_FiberHit.hh"

RES_RunManager::RES_RunManager() :
  G4RunManager()
{
  m_messenger = new RES_RunMessenger(this);
  m_dataHandler = new RES_DataHandler();
  m_eventActionGen = new RES_EventActionGeneration();
  m_eventActionRec = new RES_EventActionReconstruction();

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
  /*******************************************
   * - LOOP OVER ALL GENERATED EVENTS
   * --- READ generated event
   * --- IVOKE minimizer
   * --- STORE reconstructed event
   */
}

void RES_RunManager::SetActionsForGeneration()
{
  SetUserAction(m_eventActionGen);
}

void RES_RunManager::SetActionsForReconstruction()
{
  SetUserAction(m_eventActionRec);
}
