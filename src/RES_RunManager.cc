#include "RES_RunManager.hh"

#include "RES_DataHandler.hh"

#include "RES_FiberHit.hh"

RES_RunManager::RES_RunManager() :
  G4RunManager()
{
  m_dataHandler = new RES_DataHandler();
}

RES_RunManager::~RES_RunManager()
{
  delete m_dataHandler;
}

void RES_RunManager::TestDataHandler()
{
  // TESTING DATA HANDLER
  m_dataHandler->InitNewEvent();
  m_dataHandler->SetEventType(reconstructed);
  RES_FiberHit* hit = new RES_FiberHit();
  hit->SetPosition(G4ThreeVector(2.4,0,0));
  m_dataHandler->AddHitInformation(hit);
  delete hit;
  m_dataHandler->FinalizeEvent();
  m_dataHandler->WriteFile();
}
