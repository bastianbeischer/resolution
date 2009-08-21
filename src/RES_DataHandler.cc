#include "RES_DataHandler.hh"

#include "RES_FiberHit.hh"

#include "TFile.h"
#include "TTree.h"

RES_DataHandler::RES_DataHandler() :
  m_fileName("output.root")
{
  Initialize();
}

RES_DataHandler::RES_DataHandler(G4String fileName)
{
  m_fileName = fileName;
  Initialize();
}

RES_DataHandler::~RES_DataHandler()
{
  delete m_file;
  delete m_event;
}

void RES_DataHandler::Initialize()
{
  m_file = new TFile(m_fileName.c_str(), "UPDATE");
  m_genTree = (TTree*) m_file->Get("resolution_gen_tree");
  m_recTree = (TTree*) m_file->Get("resolution_rec_tree");

  InitNewEvent();

  if (!m_genTree) {
    m_genTree = new TTree("resolution_gen_tree", "tree for momentum resolution calculation");
    m_genTree->Branch("event", "RES_Event", &m_event);
  }
  else {
    m_genTree->SetBranchAddress("event", &m_event);
  }

  if (!m_recTree) {
    m_recTree = new TTree("resolution_rec_tree", "tree for momentum resolution calculation");
    m_recTree->Branch("event", "RES_Event", &m_event);
  }
  else {
    m_recTree->SetBranchAddress("event", &m_event);
  }

}

void RES_DataHandler::InitNewEvent()
{
  m_event = new RES_Event();
}

void RES_DataHandler::SetEventType(EventType type)
{
  m_event->SetEventType(type);
}

void RES_DataHandler::AddHitInformation(RES_FiberHit* hit)
{
  G4ThreeVector pos = hit->GetPosition();
  m_event->AddHit(pos.x(), pos.y(), pos.z());
}

void RES_DataHandler::FinalizeEvent()
{
  if (m_event->GetEventType() == generated) {
    m_event->SetID(m_genTree->GetEntries());
    m_genTree->Fill();
  }
  else if (m_event->GetEventType() == reconstructed) {
    m_event->SetID(m_recTree->GetEntries());
    m_recTree->Fill();
  }
  delete m_event;
  m_event = 0;
}

void RES_DataHandler::WriteFile()
{
  m_genTree->Write();
  m_recTree->Write();
  m_file->Write();
  m_file->Close();
}
