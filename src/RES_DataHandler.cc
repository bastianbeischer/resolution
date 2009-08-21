#include "RES_DataHandler.hh"

#include "RES_FiberHit.hh"
#include "RES_DataMessenger.hh"

#include "TFile.h"
#include "TTree.h"

RES_DataHandler::RES_DataHandler() :
  m_fileName("output.root"),
  m_file(0),
  m_genTree(0),
  m_recTree(0),
  m_event(0),
  m_initialized(false)
{
  m_messenger = new RES_DataMessenger(this);
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
  if (!m_initialized) {
    m_file = new TFile(m_fileName.c_str(), "UPDATE");
    m_genTree = (TTree*) m_file->Get("resolution_gen_tree");
    m_recTree = (TTree*) m_file->Get("resolution_rec_tree");

    m_event = new RES_Event();

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
  m_initialized = true;
}

G4int RES_DataHandler::GetNumberOfGeneratedEvents()
{
  return m_genTree->GetEntries();
}

void RES_DataHandler::LoadGeneratedEntry(G4int i)
{
  m_genTree->GetEntry(i);
}

RES_Event RES_DataHandler::GetCurrentEvent()
{
  RES_Event currentEvent = *m_event;
  return currentEvent;
}

void RES_DataHandler::AddEvent(RES_Event event)
{
  if (m_initialized) {
    m_event = new RES_Event();
    *m_event = event;
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
}

void RES_DataHandler::WriteFile()
{
  if (m_initialized) {
    m_genTree->Write();
    m_recTree->Write();
    m_file->Write();
    m_file->Close();
  }
  else {
    G4cerr << "Data Handler not initialized! (Call SetFileName or construct with string parameter)" << G4endl;
  }
}
