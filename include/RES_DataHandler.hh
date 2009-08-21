#ifndef RES_DataHandler_hh
#define RES_DataHandler_hh

#include "G4String.hh"
#include "G4ThreeVector.hh"

#include "RES_Event.hh"

class RES_FiberHit;
class RES_DataMessenger;
class TTree;
class TFile;

class RES_DataHandler
{

public:
  RES_DataHandler();
  RES_DataHandler(G4String fileName);
  ~RES_DataHandler();

public:
  void SetFileName(G4String fileName) {m_fileName = fileName;}

  void Initialize();

  void InitNewEvent();
  void SetEventType(EventType type);
  void AddHitInformation(RES_FiberHit* hit);
  void FinalizeEvent();
  void WriteFile();

private:
  RES_DataMessenger* m_messenger;

  G4String           m_fileName;
  TFile*             m_file;
  TTree*             m_genTree;
  TTree*             m_recTree;
  RES_Event*         m_event;

  G4bool             m_initialized;

};

#endif /* RES_DataHandler_hh */
