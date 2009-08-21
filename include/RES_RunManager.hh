#ifndef RES_RunManager_hh
#define RES_RunManager_hh

#include "G4RunManager.hh"

class RES_RunMessenger;
class RES_DataHandler;
class RES_EventActionGeneration;
class RES_EventActionReconstruction;

class RES_RunManager : public G4RunManager
{

public:
  RES_RunManager();
  ~RES_RunManager();

public:
  RES_DataHandler* GetDataHandler() {return m_dataHandler;}

  void SetStoreResults(G4bool value);

  void StartGenerationRun(G4int nEvents);
  void StartReconstructionRun();

  void TestDataHandler();

private:
  void SetActionsForGeneration();
  void SetActionsForReconstruction();

private:
  RES_RunMessenger*               m_messenger;
  RES_DataHandler*                m_dataHandler;

  RES_EventActionGeneration*      m_eventActionGen;
  RES_EventActionReconstruction*  m_eventActionRec;

  G4bool                          m_storeResults;

};

#endif /* RES_RunManager_hh */
