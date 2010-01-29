// $Id: RES_RunManager.hh,v 1.9 2010/01/29 12:51:39 beischer Exp $

#ifndef RES_RunManager_hh
#define RES_RunManager_hh

#include "G4RunManager.hh"

class RES_RunMessenger;
class RES_DataHandler;
class RES_EventActionGeneration;
class RES_EventActionReconstruction;
class RES_TrackFitter;

class RES_RunManager : public G4RunManager
{

public:
  static RES_RunManager* GetRunManager();

  ~RES_RunManager();

public:
  void SetStoreResults(G4bool value);
  void SetFixedDof(G4int dof) {m_fixedDof = dof;}

  RES_DataHandler* GetDataHandler()  {return m_dataHandler;}
  G4bool           GetStoreResults() {return m_storeResults;}
  G4int            GetFixedDof()     {return m_fixedDof;}

  void StartGenerationRun(G4int nEvents);
  void StartReconstructionRun();
  void StartReconstructionRunWithoutLayer(G4int layer);
  void ScanChi2Function(G4int iPar, G4int jPar, G4String filename);

private:
  
  RES_RunManager();

  void SetActionsForGeneration();
  void SetActionsForReconstruction();

private:
  static RES_RunManager*          m_instance;

  RES_RunMessenger*               m_messenger;
  RES_DataHandler*                m_dataHandler;

  RES_EventActionGeneration*      m_eventActionGen;
  RES_EventActionReconstruction*  m_eventActionRec;

  RES_TrackFitter*                m_trackFitter;

  G4bool                          m_storeResults;

  G4int                           m_fixedDof;

};

#endif /* RES_RunManager_hh */
