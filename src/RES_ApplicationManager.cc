#include "RES_ApplicationManager.hh"

#include "RES_ApplicationMessenger.hh"

#include "RES_RunManager.hh"
#include "RES_DetectorConstruction.hh"
#include "RES_PhysicsList.hh"
#include "RES_PrimaryGeneratorAction.hh"
#include "RES_FieldManager.hh"
#include "RES_DataHandler.hh"
#include "RES_AlignmentManager.hh"

#include "G4TransportationManager.hh"
#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "G4UIsession.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

RES_ApplicationManager::RES_ApplicationManager()
{
  m_messenger = new RES_ApplicationMessenger(this);

  // create a new RunManager
  m_runManager = RES_RunManager::GetRunManager();

  RES_AlignmentManager::GetInstance();

  // set user initializations
  RES_DetectorConstruction* detectorConstruction = new RES_DetectorConstruction();
  m_runManager->SetUserInitialization(detectorConstruction);
  RES_PhysicsList* physicsList = new RES_PhysicsList();
  m_runManager->SetUserInitialization(physicsList);
  RES_PrimaryGeneratorAction* generatorAction = new RES_PrimaryGeneratorAction();
  m_runManager->SetUserAction(generatorAction);

  RES_FieldManager* fieldManager = new RES_FieldManager();
  G4TransportationManager::GetTransportationManager()->SetFieldManager(fieldManager);

#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
}

RES_ApplicationManager::~RES_ApplicationManager()
{
  delete m_messenger;
  delete m_runManager;
}

// run batch job from script name
int RES_ApplicationManager::RunBatchScript(G4String scriptName)
{
  G4String command = "/control/execute ";
  G4String fileName = scriptName;
    
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand(command+fileName);

  return 1;
}

void RES_ApplicationManager::CreateSession()
{
  G4UIsession * session = new G4UIterminal(new G4UItcsh);
  session->SessionStart(); 
  delete session;
}
