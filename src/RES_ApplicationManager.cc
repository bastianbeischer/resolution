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
#include "G4UIbatch.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif
#ifdef G4UI_USE_QT
#include "G4UIQt.hh"
#endif

#include "G4UIsession.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "CLHEP/Random/Random.h"
#include <time.h>

RES_ApplicationManager::RES_ApplicationManager(int argc, char** argv)
{
  m_argc = argc;
  m_argv = argv;

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

void RES_ApplicationManager::CreateSession(SessionType type, G4String macro)
{
  G4UIsession* currentSession = G4UImanager::GetUIpointer()->GetSession();
  if (currentSession) {
    G4UImanager::GetUIpointer()->ApplyCommand("exit");
  }

  G4UIsession* session = 0;

  // setup the session
  if (type == Terminal) {
    session = new G4UIterminal(new G4UItcsh);
  }
#ifdef G4UI_USE_QT
  else if (type == Qt) {
    G4UIQt* qtsession = new G4UIQt(m_argc, m_argv);
    qtsession->AddMenu("macs", "Macros");
    qtsession->AddButton("macs", "vis.ogl.mac", "/control/execute mac/vis.ogl.mac");
    session = qtsession;
  }
#endif
  if (macro != "") {
    G4UImanager::GetUIpointer()->ExecuteMacroFile(macro);
  }

  if (session != 0) {
    session->SessionStart();
    delete session;
  }
}

void RES_ApplicationManager::SetSeedToSystemTime()
{
  long ltime = time(NULL);
  int  stime = (unsigned) ltime/2;
  CLHEP::HepRandom::setTheSeed(stime);
}
