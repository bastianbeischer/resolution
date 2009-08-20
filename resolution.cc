#include "G4RunManager.hh"

#include "RES_DetectorConstruction.hh"
#include "RES_PhysicsList.hh"
#include "RES_PrimaryGeneratorAction.hh"
#include "RES_FieldManager.hh"

#include "G4TransportationManager.hh"
#include "G4UIsession.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

int main(int argc, char** argv)
{
  // create a new RunManager
  G4RunManager* runManager = new G4RunManager();

  // set user initializations
  RES_DetectorConstruction* detectorConstruction = new RES_DetectorConstruction();
  runManager->SetUserInitialization(detectorConstruction);
  RES_PhysicsList* physicsList = new RES_PhysicsList();
  runManager->SetUserInitialization(physicsList);
  RES_PrimaryGeneratorAction* generatorAction = new RES_PrimaryGeneratorAction();
  runManager->SetUserAction(generatorAction);

  RES_FieldManager* fieldManager = new RES_FieldManager();
  G4TransportationManager::GetTransportationManager()->SetFieldManager(fieldManager);

#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = 0;
  visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // User interface
  if(argc==1){ //Interactive mode	
    G4UIsession* session = new G4UIterminal(new G4UItcsh);
    session->SessionStart(); 
    delete session;
  }
  else { // Batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    
    G4UImanager* UI = G4UImanager::GetUIpointer();
    UI->ApplyCommand(command+fileName);
  }

  delete runManager;

  return 0;
}
