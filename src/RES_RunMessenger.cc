#include "RES_RunMessenger.hh"

#include "RES_RunManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"

RES_RunMessenger::RES_RunMessenger(RES_RunManager* manager)
{
  m_manager = manager;

  m_directory = new G4UIdirectory("/RES/Run/");
  m_directory->SetGuidance("Commands for the run manager");

  m_setStoreResultsCmd = new G4UIcmdWithABool("/RES/Run/StoreResults", this);
  m_setStoreResultsCmd->SetGuidance("Store results?");
  m_setStoreResultsCmd->SetParameterName("store", true);
  m_setStoreResultsCmd->SetDefaultValue(true);
  m_setStoreResultsCmd->AvailableForStates(G4State_PreInit);

  m_generateCmd = new G4UIcmdWithAnInteger("/RES/Run/Generate", this);
  m_generateCmd->SetGuidance("Start a generation run");
  m_generateCmd->SetParameterName("nEvents", true);
  m_generateCmd->SetDefaultValue(1);
  m_generateCmd->AvailableForStates(G4State_Idle);

  m_reconstructCmd = new G4UIcmdWithoutParameter("/RES/Run/Reconstruct", this);
  m_reconstructCmd->SetGuidance("Reconstruct the events which are stored in the DataHandler");
  m_reconstructCmd->AvailableForStates(G4State_Idle);
}

RES_RunMessenger::~RES_RunMessenger()
{
  delete m_directory;
  delete m_setStoreResultsCmd;
  delete m_generateCmd;
  delete m_reconstructCmd;
}

void RES_RunMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_setStoreResultsCmd) {
    m_manager->SetStoreResults(m_setStoreResultsCmd->GetNewBoolValue(newValue));
  }
  if (command == m_generateCmd) {
    m_manager->StartGenerationRun(m_generateCmd->GetNewIntValue(newValue));
  }
  if (command == m_reconstructCmd) {
    m_manager->StartReconstructionRun();
  }
}
