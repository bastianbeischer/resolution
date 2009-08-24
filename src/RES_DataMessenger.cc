#include "RES_DataMessenger.hh"

#include "RES_DataHandler.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

RES_DataMessenger::RES_DataMessenger(RES_DataHandler* dataHandler)
{
  m_dataHandler = dataHandler;

  m_directory = new G4UIdirectory("/RES/Data/");
  m_directory->SetGuidance("Commands for the data handler");

  m_overwriteFileCmd = new G4UIcmdWithABool("/RES/Data/OverWriteFile", this);
  m_overwriteFileCmd->SetGuidance("Overwrite files upon data handler initialization? Otherwise reopen and append new data.");
  m_overwriteFileCmd->SetParameterName("overwrite", true);
  m_overwriteFileCmd->SetDefaultValue(true);
  m_overwriteFileCmd->AvailableForStates(G4State_PreInit);

  m_setFileNameCmd = new G4UIcmdWithAString("/RES/Data/SetFileName", this);
  m_setFileNameCmd->SetGuidance("Set the file name to which results will be written or appended.");
  m_setFileNameCmd->SetParameterName("filename", false);
  m_setFileNameCmd->AvailableForStates(G4State_PreInit);
}

RES_DataMessenger::~RES_DataMessenger()
{
  delete m_directory;
  delete m_overwriteFileCmd;
  delete m_setFileNameCmd;
}

void RES_DataMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_setFileNameCmd) {
    m_dataHandler->SetFileName(newValue);
  }
  if (command == m_overwriteFileCmd) {
    m_dataHandler->SetOverWriteFile(m_overwriteFileCmd->GetNewBoolValue(newValue));
  }
}
