#include "RES_DataMessenger.hh"

#include "RES_DataHandler.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

RES_DataMessenger::RES_DataMessenger(RES_DataHandler* dataHandler)
{
  m_dataHandler = dataHandler;

  m_directory = new G4UIdirectory("/RES/Data/");

  m_setFileNameCmd = new G4UIcmdWithAString("/RES/Data/SetFileName", this);
  m_setFileNameCmd->SetGuidance("Set the file name to which results will be written or appended.");
  m_setFileNameCmd->SetParameterName("filename", false);
  m_setFileNameCmd->AvailableForStates(G4State_PreInit);
}

RES_DataMessenger::~RES_DataMessenger()
{
  delete m_directory;
}

void RES_DataMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_setFileNameCmd) {
    m_dataHandler->SetFileName(newValue);
  }
}
