#include "RES_ApplicationMessenger.hh"

#include "RES_ApplicationManager.hh"
#include "G4UIcmdWithoutParameter.hh"

RES_ApplicationMessenger::RES_ApplicationMessenger(RES_ApplicationManager* manager)
{
  m_manager = manager;

  m_createSessionCmd = new G4UIcmdWithoutParameter("/RES/CreateSession", this);
  m_setSeedToSystemTimeCmd = new G4UIcmdWithoutParameter("/RES/SetSeedToSystemTime", this);
}

RES_ApplicationMessenger::~RES_ApplicationMessenger()
{
  delete m_createSessionCmd;
  delete m_setSeedToSystemTimeCmd;
}

void RES_ApplicationMessenger::SetNewValue(G4UIcommand* command, G4String /*newValue*/)
{
  if (command == m_createSessionCmd) {
    m_manager->CreateSession();
  }
  if (command == m_setSeedToSystemTimeCmd) {
    m_manager->SetSeedToSystemTime();
  }
}
