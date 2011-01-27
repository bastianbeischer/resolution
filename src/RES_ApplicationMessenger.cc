#include "RES_ApplicationMessenger.hh"

#include "RES_ApplicationManager.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"

RES_ApplicationMessenger::RES_ApplicationMessenger(RES_ApplicationManager* manager)
{
  m_manager = manager;

  m_createTerminalSessionCmd = new G4UIcmdWithAString("/RES/CreateTerminalSession", this);
  m_createTerminalSessionCmd->SetGuidance("Create an interactive terminal session)");
  m_createTerminalSessionCmd->SetParameterName("macro", true);
#ifdef G4UI_USE_QT  
  m_createQtSessionCmd = new G4UIcmdWithAString("/RES/CreateQtSession", this);
  m_createQtSessionCmd->SetGuidance("Create an interactive Qt session)");
  m_createQtSessionCmd->SetParameterName("macro", true);
#endif
  m_setSeedToSystemTimeCmd = new G4UIcmdWithoutParameter("/RES/SetSeedToSystemTime", this);
  m_setSeedToSystemTimeCmd->SetGuidance("Sets the CLHEP Random Seed to the current system time");
}

RES_ApplicationMessenger::~RES_ApplicationMessenger()
{
  delete m_createTerminalSessionCmd;
#ifdef G4UI_USE_QT
  delete m_createQtSessionCmd;
#endif
  delete m_setSeedToSystemTimeCmd;
}

void RES_ApplicationMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_createTerminalSessionCmd) {
    RES_ApplicationManager::SessionType type = RES_ApplicationManager::Terminal;
    m_manager->CreateSession(type, newValue);
  }
#ifdef G4UI_USE_QT
  if (command == m_createQtSessionCmd) {
    RES_ApplicationManager::SessionType type = RES_ApplicationManager::Qt;
    m_manager->CreateSession(type, newValue);
  }
#endif
  if (command == m_setSeedToSystemTimeCmd) {
    m_manager->SetSeedToSystemTime();
  }
}
