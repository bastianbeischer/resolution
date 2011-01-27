#ifndef RES_ApplicationMessenger_hh
#define RES_ApplicationMessenger_hh

#include "G4UImessenger.hh"

class RES_ApplicationManager;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class RES_ApplicationMessenger : G4UImessenger
{

public:
  RES_ApplicationMessenger(RES_ApplicationManager* manager);
  ~RES_ApplicationMessenger();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:
  RES_ApplicationManager* m_manager;

  G4UIcmdWithAString* m_createTerminalSessionCmd;
#ifdef G4UI_USE_QT
  G4UIcmdWithAString* m_createQtSessionCmd;
#endif
  G4UIcmdWithoutParameter* m_setSeedToSystemTimeCmd;

};

#endif /* RES_ApplicationMessenger_hh */
