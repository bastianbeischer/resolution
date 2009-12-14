#ifndef RES_ApplicationMessenger_hh
#define RES_ApplicationMessenger_hh

#include "G4UImessenger.hh"

class RES_ApplicationManager;
class G4UIcmdWithoutParameter;

class RES_ApplicationMessenger : G4UImessenger
{

public:
  RES_ApplicationMessenger(RES_ApplicationManager* manager);
  ~RES_ApplicationMessenger();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:
  RES_ApplicationManager* m_manager;

  G4UIcmdWithoutParameter* m_createSessionCmd;

};

#endif /* RES_ApplicationMessenger_hh */
