#ifndef RES_StepperMessenger_hh
#define RES_StepperMessenger_hh

#include "G4UImessenger.hh"

#include "globals.hh"

class RES_FieldManager;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithoutParameter;

class RES_StepperMessenger : G4UImessenger
{

public:
  explicit RES_StepperMessenger(RES_FieldManager* manager);
  ~RES_StepperMessenger();

public:
  void SetNewValue(G4UIcommand* command, G4String newValue);

private:
  RES_FieldManager*        m_manager;

  G4UIdirectory*           m_directory;
  G4UIcmdWithoutParameter* m_activateMyStepperCmd;
  G4UIcmdWithoutParameter* m_deactivateMyStepperCmd;

};

#endif /* RES_StepperMessenger_hh */
