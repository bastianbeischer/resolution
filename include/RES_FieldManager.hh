#ifndef RES_FieldManager_hh
#define RES_FieldManager_hh

#include "G4FieldManager.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class RES_FieldMessenger;
class RES_StepperMessenger;
class G4MagIntegratorStepper;

class RES_FieldManager : public G4FieldManager
{

public:
  RES_FieldManager();
  ~RES_FieldManager();

public:
  void SwitchOnInhomField(G4String dataFileName);
  void SwitchOnUniformField(G4ThreeVector fieldVector);

  void ActivateMyStepper();
  void DeactivateMyStepper();

private:
  RES_FieldMessenger*     m_fieldMessenger;
  RES_StepperMessenger*   m_stepperMessenger;

  G4MagIntegratorStepper* m_myStepper;

};

#endif /* RES_FieldManager_hh */
