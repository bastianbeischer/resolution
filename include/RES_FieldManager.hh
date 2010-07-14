// $Id: RES_FieldManager.hh,v 1.6 2010/07/14 12:17:35 beischer Exp $

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
  void SwitchOnAMS02Field(G4String dataFileName);
  void SwitchOnUniformField(G4ThreeVector fieldVector);

  void SetDisplacement(G4ThreeVector displacement);

  void ActivateMyStepper();
  void DeactivateMyStepper();

private:
  RES_FieldMessenger*     m_fieldMessenger;
  RES_StepperMessenger*   m_stepperMessenger;

  G4MagIntegratorStepper* m_myStepper;

};

#endif /* RES_FieldManager_hh */
