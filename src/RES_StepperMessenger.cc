// $Id: RES_StepperMessenger.cc,v 1.3 2009/10/14 09:24:23 beischer Exp $

#include "RES_StepperMessenger.hh"

#include "RES_FieldManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithoutParameter.hh"


RES_StepperMessenger::RES_StepperMessenger(RES_FieldManager* manager)
{
  m_manager = manager;

  m_directory = new G4UIdirectory("/RES/MyStepper/");
  m_directory->SetGuidance("Commands to alter the integration stepper");

  m_activateMyStepperCmd = new G4UIcmdWithoutParameter("/RES/MyStepper/Activate", this);
  m_activateMyStepperCmd->SetGuidance("Activate my stepper");
  m_activateMyStepperCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_deactivateMyStepperCmd = new G4UIcmdWithoutParameter("/RES/MyStepper/Deactivatep", this);
  m_deactivateMyStepperCmd->SetGuidance("Deactivat my stepper");
  m_deactivateMyStepperCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RES_StepperMessenger::~RES_StepperMessenger()
{
  delete m_directory;
  delete m_activateMyStepperCmd;
  delete m_deactivateMyStepperCmd;
}

void RES_StepperMessenger::SetNewValue(G4UIcommand* command, G4String)
{
  if (command == m_activateMyStepperCmd) {
    m_manager->ActivateMyStepper();
  }
  if (command == m_deactivateMyStepperCmd) {
    m_manager->DeactivateMyStepper();
  }
}



