// $Id: RES_PrimaryGeneratorMessenger.cc,v 1.6 2010/01/04 15:06:20 beischer Exp $

#include "RES_PrimaryGeneratorMessenger.hh"

#include "RES_PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

RES_PrimaryGeneratorMessenger::RES_PrimaryGeneratorMessenger(RES_PrimaryGeneratorAction* generator)
{
  m_generator = generator;

  m_directory = new G4UIdirectory("/RES/Gun/");
  m_directory->SetGuidance("Commands for the primary generator action");

  m_randomOriginCmd = new G4UIcmdWithABool("/RES/Gun/RandomOrigin", this);
  m_randomOriginCmd->SetGuidance("Set random origin for the particle?");
  m_randomOriginCmd->SetParameterName("randOrigin", true);
  m_randomOriginCmd->SetDefaultValue(true);
  m_randomOriginCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_randomDirectionCmd = new G4UIcmdWithABool("/RES/Gun/RandomDirection", this);
  m_randomDirectionCmd->SetGuidance("Set random direction for the particle?");
  m_randomDirectionCmd->SetParameterName("randDirection", true);
  m_randomDirectionCmd->SetDefaultValue(true);
  m_randomDirectionCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setStartZCmd = new G4UIcmdWithADoubleAndUnit("/RES/Gun/StartZ", this);
  m_setStartZCmd->SetGuidance("Set the height for the particle gun.");
  m_setStartZCmd->SetParameterName("z", false);
  m_setStartZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RES_PrimaryGeneratorMessenger::~RES_PrimaryGeneratorMessenger()
{
  delete m_directory;
  delete m_randomOriginCmd;
  delete m_randomDirectionCmd;
  delete m_setStartZCmd;
}

void RES_PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_randomOriginCmd) {
    m_generator->SetRandomOrigin(m_randomOriginCmd->GetNewBoolValue(newValue));
  }
  if (command == m_randomDirectionCmd) {
    m_generator->SetRandomDirection(m_randomDirectionCmd->GetNewBoolValue(newValue));
  }
  if (command == m_setStartZCmd) {
    m_generator->SetStartZ(m_setStartZCmd->GetNewDoubleValue(newValue));
  }
}
