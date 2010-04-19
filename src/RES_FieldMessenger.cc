// $Id: RES_FieldMessenger.cc,v 1.6 2010/04/19 13:40:21 beischer Exp $

#include "RES_FieldMessenger.hh"

#include "RES_FieldManager.hh"

#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

RES_FieldMessenger::RES_FieldMessenger(RES_FieldManager* manager)
{
  m_manager = manager;

  m_directory = new G4UIdirectory("/RES/Field/");
  m_directory->SetGuidance("Commands for the magnetic field");

  m_setInhomFieldFromFileCmd = new G4UIcmdWithAString("/RES/Field/SetInhomFieldFrom", this);
  m_setInhomFieldFromFileCmd->SetGuidance("Read the field map from the specified data file");
  m_setInhomFieldFromFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setAMS02FieldFromFileCmd = new G4UIcmdWithAString("/RES/Field/SetAMS02FieldFrom", this);
  m_setAMS02FieldFromFileCmd->SetGuidance("Read the field map from the specified AMS02 data file");
  m_setAMS02FieldFromFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setDummyFieldCmd = new G4UIcmdWith3VectorAndUnit("/RES/Field/SetDummyField", this);
  m_setDummyFieldCmd->SetGuidance("Set a dummy magnetic field with the value given by the specified ThreeVector");
  m_setDummyFieldCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setUniformFieldCmd = new G4UIcmdWith3VectorAndUnit("/RES/Field/SetUniformField", this);
  m_setUniformFieldCmd->SetGuidance("Set a uniform magnetic field with the value given by the specified ThreeVector");
  m_setUniformFieldCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setDisplacementCmd = new G4UIcmdWith3VectorAndUnit("/RES/Field/SetDisplacement", this);
  m_setDisplacementCmd->SetGuidance("Set displacement vector of the magnet.");
  m_setDisplacementCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RES_FieldMessenger::~RES_FieldMessenger()
{
  delete m_directory;
  delete m_setInhomFieldFromFileCmd;
  delete m_setAMS02FieldFromFileCmd;
  delete m_setDummyFieldCmd;
  delete m_setUniformFieldCmd;
  delete m_setDisplacementCmd;
}

void RES_FieldMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_setInhomFieldFromFileCmd) {
    m_manager->SwitchOnInhomField(newValue);
  }
  if (command == m_setAMS02FieldFromFileCmd) {
    m_manager->SwitchOnAMS02Field(newValue);
  }
  if (command == m_setDummyFieldCmd) {
    m_manager->SwitchOnDummyField(m_setDummyFieldCmd->GetNew3VectorValue(newValue));
  }
  if (command == m_setUniformFieldCmd) {
    m_manager->SwitchOnUniformField(m_setUniformFieldCmd->GetNew3VectorValue(newValue));
  }
  if (command == m_setDummyFieldCmd) {
    m_manager->SetDisplacement(m_setDisplacementCmd->GetNew3VectorValue(newValue));
  }
}
