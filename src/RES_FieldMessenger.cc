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

  m_setUniformFieldCmd = new G4UIcmdWith3VectorAndUnit("/RES/Field/SetUniformField", this);
  m_setUniformFieldCmd->SetGuidance("Set a uniform magnetic field with the value given by the specified ThreeVector");
  m_setUniformFieldCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RES_FieldMessenger::~RES_FieldMessenger()
{
  delete m_directory;
}

void RES_FieldMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_setInhomFieldFromFileCmd) {
    m_manager->SwitchOnInhomField(newValue);
  }
  if (command == m_setUniformFieldCmd) {
    m_manager->SwitchOnUniformField(m_setUniformFieldCmd->GetNew3VectorValue(newValue));
  }
}
