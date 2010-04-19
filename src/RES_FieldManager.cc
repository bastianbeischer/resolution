// $Id: RES_FieldManager.cc,v 1.5 2010/04/19 13:40:22 beischer Exp $

#include "RES_FieldManager.hh"

#include "RES_FieldMessenger.hh"
#include "RES_StepperMessenger.hh"
#include "RES_InhomField.hh"
#include "RES_AMS02Field.hh"
#include "RES_DummyField.hh"

#include "G4ChordFinder.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4UniformMagField.hh"

RES_FieldManager::RES_FieldManager() :
  m_myStepper(0)
{
  m_fieldMessenger = new RES_FieldMessenger(this);
  m_stepperMessenger = new RES_StepperMessenger(this);
}

RES_FieldManager::~RES_FieldManager()
{
  delete m_fieldMessenger;
  delete m_stepperMessenger;
}

void RES_FieldManager::SwitchOnInhomField(G4String dataFileName)
{
  SetDetectorField(new RES_InhomField(dataFileName));
  CreateChordFinder((G4MagneticField*) GetDetectorField());
}

void RES_FieldManager::SwitchOnAMS02Field(G4String dataFileName)
{
  SetDetectorField(new RES_AMS02Field(dataFileName));
  CreateChordFinder((G4MagneticField*) GetDetectorField());
}

void RES_FieldManager::SwitchOnDummyField(G4ThreeVector fieldVector)
{
  SetDetectorField(new RES_DummyField(fieldVector));
  CreateChordFinder((G4MagneticField*) GetDetectorField());
}

void RES_FieldManager::SwitchOnUniformField(G4ThreeVector fieldVector)
{
  SetDetectorField(new G4UniformMagField(fieldVector));
  CreateChordFinder((G4MagneticField*) GetDetectorField());
}

void RES_FieldManager::SetDisplacement(G4ThreeVector displacement)
{
  G4Field* field = (G4Field*)GetDetectorField();
  if (dynamic_cast<RES_DummyField*>(field)) {
    RES_DummyField* dummyField = (RES_DummyField*)field;
    dummyField->SetDisplacement(displacement);
  }
  else if (dynamic_cast<RES_InhomField*>(field)) {
    RES_InhomField* inhomField = (RES_InhomField*)field;
    inhomField->SetDisplacement(displacement);
  }
}

void RES_FieldManager::ActivateMyStepper()
{
  if (m_myStepper) GetChordFinder()->GetIntegrationDriver()->RenewStepperAndAdjust(m_myStepper);
  else {
    G4cerr << "ERROR: My stepper has not been set!" << G4endl;
  }
}

void RES_FieldManager::DeactivateMyStepper()
{
  CreateChordFinder((G4MagneticField*) GetDetectorField());
}
