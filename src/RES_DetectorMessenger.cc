#include "RES_DetectorMessenger.hh"

#include "RES_DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4String.hh"

RES_DetectorMessenger::RES_DetectorMessenger(RES_DetectorConstruction* detector)
{
  m_detector = detector;

  //setup commands and directories
  m_topDirectory = new G4UIdirectory("/RES/");
  m_topDirectory->SetGuidance("Commands for the resolution software.");
  m_detDirectory = new G4UIdirectory("/RES/Det/");
  m_detDirectory->SetGuidance("Commands for the geometry specification of the detector.");

  m_addModulePlacementCmd = new G4UIcmdWith3VectorAndUnit("/RES/Det/AddModule", this);
  m_addModulePlacementCmd->SetGuidance("Add a module centered on the given point (relative to the world volume)");
  m_addModulePlacementCmd->SetParameterName("x","y","z",true,true);
  m_addModulePlacementCmd->SetDefaultUnit("cm");
  m_addModulePlacementCmd->SetUnitCategory("Length");
  m_addModulePlacementCmd->AvailableForStates(G4State_PreInit);

  m_setModuleRotationCmd = new G4UIcmdWithAString("/RES/Det/SetModuleRotation", this);
  m_setModuleRotationCmd->SetGuidance("Set the rotation for the module.");
  m_setModuleRotationCmd->SetParameterName("rotationString", false);
  m_setModuleRotationCmd->AvailableForStates(G4State_PreInit);
}

RES_DetectorMessenger::~RES_DetectorMessenger()
{
  delete m_topDirectory;
  delete m_detDirectory;
  delete m_addModulePlacementCmd;
  delete m_setModuleRotationCmd;
}

void RES_DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_addModulePlacementCmd) {
    m_detector->AddModulePlacement(m_addModulePlacementCmd->GetNew3VectorValue(newValue));
  }
  if (command == m_setModuleRotationCmd) {
    G4int iModule = m_setModuleRotationCmd->ConvertToInt(newValue.strip(0));
    G4double angle = m_setModuleRotationCmd->ConvertToDouble(newValue.strip(1)) * M_PI/180.;;
    m_detector->SetModuleAngle(iModule, angle);
  }
}

