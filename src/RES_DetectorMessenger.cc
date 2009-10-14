// $Id: RES_DetectorMessenger.cc,v 1.8 2009/10/14 16:51:31 beischer Exp $

#include "RES_DetectorMessenger.hh"

#include "RES_DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
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

  m_setModuleInternalRotationCmd = new G4UIcmdWithAString("/RES/Det/SetModuleInternalRotation", this);
  m_setModuleInternalRotationCmd->SetGuidance("Set the internal rotation for the module.");
  m_setModuleInternalRotationCmd->SetParameterName("internalRotationString", false);
  m_setModuleInternalRotationCmd->AvailableForStates(G4State_PreInit);

  m_setModuleWidthCmd = new G4UIcmdWithAString("/RES/Det/SetModuleWidth", this);
  m_setModuleWidthCmd->SetGuidance("Set the width for the module");
  m_setModuleWidthCmd->SetParameterName("widthString", false);
  m_setModuleWidthCmd->AvailableForStates(G4State_PreInit);

  m_setModuleLengthCmd = new G4UIcmdWithAString("/RES/Det/SetModuleLength", this);
  m_setModuleLengthCmd->SetGuidance("Set the module length for the given module");
  m_setModuleLengthCmd->SetParameterName("lengthString", false);
  m_setModuleLengthCmd->AvailableForStates(G4State_PreInit);

  m_setModuleFiberThicknessCmd = new G4UIcmdWithADoubleAndUnit("/RES/Det/ModuleFiberThickness", this);
  m_setModuleFiberThicknessCmd->SetGuidance("Set the module fiberThickness");
  m_setModuleFiberThicknessCmd->SetParameterName("fiberThickness", false);
  m_setModuleFiberThicknessCmd->SetDefaultUnit("cm");
  m_setModuleFiberThicknessCmd->SetUnitCategory("Length");
  m_setModuleFiberThicknessCmd->AvailableForStates(G4State_PreInit);

  m_setModuleGapCmd = new G4UIcmdWithADoubleAndUnit("/RES/Det/ModuleGap", this);
  m_setModuleGapCmd->SetGuidance("Set the module gap");
  m_setModuleGapCmd->SetParameterName("gap", false);
  m_setModuleGapCmd->SetDefaultUnit("cm");
  m_setModuleGapCmd->SetUnitCategory("Length");
  m_setModuleGapCmd->AvailableForStates(G4State_PreInit);

}

RES_DetectorMessenger::~RES_DetectorMessenger()
{
  delete m_topDirectory;
  delete m_detDirectory;
  delete m_addModulePlacementCmd;
  delete m_setModuleRotationCmd;
  delete m_setModuleInternalRotationCmd;
  delete m_setModuleWidthCmd;
  delete m_setModuleLengthCmd;
  delete m_setModuleFiberThicknessCmd;
  delete m_setModuleGapCmd;
}

void RES_DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_addModulePlacementCmd) {
    m_detector->AddModulePlacement(m_addModulePlacementCmd->GetNew3VectorValue(newValue));
  }
  if (command == m_setModuleRotationCmd) {
    G4int iModule = m_setModuleRotationCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double angle = m_setModuleRotationCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * M_PI/180.;
    m_detector->SetModuleAngle(iModule, angle);
  }
  if (command == m_setModuleInternalRotationCmd) {
    G4int iModule = m_setModuleInternalRotationCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double angle = m_setModuleInternalRotationCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * M_PI/180.;
    m_detector->SetModuleInternalAngle(iModule, angle);
  }
  if (command == m_setModuleWidthCmd) {
    G4int iModule = m_setModuleWidthCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double width = m_setModuleWidthCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * cm;
    m_detector->SetModuleWidth(iModule, width);
  }
  if (command == m_setModuleLengthCmd) {
    G4int iModule = m_setModuleLengthCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double length = m_setModuleLengthCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * cm;
    m_detector->SetModuleLength(iModule, length);
  }
  if (command == m_setModuleFiberThicknessCmd) {
    m_detector->SetModuleFiberThickness(m_setModuleFiberThicknessCmd->GetNewDoubleValue(newValue));
  }
  if (command == m_setModuleGapCmd) {
    m_detector->SetModuleGap(m_setModuleGapCmd->GetNewDoubleValue(newValue));
  }
}

