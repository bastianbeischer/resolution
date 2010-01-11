// $Id: RES_DetectorMessenger.cc,v 1.13 2010/01/11 09:59:58 beischer Exp $

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

  m_setModuleTypeCmd = new G4UIcmdWithAString("/RES/Det/SetModuleType", this);
  m_setModuleTypeCmd->SetGuidance("Set the type of the module (silicon ladder/fiber module)");
  m_setModuleTypeCmd->SetParameterName("typeString", false);
  m_setModuleTypeCmd->AvailableForStates(G4State_PreInit);

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

  m_setModuleUpperSigmaUCmd = new G4UIcmdWithAString("/RES/Det/SetModuleUpperSigmaU", this);
  m_setModuleUpperSigmaUCmd->SetGuidance("Specify the resolution along the fiber (in um!).");
  m_setModuleUpperSigmaUCmd->SetParameterName("sigmaUString", false);
  m_setModuleUpperSigmaUCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setModuleUpperSigmaVCmd = new G4UIcmdWithAString("/RES/Det/SetModuleUpperSigmaV", this);
  m_setModuleUpperSigmaVCmd->SetGuidance("Specify the resolution along the fiber (in um!).");
  m_setModuleUpperSigmaVCmd->SetParameterName("sigmaVString", false);
  m_setModuleUpperSigmaVCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setModuleUpperSigmaZCmd = new G4UIcmdWithAString("/RES/Det/SetModuleUpperSigmaZ", this);
  m_setModuleUpperSigmaZCmd->SetGuidance("Specify the resolution along the fiber (in um!).");
  m_setModuleUpperSigmaZCmd->SetParameterName("sigmaZString", false);
  m_setModuleUpperSigmaZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  m_setModuleLowerSigmaUCmd = new G4UIcmdWithAString("/RES/Det/SetModuleLowerSigmaU", this);
  m_setModuleLowerSigmaUCmd->SetGuidance("Specify the resolution along the fiber (in um!).");
  m_setModuleLowerSigmaUCmd->SetParameterName("sigmaUString", false);
  m_setModuleLowerSigmaUCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setModuleLowerSigmaVCmd = new G4UIcmdWithAString("/RES/Det/SetModuleLowerSigmaV", this);
  m_setModuleLowerSigmaVCmd->SetGuidance("Specify the resolution along the fiber (in um!).");
  m_setModuleLowerSigmaVCmd->SetParameterName("sigmaVString", false);
  m_setModuleLowerSigmaVCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setModuleLowerSigmaZCmd = new G4UIcmdWithAString("/RES/Det/SetModuleLowerSigmaZ", this);
  m_setModuleLowerSigmaZCmd->SetGuidance("Specify the resolution along the fiber (in um!).");
  m_setModuleLowerSigmaZCmd->SetParameterName("sigmaZString", false);
  m_setModuleLowerSigmaZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RES_DetectorMessenger::~RES_DetectorMessenger()
{
  delete m_topDirectory;
  delete m_detDirectory;
  delete m_addModulePlacementCmd;
  delete m_setModuleTypeCmd;
  delete m_setModuleRotationCmd;
  delete m_setModuleInternalRotationCmd;
  delete m_setModuleWidthCmd;
  delete m_setModuleLengthCmd;
  delete m_setModuleUpperSigmaUCmd;
  delete m_setModuleUpperSigmaVCmd;
  delete m_setModuleUpperSigmaZCmd;
  delete m_setModuleLowerSigmaUCmd;
  delete m_setModuleLowerSigmaVCmd;
  delete m_setModuleLowerSigmaZCmd;
}

void RES_DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_addModulePlacementCmd) {
    m_detector->AddModulePlacement(m_addModulePlacementCmd->GetNew3VectorValue(newValue));
  }

  if (command == m_setModuleRotationCmd) {
    G4int iModule = m_setModuleRotationCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double angle = m_setModuleRotationCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * M_PI/180.;
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetAngle(angle);
  }

  if (command == m_setModuleTypeCmd) {
    G4int iModule = m_setModuleTypeCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4String value = newValue.substr(2,newValue.length());
    value.toLower();
    RES_Module::ModuleType type;
    if (value == "silicon")
      type = RES_Module::silicon;
    if (value == "fiber")
      type = RES_Module::fiber;
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetType(type);
  }

  if (command == m_setModuleInternalRotationCmd) {
    G4int iModule = m_setModuleInternalRotationCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double angle = m_setModuleInternalRotationCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * M_PI/180.;
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetInternalAngle(angle);
  }

  if (command == m_setModuleWidthCmd) {
    G4int iModule = m_setModuleWidthCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double width = m_setModuleWidthCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * cm;
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetWidth(width);
  }

  if (command == m_setModuleLengthCmd) {
    G4int iModule = m_setModuleLengthCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double length = m_setModuleLengthCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * cm;
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetLength(length);
  }

  if (command == m_setModuleUpperSigmaUCmd) {
    G4int iModule = m_setModuleUpperSigmaUCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double sigmaU = m_setModuleUpperSigmaUCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * um;
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetUpperSigmaU(sigmaU);
  }
  if (command == m_setModuleUpperSigmaVCmd) {
    G4int iModule = m_setModuleUpperSigmaVCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double sigmaV = m_setModuleUpperSigmaVCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * um;
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetUpperSigmaV(sigmaV);
  }
  if (command == m_setModuleUpperSigmaZCmd) {
    G4int iModule = m_setModuleUpperSigmaZCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double sigmaZ = m_setModuleUpperSigmaZCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * um;
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetUpperSigmaZ(sigmaZ);
  }
  if (command == m_setModuleLowerSigmaUCmd) {
    G4int iModule = m_setModuleLowerSigmaUCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double sigmaU = m_setModuleLowerSigmaUCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * um;
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetLowerSigmaU(sigmaU);
  }
  if (command == m_setModuleLowerSigmaVCmd) {
    G4int iModule = m_setModuleLowerSigmaVCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double sigmaV = m_setModuleLowerSigmaVCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * um;
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetLowerSigmaV(sigmaV);
  }
  if (command == m_setModuleLowerSigmaZCmd) {
    G4int iModule = m_setModuleLowerSigmaZCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double sigmaZ = m_setModuleLowerSigmaZCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * um;
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetLowerSigmaZ(sigmaZ);
  }
}

