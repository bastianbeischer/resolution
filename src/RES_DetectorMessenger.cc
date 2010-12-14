// $Id: RES_DetectorMessenger.cc,v 1.20 2010/07/02 11:13:30 beischer Exp $

#include "RES_DetectorMessenger.hh"

#include "RES_DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4String.hh"

RES_DetectorMessenger::RES_DetectorMessenger(RES_DetectorConstruction* detector)
{
  m_detector = detector;

  //setup commands and directories
  m_topDirectory = new G4UIdirectory("/RES/");
  m_topDirectory->SetGuidance("Commands for the resolution software.");
  m_detDirectory = new G4UIdirectory("/RES/Det/");
  m_detDirectory->SetGuidance("Commands for the geometry specification of the detector.");

  m_setWorldXCmd = new G4UIcmdWithADoubleAndUnit("/RES/Det/SetWorldX", this);
  m_setWorldXCmd->SetGuidance("Set the length of the world in x");
  m_setWorldXCmd->SetParameterName("wx", false);
  m_setWorldXCmd->SetUnitCategory("Length");
  m_setWorldXCmd->AvailableForStates(G4State_PreInit);

  m_setWorldYCmd = new G4UIcmdWithADoubleAndUnit("/RES/Det/SetWorldY", this);
  m_setWorldYCmd->SetGuidance("Set the length of the world in y");
  m_setWorldYCmd->SetParameterName("wy", false);
  m_setWorldYCmd->SetUnitCategory("Length");
  m_setWorldYCmd->AvailableForStates(G4State_PreInit);

  m_setWorldZCmd = new G4UIcmdWithADoubleAndUnit("/RES/Det/SetWorldZ", this);
  m_setWorldZCmd->SetGuidance("Set the length of the world in z");
  m_setWorldZCmd->SetParameterName("wz", false);
  m_setWorldZCmd->SetUnitCategory("Length");
  m_setWorldZCmd->AvailableForStates(G4State_PreInit);

  m_addModulePlacementCmd = new G4UIcmdWith3VectorAndUnit("/RES/Det/AddModule", this);
  m_addModulePlacementCmd->SetGuidance("Add a module centered on the given point (relative to the world volume)");
  m_addModulePlacementCmd->SetParameterName("x","y","z",true,true);
  m_addModulePlacementCmd->SetDefaultUnit("cm");
  m_addModulePlacementCmd->SetUnitCategory("Length");
  m_addModulePlacementCmd->AvailableForStates(G4State_PreInit);

  m_printMaterialsCmd = new G4UIcmdWithoutParameter("/RES/Det/PrintMaterials", this);
  m_printMaterialsCmd->SetGuidance("Print the materials used in this geometry");
  m_printMaterialsCmd->AvailableForStates(G4State_Idle);

  m_setModuleSubtractHolesCmd = new G4UIcmdWithAString("/RES/Det/SetModuleSubtractHoles", this);
  m_setModuleSubtractHolesCmd->SetGuidance("Choose whether to subtract holes in the module");
  m_setModuleSubtractHolesCmd->SetParameterName("subtractHolesString", false);
  m_setModuleSubtractHolesCmd->AvailableForStates(G4State_PreInit);

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
  m_setModuleUpperSigmaUCmd->SetGuidance("Specify the resolution along the fiber (in um!) - upper layer.");
  m_setModuleUpperSigmaUCmd->SetParameterName("sigmaUString", false);
  m_setModuleUpperSigmaUCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setModuleUpperSigmaVCmd = new G4UIcmdWithAString("/RES/Det/SetModuleUpperSigmaV", this);
  m_setModuleUpperSigmaVCmd->SetGuidance("Specify the resolution perpendicular to the fiber (in um!) - upper layer.");
  m_setModuleUpperSigmaVCmd->SetParameterName("sigmaVString", false);
  m_setModuleUpperSigmaVCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setModuleUpperSigmaZCmd = new G4UIcmdWithAString("/RES/Det/SetModuleUpperSigmaZ", this);
  m_setModuleUpperSigmaZCmd->SetGuidance("Specify the resolution in Z (in um!) - upper layer.");
  m_setModuleUpperSigmaZCmd->SetParameterName("sigmaZString", false);
  m_setModuleUpperSigmaZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  m_setModuleLowerSigmaUCmd = new G4UIcmdWithAString("/RES/Det/SetModuleLowerSigmaU", this);
  m_setModuleLowerSigmaUCmd->SetGuidance("Specify the resolution along the fiber (in um!) - lower layer.");
  m_setModuleLowerSigmaUCmd->SetParameterName("sigmaUString", false);
  m_setModuleLowerSigmaUCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setModuleLowerSigmaVCmd = new G4UIcmdWithAString("/RES/Det/SetModuleLowerSigmaV", this);
  m_setModuleLowerSigmaVCmd->SetGuidance("Specify the resolution perpendicular to the the fiber (in um!) - lower layer.");
  m_setModuleLowerSigmaVCmd->SetParameterName("sigmaVString", false);
  m_setModuleLowerSigmaVCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setModuleLowerSigmaZCmd = new G4UIcmdWithAString("/RES/Det/SetModuleLowerSigmaZ", this);
  m_setModuleLowerSigmaZCmd->SetGuidance("Specify the resolution in Z (in um!) - lower layer.");
  m_setModuleLowerSigmaZCmd->SetParameterName("sigmaZString", false);
  m_setModuleLowerSigmaZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setModuleUpperEfficiencyCmd = new G4UIcmdWithAString("/RES/Det/SetModuleUpperEfficiency", this);
  m_setModuleUpperEfficiencyCmd->SetGuidance("Specify the effiency for the upper layer of this module (number between 0 and 1).");
  m_setModuleUpperEfficiencyCmd->SetParameterName("upperEfficiencyString", false);
  m_setModuleUpperEfficiencyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setModuleLowerEfficiencyCmd = new G4UIcmdWithAString("/RES/Det/SetModuleLowerEfficiency", this);
  m_setModuleLowerEfficiencyCmd->SetGuidance("Specify the effiency for the lower layer of this module (number between 0 and 1).");
  m_setModuleLowerEfficiencyCmd->SetParameterName("lowerEfficiencyString", false);
  m_setModuleLowerEfficiencyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setLayerSubtractHolesCmd = new G4UIcmdWithAString("/RES/Det/SetLayerSubtractHoles", this);
  m_setLayerSubtractHolesCmd->SetGuidance("Choose whether to subtract holes in the layer");
  m_setLayerSubtractHolesCmd->SetParameterName("subtractHolesString", false);
  m_setLayerSubtractHolesCmd->AvailableForStates(G4State_PreInit);

  m_setLayerTypeCmd = new G4UIcmdWithAString("/RES/Det/SetLayerType", this);
  m_setLayerTypeCmd->SetGuidance("Set the type of the layer (silicon ladder/fiber layer)");
  m_setLayerTypeCmd->SetParameterName("typeString", false);
  m_setLayerTypeCmd->AvailableForStates(G4State_PreInit);

  m_setLayerRotationCmd = new G4UIcmdWithAString("/RES/Det/SetLayerRotation", this);
  m_setLayerRotationCmd->SetGuidance("Set the rotation for the layer.");
  m_setLayerRotationCmd->SetParameterName("rotationString", false);
  m_setLayerRotationCmd->AvailableForStates(G4State_PreInit);

  m_setLayerInternalRotationCmd = new G4UIcmdWithAString("/RES/Det/SetLayerInternalRotation", this);
  m_setLayerInternalRotationCmd->SetGuidance("Set the internal rotation for the layer.");
  m_setLayerInternalRotationCmd->SetParameterName("internalRotationString", false);
  m_setLayerInternalRotationCmd->AvailableForStates(G4State_PreInit);

  m_setLayerWidthCmd = new G4UIcmdWithAString("/RES/Det/SetLayerWidth", this);
  m_setLayerWidthCmd->SetGuidance("Set the width for the layer");
  m_setLayerWidthCmd->SetParameterName("widthString", false);
  m_setLayerWidthCmd->AvailableForStates(G4State_PreInit);

  m_setLayerLengthCmd = new G4UIcmdWithAString("/RES/Det/SetLayerLength", this);
  m_setLayerLengthCmd->SetGuidance("Set the layer length for the given layer");
  m_setLayerLengthCmd->SetParameterName("lengthString", false);
  m_setLayerLengthCmd->AvailableForStates(G4State_PreInit);

  m_setLayerUpperSigmaUCmd = new G4UIcmdWithAString("/RES/Det/SetLayerUpperSigmaU", this);
  m_setLayerUpperSigmaUCmd->SetGuidance("Specify the resolution along the fiber (in um!) - upper layer.");
  m_setLayerUpperSigmaUCmd->SetParameterName("sigmaUString", false);
  m_setLayerUpperSigmaUCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setLayerUpperSigmaVCmd = new G4UIcmdWithAString("/RES/Det/SetLayerUpperSigmaV", this);
  m_setLayerUpperSigmaVCmd->SetGuidance("Specify the resolution perpendicular to the fiber (in um!) - upper layer.");
  m_setLayerUpperSigmaVCmd->SetParameterName("sigmaVString", false);
  m_setLayerUpperSigmaVCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setLayerUpperSigmaZCmd = new G4UIcmdWithAString("/RES/Det/SetLayerUpperSigmaZ", this);
  m_setLayerUpperSigmaZCmd->SetGuidance("Specify the resolution in Z (in um!) - upper layer.");
  m_setLayerUpperSigmaZCmd->SetParameterName("sigmaZString", false);
  m_setLayerUpperSigmaZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  m_setLayerLowerSigmaUCmd = new G4UIcmdWithAString("/RES/Det/SetLayerLowerSigmaU", this);
  m_setLayerLowerSigmaUCmd->SetGuidance("Specify the resolution along the fiber (in um!) - lower layer.");
  m_setLayerLowerSigmaUCmd->SetParameterName("sigmaUString", false);
  m_setLayerLowerSigmaUCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setLayerLowerSigmaVCmd = new G4UIcmdWithAString("/RES/Det/SetLayerLowerSigmaV", this);
  m_setLayerLowerSigmaVCmd->SetGuidance("Specify the resolution perpendicular to the the fiber (in um!) - lower layer.");
  m_setLayerLowerSigmaVCmd->SetParameterName("sigmaVString", false);
  m_setLayerLowerSigmaVCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setLayerLowerSigmaZCmd = new G4UIcmdWithAString("/RES/Det/SetLayerLowerSigmaZ", this);
  m_setLayerLowerSigmaZCmd->SetGuidance("Specify the resolution in Z (in um!) - lower layer.");
  m_setLayerLowerSigmaZCmd->SetParameterName("sigmaZString", false);
  m_setLayerLowerSigmaZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setLayerUpperEfficiencyCmd = new G4UIcmdWithAString("/RES/Det/SetLayerUpperEfficiency", this);
  m_setLayerUpperEfficiencyCmd->SetGuidance("Specify the effiency for the upper layer of this layer (number between 0 and 1).");
  m_setLayerUpperEfficiencyCmd->SetParameterName("upperEfficiencyString", false);
  m_setLayerUpperEfficiencyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setLayerLowerEfficiencyCmd = new G4UIcmdWithAString("/RES/Det/SetLayerLowerEfficiency", this);
  m_setLayerLowerEfficiencyCmd->SetGuidance("Specify the effiency for the lower layer of this layer (number between 0 and 1).");
  m_setLayerLowerEfficiencyCmd->SetParameterName("lowerEfficiencyString", false);
  m_setLayerLowerEfficiencyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_addLayerCmd = new G4UIcmdWithAString("/RES/Det/AddLayer", this);
  m_addLayerCmd->SetGuidance("Add a layer centered on the given point (relative to the world volume)");
  m_addLayerCmd->SetParameterName("addLayerString",true,true);
  m_addLayerCmd->AvailableForStates(G4State_PreInit);
}

RES_DetectorMessenger::~RES_DetectorMessenger()
{
  delete m_topDirectory;
  delete m_detDirectory;
  delete m_setWorldXCmd;
  delete m_setWorldYCmd;
  delete m_setWorldZCmd;
  delete m_addModulePlacementCmd;
  delete m_printMaterialsCmd;
  delete m_setModuleSubtractHolesCmd;
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
  delete m_setModuleUpperEfficiencyCmd;
  delete m_setModuleLowerEfficiencyCmd;
  delete m_setLayerSubtractHolesCmd;
  delete m_setLayerTypeCmd;
  delete m_setLayerRotationCmd;
  delete m_setLayerInternalRotationCmd;
  delete m_setLayerWidthCmd;
  delete m_setLayerLengthCmd;
  delete m_setLayerUpperSigmaUCmd;
  delete m_setLayerUpperSigmaVCmd;
  delete m_setLayerUpperSigmaZCmd;
  delete m_setLayerLowerSigmaUCmd;
  delete m_setLayerLowerSigmaVCmd;
  delete m_setLayerLowerSigmaZCmd;
  delete m_setLayerUpperEfficiencyCmd;
  delete m_setLayerLowerEfficiencyCmd;
  delete m_addLayerCmd;
}

void RES_DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_addModulePlacementCmd) {
    m_detector->AddModulePlacement(m_addModulePlacementCmd->GetNew3VectorValue(newValue));
  }

  if (command == m_printMaterialsCmd) {
    m_detector->PrintMaterials();
  }

  if (command == m_setWorldXCmd) {
    m_detector->SetWorldX(m_setWorldXCmd->GetNewDoubleValue(newValue));
  }

  if (command == m_setWorldYCmd) {
    m_detector->SetWorldY(m_setWorldYCmd->GetNewDoubleValue(newValue));
  }

  if (command == m_setWorldZCmd) {
    m_detector->SetWorldZ(m_setWorldZCmd->GetNewDoubleValue(newValue));
  }

  if (command == m_setModuleRotationCmd) {
    G4int iModule = m_setModuleRotationCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double angle = m_setModuleRotationCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * M_PI/180.;
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetAngle(angle);
  }

  if (command == m_setModuleSubtractHolesCmd) {
    G4int iModule = m_setModuleSubtractHolesCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4bool value = m_setModuleSubtractHolesCmd->ConvertToBool(newValue.substr(2,newValue.length()).c_str());
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetSubtractHoles(value);
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

  if (command == m_setModuleUpperEfficiencyCmd) {
    G4int iModule = m_setModuleUpperEfficiencyCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double eff = m_setModuleUpperEfficiencyCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str());
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetUpperEfficiency(eff);
  }

  if (command == m_setModuleLowerEfficiencyCmd) {
    G4int iModule = m_setModuleLowerEfficiencyCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double eff = m_setModuleLowerEfficiencyCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str());
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetLowerEfficiency(eff);
  }

  if (command == m_setModuleRotationCmd) {
    G4int iModule = m_setModuleRotationCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double angle = m_setModuleRotationCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * M_PI/180.;
    RES_Module* module = m_detector->GetModule(iModule);
    module->SetAngle(angle);
  }

  if (command == m_setLayerSubtractHolesCmd) {
    G4int iLayer = m_setLayerSubtractHolesCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4bool value = m_setLayerSubtractHolesCmd->ConvertToBool(newValue.substr(2,newValue.length()).c_str());
    RES_Layer* layer = m_detector->GetLayer(iLayer);
    layer->SetSubtractHoles(value);
  }

  if (command == m_setLayerTypeCmd) {
    G4int iLayer = m_setLayerTypeCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4String value = newValue.substr(2,newValue.length());
    value.toLower();
    RES_Module::ModuleType type;
    if (value == "silicon")
      type = RES_Module::silicon;
    if (value == "fiber")
      type = RES_Module::fiber;
    RES_Layer* layer = m_detector->GetLayer(iLayer);
    layer->SetType(type);
  }

  if (command == m_setLayerInternalRotationCmd) {
    G4int iLayer = m_setLayerInternalRotationCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double angle = m_setLayerInternalRotationCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * M_PI/180.;
    RES_Layer* layer = m_detector->GetLayer(iLayer);
    layer->SetInternalAngle(angle);
  }

  if (command == m_setLayerWidthCmd) {
    G4int iLayer = m_setLayerWidthCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double width = m_setLayerWidthCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * cm;
    RES_Layer* layer = m_detector->GetLayer(iLayer);
    layer->SetWidth(width);
  }

  if (command == m_setLayerLengthCmd) {
    G4int iLayer = m_setLayerLengthCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double length = m_setLayerLengthCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * cm;
    RES_Layer* layer = m_detector->GetLayer(iLayer);
    layer->SetLength(length);
  }

  if (command == m_setLayerUpperSigmaUCmd) {
    G4int iLayer = m_setLayerUpperSigmaUCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double sigmaU = m_setLayerUpperSigmaUCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * um;
    RES_Layer* layer = m_detector->GetLayer(iLayer);
    layer->SetUpperSigmaU(sigmaU);
  }
  if (command == m_setLayerUpperSigmaVCmd) {
    G4int iLayer = m_setLayerUpperSigmaVCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double sigmaV = m_setLayerUpperSigmaVCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * um;
    RES_Layer* layer = m_detector->GetLayer(iLayer);
    layer->SetUpperSigmaV(sigmaV);
  }

  if (command == m_setLayerUpperSigmaZCmd) {
    G4int iLayer = m_setLayerUpperSigmaZCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double sigmaZ = m_setLayerUpperSigmaZCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * um;
    RES_Layer* layer = m_detector->GetLayer(iLayer);
    layer->SetUpperSigmaZ(sigmaZ);
  }

  if (command == m_setLayerLowerSigmaUCmd) {
    G4int iLayer = m_setLayerLowerSigmaUCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double sigmaU = m_setLayerLowerSigmaUCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * um;
    RES_Layer* layer = m_detector->GetLayer(iLayer);
    layer->SetLowerSigmaU(sigmaU);
  }

  if (command == m_setLayerLowerSigmaVCmd) {
    G4int iLayer = m_setLayerLowerSigmaVCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double sigmaV = m_setLayerLowerSigmaVCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * um;
    RES_Layer* layer = m_detector->GetLayer(iLayer);
    layer->SetLowerSigmaV(sigmaV);
  }

  if (command == m_setLayerLowerSigmaZCmd) {
    G4int iLayer = m_setLayerLowerSigmaZCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double sigmaZ = m_setLayerLowerSigmaZCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * um;
    RES_Layer* layer = m_detector->GetLayer(iLayer);
    layer->SetLowerSigmaZ(sigmaZ);
  }

  if (command == m_setLayerUpperEfficiencyCmd) {
    G4int iLayer = m_setLayerUpperEfficiencyCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double eff = m_setLayerUpperEfficiencyCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str());
    RES_Layer* layer = m_detector->GetLayer(iLayer);
    layer->SetUpperEfficiency(eff);
  }

  if (command == m_setLayerLowerEfficiencyCmd) {
    G4int iLayer = m_setLayerLowerEfficiencyCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double eff = m_setLayerLowerEfficiencyCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str());
    RES_Layer* layer = m_detector->GetLayer(iLayer);
    layer->SetLowerEfficiency(eff);
  }

  if (command == m_addLayerCmd) {
    G4int nModules = m_addLayerCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4double z = m_addLayerCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * cm;
    m_detector->AddLayer(nModules, z);
  }
}

