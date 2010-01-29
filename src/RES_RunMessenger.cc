// $Id: RES_RunMessenger.cc,v 1.7 2010/01/29 12:51:38 beischer Exp $

#include "RES_RunMessenger.hh"

#include "RES_RunManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"

RES_RunMessenger::RES_RunMessenger(RES_RunManager* manager)
{
  m_manager = manager;

  m_directory = new G4UIdirectory("/RES/Run/");
  m_directory->SetGuidance("Commands for the run manager");

  m_setStoreResultsCmd = new G4UIcmdWithABool("/RES/Run/StoreResults", this);
  m_setStoreResultsCmd->SetGuidance("Store results?");
  m_setStoreResultsCmd->SetParameterName("store", true);
  m_setStoreResultsCmd->SetDefaultValue(true);
  m_setStoreResultsCmd->AvailableForStates(G4State_PreInit);

  m_generateCmd = new G4UIcmdWithAnInteger("/RES/Run/Generate", this);
  m_generateCmd->SetGuidance("Start a generation run");
  m_generateCmd->SetParameterName("nEvents", true);
  m_generateCmd->SetDefaultValue(1);
  m_generateCmd->AvailableForStates(G4State_Idle);

  m_reconstructCmd = new G4UIcmdWithoutParameter("/RES/Run/Reconstruct", this);
  m_reconstructCmd->SetGuidance("Reconstruct the events which are stored in the DataHandler");
  m_reconstructCmd->AvailableForStates(G4State_Idle);

  m_reconstructWithoutLayerCmd = new G4UIcmdWithAnInteger("/RES/Run/ReconstructWithoutLayer", this);
  m_reconstructWithoutLayerCmd->SetParameterName("layer", false);
  m_reconstructWithoutLayerCmd->SetGuidance("Reconstruct the events which are stored in the DataHandler without the given layer");
  m_reconstructWithoutLayerCmd->AvailableForStates(G4State_Idle);

  m_scanChi2FuncCmd = new G4UIcmdWithAString("/RES/Run/ScanChi2", this);
  m_scanChi2FuncCmd->SetGuidance("Scan the chi2 function and write the output to the given file");
  m_scanChi2FuncCmd->SetParameterName("filename", false);
  m_scanChi2FuncCmd->AvailableForStates(G4State_Idle);

  m_fixedDofCmd = new G4UIcmdWithAnInteger("/RES/Run/SetFixedDof", this);
  m_fixedDofCmd->SetGuidance("Fix the degrees of freedom for event generation (random layers will be removed)");
  m_fixedDofCmd->SetParameterName("dof", false);
  m_fixedDofCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RES_RunMessenger::~RES_RunMessenger()
{
  delete m_directory;
  delete m_setStoreResultsCmd;
  delete m_generateCmd;
  delete m_reconstructCmd;
  delete m_reconstructWithoutLayerCmd;
  delete m_scanChi2FuncCmd;
  delete m_fixedDofCmd;
}

void RES_RunMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_setStoreResultsCmd) {
    m_manager->SetStoreResults(m_setStoreResultsCmd->GetNewBoolValue(newValue));
  }
  if (command == m_generateCmd) {
    m_manager->StartGenerationRun(m_generateCmd->GetNewIntValue(newValue));
  }
  if (command == m_reconstructCmd) {
    m_manager->StartReconstructionRun();
  }
  if (command == m_reconstructWithoutLayerCmd) {
    G4int layer = m_reconstructWithoutLayerCmd->GetNewIntValue(newValue);
    m_manager->StartReconstructionRunWithoutLayer(layer);
  }
  if (command == m_scanChi2FuncCmd) {
    G4int iPar = m_scanChi2FuncCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4int jPar = m_scanChi2FuncCmd->ConvertToInt(newValue.substr(2,3).c_str());
    G4String filename = newValue.substr(4).c_str();
    m_manager->ScanChi2Function(iPar, jPar, filename);
  }
  if (command == m_fixedDofCmd) {
    G4int value = m_fixedDofCmd->GetNewIntValue(newValue);
    m_manager->SetFixedDof(value);
  }
}
