// $Id: RES_AlignmentMessenger.cc,v 1.5 2009/11/09 09:27:21 beischer Exp $

#include "RES_AlignmentMessenger.hh"

#include "RES_AlignmentManager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"


RES_AlignmentMessenger::RES_AlignmentMessenger(RES_AlignmentManager* manager)
{
  m_manager = manager;

  m_directory = new G4UIdirectory("/RES/Alignment/");
  m_directory->SetGuidance("Commands for the alignment manager");

  m_startAlignmentCmd = new G4UIcmdWithoutParameter("/RES/Alignment/Start", this);
  m_startAlignmentCmd->SetGuidance("Try to find the alignment parameters");
  m_startAlignmentCmd->AvailableForStates(G4State_Idle);

  m_setXshiftCmd = new G4UIcmdWithAString("/RES/Alignment/SetXshiftForModule", this);
  m_setXshiftCmd->SetGuidance("Set misalignment shift in x for module");
  m_setXshiftCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setYshiftCmd = new G4UIcmdWithAString("/RES/Alignment/SetYshiftForModule", this);
  m_setYshiftCmd->SetGuidance("Set misalignment shift in y for module");
  m_setYshiftCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setAngleShiftCmd = new G4UIcmdWithAString("/RES/Alignment/SetAngleShiftForModule", this);
  m_setAngleShiftCmd->SetGuidance("Set misalignment angle for module");
  m_setAngleShiftCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_verboseCmd = new G4UIcmdWithAnInteger("/RES/Alignment/Verbose", this);
  m_verboseCmd->SetGuidance("Set verbosity of the alignment part of the software");
  m_verboseCmd->SetParameterName("verbose", false);
  m_verboseCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RES_AlignmentMessenger::~RES_AlignmentMessenger()
{
  delete m_directory;
  delete m_startAlignmentCmd;
  delete m_setXshiftCmd;
  delete m_setYshiftCmd;
  delete m_setAngleShiftCmd;
  delete m_verboseCmd;
}


void RES_AlignmentMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
  if (cmd == m_startAlignmentCmd) {
    m_manager->StartAlignment();
  }
  if (cmd == m_setXshiftCmd) {
    G4int iModule = m_setXshiftCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4float shift = m_setXshiftCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * mm;
    m_manager->SetXshift(iModule, shift);
  }
  if (cmd == m_setYshiftCmd) {
    G4int iModule = m_setYshiftCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4float shift = m_setYshiftCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * mm;
    m_manager->SetYshift(iModule, shift);
  }
  if (cmd == m_setAngleShiftCmd) {
    G4int iModule = m_setAngleShiftCmd->ConvertToInt(newValue.substr(0,1).c_str());
    G4float shift = m_setAngleShiftCmd->ConvertToDouble(newValue.substr(2,newValue.length()).c_str()) * mm;
    m_manager->SetAngleShift(iModule, shift);
  }
  if (cmd == m_verboseCmd) {
    m_manager->SetVerbose(m_verboseCmd->GetNewIntValue(newValue));
  }
}
