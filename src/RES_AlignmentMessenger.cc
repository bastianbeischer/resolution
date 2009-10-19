// $Id: RES_AlignmentMessenger.cc,v 1.3 2009/10/19 12:27:44 beischer Exp $

#include "RES_AlignmentMessenger.hh"

#include "RES_AlignmentManager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"


RES_AlignmentMessenger::RES_AlignmentMessenger(RES_AlignmentManager* manager)
{
  m_manager = manager;

  m_directory = new G4UIdirectory("/RES/Alignment/");
  m_directory->SetGuidance("Commands for the alignment manager");

  m_startAlignmentCmd = new G4UIcmdWithoutParameter("/RES/Alignment/Start", this);
  m_startAlignmentCmd->SetGuidance("Try to find the alignment parameters");
  m_startAlignmentCmd->AvailableForStates(G4State_Idle);

  m_verboseCmd = new G4UIcmdWithAnInteger("/RES/Alignment/Verbose", this);
  m_verboseCmd->SetGuidance("Set verbosity of the alignment part of the software");
  m_verboseCmd->SetParameterName("verbose", false);
  m_verboseCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RES_AlignmentMessenger::~RES_AlignmentMessenger()
{
  delete m_directory;
  delete m_startAlignmentCmd;
  delete m_verboseCmd;
}


void RES_AlignmentMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
  if (cmd == m_startAlignmentCmd) {
    m_manager->StartAlignment();
  }
  if (cmd == m_verboseCmd) {
    m_manager->SetVerbose(m_verboseCmd->GetNewIntValue(newValue));
  }
}
