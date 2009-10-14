// $Id: RES_AlignmentMessenger.cc,v 1.1 2009/10/14 16:51:31 beischer Exp $

#include "RES_AlignmentMessenger.hh"

#include "RES_AlignmentManager.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithoutParameter.hh"

RES_AlignmentMessenger::RES_AlignmentMessenger(RES_AlignmentManager* manager)
{
  m_manager = manager;

  m_startAlignmentCmd = new G4UIcmdWithoutParameter("/RES/Alignment/Start", this);
  m_startAlignmentCmd->AvailableForStates(G4State_Init);
}

RES_AlignmentMessenger::~RES_AlignmentMessenger()
{
  delete m_startAlignmentCmd;
}

void RES_AlignmentMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
  if (cmd == m_startAlignmentCmd) {
    m_manager->StartAlignment();
  }
}
