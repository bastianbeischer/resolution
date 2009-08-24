#include "RES_TrackFitMessenger.hh"

#include "RES_TrackFitter.hh"

#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"

RES_TrackFitMessenger::RES_TrackFitMessenger(RES_TrackFitter* fitter)
{
  m_fitter = fitter;
  
  m_directory = new G4UIdirectory("/RES/Fit/");
  m_directory->SetGuidance("Commands for the track fit routine");

  m_setVerboseCmd = new G4UIcmdWithAnInteger("/RES/Fit/Verbose", this);
  m_setVerboseCmd->SetGuidance("Set verbosity of fit");
  m_setVerboseCmd->SetParameterName("verbose", false);
  m_setVerboseCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RES_TrackFitMessenger::~RES_TrackFitMessenger()
{
  delete m_directory;
  delete m_setVerboseCmd;
}

void RES_TrackFitMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_setVerboseCmd) {
    m_fitter->SetVerbose(m_setVerboseCmd->GetNewIntValue(newValue));
  }
}
