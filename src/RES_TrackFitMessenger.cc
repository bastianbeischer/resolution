// $Id: RES_TrackFitMessenger.cc,v 1.6 2010/02/03 15:16:23 beischer Exp $

#include "RES_TrackFitMessenger.hh"

#include "RES_TrackFitter.hh"

#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

RES_TrackFitMessenger::RES_TrackFitMessenger(RES_TrackFitter* fitter)
{
  m_fitter = fitter;
  
  m_directory = new G4UIdirectory("/RES/Fit/");
  m_directory->SetGuidance("Commands for the track fit routine");

  m_setVerboseCmd = new G4UIcmdWithAnInteger("/RES/Fit/Verbose", this);
  m_setVerboseCmd->SetGuidance("Set verbosity of fit");
  m_setVerboseCmd->SetParameterName("verbose", false);
  m_setVerboseCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_setFitMethodCmd = new G4UIcmdWithAString("/RES/Fit/Method", this);
  m_setFitMethodCmd->SetGuidance("Set fit method for minimazation (blobel/minuit)");
  m_setFitMethodCmd->SetParameterName("method", false);
  m_setFitMethodCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RES_TrackFitMessenger::~RES_TrackFitMessenger()
{
  delete m_directory;
  delete m_setVerboseCmd;
  delete m_setFitMethodCmd;
}

void RES_TrackFitMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_setVerboseCmd) {
    m_fitter->SetVerbose(m_setVerboseCmd->GetNewIntValue(newValue));
  }
  if (command == m_setFitMethodCmd) {
    FitMethod method;
    newValue.toLower();
    if (newValue == "blobel")
      method = blobel;
    if (newValue == "minuit")
      method = minuit;
    if (newValue == "oneline")
      method = oneline;
    if (newValue == "twolines")
      method = twolines;
    if (newValue == "transverse")
      method = transverse;
    if (newValue == "fullmatrix")
      method = fullmatrix;
    if (newValue == "testbeam")
      method = testbeam;

    m_fitter->SetFitMethod(method);
  }
}
