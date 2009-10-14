// $Id: RES_TrackFitMessenger.hh,v 1.3 2009/10/14 09:24:29 beischer Exp $

#ifndef RES_TrackFitMessenger_hh
#define RES_TrackFitMessenger_hh

#include "G4UImessenger.hh"

class RES_TrackFitter;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

class RES_TrackFitMessenger : public G4UImessenger
{

public:
  RES_TrackFitMessenger(RES_TrackFitter* fitter);
  ~RES_TrackFitMessenger();

public:
  void SetNewValue(G4UIcommand* command, G4String newValue);

private:
  RES_TrackFitter*      m_fitter;

  G4UIdirectory*        m_directory;
  G4UIcmdWithAnInteger* m_setVerboseCmd;
  G4UIcmdWithAString*   m_setFitMethodCmd;

};

#endif /* RES_TrackFitMessenger_hh */
