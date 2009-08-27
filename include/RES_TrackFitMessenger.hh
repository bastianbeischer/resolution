#ifndef RES_TrackFitMessenger_hh
#define RES_TrackFitMessenger_hh

#include "G4UImessenger.hh"

class RES_TrackFitter;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAnInteger;

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

};

#endif /* RES_TrackFitMessenger_hh */