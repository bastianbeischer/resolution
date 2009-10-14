// $Id: RES_AlignmentMessenger.hh,v 1.2 2009/10/14 17:29:04 beischer Exp $

#ifndef RES_AlignmentMessenger_hh
#define RES_AlignmentMessenger_hh

#include "G4UImessenger.hh"

class RES_AlignmentManager;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;

class RES_AlignmentMessenger : public G4UImessenger
{

public:
  RES_AlignmentMessenger(RES_AlignmentManager* manager);
  ~RES_AlignmentMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String newValue);

private:
  RES_AlignmentManager*    m_manager;

  G4UIdirectory*           m_directory;

  G4UIcmdWithoutParameter* m_startAlignmentCmd;

};

#endif /* RES_AlignmentMessenger_hh */
