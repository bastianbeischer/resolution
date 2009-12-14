// $Id: RES_RunMessenger.hh,v 1.4 2009/12/14 08:52:51 beischer Exp $

#ifndef RES_RunMessenger_hh
#define RES_RunMessenger_hh

#include "G4UImessenger.hh"

class RES_RunManager;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;

class RES_RunMessenger : public G4UImessenger
{

public:
  RES_RunMessenger(RES_RunManager* manager);
  ~RES_RunMessenger();

public:
  void SetNewValue(G4UIcommand* command, G4String newValue);

private:
  RES_RunManager*          m_manager;

  G4UIdirectory*           m_directory;
  G4UIcmdWithABool*        m_setStoreResultsCmd;
  G4UIcmdWithAnInteger*    m_generateCmd;
  G4UIcmdWithoutParameter* m_reconstructCmd;
  G4UIcmdWithAnInteger*    m_reconstructWithoutLayerCmd;
  G4UIcmdWithAString*      m_scanChi2FuncCmd;

};

#endif /* RES_RunMessenger_hh */
