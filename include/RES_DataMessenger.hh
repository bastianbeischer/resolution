// $Id: RES_DataMessenger.hh,v 1.3 2009/10/14 09:24:32 beischer Exp $

#ifndef RES_DataMessenger_hh
#define RES_DataMessenger_hh

#include "G4UImessenger.hh"

class RES_DataHandler;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;

class RES_DataMessenger : public G4UImessenger
{

public:
  RES_DataMessenger(RES_DataHandler* dataHandler);
  ~RES_DataMessenger();

public:
  void SetNewValue(G4UIcommand* command, G4String newValue);

private:
  RES_DataHandler*    m_dataHandler;
  
  G4UIdirectory*      m_directory;

  G4UIcmdWithABool*   m_overwriteFileCmd;
  G4UIcmdWithAString* m_setFileNameCmd;
};

#endif /* RES_DataMessenger_hh */
