// $Id: RES_FieldMessenger.hh,v 1.6 2010/07/14 12:17:35 beischer Exp $

#ifndef RES_FieldMessenger_hh
#define RES_FieldMessenger_hh

#include "G4UImessenger.hh"

class RES_FieldManager;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWith3VectorAndUnit;
class G4UIcommand;

class RES_FieldMessenger : public G4UImessenger
{

public:
  explicit RES_FieldMessenger(RES_FieldManager* manager);
  ~RES_FieldMessenger();

public:
  void SetNewValue(G4UIcommand* command, G4String newValue);

private:
  RES_FieldManager*          m_manager;
  
  G4UIdirectory*             m_directory;
  G4UIcmdWithAString*        m_setInhomFieldFromFileCmd;
  G4UIcmdWithAString*        m_setAMS02FieldFromFileCmd;
  G4UIcmdWith3VectorAndUnit* m_setUniformFieldCmd;
  G4UIcmdWith3VectorAndUnit* m_setDisplacementCmd;

};

#endif /* RES_FieldMessenger_hh */
