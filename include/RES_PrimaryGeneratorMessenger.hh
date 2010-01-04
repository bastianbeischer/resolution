// $Id: RES_PrimaryGeneratorMessenger.hh,v 1.4 2010/01/04 15:06:24 beischer Exp $

#ifndef RES_PrimaryGeneratorMessenger_hh
#define RES_PrimaryGeneratorMessenger_hh

#include "G4UImessenger.hh"
#include "G4String.hh"

class RES_PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;

class RES_PrimaryGeneratorMessenger : public G4UImessenger
{

public:
  RES_PrimaryGeneratorMessenger(RES_PrimaryGeneratorAction* generator);
  ~RES_PrimaryGeneratorMessenger();

public:
  void SetNewValue(G4UIcommand* command, G4String newValue);

private:
  RES_PrimaryGeneratorAction* m_generator;

  G4UIdirectory*              m_directory;

  G4UIcmdWithABool*           m_randomOriginCmd;
  G4UIcmdWithABool*           m_randomDirectionCmd;
  G4UIcmdWithADoubleAndUnit*  m_setStartZCmd;
  
};

#endif /* RES_PrimaryGeneratorMessenger_hh */
