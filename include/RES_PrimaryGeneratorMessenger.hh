// $Id: RES_PrimaryGeneratorMessenger.hh,v 1.3 2009/10/14 09:24:30 beischer Exp $

#ifndef RES_PrimaryGeneratorMessenger_hh
#define RES_PrimaryGeneratorMessenger_hh

#include "G4UImessenger.hh"
#include "G4String.hh"

class RES_PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithABool;

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
  
};

#endif /* RES_PrimaryGeneratorMessenger_hh */
