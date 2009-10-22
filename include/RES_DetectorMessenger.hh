// $Id: RES_DetectorMessenger.hh,v 1.9 2009/10/22 14:45:15 beischer Exp $

#ifndef RES_DetectorMessenger_hh
#define RES_DetectorMessenger_hh

#include "globals.hh"

#include "G4UImessenger.hh"

class RES_DetectorConstruction;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class RES_DetectorMessenger : public G4UImessenger
{

public:
  explicit RES_DetectorMessenger(RES_DetectorConstruction* detector);
  ~RES_DetectorMessenger();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:
  RES_DetectorConstruction*  m_detector;
  G4UIdirectory*             m_topDirectory;
  G4UIdirectory*             m_detDirectory;
  
  G4UIcmdWith3VectorAndUnit* m_addModulePlacementCmd;
  G4UIcmdWithAString*        m_setModuleRotationCmd;
  G4UIcmdWithAString*        m_setModuleInternalRotationCmd;
  G4UIcmdWithAString*        m_setModuleWidthCmd;
  G4UIcmdWithAString*        m_setModuleLengthCmd;
  G4UIcmdWithADoubleAndUnit* m_setModuleFiberThicknessCmd;
  G4UIcmdWithADoubleAndUnit* m_setModuleGapCmd;
  G4UIcmdWithAString*        m_setModuleSigmaUCmd;
  G4UIcmdWithAString*        m_setModuleSigmaVCmd;
  G4UIcmdWithAString*        m_setModuleSigmaZCmd;

};

#endif /* RES_DetectorMessenger_hh */
