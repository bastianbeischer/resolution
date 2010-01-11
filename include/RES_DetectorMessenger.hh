// $Id: RES_DetectorMessenger.hh,v 1.14 2010/01/11 15:25:51 beischer Exp $

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
  G4UIcmdWithAString*        m_setModuleTypeCmd;
  G4UIcmdWithAString*        m_setModuleRotationCmd;
  G4UIcmdWithAString*        m_setModuleInternalRotationCmd;
  G4UIcmdWithAString*        m_setModuleWidthCmd;
  G4UIcmdWithAString*        m_setModuleLengthCmd;
  G4UIcmdWithAString*        m_setModuleUpperSigmaUCmd;
  G4UIcmdWithAString*        m_setModuleUpperSigmaVCmd;
  G4UIcmdWithAString*        m_setModuleUpperSigmaZCmd;
  G4UIcmdWithAString*        m_setModuleLowerSigmaUCmd;
  G4UIcmdWithAString*        m_setModuleLowerSigmaVCmd;
  G4UIcmdWithAString*        m_setModuleLowerSigmaZCmd;
  G4UIcmdWithAString*        m_setModuleUpperEfficiencyCmd;
  G4UIcmdWithAString*        m_setModuleLowerEfficiencyCmd;
};

#endif /* RES_DetectorMessenger_hh */
