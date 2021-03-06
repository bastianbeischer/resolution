// $Id: RES_DetectorMessenger.hh,v 1.17 2010/07/02 11:13:31 beischer Exp $

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
class G4UIcmdWithoutParameter;

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
  
  G4UIcmdWithADoubleAndUnit* m_setWorldXCmd;
  G4UIcmdWithADoubleAndUnit* m_setWorldYCmd;
  G4UIcmdWithADoubleAndUnit* m_setWorldZCmd;
  G4UIcmdWith3VectorAndUnit* m_addModulePlacementCmd;
  G4UIcmdWithoutParameter*   m_printMaterialsCmd;
  G4UIcmdWithAString*        m_setModuleSubtractHolesCmd;
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
  G4UIcmdWithAString*        m_setLayerSubtractHolesCmd;
  G4UIcmdWithAString*        m_setLayerTypeCmd;
  G4UIcmdWithAString*        m_setLayerRotationCmd;
  G4UIcmdWithAString*        m_setLayerInternalRotationCmd;
  G4UIcmdWithAString*        m_setLayerWidthCmd;
  G4UIcmdWithAString*        m_setLayerLengthCmd;
  G4UIcmdWithAString*        m_setLayerUpperSigmaUCmd;
  G4UIcmdWithAString*        m_setLayerUpperSigmaVCmd;
  G4UIcmdWithAString*        m_setLayerUpperSigmaZCmd;
  G4UIcmdWithAString*        m_setLayerLowerSigmaUCmd;
  G4UIcmdWithAString*        m_setLayerLowerSigmaVCmd;
  G4UIcmdWithAString*        m_setLayerLowerSigmaZCmd;
  G4UIcmdWithAString*        m_setLayerUpperEfficiencyCmd;
  G4UIcmdWithAString*        m_setLayerLowerEfficiencyCmd;
  G4UIcmdWithAString*        m_addLayerCmd;
};

#endif /* RES_DetectorMessenger_hh */
