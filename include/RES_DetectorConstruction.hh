#ifndef RES_DetectorConstruction_hh
#define RES_DetectorConstruction_hh

#include <vector>

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"

class G4Material;
class RES_DetectorMessenger;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;

class RES_DetectorConstruction : public G4VUserDetectorConstruction
{

public:
  RES_DetectorConstruction();
  ~RES_DetectorConstruction();

  G4VPhysicalVolume* Construct();

  inline G4double GetWorldX() {return m_world_x;}
  inline G4double GetWorldY() {return m_world_y;}
  inline G4double GetWorldZ() {return m_world_z;}

  inline void AddModulePlacement(G4ThreeVector where) {
    m_modulePlacements.push_back(where);
  }

private:
  void ComputeParameters();
  void SetVisibility();

private:
  RES_DetectorMessenger*          m_messenger;

  G4Material*                     m_world_material;
  G4double                        m_world_x;
  G4double                        m_world_y;
  G4double                        m_world_z;

  G4Material*                     m_module_material;
  G4double                        m_moduleStereoAngle;
  G4double                        m_moduleWidth;
  G4double                        m_moduleLength;
  G4double                        m_moduleLayerThickness;
  G4double                        m_moduleGap; 
  G4double                        m_moduleHeight;

  G4int                           m_numberOfModules;
  std::vector<G4ThreeVector>      m_modulePlacements;
  
  G4Box*                          m_world_solid;
  G4LogicalVolume*                m_world_log;
  G4VPhysicalVolume*              m_world_phys;
  
  G4Box*                          m_module_solid;
  G4LogicalVolume*                m_module_log;
  std::vector<G4VPhysicalVolume*> m_module_phys;
  
  G4Box*                          m_module_upper_layer_solid;
  G4LogicalVolume*                m_module_upper_layer_log;
  G4VPhysicalVolume*              m_module_upper_layer_phys;

  G4Box*                          m_module_lower_layer_solid;
  G4LogicalVolume*                m_module_lower_layer_log;
  G4VPhysicalVolume*              m_module_lower_layer_phys;

  G4Box*                          m_module_bulk_solid;
  G4LogicalVolume*                m_module_bulk_log;
  G4VPhysicalVolume*              m_module_bulk_phys;

};

#endif /* RES_DetectorConstruction_hh */
