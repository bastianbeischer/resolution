#ifndef RES_DetectorConstruction_hh
#define RES_DetectorConstruction_hh

#include <vector>
#include <assert.h>

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

  G4bool TrackInAcceptance(G4ThreeVector position, G4ThreeVector direction);

  inline G4double GetWorldX() {return m_world_x;}
  inline G4double GetWorldY() {return m_world_y;}
  inline G4double GetWorldZ() {return m_world_z;}

  inline G4double GetModuleLength()       {return m_moduleLength;}
  inline G4double GetModuleAngle(unsigned int i) {assert(i < m_moduleAngles.size()); return m_moduleAngles.at(i);}
  inline G4double GetModuleInternalAngle(unsigned int i) {assert(i < m_moduleInternalAngles.size()); return m_moduleInternalAngles.at(i);}

  inline void AddModulePlacement(G4ThreeVector where) {
    m_modulePlacements.push_back(where);
    m_moduleAngles.push_back(0);
    m_moduleInternalAngles.push_back(0);
    m_moduleWidth.push_back(m_moduleDefaultWidth);
  }
  inline void SetModuleAngle(G4int module, G4double angle) {
    m_moduleAngles[module] = angle;
  }
  inline void SetModuleInternalAngle(G4int module, G4double angle) {
    m_moduleInternalAngles[module] = angle;
  }
  inline void SetModuleWidth(G4int module, G4double width) {
    m_moduleWidth[module] = width;
  }
  inline void SetModuleLength(G4double length) {m_moduleLength = length;}
  inline void SetModuleFiberThickness(G4double fiberThickness) {m_moduleFiberThickness = fiberThickness;}
  inline void SetModuleGap(G4double gap) {m_moduleGap = gap;}

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
  G4Material*                     m_module_fiber_material;
  G4Material*                     m_module_bulk_material;
  G4double                        m_moduleDefaultWidth;
  G4double                        m_moduleLength;
  G4double                        m_moduleFiberThickness;
  G4double                        m_moduleGap; 
  G4double                        m_moduleHeight;

  std::vector<G4ThreeVector>      m_modulePlacements;
  std::vector<G4double>           m_moduleAngles;
  std::vector<G4double>           m_moduleInternalAngles;
  std::vector<G4double>           m_moduleWidth;

  G4Box*                          m_world_solid;
  G4LogicalVolume*                m_world_log;
  G4VPhysicalVolume*              m_world_phys;
  
  std::vector<G4Box*>             m_module_solid;
  std::vector<G4LogicalVolume*>   m_module_log;
  std::vector<G4VPhysicalVolume*> m_module_phys;
  
  std::vector<G4Box*>             m_module_fiber_solid;
  std::vector<G4LogicalVolume*>   m_module_fiber_log;
  std::vector<G4VPhysicalVolume*> m_module_fiber_upper_phys;
  std::vector<G4VPhysicalVolume*> m_module_fiber_lower_phys;

  std::vector<G4Box*>             m_module_bulk_solid;
  std::vector<G4LogicalVolume*>   m_module_bulk_log;
  std::vector<G4VPhysicalVolume*> m_module_bulk_phys;

};

#endif /* RES_DetectorConstruction_hh */
