// $Id: RES_DetectorConstruction.hh,v 1.11 2009/10/14 16:51:35 beischer Exp $

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

  inline G4double GetWorldX() {return m_worldX;}
  inline G4double GetWorldY() {return m_worldY;}
  inline G4double GetWorldZ() {return m_worldZ;}

  inline G4double GetModuleLength(unsigned int i) {assert(i < m_moduleLength.size()); return m_moduleLength[i];}
  inline G4double GetModuleAngle(unsigned int i) {assert(i < m_moduleAngles.size()); return m_moduleAngles.at(i);}
  inline G4double GetModuleInternalAngle(unsigned int i) {assert(i < m_moduleInternalAngles.size()); return m_moduleInternalAngles.at(i);}
  inline G4double GetModuleSigmaU(unsigned int i) {assert(i < m_moduleSigmaU.size()); return m_moduleSigmaU.at(i);}
  inline G4double GetModuleSigmaV(unsigned int i) {assert(i < m_moduleSigmaV.size()); return m_moduleSigmaV.at(i);}
  inline G4double GetModuleSigmaZ(unsigned int i) {assert(i < m_moduleSigmaZ.size()); return m_moduleSigmaZ.at(i);}

  inline void AddModulePlacement(G4ThreeVector where) {
    m_modulePlacements.push_back(where);
    m_moduleAngles.push_back(0);
    m_moduleInternalAngles.push_back(0);
    m_moduleLength.push_back(m_moduleDefaultLength);
    m_moduleWidth.push_back(m_moduleDefaultWidth);
    m_moduleSigmaU.push_back(m_moduleDefaultSigmaU);
    m_moduleSigmaV.push_back(m_moduleDefaultSigmaV);
    m_moduleSigmaZ.push_back(m_moduleDefaultSigmaZ);
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
  inline void SetModuleLength(G4int module, G4double length) {
    m_moduleLength[module] = length;
  }
  inline void SetModuleFiberThickness(G4double fiberThickness) {
    m_moduleFiberThickness = fiberThickness;
  }
  inline void SetModuleGap(G4double gap) {
    m_moduleGap = gap;
  }

private:
  void ComputeParameters();
  void SetVisibility();

private:
  RES_DetectorMessenger*          m_messenger;

  G4Material*                     m_worldMaterial;
  G4double                        m_worldX;
  G4double                        m_worldY;
  G4double                        m_worldZ;

  G4Material*                     m_moduleMaterial;
  G4Material*                     m_moduleFiberMaterial;
  G4Material*                     m_modulePlasticMaterial;
  G4Material*                     m_moduleFoamMaterial;
  G4double                        m_moduleDefaultLength;
  G4double                        m_moduleDefaultWidth;
  G4double                        m_moduleFiberThickness;
  G4double                        m_modulePlasticThickness;
  G4double                        m_moduleFoamThickness;
  G4double                        m_moduleGap; 
  G4double                        m_moduleHeight;
  G4double                        m_moduleDefaultSigmaU;
  G4double                        m_moduleDefaultSigmaV;
  G4double                        m_moduleDefaultSigmaZ;

  std::vector<G4ThreeVector>      m_modulePlacements;
  std::vector<G4double>           m_moduleAngles;
  std::vector<G4double>           m_moduleInternalAngles;
  std::vector<G4double>           m_moduleLength;
  std::vector<G4double>           m_moduleWidth;
  std::vector<G4double>           m_moduleSigmaU;
  std::vector<G4double>           m_moduleSigmaV;
  std::vector<G4double>           m_moduleSigmaZ;

  G4VPhysicalVolume*              m_world;
  std::vector<G4VPhysicalVolume*> m_modules;
  std::vector<G4VPhysicalVolume*> m_moduleUpperFiber;
  std::vector<G4VPhysicalVolume*> m_moduleLowerFiber;
  std::vector<G4VPhysicalVolume*> m_moduleUpperFoam;
  std::vector<G4VPhysicalVolume*> m_moduleLowerFoam;
  std::vector<G4VPhysicalVolume*> m_moduleUpperPlastic;
  std::vector<G4VPhysicalVolume*> m_moduleLowerPlastic;

};

#endif /* RES_DetectorConstruction_hh */
