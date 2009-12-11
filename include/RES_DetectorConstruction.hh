// $Id: RES_DetectorConstruction.hh,v 1.15 2009/12/11 12:52:22 beischer Exp $

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

enum ModuleType {fiber, silicon};

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
  inline G4double GetModuleUpperSigmaU(unsigned int i) {assert(i < m_moduleUpperSigmaU.size()); return m_moduleUpperSigmaU.at(i);}
  inline G4double GetModuleUpperSigmaV(unsigned int i) {assert(i < m_moduleUpperSigmaV.size()); return m_moduleUpperSigmaV.at(i);}
  inline G4double GetModuleUpperSigmaZ(unsigned int i) {assert(i < m_moduleUpperSigmaZ.size()); return m_moduleUpperSigmaZ.at(i);}
  inline G4double GetModuleLowerSigmaU(unsigned int i) {assert(i < m_moduleLowerSigmaU.size()); return m_moduleLowerSigmaU.at(i);}
  inline G4double GetModuleLowerSigmaV(unsigned int i) {assert(i < m_moduleLowerSigmaV.size()); return m_moduleLowerSigmaV.at(i);}
  inline G4double GetModuleLowerSigmaZ(unsigned int i) {assert(i < m_moduleLowerSigmaZ.size()); return m_moduleLowerSigmaZ.at(i);}

  inline void AddModulePlacement(G4ThreeVector where) {
    m_modulePlacements.push_back(where);
    m_moduleAngles.push_back(0);
    m_moduleInternalAngles.push_back(0);
    m_moduleLength.push_back(m_moduleDefaultLengthFiber);
    m_moduleWidth.push_back(m_moduleDefaultWidthFiber);
    m_moduleType.push_back(fiber);
    m_moduleUpperSigmaU.push_back(m_moduleDefaultSigmaUFiber);
    m_moduleUpperSigmaV.push_back(m_moduleDefaultSigmaVFiber);
    m_moduleUpperSigmaZ.push_back(m_moduleDefaultSigmaZFiber);
    m_moduleLowerSigmaU.push_back(m_moduleDefaultSigmaUFiber);
    m_moduleLowerSigmaV.push_back(m_moduleDefaultSigmaVFiber);
    m_moduleLowerSigmaZ.push_back(m_moduleDefaultSigmaZFiber);
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
  inline void SetModuleUpperSigmaU(G4int module, G4double sigmaU) {
    m_moduleUpperSigmaU[module] = sigmaU;
  }
  inline void SetModuleUpperSigmaV(G4int module, G4double sigmaV) {
    m_moduleUpperSigmaV[module] = sigmaV;
  }
  inline void SetModuleUpperSigmaZ(G4int module, G4double sigmaZ) {
    m_moduleUpperSigmaZ[module] = sigmaZ;
  }
  inline void SetModuleLowerSigmaU(G4int module, G4double sigmaU) {
    m_moduleLowerSigmaU[module] = sigmaU;
  }
  inline void SetModuleLowerSigmaV(G4int module, G4double sigmaV) {
    m_moduleLowerSigmaV[module] = sigmaV;
  }
  inline void SetModuleLowerSigmaZ(G4int module, G4double sigmaZ) {
    m_moduleLowerSigmaZ[module] = sigmaZ;
  }
  inline void SetModuleType(G4int iModule, ModuleType type) {
    m_moduleType[iModule] = type;
    if (type == silicon) {
      m_moduleInternalAngles[iModule] = M_PI/2.;
      m_moduleLength[iModule] = m_moduleDefaultLengthSilicon;
      m_moduleWidth[iModule] = m_moduleDefaultWidthSilicon;
      m_moduleUpperSigmaU[iModule] = m_moduleDefaultUpperSigmaUSilicon;
      m_moduleUpperSigmaV[iModule] = m_moduleDefaultUpperSigmaVSilicon;
      m_moduleUpperSigmaZ[iModule] = m_moduleDefaultUpperSigmaZSilicon;
      m_moduleLowerSigmaU[iModule] = m_moduleDefaultLowerSigmaUSilicon;
      m_moduleLowerSigmaV[iModule] = m_moduleDefaultLowerSigmaVSilicon;
      m_moduleLowerSigmaZ[iModule] = m_moduleDefaultLowerSigmaZSilicon;
    }
    else {
      m_moduleInternalAngles[iModule] = 0.;
      m_moduleLength[iModule] = m_moduleDefaultLengthFiber;
      m_moduleWidth[iModule] = m_moduleDefaultWidthFiber;
      m_moduleUpperSigmaU[iModule] = m_moduleDefaultSigmaUFiber;
      m_moduleUpperSigmaV[iModule] = m_moduleDefaultSigmaVFiber;
      m_moduleUpperSigmaZ[iModule] = m_moduleDefaultSigmaZFiber;
      m_moduleLowerSigmaU[iModule] = m_moduleDefaultSigmaUFiber;
      m_moduleLowerSigmaV[iModule] = m_moduleDefaultSigmaVFiber;
      m_moduleLowerSigmaZ[iModule] = m_moduleDefaultSigmaZFiber;
    }
  }
  inline void SetModuleFiberThickness(G4double fiberThickness) {
    m_moduleFiberThickness = fiberThickness;
  }
  inline void SetModuleGapFiber(G4double gap) {
    m_moduleGapFiber = gap;
  }
  inline void SetModuleKaptonThickness(G4double kaptonThickness) {
    m_moduleKaptonThickness = kaptonThickness;
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
  G4double                        m_moduleDefaultLengthFiber;
  G4double                        m_moduleDefaultWidthFiber;
  G4double                        m_moduleDefaultSigmaUFiber;
  G4double                        m_moduleDefaultSigmaVFiber;
  G4double                        m_moduleDefaultSigmaZFiber;
  G4double                        m_moduleFiberThickness;
  G4double                        m_modulePlasticThickness;
  G4double                        m_moduleFoamThickness;
  G4double                        m_moduleGapFiber;
  G4double                        m_moduleHeightFiber;

  G4Material*                     m_moduleSiliconMaterial;
  G4Material*                     m_moduleKaptonMaterial;
  G4double                        m_moduleDefaultLengthSilicon;
  G4double                        m_moduleDefaultWidthSilicon;
  G4double                        m_moduleSiliconThickness;
  G4double                        m_moduleKaptonThickness;
  G4double                        m_moduleGapSilicon;
  G4double                        m_moduleHeightSilicon;
  G4double                        m_moduleDefaultUpperSigmaUSilicon;
  G4double                        m_moduleDefaultUpperSigmaVSilicon;
  G4double                        m_moduleDefaultUpperSigmaZSilicon;
  G4double                        m_moduleDefaultLowerSigmaUSilicon;
  G4double                        m_moduleDefaultLowerSigmaVSilicon;
  G4double                        m_moduleDefaultLowerSigmaZSilicon;

  std::vector<G4ThreeVector>      m_modulePlacements;
  std::vector<ModuleType>         m_moduleType;
  std::vector<G4double>           m_moduleAngles;
  std::vector<G4double>           m_moduleInternalAngles;
  std::vector<G4double>           m_moduleLength;
  std::vector<G4double>           m_moduleWidth;
  std::vector<G4double>           m_moduleUpperSigmaU;
  std::vector<G4double>           m_moduleUpperSigmaV;
  std::vector<G4double>           m_moduleUpperSigmaZ;
  std::vector<G4double>           m_moduleLowerSigmaU;
  std::vector<G4double>           m_moduleLowerSigmaV;
  std::vector<G4double>           m_moduleLowerSigmaZ;

  G4VPhysicalVolume*              m_world;
  std::vector<G4VPhysicalVolume*> m_modules;
  std::vector<G4VPhysicalVolume*> m_moduleUpperFiber;
  std::vector<G4VPhysicalVolume*> m_moduleLowerFiber;
  std::vector<G4VPhysicalVolume*> m_moduleUpperFoam;
  std::vector<G4VPhysicalVolume*> m_moduleLowerFoam;
  std::vector<G4VPhysicalVolume*> m_moduleUpperPlastic;
  std::vector<G4VPhysicalVolume*> m_moduleLowerPlastic;
  std::vector<G4VPhysicalVolume*> m_moduleUpperSilicon;
  std::vector<G4VPhysicalVolume*> m_moduleLowerSilicon;
  std::vector<G4VPhysicalVolume*> m_moduleUpperKapton;
  std::vector<G4VPhysicalVolume*> m_moduleLowerKapton;
};

#endif /* RES_DetectorConstruction_hh */
