// $Id: RES_DetectorConstruction.hh,v 1.19 2010/01/13 15:24:30 beischer Exp $

#ifndef RES_DetectorConstruction_hh
#define RES_DetectorConstruction_hh

#include <vector>
#include <assert.h>

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "RES_Module.hh"

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

  G4double GetWorldX() {return m_worldX;}
  G4double GetWorldY() {return m_worldY;}
  G4double GetWorldZ() {return m_worldZ;}
  unsigned int GetNumberOfModules() {return m_modules.size();}
  RES_Module* GetModule(unsigned int i) {assert(i < m_modules.size()); return m_modules.at(i);}

  void AddModulePlacement(G4ThreeVector where) {m_modules.push_back(new RES_Module(where));}

  G4VPhysicalVolume* Construct();
  G4bool TrackInAcceptance(G4ThreeVector position, G4ThreeVector direction);

private:
  void SortModules();

  void SetSensitiveDetectors();
  void SetVisibility();

private:
  RES_DetectorMessenger*          m_messenger;

  G4Material*                     m_worldMaterial;
  G4double                        m_worldX;
  G4double                        m_worldY;
  G4double                        m_worldZ;

  std::vector<RES_Module*>        m_modules;

  G4VPhysicalVolume*              m_world;
  std::vector<G4PVPlacement*>     m_modulePlacements;

};

#endif /* RES_DetectorConstruction_hh */
