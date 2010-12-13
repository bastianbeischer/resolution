// $Id: RES_DetectorConstruction.hh,v 1.21.2.1 2010/09/08 17:44:59 beischer Exp $

#ifndef RES_DetectorConstruction_hh
#define RES_DetectorConstruction_hh

#include <vector>
#include <assert.h>

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "RES_Layer.hh"
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
  unsigned int GetNumberOfLayers();
  unsigned int GetNumberOfLayers(std::vector<RES_Module*>& modules);
  RES_Layer* GetLayer(unsigned int i) {assert(i < m_layers.size()); return m_layers.at(i);}
  RES_Module* GetModule(unsigned int i) {assert(i < m_modules.size()); return m_modules.at(i);}

  void SetWorldX(G4double x) {m_worldX = x;}
  void SetWorldY(G4double y) {m_worldY = y;}
  void SetWorldZ(G4double z) {m_worldZ = z;}
  void AddLayer(unsigned int nModules, G4double z) {m_layers.push_back(new RES_Layer(nModules, z));}
  void AddModulePlacement(G4ThreeVector where) {m_modules.push_back(new RES_Module(where));}

  G4VPhysicalVolume* Construct();
  G4bool TrackInAcceptance(G4ThreeVector position, G4ThreeVector direction);

  void PrintMaterials();

private:
  void SortModules();
  void SortModules(std::vector<RES_Module*>& modules);

  void SetSensitiveDetectors();
  void SetVisibility();

private:
  RES_DetectorMessenger*          m_messenger;

  G4Material*                     m_worldMaterial;
  G4double                        m_worldX;
  G4double                        m_worldY;
  G4double                        m_worldZ;

  std::vector<RES_Layer*>         m_layers;
  std::vector<RES_Module*>        m_modules;

  G4VPhysicalVolume*              m_world;
  std::vector<G4PVPlacement*>     m_modulePlacements;

};

#endif /* RES_DetectorConstruction_hh */
