// $Id: RES_DetectorConstruction.cc,v 1.32 2010/07/19 20:20:13 beischer Exp $

#include "RES_DetectorConstruction.hh"

#include "RES_MagneticField.hh"
#include "RES_DetectorMessenger.hh"
#include "RES_Module.hh"
#include "RES_SD.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4RotationMatrix.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "globals.hh"

RES_DetectorConstruction::RES_DetectorConstruction() :
  m_world(0)
{
  // create the messenger
  m_messenger = new RES_DetectorMessenger(this);
  
  // define materials
  m_worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_AIR" );

  // define world dimensions
  m_worldX = 3.0*m;
  m_worldY = 3.0*m;
  m_worldZ = 3.0*m;
}

RES_DetectorConstruction::~RES_DetectorConstruction()
{
  for(std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++)
    delete *it;
  delete m_messenger;
}

G4VPhysicalVolume* RES_DetectorConstruction::Construct()
{
  // --------------------------------------------
  // Volumes and Placements
  // --------------------------------------------
  // world
  G4Box* worldSolid         = new G4Box("world", 0.5*m_worldX, 0.5*m_worldY, 0.5*m_worldZ);
  G4LogicalVolume* worldLog = new G4LogicalVolume(worldSolid, m_worldMaterial, "world", 0, 0, 0);
  m_world                   = new G4PVPlacement(0, G4ThreeVector(), worldLog, "world", 0, false, 0);

  // detector modules
  SortModules();
  int copyNumber = 0;
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    RES_Module* module = *it;
    m_modulePlacements.push_back(module->Construct(m_world, copyNumber));
    copyNumber++;
  }

  // --------------------------------------------
  // Sensitive Detectors
  // --------------------------------------------
  SetSensitiveDetectors();

  // --------------------------------------------
  // Visualization attributes
  // --------------------------------------------
  SetVisibility();

  // return world
  return m_world;
}

void RES_DetectorConstruction::SortModules()
{
  module_less module_comp;
  std::sort(m_modules.begin(), m_modules.end(), module_comp);
}

void RES_DetectorConstruction::SetSensitiveDetectors()
{
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String fiberSDName = "fiberSD";
  RES_SD* aSD = new RES_SD( fiberSDName );
  SDman->AddNewDetector( aSD );

  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    RES_Module* module = *it;
    module->SetSD(aSD);
  }
}

void RES_DetectorConstruction::SetVisibility()
{
  m_world->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);

  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    RES_Module* module = *it;
    module->SetVisibility();
  }  
}  

G4bool RES_DetectorConstruction::TrackInAcceptance(G4ThreeVector position, G4ThreeVector direction)
{
  // count the number of distinct layers in the setup
  G4double currentZ = DBL_MAX;
  G4int nLayers = 0;
  for (unsigned int i = 0; i < m_modules.size(); i++) {
    if (currentZ != m_modules.at(i)->GetPlacement().z()) {
      nLayers++;
      currentZ = m_modules.at(i)->GetPlacement().z();
    }
  }

  // fill the vector with passed modules
  std::vector<RES_Module*> modulesPassed;
  for (unsigned int i = 0; i < m_modules.size(); i++)
    if (m_modules.at(i)->CheckIfTrackPassesThrough(position, direction))
      modulesPassed.push_back(m_modules.at(i));

  // sort passed modules by placement.z()
  module_less module_comp;
  std::sort(modulesPassed.begin(), modulesPassed.end(), module_comp);

  // count the number of distinct layers in passed modules
  G4int nPassed = 0;
  currentZ = DBL_MAX;
  for (unsigned int i = 0; i < modulesPassed.size(); i++) {
    if (currentZ != modulesPassed.at(i)->GetPlacement().z()) {
      nPassed++;
      currentZ = modulesPassed.at(i)->GetPlacement().z();
    }
  }

  // we require passage through all layers
  if (nLayers != nPassed)
    return false;

  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  RES_MagneticField* field = (RES_MagneticField*) fieldMgr->GetDetectorField();
  if (field && !field->CheckIfTrackIsInsideMagnet(position, direction))
    return false;

  return true;
}

void RES_DetectorConstruction::PrintMaterials()
{
  m_modules.at(0)->PrintMaterials();
}
