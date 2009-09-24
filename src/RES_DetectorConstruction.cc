#include "RES_DetectorConstruction.hh"

#include "RES_DetectorMessenger.hh"
#include "RES_FiberSD.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4RotationMatrix.hh"
#include "globals.hh"

RES_DetectorConstruction::RES_DetectorConstruction() :
  m_world_solid(0),
  m_world_log(0),
  m_world_phys(0)
{
  // create the messenger
  m_messenger = new RES_DetectorMessenger(this);
  
  // define world material and dimensions
  m_world_material = G4NistManager::Instance()->FindOrBuildMaterial( "G4_AIR" );
  m_world_x = 1.0*m;
  m_world_y = 1.0*m;
  m_world_z = 1.0*m;

  // define default module parameters
  m_module_material = m_world_material;
  m_moduleDefaultWidth = 6.912 * cm;
  m_moduleGap = 0.8 * cm;
  m_moduleLength = 30. * cm;
  m_moduleFiberThickness = 0.1 * cm;
}

RES_DetectorConstruction::~RES_DetectorConstruction()
{
  delete m_messenger;
}

G4VPhysicalVolume* RES_DetectorConstruction::Construct()
{
  ComputeParameters();

  // --------------------------------------------
  // Volumes and Placements
  // --------------------------------------------

  // world
  m_world_solid = new G4Box("world", 0.5*m_world_x, 0.5*m_world_y, 0.5*m_world_z);
  m_world_log   = new G4LogicalVolume(m_world_solid, m_world_material, "world", 0, 0, 0);
  m_world_phys  = new G4PVPlacement(0, G4ThreeVector(), m_world_log, "world", 0, false, 0);

  for (unsigned int i = 0; i < m_modulePlacements.size(); i++) {
    // detector modules
    G4Box* currentModuleSolid = new G4Box("module", 0.5*m_moduleLength, 0.5*m_moduleWidth[i], 0.5*m_moduleHeight);
    G4LogicalVolume* currentModuleLogic = new G4LogicalVolume(currentModuleSolid, m_module_material, "module", 0, 0, 0); // CHANGE MATERIAL HERE
    G4RotationMatrix* currentModuleRotation = new G4RotationMatrix(m_moduleAngles[i], 0., 0.);
    G4PVPlacement* currentModulePlacement = new G4PVPlacement(currentModuleRotation, m_modulePlacements[i], currentModuleLogic, "module", m_world_log, false, i);

    m_module_solid.push_back(currentModuleSolid);
    m_module_log.push_back(currentModuleLogic);
    m_module_phys.push_back(currentModulePlacement);

    // interior of modules
    
    G4Box* currentFiberSolid = new G4Box("module_fiber", 0.5*m_moduleLength, 0.5*m_moduleWidth[i], 0.5*m_moduleFiberThickness);
    G4LogicalVolume* currentFiberLogic = new G4LogicalVolume(currentFiberSolid, m_world_material, "module_fiber", 0, 0, 0);
    G4PVPlacement* currentUpperFiberPlacement = new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*m_moduleGap + 0.5*m_moduleFiberThickness), currentFiberLogic, "module_fiber", currentModuleLogic, false, 0);
    G4PVPlacement* currentLowerFiberPlacement = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*m_moduleGap - 0.5*m_moduleFiberThickness), currentFiberLogic, "module_fiber", currentModuleLogic, false, 1);
    m_module_fiber_solid.push_back(currentFiberSolid);
    m_module_fiber_log.push_back(currentFiberLogic);
    m_module_fiber_upper_phys.push_back(currentUpperFiberPlacement);
    m_module_fiber_lower_phys.push_back(currentLowerFiberPlacement);
  
    G4Box* currentBulkSolid = new G4Box("module_bulk", 0.5*m_moduleLength, 0.5*m_moduleWidth[i], 0.5*m_moduleGap);
    G4LogicalVolume* currentBulkLogic = new G4LogicalVolume(currentBulkSolid, m_world_material, "module_bulk", 0, 0, 0);
    G4PVPlacement* currentBulkPlacement = new G4PVPlacement(0, G4ThreeVector(), currentBulkLogic, "module_bulk", currentModuleLogic, false, 0);
    m_module_bulk_solid.push_back(currentBulkSolid);
    m_module_bulk_log.push_back(currentBulkLogic); 
    m_module_bulk_phys.push_back(currentBulkPlacement); 
  }


  // --------------------------------------------
  // Sensitive detectors
  // --------------------------------------------
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String fiberSDName = "fiberSD";
  RES_FiberSD* aFiberSD = new RES_FiberSD( fiberSDName );
  SDman->AddNewDetector( aFiberSD );
  for (unsigned int i = 0; i < m_module_fiber_log.size(); i++)
    m_module_fiber_log[i]->SetSensitiveDetector(aFiberSD);

  // --------------------------------------------
  // Visualization attributes
  // --------------------------------------------
  SetVisibility();

  // return world
  return m_world_phys;
}

void RES_DetectorConstruction::ComputeParameters()
{
  m_moduleHeight = 2. * m_moduleFiberThickness + m_moduleGap;
}

void RES_DetectorConstruction::SetVisibility()
{
  m_world_log->SetVisAttributes(G4VisAttributes::Invisible);
  
  for (unsigned int i = 0; i < m_module_log.size(); i++)
    m_module_log[i]->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red
  for (unsigned int i = 0; i < m_module_fiber_log.size(); i++)
    m_module_fiber_log[i]->SetVisAttributes(vis_att);

  vis_att = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // blue
  for (unsigned int i = 0; i < m_module_bulk_log.size(); i++)
    m_module_bulk_log[i]->SetVisAttributes(vis_att);
}

G4bool RES_DetectorConstruction::TrackInAcceptance(G4ThreeVector position, G4ThreeVector direction)
{
  G4bool retVal = true;

  for (unsigned int i = 0; i < m_modulePlacements.size(); i++) {
    G4double dz = position.z() - m_modulePlacements[i].z();
    G4double l = dz / direction.z();
    G4ThreeVector currentPosition = position + l*direction;
    if (fabs(currentPosition.x() - m_modulePlacements[i].x()) > 0.5*m_moduleLength)   {retVal = false;}
    if (fabs(currentPosition.y() - m_modulePlacements[i].y()) > 0.5*m_moduleWidth[i]) {retVal = false;}
  }

  return retVal;
}
