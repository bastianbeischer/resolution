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
  m_moduleWidth = 6.912 * cm;
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

  // detector modules
  m_module_solid = new G4Box("module", 0.5*m_moduleLength, 0.5*m_moduleWidth, 0.5*m_moduleHeight);
  m_module_log   = new G4LogicalVolume(m_module_solid, m_module_material, "module", 0, 0, 0); // CHANGE MATERIAL HERE
  for (unsigned int i = 0; i < m_modulePlacements.size(); i++) {
    G4RotationMatrix* rotation = new G4RotationMatrix(m_moduleAngles[i], 0., 0.);
    m_module_phys.push_back(new G4PVPlacement(rotation, m_modulePlacements[i], m_module_log, "module", m_world_log, false, i));
  }

  // interior of modules
  m_module_fiber_solid = new G4Box("module_fiber", 0.5*m_moduleLength, 0.5*m_moduleWidth, 0.5*m_moduleFiberThickness);
  m_module_fiber_log = new G4LogicalVolume(m_module_fiber_solid, m_world_material, "module_fiber", 0, 0, 0);
  m_module_fiber_phys.push_back(new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*m_moduleGap + 0.5*m_moduleFiberThickness), m_module_fiber_log, "module_fiber", m_module_log, false, 0));
  m_module_fiber_phys.push_back(new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*m_moduleGap - 0.5*m_moduleFiberThickness), m_module_fiber_log, "module_fiber", m_module_log, false, 1));
  
  m_module_bulk_solid = new G4Box("module_bulk", 0.5*m_moduleLength, 0.5*m_moduleWidth, 0.5*m_moduleGap);
  m_module_bulk_log = new G4LogicalVolume(m_module_bulk_solid, m_world_material, "module_bulk", 0, 0, 0);
  m_module_bulk_phys = new G4PVPlacement(0, G4ThreeVector(), m_module_bulk_log, "module_bulk", m_module_log, false, 0);

  // --------------------------------------------
  // Sensitive detectors
  // --------------------------------------------
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String fiberSDName = "fiberSD";
  RES_FiberSD* aFiberSD = new RES_FiberSD( fiberSDName );
  SDman->AddNewDetector( aFiberSD );
  m_module_fiber_log->SetSensitiveDetector(aFiberSD);

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
  m_module_log->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red
  m_module_fiber_log->SetVisAttributes(vis_att);

  vis_att = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // blue
  m_module_bulk_log->SetVisAttributes(vis_att);
}

G4bool RES_DetectorConstruction::TrackInAcceptance(G4ThreeVector /*position*/, G4ThreeVector /*direction*/)
{
  return true;
}
