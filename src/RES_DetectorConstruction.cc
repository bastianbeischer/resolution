// $Id: RES_DetectorConstruction.cc,v 1.24 2010/04/23 01:07:09 beischer Exp $

#include "RES_DetectorConstruction.hh"

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
#include "globals.hh"

RES_DetectorConstruction::RES_DetectorConstruction() :
  m_world(0)
{
  // create the messenger
  m_messenger = new RES_DetectorMessenger(this);
  
  // define materials
  m_worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_AIR" );

  // define world dimensions
  m_worldX = 1.6*m;
  m_worldY = 1.6*m;
  m_worldZ = 2.5*m;
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

// void RES_DetectorConstruction::ComputeParameters()
// {
//   m_moduleGapFiber = 2.*m_modulePlasticThickness + 2.*m_moduleFoamThickness;
//   m_moduleHeightFiber = 2. * m_moduleFiberThickness + m_moduleGapFiber;

//   m_moduleHeightSilicon = 2.*m_moduleGapSilicon + 2.*m_moduleSiliconThickness + 2*m_moduleKaptonThickness;
// }

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
  G4bool retVal = true;

  for (unsigned int i = 0; i < m_modules.size(); i++)
    if (!m_modules.at(i)->CheckIfTrackPassesThrough(position, direction))
      retVal = false;

  G4double dz = - position.z();
  G4double l = dz / direction.z();
  G4ThreeVector currentPosition = position + l*direction;
  if (sqrt(pow(currentPosition.x(),2.) + pow(currentPosition.x(),2.)) > 7.5*cm)   {retVal = false;}

  return retVal;
}
