// $Id: RES_DetectorConstruction.cc,v 1.22 2010/01/11 09:59:58 beischer Exp $

#include "RES_DetectorConstruction.hh"

#include "RES_DetectorMessenger.hh"
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
  ComputeParameters();

  // --------------------------------------------
  // Sensitive detectors
  // --------------------------------------------
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String fiberSDName = "fiberSD";
  RES_SD* aSD = new RES_SD( fiberSDName );
  SDman->AddNewDetector( aSD );

  // --------------------------------------------
  // Volumes and Placements (also add SDs to volumes)
  // --------------------------------------------

  // world
  G4Box* worldSolid         = new G4Box("world", 0.5*m_worldX, 0.5*m_worldY, 0.5*m_worldZ);
  G4LogicalVolume* worldLog = new G4LogicalVolume(worldSolid, m_worldMaterial, "world", 0, 0, 0);
  m_world                   = new G4PVPlacement(0, G4ThreeVector(), worldLog, "world", 0, false, 0);

  // detector modules
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    RES_Module* module = *it;

    m_modulesPhys.push_back(module->Construct(m_world));

    
  }

  // --------------------------------------------
  // Visualization attributes
  // --------------------------------------------
  SetVisibility();

  // return world
  return m_world;
}

void RES_DetectorConstruction::ComputeParameters()
{
  m_moduleGapFiber = 2.*m_modulePlasticThickness + 2.*m_moduleFoamThickness;
  m_moduleHeightFiber = 2. * m_moduleFiberThickness + m_moduleGapFiber;

  m_moduleHeightSilicon = 2.*m_moduleGapSilicon + 2.*m_moduleSiliconThickness + 2*m_moduleKaptonThickness;
}

void RES_DetectorConstruction::SetVisibility()
{
  m_world->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);
  
  for (unsigned int i = 0; i < m_modules.size(); i++)
    m_modules[i]->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red
  for (unsigned int i = 0; i < m_moduleUpperFiber.size(); i++)
    m_moduleUpperFiber[i]->GetLogicalVolume()->SetVisAttributes(visAtt);

  visAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // green
  for (unsigned int i = 0; i < m_moduleUpperFoam.size(); i++)
    m_moduleUpperFoam[i]->GetLogicalVolume()->SetVisAttributes(visAtt);

  visAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // green
  for (unsigned int i = 0; i < m_moduleUpperPlastic.size(); i++)
    m_moduleUpperPlastic[i]->GetLogicalVolume()->SetVisAttributes(visAtt);

  visAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // blue
  for (unsigned int i = 0; i < m_moduleUpperSilicon.size(); i++)
    m_moduleUpperSilicon[i]->GetLogicalVolume()->SetVisAttributes(visAtt);

  visAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // green
  for (unsigned int i = 0; i < m_moduleUpperKapton.size(); i++)
    m_moduleUpperKapton[i]->GetLogicalVolume()->SetVisAttributes(visAtt);
}

G4bool RES_DetectorConstruction::TrackInAcceptance(G4ThreeVector position, G4ThreeVector direction)
{
  G4bool retVal = true;

  for (unsigned int i = 0; i < m_modulePlacements.size(); i++) {
    G4double dz = m_modulePlacements[i].z() - position.z();
    G4double l = dz / direction.z();
    G4ThreeVector currentPosition = position + l*direction;
    if (fabs(currentPosition.x() - m_modulePlacements[i].x()) > 0.5*m_moduleLength[i])   {retVal = false;}
    if (fabs(currentPosition.y() - m_modulePlacements[i].y()) > 0.5*m_moduleWidth[i]) {retVal = false;}
  }

  return retVal;
}
