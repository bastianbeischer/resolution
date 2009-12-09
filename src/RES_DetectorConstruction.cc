// $Id: RES_DetectorConstruction.cc,v 1.17 2009/12/09 21:43:31 beischer Exp $

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
  m_world(0)
{
  // create the messenger
  m_messenger = new RES_DetectorMessenger(this);
  
  G4String symbol;
  G4double density,z,a;
  G4int components,natoms;

  // Elements
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);

  // Rohacell
  G4Material* rohacell = new G4Material( "rohacell", density = 32.*kg/m3, components = 4 );
  rohacell->AddElement(H, natoms=5); // ???
  rohacell->AddElement(C, natoms=3); // ???
  rohacell->AddElement(N, natoms=1); // ???
  rohacell->AddElement(O, natoms=1); // ???

  // define materials
  m_worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_AIR" );
  m_moduleMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_AIR" );
  m_modulePlasticMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_POLYCARBONATE" );
  m_moduleFoamMaterial = rohacell;
  m_moduleFiberMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_POLYSTYRENE" );


  // define world dimensions
  m_worldX = 1.0*m;
  m_worldY = 1.0*m;
  m_worldZ = 1.0*m;


  // define default module parameters
  m_moduleDefaultWidth = 6.912 * cm;
  m_moduleDefaultLength = 30. * cm;
  m_moduleFoamThickness = 0.3 * cm;
  m_modulePlasticThickness = 0.01 * cm;
  m_moduleFiberThickness = 0.1 * cm;
  m_moduleDefaultSigmaU = m_moduleDefaultLength/sqrt(12);
  m_moduleDefaultSigmaV = 50. * um;
  m_moduleDefaultSigmaZ = 0. * cm;
}

RES_DetectorConstruction::~RES_DetectorConstruction()
{
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
  RES_FiberSD* aFiberSD = new RES_FiberSD( fiberSDName );
  SDman->AddNewDetector( aFiberSD );

  // --------------------------------------------
  // Volumes and Placements (also add SDs to volumes)
  // --------------------------------------------

  // world
  G4Box* worldSolid         = new G4Box("world", 0.5*m_worldX, 0.5*m_worldY, 0.5*m_worldZ);
  G4LogicalVolume* worldLog = new G4LogicalVolume(worldSolid, m_worldMaterial, "world", 0, 0, 0);
  m_world                   = new G4PVPlacement(0, G4ThreeVector(), worldLog, "world", 0, false, 0);

  for (unsigned int i = 0; i < m_modulePlacements.size(); i++) {
    // detector modules
    G4Box* currentModuleSolid = new G4Box("module", 0.5*m_moduleLength[i], 0.5*m_moduleWidth[i], 0.5*m_moduleHeight);
    G4LogicalVolume* currentModuleLogic = new G4LogicalVolume(currentModuleSolid, m_moduleMaterial, "module", 0, 0, 0); // CHANGE MATERIAL HERE
    G4RotationMatrix* currentModuleRotation = new G4RotationMatrix(m_moduleAngles[i], 0., 0.);
    G4PVPlacement* currentModulePlacement = new G4PVPlacement(currentModuleRotation, m_modulePlacements[i], currentModuleLogic, "module", worldLog, false, i);
    m_modules.push_back(currentModulePlacement);

    // interior of modules
    G4Box* currentFiberSolid = new G4Box("moduleFiber", 0.5*m_moduleLength[i], 0.5*m_moduleWidth[i], 0.5*m_moduleFiberThickness);
    G4LogicalVolume* currentFiberLogic = new G4LogicalVolume(currentFiberSolid, m_moduleFiberMaterial, "moduleFiber", 0, 0, 0);
    G4RotationMatrix* currentInternalRotation = new G4RotationMatrix(m_moduleInternalAngles[i], 0., 0.);
    G4PVPlacement* currentUpperFiberPlacement = new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*m_moduleGap + 0.5*m_moduleFiberThickness),
                                                                  currentFiberLogic, "moduleFiber", currentModuleLogic, false, 0);
    G4PVPlacement* currentLowerFiberPlacement = new G4PVPlacement(currentInternalRotation, G4ThreeVector(0, 0, -0.5*m_moduleGap - 0.5*m_moduleFiberThickness),
                                                                  currentFiberLogic, "moduleFiber", currentModuleLogic, false, 1);
    m_moduleUpperFiber.push_back(currentUpperFiberPlacement);
    m_moduleLowerFiber.push_back(currentLowerFiberPlacement);

    G4Box* currentPlasticSolid = new G4Box("modulePlastic", 0.5*m_moduleLength[i], 0.5*m_moduleWidth[i], 0.5*m_modulePlasticThickness);
    G4LogicalVolume* currentPlasticLogic = new G4LogicalVolume(currentPlasticSolid, m_modulePlasticMaterial, "modulePlastic", 0, 0, 0);
    G4PVPlacement* currentUpperPlasticPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,0.5*m_moduleFoamThickness + 0.5*m_modulePlasticThickness),
                                                                    currentPlasticLogic, "modulePlastic", currentModuleLogic, false, 0);
    G4PVPlacement* currentLowerPlasticPlacement = new G4PVPlacement(currentInternalRotation, G4ThreeVector(0,0,-0.5*m_moduleFoamThickness - 0.5*m_modulePlasticThickness),
                                                                    currentPlasticLogic, "modulePlastic", currentModuleLogic, false, 1);
    m_moduleUpperPlastic.push_back(currentUpperPlasticPlacement);
    m_moduleLowerPlastic.push_back(currentLowerPlasticPlacement);

    G4Box* currentFoamSolid = new G4Box("moduleFoam", 0.5*m_moduleLength[i], 0.5*m_moduleWidth[i], 0.5*m_moduleFoamThickness);
    G4LogicalVolume* currentFoamLogic = new G4LogicalVolume(currentFoamSolid, m_moduleFoamMaterial, "moduleFoam", 0, 0, 0);
    G4PVPlacement* currentUpperFoamPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,0.5*m_moduleFoamThickness), currentFoamLogic, "moduleFoam", currentModuleLogic, false, 0);
    G4PVPlacement* currentLowerFoamPlacement = new G4PVPlacement(currentInternalRotation, G4ThreeVector(0,0,-0.5*m_moduleFoamThickness),
                                                                 currentFoamLogic, "moduleFoam", currentModuleLogic, false, 1);
    m_moduleUpperFoam.push_back(currentUpperFoamPlacement); 
    m_moduleLowerFoam.push_back(currentLowerFoamPlacement); 

    // SENSTIVE DETECTOR FOR THESE FIBERS
    currentFiberLogic->SetSensitiveDetector(aFiberSD);  
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
  m_moduleGap = 2.*m_modulePlasticThickness + 2.*m_moduleFoamThickness;
  m_moduleHeight = 2. * m_moduleFiberThickness + m_moduleGap;
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
