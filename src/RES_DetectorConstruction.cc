// $Id: RES_DetectorConstruction.cc,v 1.19 2009/12/11 12:52:24 beischer Exp $

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

  G4Material* Si = new G4Material("Silicon", z=14, a=28.09*g/mole, density=2.33*g/cm3);

  // define materials
  m_worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_AIR" );
  m_moduleMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_AIR" );
  m_modulePlasticMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_POLYCARBONATE" );
  m_moduleFoamMaterial = rohacell;
  m_moduleFiberMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_POLYSTYRENE" );
  m_moduleSiliconMaterial = Si;
  m_moduleKaptonMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_KAPTON" );

  // define world dimensions
  m_worldX = 1.0*m;
  m_worldY = 1.0*m;
  m_worldZ = 1.0*m;


  // define default module parameters
  m_moduleDefaultWidthFiber = 6.912 * cm;
  m_moduleDefaultLengthFiber = 30. * cm;
  m_moduleFoamThickness = 0.3 * cm;
  m_modulePlasticThickness = 0.01 * cm;
  m_moduleFiberThickness = 0.1 * cm;
  m_moduleDefaultSigmaUFiber = m_moduleDefaultLengthFiber/sqrt(12);
  m_moduleDefaultSigmaVFiber = 50. * um;
  m_moduleDefaultSigmaZFiber = 0. * cm;

  m_moduleDefaultWidthSilicon = 10.*cm;
  m_moduleDefaultLengthSilicon = 30.*cm;
  m_moduleKaptonThickness = 200.*um;
  m_moduleSiliconThickness = 160.*um;
  m_moduleGapSilicon = 0.964*cm; // chosen so that the 2*gap+2*kapton+2*silicon = 2cm
  m_moduleDefaultUpperSigmaUSilicon = m_moduleDefaultLengthSilicon/sqrt(12);
  m_moduleDefaultUpperSigmaVSilicon = 10. * um;
  m_moduleDefaultUpperSigmaZSilicon = 0. * cm;
  m_moduleDefaultLowerSigmaUSilicon = m_moduleDefaultLengthSilicon/sqrt(12);
  m_moduleDefaultLowerSigmaVSilicon = 30. * um;
  m_moduleDefaultLowerSigmaZSilicon = 0. * cm;
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
  RES_SD* aSD = new RES_SD( fiberSDName );
  SDman->AddNewDetector( aSD );

  // --------------------------------------------
  // Volumes and Placements (also add SDs to volumes)
  // --------------------------------------------

  // world
  G4Box* worldSolid         = new G4Box("world", 0.5*m_worldX, 0.5*m_worldY, 0.5*m_worldZ);
  G4LogicalVolume* worldLog = new G4LogicalVolume(worldSolid, m_worldMaterial, "world", 0, 0, 0);
  m_world                   = new G4PVPlacement(0, G4ThreeVector(), worldLog, "world", 0, false, 0);

  for (unsigned int i = 0; i < m_modulePlacements.size(); i++) {
    // detector modules

    G4double height = m_moduleType[i] == fiber ? m_moduleHeightFiber : m_moduleHeightSilicon;
    G4Box* currentModuleSolid = new G4Box("module", 0.5*m_moduleLength[i], 0.5*m_moduleWidth[i], 0.5*height);
    G4LogicalVolume* currentModuleLogic = new G4LogicalVolume(currentModuleSolid, m_moduleMaterial, "module", 0, 0, 0); // CHANGE MATERIAL HERE
    G4RotationMatrix* currentModuleRotation = new G4RotationMatrix(m_moduleAngles[i], 0., 0.);
    G4PVPlacement* currentModulePlacement = new G4PVPlacement(currentModuleRotation, m_modulePlacements[i], currentModuleLogic, "module", worldLog, false, i);
    m_modules.push_back(currentModulePlacement);

    // interior of modules
    if (m_moduleType[i] == fiber) {
      G4Box* currentFiberSolid = new G4Box("moduleFiber", 0.5*m_moduleLength[i], 0.5*m_moduleWidth[i], 0.5*m_moduleFiberThickness);
      G4LogicalVolume* currentFiberLogic = new G4LogicalVolume(currentFiberSolid, m_moduleFiberMaterial, "moduleFiber", 0, 0, 0);
      G4PVPlacement* currentUpperFiberPlacement = new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*m_moduleGapFiber + 0.5*m_moduleFiberThickness),
                                                                    currentFiberLogic, "moduleFiber", currentModuleLogic, false, 0);
      G4PVPlacement* currentLowerFiberPlacement = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*m_moduleGapFiber - 0.5*m_moduleFiberThickness),
                                                                    currentFiberLogic, "moduleFiber", currentModuleLogic, false, 1);
      m_moduleUpperFiber.push_back(currentUpperFiberPlacement);
      m_moduleLowerFiber.push_back(currentLowerFiberPlacement);

      G4Box* currentPlasticSolid = new G4Box("modulePlastic", 0.5*m_moduleLength[i], 0.5*m_moduleWidth[i], 0.5*m_modulePlasticThickness);
      G4LogicalVolume* currentPlasticLogic = new G4LogicalVolume(currentPlasticSolid, m_modulePlasticMaterial, "modulePlastic", 0, 0, 0);
      G4PVPlacement* currentUpperPlasticPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,m_moduleFoamThickness + 0.5*m_modulePlasticThickness),
                                                                      currentPlasticLogic, "modulePlastic", currentModuleLogic, false, 0);
      G4PVPlacement* currentLowerPlasticPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,m_moduleFoamThickness - 0.5*m_modulePlasticThickness),
                                                                      currentPlasticLogic, "modulePlastic", currentModuleLogic, false, 1);
      m_moduleUpperPlastic.push_back(currentUpperPlasticPlacement);
      m_moduleLowerPlastic.push_back(currentLowerPlasticPlacement);

      G4Box* currentFoamSolid = new G4Box("moduleFoam", 0.5*m_moduleLength[i], 0.5*m_moduleWidth[i], 0.5*m_moduleFoamThickness);
      G4LogicalVolume* currentFoamLogic = new G4LogicalVolume(currentFoamSolid, m_moduleFoamMaterial, "moduleFoam", 0, 0, 0);
      G4PVPlacement* currentUpperFoamPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,0.5*m_moduleFoamThickness), currentFoamLogic, "moduleFoam", currentModuleLogic, false, 0);
      G4PVPlacement* currentLowerFoamPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,-0.5*m_moduleFoamThickness),currentFoamLogic, "moduleFoam", currentModuleLogic, false, 1);
      m_moduleUpperFoam.push_back(currentUpperFoamPlacement); 
      m_moduleLowerFoam.push_back(currentLowerFoamPlacement); 

      // SENSTIVE DETECTOR FOR THESE FIBERS
      currentFiberLogic->SetSensitiveDetector(aSD);  
    }
    else if (m_moduleType[i] == silicon) {
      G4Box* currentSiliconSolid = new G4Box("moduleSilicon", 0.5*m_moduleLength[i], 0.5*m_moduleWidth[i], 0.5*m_moduleSiliconThickness);
      G4LogicalVolume* currentSiliconLogic = new G4LogicalVolume(currentSiliconSolid, m_moduleSiliconMaterial, "moduleSilicon", 0, 0, 0);
      G4PVPlacement* currentUpperSiliconPlacement = new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*m_moduleSiliconThickness),
                                                                    currentSiliconLogic, "moduleSilicon", currentModuleLogic, false, 0);
      G4PVPlacement* currentLowerSiliconPlacement = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*m_moduleSiliconThickness),
                                                                    currentSiliconLogic, "moduleSilicon", currentModuleLogic, false, 1);
      m_moduleUpperSilicon.push_back(currentUpperSiliconPlacement);
      m_moduleLowerSilicon.push_back(currentLowerSiliconPlacement);

      G4Box* currentKaptonSolid = new G4Box("moduleKapton", 0.5*m_moduleLength[i], 0.5*m_moduleWidth[i], 0.5*m_moduleKaptonThickness);
      G4LogicalVolume* currentKaptonLogic = new G4LogicalVolume(currentKaptonSolid, m_moduleKaptonMaterial, "moduleKapton", 0, 0, 0);
      G4PVPlacement* currentUpperKaptonPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,m_moduleSiliconThickness + m_moduleGapSilicon + 0.5*m_moduleKaptonThickness),
                                                                      currentKaptonLogic, "moduleKapton", currentModuleLogic, false, 0);
      G4PVPlacement* currentLowerKaptonPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,-m_moduleSiliconThickness - m_moduleGapSilicon - 0.5*m_moduleKaptonThickness),
                                                                      currentKaptonLogic, "moduleKapton", currentModuleLogic, false, 1);
      m_moduleUpperKapton.push_back(currentUpperKaptonPlacement);
      m_moduleLowerKapton.push_back(currentLowerKaptonPlacement);

      // SENSTIVE DETECTOR FOR THESE SILICON LADDERS
      currentSiliconLogic->SetSensitiveDetector(aSD);  
    }
    
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
