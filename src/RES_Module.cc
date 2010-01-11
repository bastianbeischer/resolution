#include "RES_Module.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4VPhysicalVolume.hh"

#include <cmath>

RES_Module::RES_Module()
{
  m_placement = G4ThreeVector(0,0,0);
  InitializeCommonValues();
  SetDefaultValuesForFiber();
}

RES_Module::RES_Module(G4ThreeVector placement)
{
  m_placement = placement;
  InitializeCommonValues();
  SetDefaultValuesForFiber();
}

RES_Module::~RES_Module()
{
}

void RES_Module::InitializeCommonValues()
{
  // Thickness
  m_foamThickness = 0.3*cm;
  m_plasticThickness = 0.01*cm;
  m_fiberThickness = 0.1*cm;

  m_kaptonThickness = 200.*um;
  m_siliconThickness = 160.*um;
  m_gapSiliconThickness = 0.964*cm; // chosen so that the 2*gap+2*kapton+2*silicon = 2cm

  // Placeholders
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

  // Define materials
  m_moduleMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_AIR" );
  m_plasticMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_POLYCARBONATE" );
  m_foamMaterial = rohacell;
  m_fiberMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_POLYSTYRENE" );
  m_siliconMaterial = Si;
  m_kaptonMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_KAPTON" );
}

void RES_Module::SetDefaultValuesForFiber()
{
  m_type          = fiber;
  m_angle         = 0.;
  m_internalAngle = 0.;
  m_length        = 30.*cm;
  m_width         = 6.912*cm;
  m_upperSigmaU   = m_length/sqrt(12);
  m_upperSigmaV   = 50.*um;
  m_upperSigmaZ   = 0.*cm;
  m_lowerSigmaU   = m_length/sqrt(12);
  m_lowerSigmaV   = 50.*um;
  m_lowerSigmaZ   = 0.*cm;
  m_efficiency    = 1.;
}

void RES_Module::SetDefaultValuesForSilicon()
{
  m_type          = silicon;
  m_angle         = 0.;
  m_internalAngle = M_PI/2.;
  m_length        = 30.*cm;
  m_width         = 10.*cm;
  m_upperSigmaU   = m_length/sqrt(12);
  m_upperSigmaV   = 10.*um;
  m_upperSigmaZ   = 0.*cm;
  m_lowerSigmaU   = m_length/sqrt(12);
  m_lowerSigmaV   = 30.*um;
  m_lowerSigmaZ   = 0.*cm;
  m_efficiency    = 1.;
}

G4VPhysicalVolume* RES_Module::Construct(G4VPhysicalVolume* mother)
{
  G4double height = m_type == fiber ? m_moduleHeightFiber : m_moduleHeightSilicon;
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
    G4PVPlacement* currentLowerPlasticPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,-m_moduleFoamThickness - 0.5*m_modulePlasticThickness),
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
