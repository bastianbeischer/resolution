#include "RES_Module.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4VSensitiveDetector.hh"

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

void RES_Module::ComputeParameters()
{
  if (m_type == fiber) {
    m_height = 2*(m_foamThickness + m_plasticThickness + m_fiberThickness);
    G4double zInModule = m_foamThickness + m_plasticThickness + 0.5*m_fiberThickness;
    m_upperZ = m_placement.z() + zInModule + 0.5*m_fiberThickness;
    m_lowerZ = m_placement.z() - zInModule + 0.5*m_fiberThickness;
  }
  else if (m_type == silicon) {
    m_height = 2*(m_kaptonThickness + m_siliconThickness + m_gapSiliconThickness);
    G4double zInModule = 0.5*m_siliconThickness;
    m_upperZ = m_placement.z() + zInModule + 0.5*m_siliconThickness;
    m_lowerZ = m_placement.z() - zInModule + 0.5*m_siliconThickness;
  }

}

void RES_Module::SetDefaultValuesForFiber()
{
  m_type            = fiber;
  m_angle           = 0.;
  m_internalAngle   = 0.;
  m_length          = 30.*cm;
  m_width           = 6.912*cm;
  m_upperSigmaU     = m_length/sqrt(12);
  m_upperSigmaV     = 50.*um;
  m_upperSigmaZ     = 0.*cm;
  m_lowerSigmaU     = m_length/sqrt(12);
  m_lowerSigmaV     = 50.*um;
  m_lowerSigmaZ     = 0.*cm;
  m_upperEfficiency = 1.;
  m_lowerEfficiency = 1.;
}

void RES_Module::SetDefaultValuesForSilicon()
{
  m_type            = silicon;
  m_angle           = 0.;
  m_internalAngle   = M_PI/2.;
  m_length          = 30.*cm;
  m_width           = 10.*cm;
  m_upperSigmaU     = m_length/sqrt(12);
  m_upperSigmaV     = 10.*um;
  m_upperSigmaZ     = 0.*cm;
  m_lowerSigmaU     = m_length/sqrt(12);
  m_lowerSigmaV     = 30.*um;
  m_lowerSigmaZ     = 0.*cm;
  m_upperEfficiency = 1.;
  m_lowerEfficiency = 1.;
}

G4PVPlacement* RES_Module::Construct(G4VPhysicalVolume* mother, G4int copyNumber)
{
  // Compute Module Parameters
  ComputeParameters();

  G4Box* moduleSolid = new G4Box("module", 0.5*m_length, 0.5*m_width, 0.5*m_height);
  G4LogicalVolume* moduleLogic = new G4LogicalVolume(moduleSolid, m_moduleMaterial, "module", 0, 0, 0); // CHANGE MATERIAL HERE
  G4RotationMatrix* moduleRotation = new G4RotationMatrix(m_angle, 0., 0.);
  m_modulePlacement = new G4PVPlacement(moduleRotation, m_placement, moduleLogic, "module", mother->GetLogicalVolume(), false, copyNumber);

  // interior of modules
  if (m_type == fiber) {
    G4Box* fiberSolid = new G4Box("moduleFiber", 0.5*m_length, 0.5*m_width, 0.5*m_fiberThickness);
    G4LogicalVolume* fiberLogic = new G4LogicalVolume(fiberSolid, m_fiberMaterial, "moduleFiber", 0, 0, 0);
    m_upperFiberPlacement = new G4PVPlacement(0, G4ThreeVector(0, 0, m_foamThickness + m_plasticThickness + 0.5*m_fiberThickness),
                                                                  fiberLogic, "moduleFiber", moduleLogic, false, 0);
    m_lowerFiberPlacement = new G4PVPlacement(0, G4ThreeVector(0, 0, -m_foamThickness - m_plasticThickness - 0.5*m_fiberThickness),
                                                                  fiberLogic, "moduleFiber", moduleLogic, false, 1);

    G4Box* plasticSolid = new G4Box("modulePlastic", 0.5*m_length, 0.5*m_width, 0.5*m_plasticThickness);
    G4LogicalVolume* plasticLogic = new G4LogicalVolume(plasticSolid, m_plasticMaterial, "modulePlastic", 0, 0, 0);
    m_upperPlasticPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,m_foamThickness + 0.5*m_plasticThickness),
                                                                    plasticLogic, "modulePlastic", moduleLogic, false, 0);
    m_lowerPlasticPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,-m_foamThickness - 0.5*m_plasticThickness),
                                                                    plasticLogic, "modulePlastic", moduleLogic, false, 1);

    G4Box* foamSolid = new G4Box("moduleFoam", 0.5*m_length, 0.5*m_width, 0.5*m_foamThickness);
    G4LogicalVolume* foamLogic = new G4LogicalVolume(foamSolid, m_foamMaterial, "moduleFoam", 0, 0, 0);
    m_upperFoamPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,0.5*m_foamThickness), foamLogic, "moduleFoam", moduleLogic, false, 0);
    m_lowerFoamPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,-0.5*m_foamThickness),foamLogic, "moduleFoam", moduleLogic, false, 1);
  }
  else if (m_type == silicon) {
    G4Box* siliconSolid = new G4Box("moduleSilicon", 0.5*m_length, 0.5*m_width, 0.5*m_siliconThickness);
    G4LogicalVolume* siliconLogic = new G4LogicalVolume(siliconSolid, m_siliconMaterial, "moduleSilicon", 0, 0, 0);
    m_upperSiliconPlacement = new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*m_siliconThickness),
                                                                    siliconLogic, "moduleSilicon", moduleLogic, false, 0);
    m_lowerSiliconPlacement = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*m_siliconThickness),
                                                                    siliconLogic, "moduleSilicon", moduleLogic, false, 1);

    G4Box* kaptonSolid = new G4Box("moduleKapton", 0.5*m_length, 0.5*m_width, 0.5*m_kaptonThickness);
    G4LogicalVolume* kaptonLogic = new G4LogicalVolume(kaptonSolid, m_kaptonMaterial, "moduleKapton", 0, 0, 0);
    m_upperKaptonPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,m_siliconThickness + m_gapSiliconThickness + 0.5*m_kaptonThickness),
                                                                   kaptonLogic, "moduleKapton", moduleLogic, false, 0);
    m_lowerKaptonPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,-m_siliconThickness - m_gapSiliconThickness - 0.5*m_kaptonThickness),
                                                                   kaptonLogic, "moduleKapton", moduleLogic, false, 1);
  }

  return m_modulePlacement;
}

G4bool RES_Module::CheckIfTrackPassesThrough(G4ThreeVector position, G4ThreeVector direction)
{
  G4bool retVal = true;

  G4double dz = m_placement.z() - position.z();
  G4double l = dz / direction.z();
  G4ThreeVector currentPosition = position + l*direction;
  if (fabs(currentPosition.x() - m_placement.x()) > 0.5*m_length)   {retVal = false;}
  if (fabs(currentPosition.y() - m_placement.y()) > 0.5*m_width) {retVal = false;}

  return retVal;
}

void RES_Module::SetSD(G4VSensitiveDetector* aSD)
{
  // UPPER AND LOWER HAVE THE SAME POINTER TO THEIR LOGICAL VOLUME!
  if(m_type == fiber)
    m_upperFiberPlacement->GetLogicalVolume()->SetSensitiveDetector(aSD);
  else if(m_type == silicon)
    m_upperSiliconPlacement->GetLogicalVolume()->SetSensitiveDetector(aSD);
}

void RES_Module::SetVisibility()
{
  m_modulePlacement->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* visAtt;
  // UPPER AND LOWER HAVE THE SAME POINTER TO THEIR LOGICAL VOLUME!
  if (m_type == fiber) {
    visAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red
    m_upperFiberPlacement->GetLogicalVolume()->SetVisAttributes(visAtt);

    //    visAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // green
    visAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0)); // black
    m_upperFoamPlacement->GetLogicalVolume()->SetVisAttributes(visAtt);

    //    visAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // green
    visAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0)); // black
    m_upperPlasticPlacement->GetLogicalVolume()->SetVisAttributes(visAtt);

  }
  else if (m_type == silicon) {
    visAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // blue
    m_upperSiliconPlacement->GetLogicalVolume()->SetVisAttributes(visAtt);

    visAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // green
    m_upperKaptonPlacement->GetLogicalVolume()->SetVisAttributes(visAtt);
  }
}
