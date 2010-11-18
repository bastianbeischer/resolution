#include "RES_Module.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
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
  m_fiberThickness = 1.15*mm;
  m_carbonFiberThickness = 0.25*mm;
  m_epoxyThickness = 0.5*mm;
  m_foamThickness = 7.6*mm - 2*m_epoxyThickness;

  m_subtractHoles = false;

  m_kaptonThickness = 200.*um;
  m_siliconThickness = 160.*um;
  m_gapSiliconThickness = 0.964*cm; // chosen so that the 2*gap+2*kapton+2*silicon = 2cm

  // Placeholders
  G4String symbol;
  G4double density,z,a;
  G4double fraction;
  G4int components,natoms;

  // Elements
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);

  // Epoxy
  G4Material* Epoxy = new G4Material("Epoxy", 1.3*g/cm3, components=3);
  Epoxy->AddElement(H,natoms=44);
  Epoxy->AddElement(C,natoms=15);
  Epoxy->AddElement(O,natoms=7);

  // Carbon Fiber Layer
  G4Material* CarbonFiber = new G4Material("Carbon Fibre", density=1.8*g/cm3, components = 2);
  CarbonFiber->AddElement(C, fraction=0.6);
  CarbonFiber->AddMaterial(Epoxy, fraction=0.4);

  // Polystyrene
  G4Material* Polystyrene = G4NistManager::Instance()->FindOrBuildMaterial( "G4_POLYSTYRENE" );
  
  // Fiber with Epoxy (80% Fiber, 20% Epoxy) (amount of glue derived from measured mass of a mat)
  G4Material* Fiber = new G4Material("Fiber", density=0.8*Polystyrene->GetDensity() + 0.2*Epoxy->GetDensity(), components=2);
  Fiber->AddMaterial(Polystyrene,fraction=0.8);
  Fiber->AddMaterial(Epoxy,fraction=0.2);

  // Rohacell
  G4Material* Rohacell = new G4Material( "rohacell", density=32.*kg/m3, components = 4 );
  Rohacell->AddElement(H, natoms=5); // ???
  Rohacell->AddElement(C, natoms=3); // ???
  Rohacell->AddElement(N, natoms=1); // ???
  Rohacell->AddElement(O, natoms=1); // ???

  // Silicon
  G4Material* Si = new G4Material("Silicon", z=14, a=28.09*g/mole, density=2.33*g/cm3);

  // Define materials
  m_moduleMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_AIR" );

  m_carbonFiberMaterial = CarbonFiber;
  m_epoxyMaterial = Epoxy;
  m_foamMaterial = Rohacell;
  // m_carbonFiberMaterial = m_moduleMaterial;
  // m_epoxyMaterial = m_moduleMaterial;
  // m_foamMaterial = m_moduleMaterial;

  m_fiberMaterial = Fiber;

  m_siliconMaterial = Si;
  m_kaptonMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_KAPTON" );
}

void RES_Module::ComputeParameters()
{
  if (m_type == fiber) {
    m_height = 2*(m_foamThickness + 2*m_carbonFiberThickness + m_fiberThickness + 2*m_epoxyThickness);
    G4double zInModule = m_epoxyThickness + m_carbonFiberThickness + m_foamThickness + m_carbonFiberThickness + m_epoxyThickness + 0.5*m_fiberThickness;
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
  m_angle           = -0.5 * M_PI/180.;
  m_internalAngle   = 1.0 * M_PI/180.;
  m_length          = 40.*cm;
  m_width           = 6.5*cm;
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
  m_length          = 40.*cm;
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
  //  G4RotationMatrix* moduleRotation = new G4RotationMatrix(m_angle, 0., 0.);
  G4RotationMatrix* moduleRotation = new G4RotationMatrix(0., 0., 0.);
  m_modulePlacement = new G4PVPlacement(moduleRotation, m_placement, moduleLogic, "module", mother->GetLogicalVolume(), false, copyNumber);

  // interior of modules
  if (m_type == fiber) {
    G4VSolid*           fiberSolid       = new G4Box("moduleFiber", 0.5*m_length, 0.5*m_width, 0.5*m_fiberThickness);
    G4VSolid*           carbonFiberSolid = new G4Box("moduleCarbonFiber", 0.5*m_length, 0.5*m_width, 0.5*m_carbonFiberThickness);
    G4VSolid*           epoxySolid       = new G4Box("moduleFoam", 0.5*m_length, 0.5*m_width, 0.5*m_epoxyThickness);
    G4VSolid*           foamSolid        = new G4Box("moduleFoam", 0.5*m_length, 0.5*m_width, 0.5*m_foamThickness);

    if (m_subtractHoles) {
      G4double holePositions[7] = {-16.5*cm, -11.*cm, -5.5*cm, 0.*cm, 5.5*cm, 11.*cm, 16.5*cm};
      G4double holeSideLength = 4.5*cm;
      
      G4Box* hole = new G4Box("hole", 0.5*holeSideLength, 0.5*holeSideLength, 0.51*m_carbonFiberThickness);
      G4SubtractionSolid* carbonFiberSubSolid = new G4SubtractionSolid("carbonFiberSubSolid", carbonFiberSolid, hole, 0, G4ThreeVector(holePositions[0], 0, 0));
      for (int i = 1; i < 7; i++)
        carbonFiberSubSolid = new G4SubtractionSolid("carbonFiberSubSolid", carbonFiberSubSolid, hole, 0, G4ThreeVector(holePositions[i], 0, 0));
      carbonFiberSolid = carbonFiberSubSolid;

      hole = new G4Box("hole", 0.5*holeSideLength, 0.5*holeSideLength, 0.51*m_epoxyThickness);
      G4SubtractionSolid* epoxySubSolid = new G4SubtractionSolid("epoxySubSolid", epoxySolid, hole, 0, G4ThreeVector(holePositions[0], 0, 0));
      for (int i = 1; i < 7; i++)
        epoxySubSolid = new G4SubtractionSolid("epoxySubSolid", epoxySubSolid, hole, 0, G4ThreeVector(holePositions[i], 0, 0));
      epoxySolid = epoxySubSolid;

      hole = new G4Box("hole", 0.5*holeSideLength, 0.5*holeSideLength, 0.51*m_foamThickness);
      G4SubtractionSolid* foamSubSolid = new G4SubtractionSolid("foamSubSolid", foamSolid, hole, 0, G4ThreeVector(holePositions[0], 0, 0));
      for (int i = 1; i < 7; i++)
        foamSubSolid = new G4SubtractionSolid("foamSubSolid", foamSubSolid, hole, 0, G4ThreeVector(holePositions[i], 0, 0));
      foamSolid = foamSubSolid;
    }

    G4LogicalVolume* fiberLogic       = new G4LogicalVolume(fiberSolid, m_fiberMaterial, "moduleFiber", 0, 0, 0);
    G4LogicalVolume* carbonFiberLogic = new G4LogicalVolume(carbonFiberSolid, m_carbonFiberMaterial, "moduleCarbonFiber", 0, 0, 0);
    G4LogicalVolume* epoxyLogic       = new G4LogicalVolume(epoxySolid, m_epoxyMaterial, "moduleEpoxy", 0, 0, 0);
    G4LogicalVolume* foamLogic        = new G4LogicalVolume(foamSolid, m_foamMaterial, "moduleFoam", 0, 0, 0);

    m_upperFiberPlacement        = new G4PVPlacement(0, G4ThreeVector(0,0,m_epoxyThickness + m_carbonFiberThickness + m_foamThickness + m_carbonFiberThickness + m_epoxyThickness + 0.5*m_fiberThickness),
                                                     fiberLogic, "moduleFiber", moduleLogic, false, 0);
    m_firstEpoxyPlacement        = new G4PVPlacement(0, G4ThreeVector(0,0,m_epoxyThickness + m_carbonFiberThickness + m_foamThickness + m_carbonFiberThickness + 0.5*m_epoxyThickness),
                                                     epoxyLogic, "moduleEpoxy", moduleLogic, false, 0);
    m_firstCarbonFiberPlacement  = new G4PVPlacement(0, G4ThreeVector(0,0,m_epoxyThickness + m_carbonFiberThickness + m_foamThickness + 0.5*m_carbonFiberThickness),
                                                     carbonFiberLogic, "moduleCarbonFiber", moduleLogic, false, 0);
    m_upperFoamPlacement         = new G4PVPlacement(0, G4ThreeVector(0,0,m_epoxyThickness + m_carbonFiberThickness + 0.5*m_foamThickness),
                                                     foamLogic, "moduleFoam", moduleLogic, false, 0);
    m_secondCarbonFiberPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,m_epoxyThickness + 0.5*m_carbonFiberThickness),
                                                     carbonFiberLogic, "moduleCarbonFiber", moduleLogic, false, 1);
    m_secondEpoxyPlacement       = new G4PVPlacement(0, G4ThreeVector(0,0,0.5*m_epoxyThickness),
                                                     epoxyLogic, "moduleEpoxy",moduleLogic, false, 1);
    m_thirdEpoxyPlacement        = new G4PVPlacement(0, G4ThreeVector(0,0,-0.5*m_epoxyThickness),
                                                     epoxyLogic, "moduleEpoxy", moduleLogic, false, 2);
    m_thirdCarbonFiberPlacement  = new G4PVPlacement(0, G4ThreeVector(0,0,-m_epoxyThickness - 0.5*m_carbonFiberThickness),
                                                     carbonFiberLogic, "moduleCarbonFiber", moduleLogic, false, 2);
    m_lowerFoamPlacement         = new G4PVPlacement(0, G4ThreeVector(0,0,-m_epoxyThickness - m_carbonFiberThickness - 0.5*m_foamThickness),
                                                     foamLogic, "moduleFoam", moduleLogic, false, 1);
    m_fourthCarbonFiberPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,-m_epoxyThickness - m_carbonFiberThickness - m_foamThickness - 0.5*m_carbonFiberThickness),
                                                     carbonFiberLogic, "moduleCarbonFiber", moduleLogic, false, 3);
    m_fourthEpoxyPlacement       = new G4PVPlacement(0, G4ThreeVector(0,0,-m_epoxyThickness - m_carbonFiberThickness - m_foamThickness - m_carbonFiberThickness - 0.5*m_epoxyThickness),
                                                     epoxyLogic, "moduleEpoxy", moduleLogic, false, 3);
    m_lowerFiberPlacement        = new G4PVPlacement(0, G4ThreeVector(0,0,-m_epoxyThickness - m_carbonFiberThickness - m_foamThickness - m_carbonFiberThickness - m_epoxyThickness - 0.5*m_fiberThickness),
                                                     fiberLogic, "moduleFiber", moduleLogic, false, 1);
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
  G4double dz = m_placement.z() - position.z();
  G4double l = dz / direction.z();
  G4ThreeVector currentPosition = position + l*direction;
  if (fabs(currentPosition.x() - m_placement.x()) > 0.5*m_length)   {return false;}
  if (fabs(currentPosition.y() - m_placement.y()) > 0.5*m_width) {return false;}

  return true;
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
    //m_upperFiberPlacement->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);

    visAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // green
    //    visAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0)); // black
    m_upperFoamPlacement->GetLogicalVolume()->SetVisAttributes(visAtt);

    visAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0)); // yellow
    //visAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0)); // black
    m_firstCarbonFiberPlacement->GetLogicalVolume()->SetVisAttributes(visAtt);

    visAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // blue
    //visAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0)); // black
    m_firstEpoxyPlacement->GetLogicalVolume()->SetVisAttributes(visAtt);
  }
  else if (m_type == silicon) {
    visAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // blue
    m_upperSiliconPlacement->GetLogicalVolume()->SetVisAttributes(visAtt);

    visAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // green
    m_upperKaptonPlacement->GetLogicalVolume()->SetVisAttributes(visAtt);
  }
}

void RES_Module::PrintMaterials()
{
  G4cout << " <------------------------ MATERIALS --------------------------> " << G4endl;
  G4cout << "  Module: " << G4endl;
  G4cout << *m_moduleMaterial << G4endl;
  G4cout << "  Carbon Fiber: " << G4endl;
  G4cout << *m_carbonFiberMaterial << G4endl;
  G4cout << "  Epoxy: " << G4endl;
  G4cout << *m_epoxyMaterial << G4endl;
  G4cout << "  Foam: " << G4endl;
  G4cout << *m_foamMaterial << G4endl;
  G4cout << "  Fiber: " << G4endl;
  G4cout << *m_fiberMaterial << G4endl;

  // // Polystyrene
  // G4Material* Polystyrene = G4NistManager::Instance()->FindOrBuildMaterial( "G4_POLYSTYRENE" );
  // G4cout << *Polystyrene << G4endl;

  G4cout << m_modulePlacement->GetLogicalVolume()->GetMass() / g << " g" << G4endl;
  //  G4cout << (1./3.)* m_upperFiberPlacement->GetLogicalVolume()->GetMass() / g << " g" << G4endl;

}
