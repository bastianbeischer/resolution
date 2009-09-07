#include "RES_PrimaryGeneratorAction.hh"

#include "RES_PrimaryGeneratorMessenger.hh"
#include "RES_DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4ParticleTypes.hh"

#include "CLHEP/Random/RandFlat.h"

RES_PrimaryGeneratorAction::RES_PrimaryGeneratorAction() :
  m_randomOrigin(false),
  m_randomDirection(false),
  m_z_start(0.5*m)
{
  m_messenger = new RES_PrimaryGeneratorMessenger(this);

  G4int nParticle = 1;
  m_particleGun = new G4ParticleGun(nParticle);
  
  m_particleGun->SetParticleDefinition(G4ChargedGeantino::ChargedGeantino());
  m_particleGun->SetParticleEnergy(1.0*GeV);
  m_particleGun->SetParticlePosition(G4ThreeVector(0.0*cm, 0.0*cm, m_z_start));
  m_particleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, -1.0));
}

RES_PrimaryGeneratorAction::~RES_PrimaryGeneratorAction()
{
  delete m_particleGun;
}

void RES_PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  if (m_randomOrigin) {
    RES_DetectorConstruction* detector = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction(); 
    G4double world_x = detector->GetWorldX();
    G4double world_y = detector->GetWorldY();
    G4double x_start = CLHEP::RandFlat::shoot(-0.5*world_x, 0.5*world_x);
    G4double y_start = CLHEP::RandFlat::shoot(-0.5*world_y, 0.5*world_y);
    m_particleGun->SetParticlePosition(G4ThreeVector(x_start, y_start, m_z_start));
  }
  
  if (m_randomDirection) {
    G4double cos_theta_square = CLHEP::RandFlat::shoot(0.0, 1.0);
    G4double theta = acos(sqrt(cos_theta_square));
    G4double phi = CLHEP::RandFlat::shoot(0.0, 2*M_PI);
    G4double direction_x = cos(phi) * sin(theta);
    G4double direction_y = sin(phi) * sin(theta);
    G4double direction_z = -cos(theta);
    m_particleGun->SetParticleMomentumDirection(G4ThreeVector(direction_x, direction_y, direction_z));
  }

  m_particleGun->GeneratePrimaryVertex(event);
}
