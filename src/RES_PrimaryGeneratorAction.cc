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
  if (m_randomOrigin || m_randomDirection) {
    RES_DetectorConstruction* detector = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction(); 
    G4double x_start = 0.;
    G4double y_start = 0.;
    G4double direction_x = 0.;
    G4double direction_y = 0.;
    G4double direction_z = -1.;
    G4ThreeVector position,direction;

    G4bool trackInAcceptance = false;
    while (!trackInAcceptance) {
      if (m_randomOrigin) {
        G4double world_x = detector->GetWorldX();
        G4double world_y = detector->GetWorldY();
        x_start = CLHEP::RandFlat::shoot(-0.5*world_x, 0.5*world_x);
        y_start = CLHEP::RandFlat::shoot(-0.5*world_y, 0.5*world_y);
      }
  
      if (m_randomDirection) {
        G4double cos_theta_square = CLHEP::RandFlat::shoot(0.0, 1.0);
        G4double theta = acos(sqrt(cos_theta_square));
        G4double phi = CLHEP::RandFlat::shoot(0.0, 2*M_PI);
        direction_x = cos(phi) * sin(theta);
        direction_y = sin(phi) * sin(theta);
        direction_z = -cos(theta);
      }

      position = G4ThreeVector(x_start, y_start, m_z_start);
      direction = G4ThreeVector(direction_x, direction_y, direction_z);

      G4cout << "new test particle" << G4endl;
      trackInAcceptance = detector->TrackInAcceptance(position,direction);
    }

    m_particleGun->SetParticlePosition(position);
    m_particleGun->SetParticleMomentumDirection(direction);
  }

  m_particleGun->GeneratePrimaryVertex(event);
}
