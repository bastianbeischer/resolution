#include "RES_TrackFitter.hh"

#include "RES_TrackFitMessenger.hh"
#include "RES_Event.hh"
#include "RES_RunManager.hh"
#include "RES_PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4ParticleGun.hh"

#include "globals.hh"

#include "blobel.h"

RES_TrackFitter::RES_TrackFitter() :
  m_verbose(0)
{
  m_messenger = new RES_TrackFitMessenger(this);
}

RES_TrackFitter::~RES_TrackFitter()
{
}

RES_Event RES_TrackFitter::Fit(RES_Event genEvent)
{
  G4double sigmaX = 1.*mm;
  G4double sigmaY = 50.*um;
  G4double sigmaZ = 0.*um;

  // G4double sigmaX = 0.;
  // G4double sigmaY = 0.;
  // G4double sigmaZ = 0.;

  G4int nHits = genEvent.GetNbOfHits();

  G4ThreeVector* pos = new G4ThreeVector[nHits];
  for (G4int i = 0; i < nHits; i++) {
    G4double x = CLHEP::RandGauss::shoot(genEvent.GetHit(i).x(), sigmaX);
    G4double y = CLHEP::RandGauss::shoot(genEvent.GetHit(i).y(), sigmaY);
    G4double z = CLHEP::RandGauss::shoot(genEvent.GetHit(i).z(), sigmaZ);
    pos[i] = G4ThreeVector(x,y,z);
  }

  G4double parameter[5];
  parameter[0] = pos[0].x();
  parameter[1] = pos[0].y();
  parameter[2] = (pos[1] - pos[0]).theta();
  parameter[3] = (pos[1] - pos[0]).phi();
  //  parameter[4] = genEvent.GetMomentum();
  parameter[4] = 100.*GeV;

  G4double step[5];
  step[0] = 0.01*parameter[0];
  step[1] = 0.01*parameter[1];
  step[2] = 0.01*parameter[2];
  step[3] = 0.01*parameter[3];
  step[4] = 0.01*parameter[4];

  G4int npar = 5; // negative to suppress printout
  G4int nflim = 3*abs(npar)*(abs(npar)+10); // dito
  if( m_verbose == 0 ) {
    npar = -npar;
    nflim = -nflim;
  }
  DVALLIN(npar,step,nflim);

  G4int conv = -1;
  G4int iter = 0;
  G4double chi2 = 0.;

  G4RunManager* runManager = G4RunManager::GetRunManager();
  const RES_PrimaryGeneratorAction* genAction = (RES_PrimaryGeneratorAction*) runManager->GetUserPrimaryGeneratorAction();
  G4ParticleGun* gun = genAction->GetParticleGun();

  while( conv <= 0 ){
    ++iter;
    chi2 = 0.;
    
    G4ThreeVector direction(cos(parameter[3])*sin(parameter[2]), sin(parameter[3])*sin(parameter[2]), cos(parameter[2]));
    G4ThreeVector position(parameter[0], parameter[1], pos[0].z());
    position -= 1.*cm * direction;
      
    gun->SetParticlePosition(position);
    gun->SetParticleMomentumDirection(direction);
    gun->SetParticleEnergy(parameter[4]);
    runManager->BeamOn(1);

    assert(nHits == m_currentRecEvent.GetNbOfHits());

    for( G4int i = 0 ; i < nHits ; i++ ) {
      chi2 += pow( genEvent.GetHit(i).x() - m_currentRecEvent.GetHit(i).x(), 2.) / (sigmaX*sigmaX);
      chi2 += pow( genEvent.GetHit(i).y() - m_currentRecEvent.GetHit(i).y(), 2.) / (sigmaY*sigmaY);
    }

    DVALLEY(chi2,parameter,conv);
  }

  if (conv) return m_currentRecEvent;
  else      return RES_Event();
}
