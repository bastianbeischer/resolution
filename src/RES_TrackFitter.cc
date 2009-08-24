#include "RES_TrackFitter.hh"

#include "RES_TrackFitMessenger.hh"
#include "RES_Event.hh"
#include "RES_RunManager.hh"
#include "RES_DetectorConstruction.hh"
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

  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();

  G4double moduleLength = det->GetModuleLength();
  G4double stereoAngle  = det->GetStereoAngle();
  //G4double stereoAngle  = M_PI/2.;

  G4double sigmaX = moduleLength/sqrt(12.);
  G4double sigmaY = 50.*um;
  G4double sigmaZ = 0.*um;

  G4int nHits = genEvent.GetNbOfHits();

  G4ThreeVector* pos = new G4ThreeVector[nHits];
  for (G4int i = 0; i < nHits; i++) {
    G4double angle;
    if (i % 2) angle = 0.;
    else       angle = stereoAngle;
    G4double s = sin(angle);
    G4double c = cos(angle);
      
    double localSigmaX = sqrt(c*c*sigmaX*sigmaX + s*s*sigmaY*sigmaY);
    double localSigmaY = sqrt(s*s*sigmaX*sigmaX + c*c*sigmaY*sigmaY);

    G4double x = CLHEP::RandGauss::shoot(genEvent.GetHit(i).x(), localSigmaX);
    G4double y = CLHEP::RandGauss::shoot(genEvent.GetHit(i).y(), localSigmaY);
    G4double z = CLHEP::RandGauss::shoot(genEvent.GetHit(i).z(), sigmaZ);

    pos[i] = G4ThreeVector(x,y,z);
  }

  // G4double theta = (pos[6].y() - pos[4].y()) / (pos[6].z() - pos[4].z()) - (pos[2].y() - pos[0].y()) / (pos[2].z() - pos[0].z());
  // G4double pSagitta = 0.3 * 0.3 * sqrt(pow(pos[4].y() - pos[2].y(),2.) + pow(pos[4].z() - pos[2].z(), 2.))/m / theta;
  // G4cout << " p: " << pSagitta << G4endl;
  
  G4double parameter[5];
  parameter[0] = pos[0].x();
  parameter[1] = pos[0].y();
  parameter[2] = (pos[2] - pos[0]).theta();
  parameter[3] = (pos[2] - pos[0]).phi();
  // parameter[4] = pSagitta;
  parameter[4] = 2.0 * GeV;

  G4double step[5];
  step[0] = 0.01*parameter[0];
  step[1] = 0.01*parameter[1];
  step[2] = 0.01*parameter[2];
  step[3] = 0.01*parameter[3];
  step[4] = 0.1*parameter[4];

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

    G4int nRecHits = m_currentRecEvent.GetNbOfHits();
    assert(nHits == nRecHits);

    // TODO: exception handling for sigmaX = sigmaY = 0.
    for( G4int i = 0 ; i < nHits ; i++ ) {
      G4double angle;
      if (i % 2) angle = 0.;
      else       angle = stereoAngle;
      G4double s = sin(angle);
      G4double c = cos(angle);
      G4double dx = pos[i].x() - m_currentRecEvent.GetHit(i).x();
      G4double dy = pos[i].y() - m_currentRecEvent.GetHit(i).y();

      chi2 += pow(dx*sigmaX*s, 2.);
      chi2 += pow(dy*sigmaY*s, 2.);
      chi2 += pow(dx*sigmaY*c, 2.);
      chi2 += pow(dy*sigmaX*c, 2.);
      chi2 += 2.*dx*dy*s*c*(pow(sigmaY,2.) - pow(sigmaX, 2.));

      chi2 /= pow(sigmaX*sigmaY, 2.);
    }

    DVALLEY(chi2,parameter,conv);
  }

  if (conv) return m_currentRecEvent;
  else      return RES_Event();
}
