#include "RES_TrackFitter.hh"

#include "RES_TrackFitMessenger.hh"
#include "RES_Event.hh"
#include "RES_RunManager.hh"
#include "RES_DetectorConstruction.hh"
#include "RES_PrimaryGeneratorAction.hh"
#include "RES_FiberHit.hh"

#include "G4RunManager.hh"
#include "G4ParticleGun.hh"

#include "globals.hh"

#include "blobel.h"

RES_TrackFitter::RES_TrackFitter() :
  m_verbose(0),
  m_smearedHits(0),
  m_fitMethod(blobel)
{
  m_messenger = new RES_TrackFitMessenger(this);
}

RES_TrackFitter::~RES_TrackFitter()
{
  delete m_smearedHits;
  delete m_parameter;
  delete m_messenger;
}

RES_Event RES_TrackFitter::Fit()
{
  SetSpatialResolutions();
  SmearHits();
  CalculateStartParameters();

  G4int conv = -1;

  switch (m_fitMethod) {
  case blobel:
    conv = DoBlobelFit();
    break;
  case minuit:
    conv = DoMinuitFit();
    break;
  default:
    break;
  }
    
  if (conv) return m_currentRecEvent;
  else      return RES_Event();
}

void RES_TrackFitter::SetSpatialResolutions()
{
  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  G4double moduleLength = det->GetModuleLength();

  m_sigmaX = moduleLength/sqrt(12.);
  m_sigmaY = 50.*um;
  m_sigmaZ = 0.*um;
}

void RES_TrackFitter::SmearHits()
{
  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();

  G4int nHits = m_currentGenEvent.GetNbOfHits();

  delete m_smearedHits;
  m_smearedHits = new G4ThreeVector[nHits];

  for (G4int i = 0; i < nHits; i++) {
    G4int iModule = m_currentGenEvent.GetModuleID(i);
    G4double angle = det->GetModuleAngle(iModule);
    G4double s = sin(angle);
    G4double c = cos(angle);
      
    double currentSigmaX = sqrt(c*c*m_sigmaX*m_sigmaX + s*s*m_sigmaY*m_sigmaY);
    double currentSigmaY = sqrt(s*s*m_sigmaX*m_sigmaX + c*c*m_sigmaY*m_sigmaY);

    G4double x = CLHEP::RandGauss::shoot(m_currentGenEvent.GetHitPosition(i).x(), currentSigmaX);
    G4double y = CLHEP::RandGauss::shoot(m_currentGenEvent.GetHitPosition(i).y(), currentSigmaY);
    G4double z = CLHEP::RandGauss::shoot(m_currentGenEvent.GetHitPosition(i).z(), m_sigmaZ);

    m_smearedHits[i] = G4ThreeVector(x,y,z);
  }
}

void RES_TrackFitter::CalculateStartParameters()
{
  G4double theta =   (m_smearedHits[6].y() - m_smearedHits[4].y()) / (m_smearedHits[6].z() - m_smearedHits[4].z())
                   - (m_smearedHits[2].y() - m_smearedHits[0].y()) / (m_smearedHits[2].z() - m_smearedHits[0].z());
  G4double L = sqrt(pow(m_smearedHits[4].y() - m_smearedHits[3].y(),2.) + pow(m_smearedHits[4].z() - m_smearedHits[3].z(), 2.))/m;
  G4double B = 0.3;
  G4double pStart = 0.3 * B * L / theta * GeV;
  
  delete m_parameter;
  m_parameter = new G4double[5];
  m_parameter[0] = m_smearedHits[0].x();
  m_parameter[1] = m_smearedHits[0].y();
  m_parameter[2] = (m_smearedHits[2] - m_smearedHits[0]).theta();
  m_parameter[3] = (m_smearedHits[2] - m_smearedHits[0]).phi();
  m_parameter[4] = pStart;
}

G4int RES_TrackFitter::DoBlobelFit()
{
  G4double step[5];
  step[0] = 0.01*m_parameter[0];
  step[1] = 0.01*m_parameter[1];
  step[2] = 0.01*m_parameter[2];
  step[3] = 0.01*m_parameter[3];
  step[4] = 0.01*m_parameter[4];

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

  while( conv <= 0 ){
    ++iter;
    chi2 = Chi2();
    DVALLEY(chi2,m_parameter,conv);
  }

  int nHits = m_currentGenEvent.GetNbOfHits();
  double dof = nHits - npar;
  m_currentRecEvent.SetChi2OverDof(chi2/dof);

  return conv;
}

G4int RES_TrackFitter::DoMinuitFit()
{
  return -1;
}

G4double RES_TrackFitter::Chi2()
{
  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();

  G4RunManager* runManager = G4RunManager::GetRunManager();
  const RES_PrimaryGeneratorAction* genAction = (RES_PrimaryGeneratorAction*) runManager->GetUserPrimaryGeneratorAction();
  G4ParticleGun* gun = genAction->GetParticleGun();

  G4double chi2 = 0.;
    
  G4ThreeVector direction(cos(m_parameter[3])*sin(m_parameter[2]), sin(m_parameter[3])*sin(m_parameter[2]), cos(m_parameter[2]));
  G4ThreeVector position(m_parameter[0], m_parameter[1], m_smearedHits[0].z());
  position -= 1.*cm * direction;
      
  gun->SetParticlePosition(position);
  gun->SetParticleMomentumDirection(direction);
  gun->SetParticleEnergy(m_parameter[4]);
  runManager->BeamOn(1);

  G4int nHits = m_currentGenEvent.GetNbOfHits();
  G4int nRecHits = m_currentRecEvent.GetNbOfHits();
  assert(nHits == nRecHits);

  for( G4int i = 0 ; i < nHits ; i++ ) {
    G4int iModule = m_currentRecEvent.GetModuleID(i);
    G4double angle = det->GetModuleAngle(iModule);
    G4double s = sin(angle);
    G4double c = cos(angle);
    G4double dx = m_smearedHits[i].x() - m_currentRecEvent.GetHitPosition(i).x();
    G4double dy = m_smearedHits[i].y() - m_currentRecEvent.GetHitPosition(i).y();

    chi2 += pow(dx*m_sigmaX*s, 2.);
    chi2 += pow(dy*m_sigmaY*s, 2.);
    chi2 += pow(dx*m_sigmaY*c, 2.);
    chi2 += pow(dy*m_sigmaX*c, 2.);
    chi2 += 2.*dx*dy*s*c*(pow(m_sigmaY,2.) - pow(m_sigmaX, 2.));

    chi2 /= pow(m_sigmaX*m_sigmaY, 2.);
  }
    
  return chi2;
}
