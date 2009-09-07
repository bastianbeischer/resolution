#include <cmath>

#include "RES_TrackFitter.hh"

#include "RES_TrackFitMessenger.hh"
#include "RES_Event.hh"
#include "RES_RunManager.hh"
#include "RES_DetectorConstruction.hh"
#include "RES_PrimaryGeneratorAction.hh"
#include "RES_FiberHit.hh"
#include "blobel.h"

#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

#include "TMinuit.h"

RES_TrackFitter* RES_TrackFitter::m_instance = 0;

void MinuitChi2Wrapper(int& npar, double* /*gin*/, double& f, double* par, int /*iflag*/)
{
  RES_TrackFitter* fitter = RES_TrackFitter::GetInstance();

  for (int i = 0; i < npar; i++)
    fitter->m_parameter[i] = par[i];

  double chi2 = fitter->Chi2();

  f = chi2;
}

RES_TrackFitter::RES_TrackFitter() :
  m_verbose(0),
  m_smearedHits(0),
  m_fitMethod(blobel),
  m_sigmaU(0.),
  m_sigmaV(0.),
  m_sigmaZ(0.)
{
  m_messenger = new RES_TrackFitMessenger(this);
  m_parameter = new double[5];
  m_step = new double[5];
  m_lowerBound = new double[5];
  m_upperBound = new double[5];
}

RES_TrackFitter::~RES_TrackFitter()
{
  delete[] m_smearedHits;
  delete[] m_parameter;
  delete[] m_step;
  delete[] m_lowerBound;
  delete[] m_upperBound;
  delete m_messenger;
}

RES_TrackFitter* RES_TrackFitter::GetInstance()
{
  if (!m_instance) m_instance = new RES_TrackFitter();
  return m_instance;
}

RES_Event RES_TrackFitter::Fit()
{
  SetSpatialResolutions();
  SmearHits();
  CalculateStartParameters();

  G4int conv = -1;

  switch (m_fitMethod) {
  case blobel:
    conv = DoBlobelFit(3);
    break;
  case minuit:
    conv = DoMinuitFit(3);
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
  m_sigmaU = moduleLength/sqrt(12.);
  m_sigmaV = 50.*um;
  m_sigmaZ = 0.*um;
}

void RES_TrackFitter::SmearHits()
{
  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();

  G4int nHits = m_currentGenEvent.GetNbOfHits();

  delete[] m_smearedHits;
  m_smearedHits = new G4ThreeVector[nHits];

  for (G4int i = 0; i < nHits; i++) {
    G4int iModule = m_currentGenEvent.GetModuleID(i);
    G4double angle = det->GetModuleAngle(iModule);
    
    // collect hit information
    G4double x = m_currentGenEvent.GetHitPosition(i).x();
    G4double y = m_currentGenEvent.GetHitPosition(i).y();
    G4double z = m_currentGenEvent.GetHitPosition(i).z();
    G4ThreeVector hit(x,y,z);

    // transform from detector reference frame (x,y,z) to module frame (u,v,z), smear the hits, transform back to detector frame
    G4RotationMatrix forwardRotation(angle, 0., 0.);
    G4RotationMatrix backwardRotation(-angle, 0., 0.);
    hit = forwardRotation*hit;
    //hit.setX(CLHEP::RandGauss::shoot(hit.x(), m_sigmaU));
    hit.setX(CLHEP::RandGauss::shoot(hit.x(), 0.));
    hit.setY(CLHEP::RandGauss::shoot(hit.y(), m_sigmaV));
    hit.setZ(CLHEP::RandGauss::shoot(hit.z(), m_sigmaZ));
    hit = backwardRotation*hit;

    m_smearedHits[i] = hit;
  }
}

void RES_TrackFitter::CalculateStartParameters()
{
  // G4int nHits = m_currentGenEvent.GetNbOfHits();
  // G4double theta =   (m_smearedHits[nHits-1].y() - m_smearedHits[nHits-2].y()) / (m_smearedHits[nHits-1].z() - m_smearedHits[nHits-2].z())
  //                  - (m_smearedHits[1].y()       - m_smearedHits[0].y())       / (m_smearedHits[1].z()       - m_smearedHits[0].z());
  // G4double L = sqrt(pow(m_smearedHits[nHits-1].y() - m_smearedHits[0].y(),2.) + pow(m_smearedHits[nHits-1].z() - m_smearedHits[0].z(), 2.))/m;
  G4double theta =   (m_smearedHits[7].y() - m_smearedHits[4].y()) / (m_smearedHits[7].z() - m_smearedHits[4].z())
                   - (m_smearedHits[3].y() - m_smearedHits[0].y()) / (m_smearedHits[3].z() - m_smearedHits[0].z());
  G4double L = sqrt(pow(m_smearedHits[4].y() - m_smearedHits[3].y(),2.) + pow(m_smearedHits[4].z() - m_smearedHits[3].z(), 2.))/m;
  G4double B = 0.3;
  G4double pStart = 0.3 * B * L / theta * GeV;
  //  G4double pStart = m_currentGenEvent.GetMomentum();
  
  G4double x0 = m_smearedHits[0].x();
  G4double x1 = m_smearedHits[1].x();
  G4double y0 = m_smearedHits[0].y();
  G4double y1 = m_smearedHits[1].y();
  G4double z0 = m_smearedHits[0].z();
  G4double z1 = m_smearedHits[1].z();

  m_parameter[0] = pStart;
  m_parameter[1] = y0;
  m_parameter[2] = atan((y1-y0)/(z1-z0));
  m_parameter[3] = x0;
  m_parameter[4] = atan((x1-x0)/(z1-z0));

  G4double sigmaPhi = sqrt(2) * m_sigmaV / ((z1-z0)*(1 + pow((y1-y0)/(z1-z0),2.0)));
  G4double sigmaTheta = sqrt(2) * m_sigmaU / ((z1-z0)*(1 + pow((x1-x0)/(z1-z0),2.0)));

  // for (int i = 0; i < 5; i++)
  //   m_step[i] = 0.1*m_parameter[i];
  m_step[0] = 0.1*m_parameter[0];
  m_step[1] = m_sigmaV;
  m_step[2] = sigmaPhi;
  m_step[3] = m_sigmaU;
  m_step[4] = sigmaTheta;

  for (int i = 0; i < 5; i++) {
    m_lowerBound[i] = m_parameter[i] - 5.*fabs(m_step[i]);
    m_upperBound[i] = m_parameter[i] + 5.*fabs(m_step[i]);
  }

  if (m_verbose > 0) {
    G4cout << "---------------------------------------------" << G4endl;
    G4cout << "Start values: " << G4endl;
    for (int i = 0; i < 5; i++)
      G4cout << "m_parameter["<<i<<"] = " << m_parameter[i] << ", m_step["<<i<<"] = " << m_step[i] << ", m_lowerBound["<<i<<"] = " << m_lowerBound[i] << ", m_upperBound["<<i<<"] = " << m_upperBound[i] << G4endl;
  }

}

G4int RES_TrackFitter::DoBlobelFit(G4int npar)
{
  G4int nflim = 3*abs(npar)*(abs(npar)+10); // dito
  if( m_verbose == 0 ) {
    npar = -npar;
    nflim = -nflim;
  }
  DVALLIN(npar,m_step,nflim);

  G4int conv = -1;
  G4int iter = 0;
  G4double chi2 = 0.;

  while( conv <= 0 ){
    ++iter;
    chi2 = Chi2();
    DVALLEY(chi2,m_parameter,conv);
  }

  G4int nHits = m_currentGenEvent.GetNbOfHits();
  G4double dof = nHits - npar;
  m_currentRecEvent.SetChi2OverDof(chi2/dof);

  return conv;
}

G4int RES_TrackFitter::DoMinuitFit(G4int npar)
{
  TMinuit* gMinuit = new TMinuit(npar);
  gMinuit->SetFCN(MinuitChi2Wrapper);
  gMinuit->SetPrintLevel(-1);
  
  double arglist[10];
  int ierflg = 0;
  arglist[0] = 1.;
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
  
  G4String names[5] = {"pt", "y0", "phi", "x0", "theta"};
  for (int i = 0; i < npar; i++)
    gMinuit->mnparm(i, names[i].c_str(), m_parameter[i], m_step[i], m_lowerBound[i], m_upperBound[i], ierflg);

  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1.;
  gMinuit->mnexcm("SIMPLEX", arglist, 2, ierflg);
  // gMinuit->mnexcm("MINOS", arglist, 2, ierflg);

  // // Print results
  // Double_t amin,edm,errdef;
  // Int_t nvpar,nparx,icstat;
  // gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  // gMinuit->mnprin(3,amin);

  // in this step a new rec event will be created with (hopefully) identical properties than the result of the minimization. this should be done properly in the future...
  G4double chi2 = Chi2();
  G4int nHits = m_currentGenEvent.GetNbOfHits();
  G4double dof = nHits - npar;
  m_currentRecEvent.SetChi2OverDof(chi2/dof);

  return 1;
}

G4double RES_TrackFitter::Chi2()
{
  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();

  G4RunManager* runManager = G4RunManager::GetRunManager();
  const RES_PrimaryGeneratorAction* genAction = (RES_PrimaryGeneratorAction*) runManager->GetUserPrimaryGeneratorAction();
  G4ParticleGun* gun = genAction->GetParticleGun();

  G4double chi2 = 0.;
    
  double pt    = m_parameter[0];
  double y0    = m_parameter[1];
  double phi   = m_parameter[2];
  double x0    = m_parameter[3];
  double theta = m_parameter[4];
  
  G4ThreeVector direction(sin(M_PI - theta), cos(M_PI - theta)*sin(phi), cos(M_PI - theta)*cos(phi));
  G4ThreeVector position(x0, y0, m_smearedHits[0].z());
  position -= 1.*cm * direction;
      
  gun->SetParticlePosition(position);
  gun->SetParticleMomentumDirection(direction);
  // gun->SetParticlePosition(G4ThreeVector(0,0,500));
  // gun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1));
  gun->SetParticleEnergy(pt/cos(theta));
  runManager->BeamOn(1);

  G4int nHits = m_currentGenEvent.GetNbOfHits();
  G4int nRecHits = m_currentRecEvent.GetNbOfHits();
  //assert(nHits == nRecHits);
  if (nHits != nRecHits) {
    chi2 = DBL_MAX;
    return chi2;
  }

  for( G4int i = 0 ; i < nHits ; i++ ) {
    G4int iModule = m_currentRecEvent.GetModuleID(i);
    G4double angle = det->GetModuleAngle(iModule);
    G4double s = sin(angle);
    G4double c = cos(angle);
    G4double dx = m_smearedHits[i].x() - m_currentRecEvent.GetHitPosition(i).x();
    G4double dy = m_smearedHits[i].y() - m_currentRecEvent.GetHitPosition(i).y();
    chi2 += pow(dx*m_sigmaU*s, 2.);
    chi2 += pow(dy*m_sigmaV*s, 2.);
    chi2 += pow(dx*m_sigmaV*c, 2.);
    chi2 += pow(dy*m_sigmaU*c, 2.);
    chi2 += 2.*dx*dy*s*c*(pow(m_sigmaV,2.) - pow(m_sigmaU, 2.));
  }
  
  chi2 /= pow(m_sigmaU*m_sigmaV, 2.);
    
  return chi2;
}
