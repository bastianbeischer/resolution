#include <cmath>
#include <fstream>

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
#include "TMatrixD.h"
#include "TVectorD.h"

RES_TrackFitter* RES_TrackFitter::m_instance = 0;

void MinuitChi2Wrapper(int& npar, double* /*gin*/, double& f, double* par, int /*iflag*/)
{
  RES_TrackFitter* fitter = RES_TrackFitter::GetInstance();

  for (int i = 0; i < npar; i++)
    fitter->m_parameter[i] = par[i];

  double chi2 = fitter->Chi2InModuleFrame();

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
  //SetStartParametesToGeneratedParticle();

  G4int conv = -1;

  switch (m_fitMethod) {
  case blobel:
    conv = DoBlobelFit(5);
    break;
  case minuit:
    conv = DoMinuitFit(5);
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

  // CLHEP::HepRandom::setTheSeed(1234);
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

    if (m_verbose > 0)
      G4cout << "original hit: " << i << " --> " << hit << G4endl;

    hit = forwardRotation*hit;
    //    hit.setX(CLHEP::RandGauss::shoot(0., m_sigmaU));
    hit.setX(0.);
    hit.setY(CLHEP::RandGauss::shoot(hit.y(), m_sigmaV));
    hit.setZ(CLHEP::RandGauss::shoot(hit.z(), m_sigmaZ));
    hit = backwardRotation*hit;

    m_smearedHits[i] = hit;
  }
  if (m_verbose >0)
    for (int i = 0; i < nHits; i++)
      G4cout << "smeared hit: " << i << " --> " << m_smearedHits[i] << G4endl;

}

void RES_TrackFitter::SetStartParametesToGeneratedParticle()
{
  G4double p  = m_currentGenEvent.GetMomentum();
  G4double x0 = m_currentGenEvent.GetHitPosition(0).x();
  G4double x1 = m_currentGenEvent.GetHitPosition(1).x();
  G4double y0 = m_currentGenEvent.GetHitPosition(0).y();
  G4double y1 = m_currentGenEvent.GetHitPosition(1).y();
  G4double z0 = m_currentGenEvent.GetHitPosition(0).z();
  G4double z1 = m_currentGenEvent.GetHitPosition(1).z();

  m_parameter[0] = p;
  m_parameter[1] = y0;
  m_parameter[2] = atan((y1-y0)/(z1-z0));
  m_parameter[3] = x0;
  m_parameter[4] = atan((x1-x0)/(z1-z0));

  for (int i = 0; i < 5; i++)
    m_step[i] = 0.1*m_parameter[i];

}

void RES_TrackFitter::CalculateStartParameters()
{
  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();

  G4int nHits = m_currentGenEvent.GetNbOfHits();

  G4double dx_over_dz, dy_over_dz;
  G4double* x = new G4double[nHits];
  G4double* y = new G4double[nHits];
  G4double* z = new G4double[nHits];
  G4double* k = new G4double[nHits];
  G4double* l = new G4double[nHits];
  G4double** f = new G4double*[nHits];
  for (int i = 0; i < nHits; i++)
    f[i] = new G4double[3];
  G4double** alpha = new G4double*[nHits];
  for (int i = 0; i < nHits; i++)
    alpha[i] = new G4double[2];

  for (int i = 0; i < nHits; i++) {
    G4int iModule = m_currentGenEvent.GetModuleID(i);
    G4double angle = det->GetModuleAngle(iModule);
    alpha[i][0] = cos(angle);
    alpha[i][1] = sin(angle);
    f[i][0] = m_smearedHits[i].x();
    f[i][1] = m_smearedHits[i].y();
    f[i][2] = m_smearedHits[i].z();
    k[i] = f[i][2] - f[0][2];
  }

  int nRow = 2*nHits;
  int nCol = 4+nHits;

  TMatrixD A(nRow,nCol);
  for (int i = 0; i < nRow; i++)
    for (int j = 0; j < nCol; j++)
      A(i,j) = 0.;

  for (int i = 0; i < nRow; i++) {
    A(i,i%2)   = 1.;
    A(i,i%2+2) = k[i/2];
    A(i,4+i/2) = -alpha[i/2][i%2];
  }

  TVectorD b(nRow);
  for (int i = 0; i < nRow; i++)
    b(i) = f[i/2][i%2];

  TMatrixD U(nRow,nRow);
  for (int i = 0; i < nRow; i++)
    for (int j = 0; j < nRow; j++)
      U(i,j) = 0.;

  for (int i = 0; i < nHits; i++) {
    G4int iModule = m_currentGenEvent.GetModuleID(i);
    G4double angle = det->GetModuleAngle(iModule);

    TMatrixD Rot(2,2); 
    Rot(0,0) = cos(angle);
    Rot(0,1) = -sin(angle);
    Rot(1,0) = sin(angle);
    Rot(1,1) = cos(angle);
    TMatrixD V_prime(2,2);
    V_prime(0,0) = m_sigmaU;
    V_prime(0,1) = 0.;
    V_prime(1,0) = 0.;
    V_prime(1,1) = m_sigmaV;
    TMatrixD V(2,2);
    V = Rot * V_prime * Rot.T();
    
    U(2*i,2*i)     = V(0,0);
    U(2*i,2*i+1)   = V(0,1);
    U(2*i+1,2*i)   = V(1,0);
    U(2*i+1,2*i+1) = V(1,1);
  }

  U.Invert();
  TMatrixD ATranspose(nCol,nRow);
  ATranspose.Transpose(A);
  TMatrixD M = ATranspose * U * A;
  TVectorD c = ATranspose * U * b;
  M.Invert();
  TVectorD solution(nCol);
  solution = M * c;

  dx_over_dz = solution(2);
  dy_over_dz = solution(3);
  for (int i = 0; i < nHits; i++) {
    l[i] = solution(i+4);
    x[i] = f[i][0] + l[i] * alpha[i][0];
    y[i] = f[i][1] + l[i] * alpha[i][1];
    z[i] = f[i][2];
  }

  if (m_verbose > 0) {
    G4cout << "straight line fit:" << G4endl;
    for (int i = 0; i < nHits; i++) {
      G4cout << "i: " << i << " x: " << x[i] << " y: " << y[i] << " z: " << z[i] << G4endl;
    }
  }

  G4double phi = atan(dy_over_dz);
  G4double theta = atan(-dx_over_dz*cos(phi));

  // restore u components of hits
  for (G4int i = 0; i < nHits; i++) {
    G4int iModule = m_currentGenEvent.GetModuleID(i);
    G4double angle = det->GetModuleAngle(iModule);

    G4RotationMatrix forwardRotation(angle, 0., 0.);
    G4RotationMatrix backwardRotation(-angle, 0., 0.);

    // G4double dz = m_smearedHits[i].z() - m_smearedHits[0].z();
    // G4ThreeVector start(x[0],y[0],z[0]);
    // G4ThreeVector direction(sin(theta), -cos(theta)*sin(phi), -cos(theta)*cos(phi));
    // G4double l = dz/direction.z();
    // G4ThreeVector straightLine = start + l*direction;
    G4ThreeVector straightLine(x[i],y[i],z[i]);

    m_smearedHits[i] = forwardRotation*m_smearedHits[i];
    straightLine     = forwardRotation*straightLine;

    m_smearedHits[i].setX(straightLine.x());

    m_smearedHits[i] = backwardRotation*m_smearedHits[i];

    if (m_verbose > 0)
      G4cout << "restored hit: " << i << " --> " << m_smearedHits[i] << G4endl;
  }

  G4double deltaTheta = fabs(  (m_smearedHits[7].y()-m_smearedHits[4].y())/(m_smearedHits[7].z()-m_smearedHits[4].z())
                             - (m_smearedHits[3].y()-m_smearedHits[0].y())/(m_smearedHits[3].z()-m_smearedHits[0].z()));
  G4double B = 0.3;
  G4double L = sqrt(pow(m_smearedHits[4].y()-m_smearedHits[3].y(),2.) + pow(m_smearedHits[4].z()-m_smearedHits[3].z(),2.))/m;
  G4double p = 0.3*B*L/deltaTheta*GeV;
  //  G4double p = m_currentGenEvent.GetMomentum();

  m_parameter[0] = p;
  m_parameter[1] = y[0];
  m_parameter[2] = phi;
  m_parameter[3] = x[0];
  m_parameter[4] = theta;

  G4double sigmaEllipsis = 5.*mm;
  G4double sigmaPhi   = sqrt(2) * m_sigmaV      / ((z[1]-z[0])*(1 + pow((y[1]-y[0])/(z[1]-z[0]),2.0)));
  G4double sigmaTheta = sqrt(2) * sigmaEllipsis / ((z[1]-z[0])*(1 + pow((x[1]-x[0])/(z[1]-z[0]),2.0)));

  m_step[0] = 0.2*p;
  m_step[1] = m_sigmaV;
  m_step[2] = sigmaPhi;
  m_step[3] = sigmaEllipsis;
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

  delete[] x;
  delete[] y;
  delete[] z;
  delete[] k;
  delete[] l;
  for(int i = 0; i < nHits; i++)
    delete[] f[i];
  delete[] f;
  for(int i = 0; i < nHits; i++)
    delete[] alpha[i];
  delete[] alpha;
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
    chi2 = Chi2InModuleFrame();
    DVALLEY(chi2,m_parameter,conv);
  }

  G4int nHits = m_currentGenEvent.GetNbOfHits();
  G4double dof = nHits - npar;
  m_currentRecEvent.SetChi2(chi2);
  m_currentRecEvent.SetDof(dof);

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
  // gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  // // Print results
  // Double_t amin,edm,errdef;
  // Int_t nvpar,nparx,icstat;
  // gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  // gMinuit->mnprin(3,amin);

  // in this step a new rec event will be created with (hopefully) identical properties than the result of the minimization. this should be done properly in the future...
  G4double chi2 = Chi2InModuleFrame();
  G4int nHits = m_currentRecEvent.GetNbOfHits();
  G4double dof = nHits - npar;
  m_currentRecEvent.SetChi2(chi2);
  m_currentRecEvent.SetDof(dof);

  return 1;
}

G4double RES_TrackFitter::Chi2InDetFrame()
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
  
  G4ThreeVector direction(sin(theta), -cos(theta)*sin(phi), -cos(theta)*cos(phi));
  G4ThreeVector position(x0, y0, m_smearedHits[0].z());
  position -= 1.*cm * direction;
      
  gun->SetParticlePosition(position);
  gun->SetParticleMomentumDirection(direction);
  gun->SetParticleEnergy(pt/cos(theta));
  runManager->BeamOn(1);

  G4int nHits = m_currentGenEvent.GetNbOfHits();
  G4int nRecHits = m_currentRecEvent.GetNbOfHits();
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

G4double RES_TrackFitter::Chi2InModuleFrame()
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
  gun->SetParticleEnergy(pt/cos(theta));
  runManager->BeamOn(1);

  G4int nHits = m_currentGenEvent.GetNbOfHits();
  G4int nRecHits = m_currentRecEvent.GetNbOfHits();
  if (nHits != nRecHits) {
    chi2 = DBL_MAX;
    return chi2;
  }

  for( G4int i = 0 ; i < nHits ; i++ ) {
    G4int iModule = m_currentRecEvent.GetModuleID(i);
    G4double angle = det->GetModuleAngle(iModule);
    G4RotationMatrix forwardRotation(angle, 0, 0);
    G4RotationMatrix backwardRotation(-angle, 0, 0);
    G4ThreeVector hit(m_currentRecEvent.GetHitPosition(i).x(),m_currentRecEvent.GetHitPosition(i).y(),m_currentRecEvent.GetHitPosition(i).z());

    hit = forwardRotation*hit;
    m_smearedHits[i] = forwardRotation*m_smearedHits[i];
    
    G4double dv = m_smearedHits[i].y() - hit.y();
    chi2 += pow(dv/m_sigmaV,2.);

    m_smearedHits[i] = backwardRotation*m_smearedHits[i];
  }
    
  return chi2;
}

void RES_TrackFitter::ScanChi2Function(G4int i, G4int j, G4String filename)
{
  SetSpatialResolutions();
  SmearHits();
  //CalculateStartParameters();
  SetStartParametesToGeneratedParticle();
  
  // hardcoded at the moment
  m_lowerBound[0] = 0.2*GeV;        m_upperBound[0] = 1.8*GeV;
  m_lowerBound[1] = 1.5*cm;         m_upperBound[1] = 2.5*cm;
  m_lowerBound[2] = -1.*M_PI/180.;  m_upperBound[2] = 1.*M_PI/180.;
  m_lowerBound[3] = 0.0*cm;         m_upperBound[3] = 4.0*cm;
  m_lowerBound[4] = -10.*M_PI/180.; m_upperBound[4] = 10.*M_PI/180.;

  G4double nSteps = 200.;
  for (int k = 0; k < 5; k++)
    m_step[k] = (m_upperBound[k] - m_lowerBound[k]) / nSteps;

  if (m_verbose > 0) {
    G4cout << "---------------------------------------------" << G4endl;
    G4cout << " Scanning chi2 around the following fixed values:" << G4endl;
    for (int k = 0; k < 5; k++) {
      if (k != i && k != j)
        G4cout << "  m_parameter["<<k<<"] = " << m_parameter[k] << G4endl;
    }
    G4cout << " Varying parameters " << i << " and " << j << ":" << G4endl;    
    G4cout << "  " << i <<" --> from " << m_lowerBound[i] << " to " << m_upperBound[i] << " with step size: " << m_step[i] << G4endl;
    G4cout << "  " << j <<" --> from " << m_lowerBound[j] << " to " << m_upperBound[j] << " with step size: " << m_step[j] << G4endl;
  }

  std::ofstream file(filename);
  if (!file.is_open()) {
    G4cerr << "Error opening file" << G4endl;
    return;
  }

  file << i << "\t" << j << "\t" << nSteps << "\t" << m_lowerBound[i] << "\t" << m_upperBound[i] << "\t" << m_lowerBound[j] << "\t" << m_upperBound[j] << std::endl;

  m_parameter[i] = m_lowerBound[i];
  for (int ctr1 = 0; ctr1 <= nSteps; ctr1++) {
    m_parameter[j] = m_lowerBound[j];
    for (int ctr2 = 0; ctr2 <= nSteps; ctr2++) {
      G4double chi2 = Chi2InModuleFrame();
      file << m_parameter[i] << "\t" << m_parameter[j] << "\t" << chi2 << std::endl;;
      m_parameter[j] += m_step[j];
    }
    m_parameter[i] += m_step[i];
  }

  file.close();
}
