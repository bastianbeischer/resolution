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
  G4int nHits = m_currentGenEvent.GetNbOfHits();
  if (nHits < 8) return m_currentRecEvent;

  SetSpatialResolutions();
  SmearHits();
  CalculateStartParameters();
  //SetStartParametesToGeneratedParticle();

  G4int conv = -1;

  // G4double chi2 = Chi2InModuleFrame();
  // G4int dof = nHits - 4;
  // m_currentRecEvent.SetChi2(chi2);
  // m_currentRecEvent.SetDof(dof);

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

  return m_currentRecEvent;

  // if (conv) return m_currentRecEvent;
  // else      return RES_Event();
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
    G4int iFiber  = m_currentGenEvent.GetFiberID(i);
    G4double angle = det->GetModuleAngle(iModule);
    if (iFiber > 0) angle += det->GetModuleInternalAngle(iModule);
    
    // collect hit information
    G4double x = m_currentGenEvent.GetHitPosition(i).x();
    G4double y = m_currentGenEvent.GetHitPosition(i).y();
    G4double z = m_currentGenEvent.GetHitPosition(i).z();
    G4ThreeVector hit(x,y,z);

    // transform from detector reference frame (x,y,z) to module frame (u,v,z), smear the hits, transform back to detector frame
    G4RotationMatrix forwardRotation(angle, 0., 0.);
    G4RotationMatrix backwardRotation(-angle, 0., 0.);

    if (m_verbose > 1)
      G4cout << "original hit: " << i << " --> " << hit << G4endl;

    hit = forwardRotation*hit;
    //    hit.setX(CLHEP::RandGauss::shoot(0., m_sigmaU));
    hit.setX(0.);
    hit.setY(CLHEP::RandGauss::shoot(hit.y(), m_sigmaV));
    hit.setZ(CLHEP::RandGauss::shoot(hit.z(), m_sigmaZ));
    hit = backwardRotation*hit;

    m_smearedHits[i] = hit;
  }
  if (m_verbose > 1)
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

  // there is something wrong here: it should be pt instead of p (left unfixed because this method is not used at the moment)
  m_parameter[0] = 1./p;
  m_parameter[1] = y0;
  m_parameter[2] = atan((y1-y0)/(z1-z0));
  m_parameter[3] = x0;
  m_parameter[4] = atan((x1-x0)/(z1-z0));

  for (int i = 0; i < 5; i++)
    m_step[i] = 0.1*m_parameter[i];

}

void RES_TrackFitter::FitStraightLine(G4int n0, G4int n1, G4double &x0, G4double &y0, G4double &dxdz, G4double &dydz)
{
  G4int nHits = n1 - n0;
  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  
  G4double z0 = 0.;
  G4double* k = new G4double[nHits];
  G4double** f = new G4double*[nHits];
  for (int i = 0; i < nHits; i++)
    f[i] = new G4double[3];
  G4double** alpha = new G4double*[nHits];
  for (int i = 0; i < nHits; i++)
    alpha[i] = new G4double[2];

  for (int i = 0; i < nHits; i++) {
    G4int iModule = m_currentGenEvent.GetModuleID(i);
    G4int iFiber  = m_currentGenEvent.GetFiberID(i);
    G4double angle = det->GetModuleAngle(iModule);
    if (iFiber > 0) angle += det->GetModuleInternalAngle(iModule);
    alpha[i][0] = cos(angle);
    alpha[i][1] = sin(angle);
    f[i][0] = m_smearedHits[i+n0].x();
    f[i][1] = m_smearedHits[i+n0].y();
    f[i][2] = m_smearedHits[i+n0].z();
    k[i]    = f[i][2] - z0;
  }

  int nRow = nHits;
  int nCol = 4;

  TMatrixD A(nRow,nCol);
  for (int i = 0; i < nRow; i++)
    for (int j = 0; j < nCol; j++)
      A(i,j) = 0.;

  for (int i = 0; i < nRow; i++) {
    A(i,0) = 1.;
    A(i,1) = -alpha[i][0]/alpha[i][1];
    A(i,2) = k[i];
    A(i,3) = -k[i]*alpha[i][0]/alpha[i][1];
  }

  TVectorD b(nRow);
  for (int i = 0; i < nRow; i++)
    b(i) = f[i][0] - f[i][1]*alpha[i][0]/alpha[i][1];

  TMatrixD U(nRow,nRow);
  for (int i = 0; i < nRow; i++)
    for (int j = 0; j < nRow; j++)
      U(i,j) = 0.;

  for (int i = 0; i < nHits; i++) {
    G4int iModule = m_currentGenEvent.GetModuleID(i);
    G4int iFiber  = m_currentGenEvent.GetFiberID(i);
    G4double angle = det->GetModuleAngle(iModule);
    if (iFiber > 0) angle += det->GetModuleInternalAngle(iModule);

    // Rot is the matrix that maps u,v, to x,y (i.e. the backward rotation)
    TMatrixD Rot(2,2); 
    Rot(0,0) = cos(angle);
    Rot(0,1) = sin(angle);
    Rot(1,0) = -sin(angle);
    Rot(1,1) = cos(angle);
    TMatrixD V1(2,2);
    V1(0,0) = 0.;
    V1(0,1) = 0.;
    V1(1,0) = 0.;
    V1(1,1) = m_sigmaV*m_sigmaV;
    TMatrixD RotTrans(2,2);
    RotTrans.Transpose(Rot);
    TMatrixD V2(2,2);
    V2 = Rot * V1 * RotTrans;
    
    TMatrixD Lin(1,2);
    Lin(0,0) = 1.;
    Lin(0,1) = -alpha[i][0]/alpha[i][1];
    TMatrixD LinTrans(2,1);
    LinTrans.Transpose(Lin);
    TMatrixD V3 = TMatrixD(1,1);
    V3 = Lin * V2 * LinTrans;
    
    U(i,i) = V3(0,0);
  }

  U.Invert();
  TMatrixD Uinv = U;
  TMatrixD ATranspose(nCol,nRow);
  ATranspose.Transpose(A);
  TMatrixD M = ATranspose * Uinv * A;
  TVectorD c = ATranspose * Uinv * b;
  M.Invert();
  TMatrixD Minv = M;
  TVectorD solution(nCol);
  solution = Minv * c;

  x0   = solution(0);
  y0   = solution(1);
  dxdz = solution(2);
  dydz = solution(3);

  if (m_verbose > 1) {
    G4cout << "covariance matrix for this fit:" << G4endl;
    TMatrixD Lin(2*nRow,nCol);
    for (int i = 0; i < 2*nRow; i++)
      for (int j = 0; j < nCol; j++)
        Lin(i,j) = 0.;
    for (int i = 0; i < 2*nRow; i++){
      Lin(i,i%2)     = 1.;
      Lin(i,(i%2)+2) = k[i/2];
    }
    TMatrixD LinTrans(nCol, 2*nRow);
    LinTrans.Transpose(Lin);
    TMatrixD Cov(2*nRow, 2*nRow);
    Cov = Lin * Minv * LinTrans;

    for (int i = 0; i < nHits; i++) {
      G4cout << "resolution in x" << i << " --> " << sqrt(Cov(2*i,2*i)) << " mm" << G4endl;
      G4cout << "resolution in y" << i << " --> " << sqrt(Cov(2*i+1,2*i+1)) << " mm" << G4endl;
    }

    if (m_verbose > 2)
      Cov.Print();
    

  }

  delete[] k;
  for(int i = 0; i < nHits; i++)
    delete[] f[i];
  delete[] f;
  for(int i = 0; i < nHits; i++)
    delete[] alpha[i];
  delete[] alpha;
}

void RES_TrackFitter::CalculateStartParameters()
{
  G4int nHits = m_currentGenEvent.GetNbOfHits();

  G4double z0 = 0.;

  G4double* k = new G4double[nHits];
  G4double* x = new G4double[nHits];
  G4double* y = new G4double[nHits];
  G4double* z = new G4double[nHits];

  for (int i = 0; i < nHits; i++)
    k[i] = m_smearedHits[i].z() - z0;

  G4double x0_top,y0_top,dx_over_dz_top,dy_over_dz_top;
  FitStraightLine(0, nHits/2, x0_top, y0_top, dx_over_dz_top, dy_over_dz_top);
  G4double x0_bottom,y0_bottom,dx_over_dz_bottom,dy_over_dz_bottom;
  FitStraightLine(nHits/2, nHits, x0_bottom, y0_bottom, dx_over_dz_bottom, dy_over_dz_bottom);

  for (int i = 0; i < nHits; i++) {
    if (i < nHits/2) {
      x[i] = x0_top + k[i] * dx_over_dz_top;
      y[i] = y0_top + k[i] * dy_over_dz_top;
      z[i] = z0 + k[i];
    }
    else {
      x[i] = x0_bottom + k[i] * dx_over_dz_bottom;
      y[i] = y0_bottom + k[i] * dy_over_dz_bottom;
      z[i] = z0 + k[i];
    }
  }
  G4double phi = atan(dy_over_dz_top);
  G4double theta = atan(-dx_over_dz_top*cos(phi));
 

  // G4double x0,y0,dx_over_dz,dy_over_dz;
  // FitStraightLine(0, nHits, x0, y0, dx_over_dz, dy_over_dz);

  // for (int i = 0; i < nHits; i++) {
  //   x[i] = x0 + k[i] * dx_over_dz;
  //   y[i] = y0 + k[i] * dy_over_dz;
  //   z[i] = z0 + k[i];
  //   // x[i] += xCorrection[i];
  //   // y[i] += yCorrection[i];
  // }
  // G4double phi = atan(dy_over_dz);
  // G4double theta = atan(-dx_over_dz*cos(phi));

  // G4double dy_over_dz_top = dy_over_dz;
  // G4double dy_over_dz_bottom = dy_over_dz;
 

  if (m_verbose > 1) {
    G4cout << "straight line fit:" << G4endl;
    for (int i = 0; i < nHits; i++) {
      G4cout << "i: " << i << " x: " << x[i] << " y: " << y[i] << " z: " << z[i] << G4endl;
    }
  }

  //restore u components of hits
  // for (G4int i = 0; i < nHits; i++) {
  // RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    // G4int iModule = m_currentGenEvent.GetModuleID(i);
    // G4int iFiber  = m_currentGenEvent.GetFiberID(i);
    // G4double angle = det->GetModuleAngle(iModule);
    // if (iFiber > 0) angle += det->GetModuleInternalAngle(iModule);

  //   G4RotationMatrix forwardRotation(angle, 0., 0.);
  //   G4RotationMatrix backwardRotation(-angle, 0., 0.);
  //   G4ThreeVector straightLine(x[i],y[i],z[i]);

  //   m_smearedHits[i] = forwardRotation*m_smearedHits[i];
  //   straightLine     = forwardRotation*straightLine;

  //   // m_smearedHits[i].setX(straightLine.x());

  //   m_smearedHits[i] = backwardRotation*m_smearedHits[i];

  //   if (m_verbose > 1)
  //     G4cout << "restored hit: " << i << " --> " << m_smearedHits[i] << G4endl;
  // }

  G4double deltaTheta = fabs(dy_over_dz_bottom - dy_over_dz_top);

  G4double B = 0.27;
  G4double L = sqrt(pow(y[4]-y[3],2.) + pow(z[4]-z[3],2.))/m;
  G4double p = 0.3*B*L/deltaTheta*GeV;

  m_parameter[0] = 1./p;
  m_parameter[1] = y[0];
  m_parameter[2] = phi;
  m_parameter[3] = x[0];
  m_parameter[4] = theta;

  G4double sigmaEllipsis = 1.*mm;
  G4double sigmaPhi   = sqrt(2) * m_sigmaV      / ((z[1]-z[0])*(1 + pow((y[1]-y[0])/(z[1]-z[0]),2.0)));
  G4double sigmaTheta = sqrt(2) * sigmaEllipsis / ((z[1]-z[0])*(1 + pow((x[1]-x[0])/(z[1]-z[0]),2.0)));

  m_step[0] = 0.1*m_parameter[0];
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
  // npar has been put to -npar before
  G4double dof = nHits - (-npar);
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
  
  G4String names[5] = {"1/pt", "y0", "phi", "x0", "theta"};
  for (int i = 0; i < npar; i++)
    gMinuit->mnparm(i, names[i].c_str(), m_parameter[i], m_step[i], m_lowerBound[i], m_upperBound[i], ierflg);

  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1.;
  //gMinuit->mnexcm("SIMPLEX", arglist, 2, ierflg);
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  if (m_verbose > 0) {
    // Print results
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    gMinuit->mnprin(3,amin);
  }

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
    
  G4double pt    = 1./m_parameter[0];
  G4double y0    = m_parameter[1];
  G4double phi   = m_parameter[2];
  G4double x0    = m_parameter[3];
  G4double theta = m_parameter[4];
  
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
    G4int iModule = m_currentGenEvent.GetModuleID(i);
    G4int iFiber  = m_currentGenEvent.GetFiberID(i);
    G4double angle = det->GetModuleAngle(iModule);
    if (iFiber > 0) angle += det->GetModuleInternalAngle(iModule);

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
    
  G4double pt    = 1./m_parameter[0];
  G4double y0    = m_parameter[1];
  G4double phi   = m_parameter[2];
  G4double x0    = m_parameter[3];
  G4double theta = m_parameter[4];

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
    G4int iModule = m_currentGenEvent.GetModuleID(i);
    G4int iFiber  = m_currentGenEvent.GetFiberID(i);
    G4double angle = det->GetModuleAngle(iModule);
    if (iFiber > 0) angle += det->GetModuleInternalAngle(iModule);
    G4RotationMatrix forwardRotation(angle, 0, 0);
    G4RotationMatrix backwardRotation(-angle, 0, 0);
    G4ThreeVector hit(m_currentRecEvent.GetHitPosition(i).x(),m_currentRecEvent.GetHitPosition(i).y(),m_currentRecEvent.GetHitPosition(i).z());

    hit = forwardRotation*hit;
    m_smearedHits[i] = forwardRotation*m_smearedHits[i];
    
    G4double du = m_smearedHits[i].x() - hit.x();
    G4double dv = m_smearedHits[i].y() - hit.y();
    //chi2 += pow(du/m_sigmaU,2.);
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
