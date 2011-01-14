// $Id: RES_TrackFitter.cc,v 1.64.2.1 2010/09/08 17:45:00 beischer Exp $

#include <cmath>
#include <fstream>

#include "RES_TrackFitter.hh"

#include "RES_TrackFitMessenger.hh"
#include "RES_Event.hh"
#include "RES_MagneticField.hh"
#include "RES_RunManager.hh"
#include "RES_DetectorConstruction.hh"
#include "RES_PrimaryGeneratorAction.hh"
#include "blobel.h"

#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
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
  m_fitMethod(blobel)
{
  m_messenger = new RES_TrackFitMessenger(this);
  m_initialParameter = new double[5];
  m_parameter = new double[5];
  m_step = new double[5];
  m_lowerBound = new double[5];
  m_upperBound = new double[5];
}

RES_TrackFitter::~RES_TrackFitter()
{
  delete[] m_smearedHits;
  delete[] m_initialParameter;
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
  // this has to be redone (use efficiencies, don't require all hits, remove the hardcoded 8)...
  G4int nMinimumHits = 8;
  if (m_fitMethod == oneline) nMinimumHits = 4;

  G4int nHits = m_currentGenEvent.GetNbOfHits();
  if (nHits < nMinimumHits) return m_currentRecEvent;

  G4RunManager* runManager = G4RunManager::GetRunManager();
  const RES_PrimaryGeneratorAction* genAction = (RES_PrimaryGeneratorAction*) runManager->GetUserPrimaryGeneratorAction();
  G4ParticleGun* gun = genAction->GetParticleGun();
  m_initialCharge = gun->GetParticleCharge();

  //SetStartParametesToGeneratedParticle();

  switch (m_fitMethod) {
  case oneline:
    double x0,y0,lambda_x,lambda_y;
    FitStraightLine(0,nHits,x0,y0,lambda_x,lambda_y);
    break;
 case blobel:
    CalculateStartParameters();
    DoBlobelFit(5);
    break;
  case minuit:
    CalculateStartParameters();
    DoMinuitFit(5);
    break;
  case transverse:
    CalculateStartParameters();
    DoBlobelFit(3);
    break;
  case fullmatrix:
    DoFullFit();
    break;
  case testbeam:
    // AddLayerToBeSkipped(0);
    // AddLayerToBeSkipped(10);    
    CalculateStartParameters();
    DoBlobelFit(3);
    m_layersToBeSkipped.clear();
    break;
  default:
    break;
  }

  m_currentRecEvent.SetInitialParameters(5, m_initialParameter);
  m_currentRecEvent.SetFinalParameters(5, m_parameter);
    
  if (gun->GetParticleCharge() != m_initialCharge) {
    m_currentRecEvent.SetMomentum(-m_currentRecEvent.GetMomentum());
    m_currentRecEvent.SetTransverseMomentum(-m_currentRecEvent.GetTransverseMomentum());
  }
  gun->SetParticleCharge(m_initialCharge);

  //std::cout << "p: " << m_currentRecEvent.GetMomentum()/GeV << std::endl;

  return m_currentRecEvent;
}

void RES_TrackFitter::CopyHits()
{
  G4int nHits = m_currentGenEvent.GetNbOfHits();

  delete[] m_smearedHits;
  m_smearedHits = new G4ThreeVector[nHits];

  for (int i = 0; i < nHits; i++) {
    m_smearedHits[i] = G4ThreeVector(m_currentGenEvent.GetSmearedHitPosition(i). x(), m_currentGenEvent.GetSmearedHitPosition(i). y(), m_currentGenEvent.GetSmearedHitPosition(i). z());
  }

  if (m_verbose > 2)
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

  double phi = atan((y1-y0)/(z1-z0));
  double theta = atan((x1-x0)/(z1-z0));

  m_initialParameter[0] = 1./(cos(theta)*p);
  m_initialParameter[1] = y0;
  m_initialParameter[2] = phi;
  m_initialParameter[3] = x0;
  m_initialParameter[4] = theta;

  m_initialParameter[0] = 1./GeV;
  m_initialParameter[1] = 0;
  m_initialParameter[2] = 0;
  m_initialParameter[3] = 0;
  m_initialParameter[4] = 0;

  for (int i = 0; i < 5; i++)
    m_parameter[i] = m_initialParameter[i];

  for (int i = 0; i < 5; i++)
    m_step[i] = 0.1*m_parameter[i];
}

void RES_TrackFitter::FitStraightLine(G4int n0, G4int n1, G4double &x0, G4double &y0, G4double &dxdz, G4double &dydz)
{
  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();

  // only fit from n0 to n1! (useful if segments of the track are to be fitted)
  unsigned int nHits = n1 - n0;

  // this parameter is arbitrary. z0 = 0 should minimize correlations...
  float z0 = 0.0;

  unsigned int numberOfLayersToBeSkipped = 0;
  for (unsigned int i = 0 ; i < nHits; i++) {
    G4int iModule = m_currentGenEvent.GetModuleID(i);
    G4int iLayer  = m_currentGenEvent.GetLayerID(i);
    G4int uniqueLayer = 2*iModule + iLayer;
    std::vector<int>::iterator findIt = std::find(m_layersToBeSkipped.begin(), m_layersToBeSkipped.end(), uniqueLayer);
    if (findIt != m_layersToBeSkipped.end())
      numberOfLayersToBeSkipped++;
  }

  // basic dimensions of matrices
  unsigned int nRow = nHits - numberOfLayersToBeSkipped;
  unsigned int nCol = 4;
  unsigned int nModules = det->GetNumberOfModules();

  // declare matrices for the calculation
  TMatrixD A(nRow,nCol);
  TVectorD b(nRow);
  TMatrixD U(nRow,nRow);
  TMatrixD CombineXandY(1,2);
  TMatrixD SolutionToPositions(4*nModules,nCol);

  unsigned int counter = 0;
  for (unsigned int i = 0; i < nHits; i++) {

    G4int iModule = m_currentGenEvent.GetModuleID(i);
    G4int iLayer  = m_currentGenEvent.GetLayerID(i);
    G4int uniqueLayer = 2*iModule + iLayer;

    std::vector<int>::iterator it = std::find(m_layersToBeSkipped.begin(), m_layersToBeSkipped.end(), uniqueLayer);
    if (it != m_layersToBeSkipped.end())
      continue;

    // get information from detector...
    RES_Module* module = det->GetModule(iModule);
    G4double angle = module->GetAngle();
    if (iLayer > 0) angle += module->GetInternalAngle();

    G4ThreeVector pos = m_smearedHits[i+n0];

    // fill the matrices
    G4float k = pos.z() - z0;
    G4bool useTangens = fabs(angle) < M_PI/4. ? true : false;
    if (useTangens) {
      float tangens     = sin(angle)/cos(angle);
      CombineXandY(0,0) = -tangens;
      CombineXandY(0,1) = 1.;
      A(counter,0)      = -tangens;
      A(counter,1)      = 1.;
      A(counter,2)      = -k*tangens;
      A(counter,3)      = k;
      b(counter)        = -tangens*pos.x() + pos.y();
    }
    else {
      float cotangens   = cos(angle)/sin(angle);
      CombineXandY(0,0) = 1.;
      CombineXandY(0,1) = -cotangens;
      A(counter,0)      = 1.;
      A(counter,1)      = -cotangens;
      A(counter,2)      = k;
      A(counter,3)      = -k*cotangens;
      b(counter)        = pos.x() - cotangens*pos.y();
    }

    // calculate covariance matrix
    G4double sigmaV = iLayer==0? module->GetUpperSigmaV() : module->GetLowerSigmaV();

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
    V1(1,1) = sigmaV*sigmaV;
    TMatrixD RotTrans(2,2);
    RotTrans.Transpose(Rot);
    TMatrixD V2(2,2);
    V2 = Rot * V1 * RotTrans;
    TMatrixD CombineXandYTrans(2,1);
    CombineXandYTrans.Transpose(CombineXandY);
    TMatrixD V3 = TMatrixD(1,1);
    V3 = CombineXandY * V2 * CombineXandYTrans;
    U(counter,counter) = V3(0,0); // this is the sigma for the i'th measurement
    counter++;
  }
  
  // calculate solution
  TMatrixD Uinv = U;
  Uinv.Invert();
  TMatrixD ATranspose(nCol,nRow);
  ATranspose.Transpose(A);
  TMatrixD M = ATranspose * Uinv * A;
  TVectorD c = ATranspose * Uinv * b;
  TMatrixD Minv = M;
  Minv.Invert();
  TVectorD solution(nCol);
  solution = Minv * c;

  // calculate chi2 and track positions from fit parameters
  TMatrixD vec(nRow,1);
  for (unsigned int i = 0; i < nRow; i++)
    vec(i,0) = (A*solution - b)(i);
  TMatrixD vecTrans(1,nRow);
  vecTrans.Transpose(vec);
  G4double chi2 = (vecTrans * Uinv * vec)(0,0);

  // SolutionToPositions is the linear transformation that maps the solution to positions
  for (unsigned int i = 0; i < 4*nModules; i++){
    G4double z;
    if ((i%4) < 2) z = det->GetModule(i/4)->GetUpperZ();
    else           z = det->GetModule(i/4)->GetLowerZ();
    SolutionToPositions(i,i%2)     = 1.;
    SolutionToPositions(i,(i%2)+2) = z - z0;
  }
  TVectorD positions = SolutionToPositions*solution;

  // Fill information in m_currentRecEvent (usually done in RES_EventActionReconstruction for other fit methods)
  if (m_fitMethod == oneline) {
    m_currentRecEvent = RES_Event();
    for (unsigned int i = 0; i < 2*nModules; i++) {
      G4int iModule = i/2;
      G4int iLayer  = i%2;
      G4double z;
      if (iLayer == 0) z = det->GetModule(iModule)->GetUpperZ();
      else             z = det->GetModule(iModule)->GetLowerZ();
      m_currentRecEvent.AddHit(iModule, iLayer, positions(2*i), positions(2*i+1), z);
    }
    m_currentRecEvent.SetChi2(chi2);
    m_currentRecEvent.SetDof(nRow - nCol);
    m_currentRecEvent.SetEventType(reconstructed);
    m_currentRecEvent.SetMomentum(DBL_MAX);
    m_currentRecEvent.SetTransverseMomentum(DBL_MAX);
  }

  // print covariance matrix if the user wants to
  if (m_verbose > 1) {
    TMatrixD SolutionToPositionsTrans(nCol, 4*nModules);
    SolutionToPositionsTrans.Transpose(SolutionToPositions);
    TMatrixD Cov(4*nModules, 4*nModules);
    Cov = SolutionToPositions * Minv * SolutionToPositionsTrans;

    for (unsigned int i = 0; i < 2*nModules; i++)
      G4cout << "resolution in x" << i << " --> " << sqrt(Cov(2*i,2*i)) << " mm" << G4endl;
    for (unsigned int i = 0; i < 2*nModules; i++)
      G4cout << "resolution in y" << i << " --> " << sqrt(Cov(2*i+1,2*i+1)) << " mm" << G4endl;

    if (m_verbose > 2) {
      G4cout << "covariance matrix for this fit:" << G4endl;
      Cov.Print();
    }
  }

  // return information from the fit.
  x0   = solution(0);
  y0   = solution(1);
  dxdz = solution(2);
  dydz = solution(3);
}

void RES_TrackFitter::CalculateStartParameters()
{
  G4double lambda_x_top = 0., lambda_x_bottom = 0., lambda_y_top = 0., lambda_y_bottom = 0.;

  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  RES_MagneticField* field = (RES_MagneticField*) fieldMgr->GetDetectorField();
  G4double B_estimate  = field->GetFieldEstimate();
  G4double magnetHeight = fabs(field->GetZ1() - field->GetZ0());
  G4double z0_magnet = -magnetHeight/2. + field->GetDisplacement().z();
  G4double z1_magnet =  magnetHeight/2. + field->GetDisplacement().z();

  if (m_fitMethod == transverse) {

    // hardcoded array bounds (urgh)
    // not too important because fit method is deprecated
    lambda_y_top = (m_smearedHits[3].y() - m_smearedHits[0].y()) / (m_smearedHits[3].z() - m_smearedHits[0].z());
    lambda_y_bottom = (m_smearedHits[7].y() - m_smearedHits[4].y()) / (m_smearedHits[7].z() - m_smearedHits[4].z());

    G4double deltaTheta = lambda_y_top - lambda_y_bottom;

    G4double L  = sqrt(pow(m_smearedHits[4].y()-m_smearedHits[3].y(),2.) + pow(m_smearedHits[4].z()-m_smearedHits[3].z(),2.));
    G4double pt = 0.3*(B_estimate/tesla)*(L/m)/deltaTheta*GeV;

    m_initialParameter[0] = 1./pt;
    m_initialParameter[1] = m_smearedHits[0].y();
    m_initialParameter[2] = atan(lambda_y_top);
    m_initialParameter[3] = 0.;
    m_initialParameter[4] = 0.;

  } // if (method == transverse)

  if (m_fitMethod == blobel || m_fitMethod == minuit || m_fitMethod == oneline || m_fitMethod == twolines) {
    G4int nHits = m_currentGenEvent.GetNbOfHits();

    G4double x0_top,y0_top;
    FitStraightLine(0, nHits/2, x0_top, y0_top, lambda_x_top, lambda_y_top);
    G4double x0_bottom,y0_bottom;
    FitStraightLine(nHits/2, nHits, x0_bottom, y0_bottom, lambda_x_bottom, lambda_y_bottom);

    G4double deltaTheta = atan(lambda_y_top) - atan(lambda_y_bottom);
    G4double y0_magnet = y0_bottom + z0_magnet*lambda_y_bottom;
    G4double y1_magnet = y0_top    + z1_magnet*lambda_y_top;
    G4double L  = sqrt(pow(y1_magnet - y0_magnet, 2.) + pow(z1_magnet - z0_magnet,2.));
    G4double pt = 0.3*(B_estimate/tesla)*(L/m)/deltaTheta*GeV;

    // pt is by definition positive for electrons
    if (m_initialCharge > 0)
      pt = -pt;

    G4double phi = atan(lambda_y_top);
    G4double theta = atan(-lambda_x_top*cos(phi));
    G4double z0 = 0.;
    G4double k0 = m_smearedHits[0].z() - z0;
    G4double x0 = x0_top + k0*lambda_x_top;
    G4double y0 = y0_top + k0*lambda_y_top;

    m_initialParameter[0] = 1./pt;
    m_initialParameter[1] = y0;
    m_initialParameter[2] = phi;
    m_initialParameter[3] = x0;
    m_initialParameter[4] = theta;

  } // if (method == blobel,minuit,oneline,twolines)

  if (m_fitMethod == testbeam) {

    G4int nHits = m_currentGenEvent.GetNbOfHits();
    double x0,y0_top,y0_bottom,lambda_x,lambda_y_top,lambda_y_bottom;

    G4double dummy1,dummy2;

    std::vector<int> alreadyInToBeSkipped;
    for (std::vector<int>::iterator it = m_layersToBeSkipped.begin(); it != m_layersToBeSkipped.end(); it++)
      alreadyInToBeSkipped.push_back(*it);

    for (int i = 6; i < 12; i++)
      if (i != 11)
        AddLayerToBeSkipped(i);
    FitStraightLine(0,nHits,x0,y0_top,lambda_x,lambda_y_top);
    for (int i = 6; i < 12; i++)
      if ((i != 11) && (std::find(alreadyInToBeSkipped.begin(), alreadyInToBeSkipped.end(), i) == alreadyInToBeSkipped.end()))
        RemoveLayerToBeSkipped(i);

    for (int i = 0; i < 6; i++)
      if (i != 1)
        AddLayerToBeSkipped(i);
    FitStraightLine(0,nHits,dummy1,y0_bottom,dummy2,lambda_y_bottom);
    for (int i = 0; i < 6; i++)
      if ((i != 1) && (std::find(alreadyInToBeSkipped.begin(), alreadyInToBeSkipped.end(), i) == alreadyInToBeSkipped.end()))
        RemoveLayerToBeSkipped(i);

    G4double deltaTheta = lambda_y_top - lambda_y_bottom;
    G4double y0_magnet = y0_bottom + z0_magnet*lambda_y_bottom;
    G4double y1_magnet = y0_top    + z1_magnet*lambda_y_top;
    G4double L  = sqrt(pow(y1_magnet - y0_magnet, 2.) + pow(z1_magnet - z0_magnet,2.));
    G4double pt = (0.3*(B_estimate/tesla)*(L/m)/deltaTheta)*GeV;

    G4double phi   = atan(lambda_y_top);
    G4double theta = atan(-lambda_x*cos(phi));
    G4double z0 = 0;
    G4double k0 = m_smearedHits[0].z() - z0;

    m_initialParameter[0] = 1./pt;
    m_initialParameter[1] = y0_top + k0*lambda_y_top;
    m_initialParameter[2] = phi;
    m_initialParameter[3] = x0 + k0*lambda_x;
    m_initialParameter[4] = theta;

  } // if (method == testbeam)
  
  for (int i = 0; i < 5; i++) {
    m_parameter[i] = m_initialParameter[i];
    m_step[i] = 0.1*m_parameter[i];
    m_lowerBound[i] = 0.;
    m_upperBound[i] = 0.;
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
  G4int nflim = 3*abs(npar)*(abs(npar)+10);
  if( m_verbose == 0 ) {
    npar = -npar;
    nflim = -nflim;
  }
  DVALLIN(npar,m_step,nflim);

  G4int nHits = m_currentGenEvent.GetNbOfHits();
  G4int dof =  nHits - abs(npar);
  dof += m_layersToBeSkipped.size();
  m_currentRecEvent.SetDof(dof);

  G4int conv = -1;
  G4int iter = 0;
  G4double chi2 = 0.;
  while( conv <= 0 ) {
    ++iter;
    chi2 = Chi2InModuleFrame();
    DVALLEY(chi2,m_parameter,conv);
    if (chi2 == DBL_MAX)
      break;
  }

  m_currentRecEvent.SetChi2(chi2);

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
  arglist[1] = 0.1;
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
  G4int nHits = m_currentGenEvent.GetNbOfHits();
  G4int dof = nHits - npar;
  dof = dof + m_layersToBeSkipped.size();
  
  m_currentRecEvent.SetChi2(chi2);
  m_currentRecEvent.SetDof(dof);

  return 1;
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

  if (pt < 0) {
    gun->SetParticleCharge(-m_initialCharge);
    pt = -pt;
  }
  else {
    gun->SetParticleCharge(m_initialCharge);
  }

  G4double mass = gun->GetParticleDefinition()->GetPDGMass();
  G4double momentum = pt/cos(theta);
  G4double energy = sqrt( pow(momentum, 2.) + pow(mass, 2.) ) - mass;

  G4ThreeVector direction(sin(theta), -cos(theta)*sin(phi), -cos(theta)*cos(phi));
  G4ThreeVector position(x0, y0, m_smearedHits[0].z());
  position -= 1.*cm * direction;
      
  gun->SetParticlePosition(position);
  gun->SetParticleMomentumDirection(direction);
  gun->SetParticleEnergy(energy);
  runManager->BeamOn(1);

  if (m_currentRecEvent.GetNbOfHits() != 2*det->GetNumberOfLayers())
    return DBL_MAX;

  G4int nHits = m_currentGenEvent.GetNbOfHits();

  for( G4int i = 0 ; i < nHits ; i++ ) {

    G4int iModule = m_currentGenEvent.GetModuleID(i);
    G4int iLayer  = m_currentGenEvent.GetLayerID(i);
    unsigned int uniqueLayer = 2*iModule + iLayer;

    std::vector<int>::iterator it = std::find(m_layersToBeSkipped.begin(), m_layersToBeSkipped.end(), uniqueLayer);
    if (it != m_layersToBeSkipped.end())
      continue;

    RES_Module* module = det->GetModule(iModule);
    G4double angle = module->GetAngle();
    if (iLayer > 0) angle += module->GetInternalAngle();
    G4RotationMatrix forwardRotation(angle, 0, 0);
    G4RotationMatrix backwardRotation(-angle, 0, 0);

    unsigned int j = 0;
    for (; j < m_currentRecEvent.GetNbOfHits(); j++)
      if (m_currentRecEvent.GetHitPosition(j).z() == m_smearedHits[i].z())
        break;
    if (j == m_currentRecEvent.GetNbOfHits()) // this should not be necessary, it is however?!...
      continue;

    G4ThreeVector hit(m_currentRecEvent.GetHitPosition(j).x(),m_currentRecEvent.GetHitPosition(j).y(),m_currentRecEvent.GetHitPosition(j).z());

    hit = forwardRotation*hit;
    m_smearedHits[i] = forwardRotation*m_smearedHits[i];
    
    G4double sigmaU, sigmaV;
    if (iLayer == 0) {
      sigmaU = module->GetUpperSigmaU();
      sigmaV = module->GetUpperSigmaV();
    }
    else {
      sigmaU = module->GetLowerSigmaU();
      sigmaV = module->GetLowerSigmaV();
    }

    //G4double du = m_smearedHits[i].x() - hit.x();
    //chi2 += pow(du/sigmaU,2.);

    G4double dv = m_smearedHits[i].y() - hit.y();
    chi2 += pow(dv/sigmaV,2.);

    m_smearedHits[i] = backwardRotation*m_smearedHits[i];

  }
    
  return chi2;
}

void RES_TrackFitter::ScanChi2Function(G4int i, G4int j, G4String filename)
{
  G4RunManager* runManager = G4RunManager::GetRunManager();
  const RES_PrimaryGeneratorAction* genAction = (RES_PrimaryGeneratorAction*) runManager->GetUserPrimaryGeneratorAction();
  G4ParticleGun* gun = genAction->GetParticleGun();
  m_initialCharge = gun->GetParticleCharge();

  CalculateStartParameters();
  //SetStartParametesToGeneratedParticle();
  DoBlobelFit(5);

  // hardcoded at the moment
  m_lowerBound[0] = 0.2*GeV;        m_upperBound[0] = 1.8*GeV;
  m_lowerBound[1] = -0.5*cm;         m_upperBound[1] = 0.5*cm;
  m_lowerBound[2] = -1.*M_PI/180.;  m_upperBound[2] = 1.*M_PI/180.;
  m_lowerBound[3] = -0.5*cm;         m_upperBound[3] = 0.5*cm;
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
      m_parameter[0] = 1./m_parameter[0];
      G4double chi2 = Chi2InModuleFrame();
      m_parameter[0] = 1./m_parameter[0];
      file << m_parameter[i] << "\t" << m_parameter[j] << "\t" << chi2 << std::endl;;
      m_parameter[j] += m_step[j];
      //G4cout << i << j<< " --> " << m_parameter[i] << "  " << m_parameter[j] << " ---> " << chi2 << G4endl;
    }
    m_parameter[i] += m_step[i];
  }

  file.close();
}

void RES_TrackFitter::AddLayerToBeSkipped(G4int layer)
{
  std::vector<G4int>::iterator it = std::find(m_layersToBeSkipped.begin(), m_layersToBeSkipped.end(), layer);
  if (it == m_layersToBeSkipped.end())
    m_layersToBeSkipped.push_back(layer);
}

void RES_TrackFitter::RemoveLayerToBeSkipped(G4int layer)
{
  std::vector<G4int>::iterator it = std::find(m_layersToBeSkipped.begin(), m_layersToBeSkipped.end(), layer);
  if (it != m_layersToBeSkipped.end())
    m_layersToBeSkipped.erase(it);
}

void RES_TrackFitter::DoFullFit()
{
  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();

  unsigned int nHits = m_currentGenEvent.GetNbOfHits();

  // basic dimensions of matrices
  unsigned int nRow = nHits;
  unsigned int nCol = 2;

  // declare matrices for the calculation
  TMatrixD A(nRow,nCol);
  TVectorD b(nRow);
  TMatrixD U(nRow,nRow);

  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  RES_MagneticField* field = (RES_MagneticField*) fieldMgr->GetDetectorField();
  G4double B_estimate  = field->GetFieldEstimate();
  G4double magnetHeight = fabs(field->GetZ1() - field->GetZ0());
  G4double z0_magnet = -magnetHeight/2. + field->GetDisplacement().z();
  G4double z1_magnet =  magnetHeight/2. + field->GetDisplacement().z();
  G4double L = fabs(z1_magnet - z0_magnet);
  
  std::cout << "--------------------------------------------------------" << std::endl;
  for (unsigned int i = 0; i < nHits; i++) {
    G4ThreeVector pos = m_smearedHits[i];

    // get module
    G4int iModule = m_currentGenEvent.GetModuleID(i);
    G4int iLayer  = m_currentGenEvent.GetLayerID(i);
    RES_Module* module = det->GetModule(iModule);
    G4double sigmaV = iLayer==0? module->GetUpperSigmaV() : module->GetLowerSigmaV();
    G4float k = z1_magnet - pos.z();

    // fill the matrices
    A(i,0) = 1;
    if (k <= 0)          A(i,1) = 0.;
    if (k > 0 && k <= L) A(i,1) = 0.5*0.3*(B_estimate/tesla)/m * pow(k,2.) * GeV;
    if (k > L)           A(i,1) = ( 0.5*0.3*(B_estimate/tesla)/m *pow(L,2.) + 0.3*(B_estimate/tesla)/m * L*(k-L)) * GeV;
    b(i) = pos.y();
    // G4double zahl1 = A(i,0)*1.0*cm;
    // G4double zahl2 = A(i,1)*1./(1.0*GeV);
    //std::cout << k/cm  << " ----> " << pos.y()/cm << "   " << zahl1/cm << " " << zahl2/cm << "   " << (zahl1+zahl2)/cm << std::endl;
    U(i,i) = sigmaV*sigmaV;
  }
  
  // calculate solution
  TMatrixD Uinv = U;
  Uinv.Invert();
  TMatrixD ATranspose(nCol,nRow);
  ATranspose.Transpose(A);
  TMatrixD M = ATranspose * Uinv * A;
  TVectorD c = ATranspose * Uinv * b;
  TMatrixD Minv = M;
  Minv.Invert();
  TVectorD solution(nCol);
  solution = Minv * c;

  // calculate chi2 and track positions from fit parameters
  TMatrixD vec(nRow,1);
  for (unsigned int i = 0; i < nRow; i++)
    vec(i,0) = (A*solution - b)(i);
  TMatrixD vecTrans(1,nRow);
  vecTrans.Transpose(vec);
  G4double chi2 = (vecTrans * Uinv * vec)(0,0);
  //std::cout << solution(0)/cm << " " << (1./solution(1)) / GeV << " " << chi2 << std::endl;
}
