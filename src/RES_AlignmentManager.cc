// $Id: RES_AlignmentManager.cc,v 1.5 2009/10/22 14:45:15 beischer Exp $

#include "RES_AlignmentManager.hh"

#include "RES_AlignmentMessenger.hh"
#include "RES_DetectorConstruction.hh"
#include "RES_DataHandler.hh"
#include "RES_RunManager.hh"
#include "RES_Event.hh"

#include "TMatrixD.h"

#include "millepede1.h"

#include <iostream>

RES_AlignmentManager::RES_AlignmentManager() :
  m_parameters(0),
  m_verbose(0)
{
  RES_RunManager* runMgr = (RES_RunManager*) G4RunManager::GetRunManager();

  m_messenger = new RES_AlignmentMessenger(this);
  m_dataHandler = runMgr->GetDataHandler();
}

RES_AlignmentManager::~RES_AlignmentManager()
{
  delete m_messenger;
  delete[] m_parameters;
}

void RES_AlignmentManager::StartAlignment()
{
  m_dataHandler->Initialize();

  G4int nEvents = m_dataHandler->GetNumberOfGeneratedEvents();

  m_dataHandler->LoadGeneratedEntry(0);
  RES_Event event = m_dataHandler->GetCurrentEvent();

  unsigned int nHits;
  G4int nGlobal, nLocal, nStdDev, iPar, nIter, iPrLim;
  G4float sigma, rhs, cutvalue;
  
  nHits = event.GetNbOfHits();
  nGlobal = nHits;

  delete[] m_parameters;
  m_parameters = new G4float[nGlobal];
  for (G4int i = 0 ; i < nGlobal; i++)
    m_parameters[i] = 0.;
  
  G4float* x_constraints = new G4float[nGlobal];
  G4float* y_constraints = new G4float[nGlobal];
  for(G4int i = 0; i < nGlobal; i++) {
    if ( i % 2 == 0)  {
      x_constraints[i] = 1.;
      y_constraints[i] = 0.;
    }
    else {
      x_constraints[i] = 0.;
      y_constraints[i] = 1.;
    }
  }

  INITGL(nGlobal, nLocal = 4, nStdDev = 0, iPrLim = m_verbose);
  //PARGLO(m_parameters);
  //CONSTF(x_constraints, rhs = 0.);
  //CONSTF(y_constraints, rhs = 0.);
  PARSIG(iPar = 1, sigma = 0.);
  PARSIG(iPar = 2, sigma = 0.);
  INITUN(nIter = 11, cutvalue = 100.);

  G4float* dergb = new G4float[nGlobal];
  G4float* derlc = new G4float[nLocal];
  for (G4int i = 0 ; i < nGlobal; i++)
    dergb[i] = 0.;
  for (G4int i = 0 ; i < nLocal; i++)
    derlc[i] = 0.;
  
  for(G4int iEvent = 0; iEvent < nEvents; iEvent++) {

    m_dataHandler->LoadGeneratedEntry(iEvent);
    event = m_dataHandler->GetCurrentEvent();
    
    if (event.GetNbOfHits() != nHits) continue;

    for (unsigned int iHit = 0; iHit < nHits; iHit++) {

      RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
      G4int iModule = event.GetModuleID(iHit);
      G4int iFiber  = event.GetFiberID(iHit);
      G4double angle = det->GetModuleAngle(iModule);
      if (iFiber > 0) angle += det->GetModuleInternalAngle(iModule);

      G4ThreeVector smearedHit(event.GetSmearedHitPosition(iHit).x(), event.GetSmearedHitPosition(iHit).y(), event.GetSmearedHitPosition(iHit).z());
      G4float z0 = 0.;
      G4float cotan = cos(angle)/sin(angle);
      G4float fx = smearedHit.x();
      G4float fy = smearedHit.y();
      G4float fz = smearedHit.z();
      G4float k = fz - z0;

      G4double sigmaV = det->GetModuleSigmaV(iModule);

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
    
      TMatrixD Lin(1,2);
      Lin(0,0) = 1.;
      Lin(0,1) = -cotan;
      TMatrixD LinTrans(2,1);
      LinTrans.Transpose(Lin);
      TMatrixD V3 = TMatrixD(1,1);
      V3 = Lin * V2 * LinTrans;

      dergb[2*(iHit/2)  ] = 1.;
      dergb[2*(iHit/2)+1] = -cotan;
      derlc[0]            = 1.;
      derlc[1]            = -cotan;
      derlc[2]            = k;
      derlc[3]            = -k*cotan;
      G4float y           = fx - fy*cotan;
      G4float sigma       = sqrt(V3(0,0));

      EQULOC(dergb, derlc, y, sigma);

    }

    FITLOC();

  }

  FITGLO(m_parameters);

  if (m_verbose > 0) {
    std::cout << " >>> Result of Alignment: " << std::endl;
    for(G4int i = 0 ; i < nGlobal; i++) {
      std::cout << "       parameter " << i << " --> " << m_parameters[i] << std::endl;
    }
  }

  delete[] dergb;
  delete[] derlc;
}
