// $Id: RES_AlignmentManager.cc,v 1.3 2009/10/19 09:24:18 beischer Exp $

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
  m_parameters(0)
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

  int nEvents = m_dataHandler->GetNumberOfGeneratedEvents();

  m_dataHandler->LoadGeneratedEntry(0);
  RES_Event event = m_dataHandler->GetCurrentEvent();

  unsigned int nHits;
  int nGlobal, nLocal, nStdDev, iPar, nIter, iPrLim;
  float fixvalue, rhs, cutvalue;
  
  nHits = event.GetNbOfHits();
  nGlobal = 2*nHits;

  delete[] m_parameters;
  m_parameters = new float[nGlobal];
  for (int i = 0 ; i < nGlobal; i++)
    m_parameters[i] = 0.;
  
  float* constraints = new float[nGlobal];
  for (int i = 0; i < nGlobal; i++)
    constraints[i] = 0.;

  for(int i = 1; i < nHits; i++) {
    constraints[2*i+1] = 1.0/(event.GetHitPosition(i).z() - event.GetHitPosition(0).z()) ;
  }

  INITGL(nGlobal, nLocal = 4, nStdDev = 0, iPrLim = 1);
  // PARSIG(iPar = 1, fixvalue = 0.);
  // PARSIG(iPar = 2, fixvalue = 0.);
  // PARGLO(m_parameters);
  CONSTF(constraints, rhs = 0.);
  INITUN(nIter = 11, cutvalue = 100.);

  float* dergb = new float[nGlobal];
  float* derlc = new float[nLocal];
  for (int i = 0 ; i < nGlobal; i++)
    dergb[i] = 0.;
  for (int i = 0 ; i < nLocal; i++)
    derlc[i] = 0.;
  
  for(int iEvent = 0; iEvent < nEvents; iEvent++) {

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
      float z0 = 0.;
      float cotan = cos(angle)/sin(angle);
      float fx = smearedHit.x();
      float fy = smearedHit.y();
      float fz = smearedHit.z();
      float k = fz - z0;

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

      dergb[2*iHit]   = 1.;
      dergb[2*iHit+1] = -cotan;
      derlc[0]        = 1.;
      derlc[1]        = -cotan;
      derlc[2]        = k;
      derlc[3]        = -k*cotan;
      float y         = fx - fy*cotan;
      float sigma     = sqrt(V3(0,0));

      EQULOC(dergb, derlc, y, sigma);

    }

    FITLOC();

  }

  FITGLO(m_parameters);

  for(int i = 0 ; i < nGlobal; i++) {
    std::cout << " >>> Result of Alignment: " << std::endl;
    std::cout << "       parameter " << i << " --> " << m_parameters[i] << std::endl;
  }

  delete[] dergb;
  delete[] derlc;
}
