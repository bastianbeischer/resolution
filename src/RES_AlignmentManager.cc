// $Id: RES_AlignmentManager.cc,v 1.10 2009/12/11 12:52:24 beischer Exp $

#include "RES_AlignmentManager.hh"

#include "RES_AlignmentMessenger.hh"
#include "RES_DetectorConstruction.hh"
#include "RES_TrackFitter.hh"
#include "RES_DataHandler.hh"
#include "RES_RunManager.hh"
#include "RES_Event.hh"

#include "TMatrixD.h"

#include "millepede1.h"

#include <iostream>

RES_AlignmentManager* RES_AlignmentManager::m_instance = 0;

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

RES_AlignmentManager* RES_AlignmentManager::GetInstance()
{
  if (!m_instance) m_instance = new RES_AlignmentManager();
  return m_instance;
}

void RES_AlignmentManager::StartAlignment()
{
  m_dataHandler->Initialize();

  G4int nEvents = m_dataHandler->GetNumberOfGeneratedEvents();

  m_dataHandler->LoadGeneratedEntry(0);
  RES_Event event = m_dataHandler->GetCurrentEvent();

  unsigned int nHits;
  G4int nGlobal, nLocal, nStdDev, iPar, nIter, iPrLim;
  G4float sigma, /*rhs,*/ cutvalue;
  
  nHits = event.GetNbOfHits();
  int nModules = nHits/2;
  nGlobal = 3*nModules;

  delete[] m_parameters;
  m_parameters = new G4float[nGlobal];
  for (G4int i = 0 ; i < nGlobal; i++)
    m_parameters[i] = 0.;
  
  G4float* dergb = new G4float[nGlobal];
  G4float* derlc = new G4float[nLocal];
  for (G4int i = 0 ; i < nGlobal; i++)
    dergb[i] = 0.;
  for (G4int i = 0 ; i < nLocal; i++)
    derlc[i] = 0.;

  G4float* x_constraints = new G4float[nGlobal];
  G4float* y_constraints = new G4float[nGlobal];
  for(G4int i = 0; i < nGlobal/2; i++) {
    if ( i % 2 == 0)  {
      x_constraints[i] = 1.;
      y_constraints[i] = 0.;
    }
    else {
      x_constraints[i] = 0.;
      y_constraints[i] = 1.;
    }
  }

  int nGlobalIter = 1;
  for (int iGlobalIter = 0; iGlobalIter < nGlobalIter; iGlobalIter++) {

    INITGL(nGlobal, nLocal = 4, nStdDev = 0, iPrLim = m_verbose);
    PARGLO(m_parameters);
    //CONSTF(x_constraints, rhs = 0.);
    //CONSTF(y_constraints, rhs = 0.);
    PARSIG(iPar = 1, sigma = 0.);
    PARSIG(iPar = nModules+1, sigma = 0.);
    PARSIG(iPar = 2*nModules+1, sigma = 0.);
    for (int i = 2; i <= nModules; i++)
      PARSIG(iPar = 2*nModules+i, sigma = 0.);
    INITUN(nIter = 11, cutvalue = 100.);
  
    for(G4int iEvent = 0; iEvent < nEvents; iEvent++) {

      m_dataHandler->LoadGeneratedEntry(iEvent);
      event = m_dataHandler->GetCurrentEvent();
    
      if (event.GetNbOfHits() != nHits) continue;

      RES_TrackFitter* fitter = RES_TrackFitter::GetInstance();
      fitter->SetCurrentGenEvent(event);
      G4double x0,y0,lambda_x,lambda_y;
      fitter->FitStraightLine(0,nHits,x0,y0,lambda_x,lambda_y);

      for (unsigned int iHit = 0; iHit < nHits; iHit++) {

        RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
        G4int iModule = event.GetModuleID(iHit);
        G4int iLayer  = event.GetLayerID(iHit);
        G4double angle = det->GetModuleAngle(iModule);
        if (iLayer > 0) angle += det->GetModuleInternalAngle(iModule);

        G4ThreeVector smearedHit(event.GetSmearedHitPosition(iHit).x(), event.GetSmearedHitPosition(iHit).y(), event.GetSmearedHitPosition(iHit).z());
        G4float z0 = 0.;
        G4float cotan = cos(angle)/sin(angle);
        G4float fx = smearedHit.x();
        G4float fy = smearedHit.y();
        G4float fz = smearedHit.z();
        G4float k = fz - z0;
        G4float y = fx - fy*cotan;

        G4double sigmaV = 0.;
        if (iLayer == 0)
          sigmaV = det->GetModuleUpperSigmaV(iModule);
        else
          sigmaV = det->GetModuleUpperSigmaV(iModule);
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
        G4float sigma = sqrt(V3(0,0));

        int module = iHit/2;
      
        dergb[module           ] = 1.;
        dergb[nModules+module  ] = -(cotan + m_parameters[2*nModules+module]);
        dergb[2*nModules+module] = m_parameters[nModules+module] + fy + y0 - k*lambda_y;
        derlc[0]                 = 1.;
        derlc[1]                 = -(cotan + m_parameters[2*nModules+module]);
        derlc[2]                 = k;
        derlc[3]                 = -k*(cotan + m_parameters[2*nModules+module]);


        EQULOC(dergb, derlc, y, sigma);

      }

      FITLOC();

    }

    FITGLO(m_parameters);

    if (m_verbose > 0) {
      std::cout << " >>> Result of Alignment Iteration " << iGlobalIter << ": " << std::endl;
      for(G4int i = 0 ; i < nGlobal; i++) {
        std::cout << "       parameter " << i << " --> " << m_parameters[i] << std::endl;
      }
    }

  }

  delete[] dergb;
  delete[] derlc;
}
