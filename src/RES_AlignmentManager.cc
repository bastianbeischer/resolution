// $Id: RES_AlignmentManager.cc,v 1.2 2009/10/14 17:29:04 beischer Exp $

#include "RES_AlignmentManager.hh"

#include "RES_AlignmentMessenger.hh"
#include "RES_DetectorConstruction.hh"
#include "RES_DataHandler.hh"
#include "RES_RunManager.hh"
#include "RES_Event.hh"

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

  int nHits, nGlobal, nLocal, nStdDev, iPrLim, iPar, nIter;
  float fixvalue, cutvalue;
  
  nHits = event.GetNbOfHits();
  nGlobal = 2*nHits;

  delete[] m_parameters;
  m_parameters = new float[nGlobal];

  INITGL(nGlobal, nLocal = 4, nStdDev = 3, iPrLim = 1);
  PARSIG(iPar = 1, fixvalue = 0.);
  PARSIG(iPar = 2, fixvalue = 0.);
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

    for (int iHit = 0; iHit < nHits; iHit++) {

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

      dergb[2*iHit]   = 1;
      dergb[2*iHit+1] = -cotan;
      derlc[0]     = 1.;
      derlc[1]     = -cotan;
      derlc[2]     = k;
      derlc[3]     = -k*cotan;
      float y      = fx - fy*cotan;
      float sigma  = 1.;
      EQULOC(dergb, derlc, y, sigma);

    }

    FITLOC();

  }

  FITGLO(m_parameters);

  delete[] dergb;
  delete[] derlc;
}
