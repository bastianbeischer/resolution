// $Id: RES_AlignmentManager.cc,v 1.1 2009/10/14 16:51:32 beischer Exp $

#include "RES_AlignmentManager.hh"

#include "RES_AlignmentMessenger.hh"
#include "RES_DetectorConstruction.hh"
#include "RES_DataHandler.hh"
#include "RES_RunManager.hh"
#include "RES_Event.hh"

#include "millepede1.h"

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
  INITUN(nIter = 11, cutvalue = 0.);

  float* dergb = new float[nGlobal];
  float* derlc = new float[nLocal];
  ZERLOC(dergb, derlc);
  
  for(int i = 0; i < nEvents; i++) {

    m_dataHandler->LoadGeneratedEntry(i);
    event = m_dataHandler->GetCurrentEvent();
    
    for (int j = 0; j < nHits; j++) {

      RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
      G4int iModule = event.GetModuleID(i);
      G4int iFiber  = event.GetFiberID(i);
      G4double angle = det->GetModuleAngle(iModule);
      if (iFiber > 0) angle += det->GetModuleInternalAngle(iModule);

      G4ThreeVector smearedHit(event.GetSmearedHitPosition(i).x(), event.GetSmearedHitPosition(i).y(), event.GetSmearedHitPosition(i).z());
      float z0 = 0.;
      float cotan = cos(angle)/sin(angle);
      float fx = smearedHit.x();
      float fy = smearedHit.y();
      float fz = smearedHit.z();
      float k = fz - z0;

      dergb[2*j]   = 1;
      dergb[2*j+1] = -cotan;
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
