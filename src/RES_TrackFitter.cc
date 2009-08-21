#include "RES_TrackFitter.hh"

#include "RES_Event.hh"
#include "RES_RunManager.hh"

#include "globals.hh"

#include "blobel.h"

RES_TrackFitter::RES_TrackFitter()
{
}

RES_TrackFitter::~RES_TrackFitter()
{
}

RES_Event RES_TrackFitter::Fit(RES_Event genEvent)
{
  G4double parameter[3];
  parameter[0] = genEvent.GetHit(0).x();
  parameter[1] = genEvent.GetHit(0).y();
  parameter[2] = genEvent.GetHit(0).z();

  G4double step[3];
  step[0] = 0.1*parameter[0];
  step[1] = 0.1*parameter[1];
  step[2] = 0.1*parameter[2];

  G4int npar = 3; // negative to suppress printout
  G4int nflim = 3*abs(npar)*(abs(npar)+10); // dito
  // if( fVerboseLevel == 0 ){
  //   npar = -npar;
  //   nflim = -nflim;
  // }
  DVALLIN(npar,step,nflim);

  G4int conv = -1;
  G4int iter = 0;
  G4double chi2 = 0.;
  while( conv <= 0 ){
    ++iter;
    chi2 = 0.;

    RES_RunManager* runManager = (RES_RunManager*) G4RunManager::GetRunManager();
    runManager->BeamOn(1);

    G4int nHits = genEvent.GetNbOfHits();

    G4double sigmaX = 50*um;
    G4double sigmaY = 300*um;

    for( G4int i=0 ; i<nHits ; i++ ){
      chi2 += pow( genEvent.GetHit(i).x()*cm - m_currentRecEvent.GetHit(i).x()*cm, 2.) / sigmaX;
      chi2 += pow( genEvent.GetHit(i).y()*cm - m_currentRecEvent.GetHit(i).y()*cm, 2.) / sigmaY;
    }

    DVALLEY(chi2,parameter,conv);
  }

  if (conv) return m_currentRecEvent;
  else      return RES_Event();
}
