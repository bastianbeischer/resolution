#ifndef RES_TrackFitter_hh
#define RES_TrackFitter_hh

#include "RES_Event.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"

enum FitMethod
{
  blobel, minuit
};

class RES_TrackFitMessenger;

class RES_TrackFitter
{

public:
  RES_TrackFitter();
  ~RES_TrackFitter();

public:
  inline void SetCurrentGenEvent(RES_Event event) {m_currentGenEvent = event;}
  inline void SetCurrentRecEvent(RES_Event event) {m_currentRecEvent = event;}
  inline void SetVerbose(G4int verbose) {m_verbose = verbose;}
  inline void SetFitMethod(FitMethod method) {m_fitMethod = method;}

  RES_Event Fit();

private:
  void     SetSpatialResolutions();
  void     SmearHits();
  void     CalculateStartParameters();
  G4int    DoBlobelFit();
  G4int    DoBlobelFitInPlane();
  G4int    DoMinuitFit();
  G4double Chi2();
  G4double Chi2InPlane();

private:
  RES_TrackFitMessenger* m_messenger;

  RES_Event              m_currentGenEvent;
  RES_Event              m_currentRecEvent;
  G4int                  m_verbose;

  G4ThreeVector*         m_smearedHits;
  G4double*              m_parameter;
  FitMethod              m_fitMethod;

  G4double               m_sigmaX;
  G4double               m_sigmaY;
  G4double               m_sigmaZ;

};

#endif /* RES_TrackFitter_hh */
