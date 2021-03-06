// $Id: RES_TrackFitter.hh,v 1.22 2010/07/14 12:17:36 beischer Exp $

#ifndef RES_TrackFitter_hh
#define RES_TrackFitter_hh

#include "RES_Event.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>

enum FitMethod
{
  blobel, minuit, oneline, twolines, transverse, brokenline
};

class RES_TrackFitMessenger;

class RES_TrackFitter
{

public:
  ~RES_TrackFitter();

public:
  static RES_TrackFitter* GetInstance();

  inline void SetCurrentGenEvent(RES_Event event) {m_currentGenEvent = event; CopyHits();}
  inline void SetCurrentRecEvent(RES_Event event) {m_currentRecEvent = event;}
  inline void SetVerbose(G4int verbose) {m_verbose = verbose;}
  inline void SetFitMethod(FitMethod method) {m_fitMethod = method;}

  RES_Event Fit();

  void ScanChi2Function(G4int i, G4int j, G4String filename);

public:
  void CopyHits();
  void FitStraightLine(G4int n0, G4int n1, G4double &x0, G4double &y0, G4double &dxdz, G4double &dydz);
  void BrokenLineFit();

  void ClearLayersToBeSkipped() {m_layersToBeSkipped.clear();}
  void AddLayerToBeSkipped(G4int layer);
  void RemoveLayerToBeSkipped(G4int layer);

private:
  RES_TrackFitter();

  void     SetStartParametesToGeneratedParticle();
  void     CalculateStartParameters();
  G4int    DoBlobelFit(G4int npar);
  G4int    DoMinuitFit(G4int npar);
  G4double Chi2InModuleFrame();

private:
  static RES_TrackFitter* m_instance;

  RES_TrackFitMessenger*  m_messenger;

  RES_Event               m_currentGenEvent;
  RES_Event               m_currentRecEvent;
  G4int                   m_verbose;

  G4ThreeVector*          m_smearedHits;
  G4double*               m_initialParameter;
  G4double*               m_parameter;
  G4double*               m_step;
  G4double*               m_lowerBound;
  G4double*               m_upperBound;
  FitMethod               m_fitMethod;

  G4double                m_initialCharge;

  std::vector<int>        m_layersToBeSkipped;

  friend void MinuitChi2Wrapper(int& npar, double* /*gin*/, double& f, double* par, int /*iflag*/);

};

#endif /* RES_TrackFitter_hh */
