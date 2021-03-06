// $Id: RES_Event.hh,v 1.10.2.1 2010/09/08 17:44:58 beischer Exp $

#ifndef RES_Event_hh
#define RES_Event_hh

#include "TObject.h"

#include <vector>
#include <assert.h>
#include <TVector3.h>

enum EventType {
  dummy, generated, reconstructed
};

class RES_Event : public TObject
{

public:
  RES_Event();
  RES_Event(const RES_Event& other);
  ~RES_Event();
  const RES_Event& operator=(const RES_Event& right);

public:
  inline void SetID(int ID) {m_ID = ID;}
  inline void SetEventType(EventType type) {m_eventType = type;}
  inline void SetTransverseMomentum(double pt) {m_transverseMomentum = pt;}
  inline void SetMomentum(double momentum) {m_momentum = momentum;}
  inline void SetChi2(double chi2) {m_chi2 = chi2;}
  inline void SetDof(int dof) {m_dof = dof;}
  void SetInitialParameters(unsigned int nPar, double* par);
  void SetFinalParameters(unsigned int nPar, double* par);

  inline int          GetID()                               {return m_ID;}
  inline int          GetModuleID(unsigned int i)           {assert(i <m_moduleID.size()); return m_moduleID.at(i);}
  inline int          GetLayerID(unsigned int i)            {assert(i <m_layerID.size()); return m_layerID.at(i);}
  inline EventType    GetEventType()                        {return m_eventType;}
  inline double       GetTransverseMomentum()               {return m_transverseMomentum;}
  inline double       GetMomentum()                         {return m_momentum;}
  inline double       GetChi2()                             {return m_chi2;}
  inline int          GetDof()                              {return m_dof;}
  inline unsigned int GetNbOfHits()                         {return m_hits.size();}
  inline TVector3     GetHitPosition(unsigned int i)        {assert(i < m_hits.size()); return m_hits.at(i);}
  inline TVector3     GetSmearedHitPosition(unsigned int i) {assert(i < m_smearedHits.size()); return m_smearedHits.at(i);}
  inline double       GetInitialParameter(unsigned int i)   {assert(i < m_initialParameters.size()); return m_initialParameters.at(i);}
  inline double       GetFinalParameter(unsigned int i)     {assert(i < m_finalParameters.size()); return m_finalParameters.at(i);}

  void AddHit(int module_ID, int layer_ID, double x, double y, double z);
  void AddSmearedHit(double x, double y, double z);
  TVector3 GetHitAtZ(double z);
  TVector3 GetSmearedHitAtZ(double z);

private:
  int                   m_ID;
  std::vector<int>      m_moduleID;
  std::vector<int>      m_layerID;
  std::vector<TVector3> m_hits;
  std::vector<TVector3> m_smearedHits;
  double                m_transverseMomentum;
  double                m_momentum;
  EventType             m_eventType;
  double                m_chi2;
  unsigned int          m_dof;
  
  std::vector<double>   m_initialParameters;
  std::vector<double>   m_finalParameters;

  ClassDef( RES_Event, 2 );

};

#endif /* RES_Event_hh */
