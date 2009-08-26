#ifndef RES_Event_hh
#define RES_Event_hh

#include "TObject.h"

#include <vector>
#include <assert.h>
#include <TVector3.h>

enum EventType {
  generated, reconstructed
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
  inline void SetMomentum(double momentum) {m_momentum = momentum;}
  inline void SetChi2OverDof(double chi2_over_dof) {m_chi2_over_dof = chi2_over_dof;}

  inline int          GetID()                        {return m_ID;}
  inline int          GetModuleID(unsigned int i)    {assert(i <m_moduleID.size()); return m_moduleID.at(i);}
  inline int          GetFiberID(unsigned int i)     {assert(i <m_fiberID.size()); return m_fiberID.at(i);}
  inline EventType    GetEventType()                 {return m_eventType;}
  inline double       GetMomentum()                  {return m_momentum;}
  inline double       GetChi2OverDof()               {return m_chi2_over_dof;}
  inline unsigned int GetNbOfHits()                  {return m_hits.size();}
  inline TVector3     GetHitPosition(unsigned int i) {assert(i < m_hits.size()); return m_hits.at(i);}

  void AddHit(int module_ID, int fiber_ID, double x, double y, double z);

private:
  int                   m_ID;
  std::vector<int>      m_moduleID;
  std::vector<int>      m_fiberID;
  std::vector<TVector3> m_hits;
  double                m_momentum;
  EventType             m_eventType;
  double                m_chi2_over_dof;
  
  ClassDef( RES_Event, 1 );

};

#endif /* RES_Event_hh */
