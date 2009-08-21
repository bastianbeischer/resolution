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

  inline int          GetID()                {return m_ID;}
  inline EventType    GetEventType()         {return m_eventType;}
  inline double       GetMomentum()          {return m_momentum;}
  inline unsigned int GetNbOfHits()          {return m_hits.size();}
  inline TVector3     GetHit(unsigned int i) {assert(i < m_hits.size()); return m_hits.at(i);}

  void AddHit(double x, double y, double z);

private:
  int                   m_ID;
  std::vector<TVector3> m_hits;
  double                m_momentum;
  EventType             m_eventType;

  ClassDef( RES_Event, 1 );

};

#endif /* RES_Event_hh */
