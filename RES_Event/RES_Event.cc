#include "RES_Event.hh"

ClassImp( RES_Event );

RES_Event::RES_Event() :
  m_ID(-1),
  m_momentum(0.),
  m_eventType(generated)
{
}

RES_Event::~RES_Event()
{
}

void RES_Event::AddHit(double x, double y, double z)
{
  m_hits.push_back(TVector3(x,y,z));
}
