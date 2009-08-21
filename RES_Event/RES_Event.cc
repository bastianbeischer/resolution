#include "RES_Event.hh"

ClassImp( RES_Event );

RES_Event::RES_Event() :
  m_ID(-1),
  m_momentum(0.),
  m_eventType(generated)
{
}

RES_Event::RES_Event(const RES_Event& other)
{
  m_ID = other.m_ID;
  m_hits.clear();
  for (unsigned int i = 0; i < other.m_hits.size(); i++) {
    m_hits.push_back(other.m_hits.at(i));
  }
  m_momentum = other.m_momentum;
  m_eventType = other.m_eventType;
}

RES_Event::~RES_Event()
{
}

const RES_Event& RES_Event::operator=(const RES_Event& right)
{
  m_ID = right.m_ID;
  m_hits.clear();
  for (unsigned int i = 0; i < right.m_hits.size(); i++) {
    m_hits.push_back(right.m_hits.at(i));
  }
  m_momentum = right.m_momentum;
  m_eventType = right.m_eventType;
  return *this;
}

void RES_Event::AddHit(double x, double y, double z)
{
  m_hits.push_back(TVector3(x,y,z));
}
