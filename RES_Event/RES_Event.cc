// $Id: RES_Event.cc,v 1.7 2009/10/14 09:24:35 beischer Exp $

#include "RES_Event.hh"

ClassImp( RES_Event );

RES_Event::RES_Event() :
  m_ID(-1),
  m_moduleID(0),
  m_fiberID(0),
  m_hits(0),
  m_transverseMomentum(0.),
  m_momentum(0.),
  m_eventType(generated),
  m_chi2(0.),
  m_dof(0)
{
}

RES_Event::RES_Event(const RES_Event& other)
{
  m_ID = other.m_ID;
  m_moduleID = other.m_moduleID;
  m_fiberID = other.m_fiberID;
  m_hits.clear();
  for (unsigned int i = 0; i < other.m_hits.size(); i++) {
    m_hits.push_back(other.m_hits.at(i));
  }
  m_transverseMomentum = other.m_transverseMomentum;
  m_momentum = other.m_momentum;
  m_eventType = other.m_eventType;
  m_chi2 = other.m_chi2;
  m_dof = other.m_dof;
}

RES_Event::~RES_Event()
{
}

const RES_Event& RES_Event::operator=(const RES_Event& right)
{
  m_ID = right.m_ID;
  m_moduleID = right.m_moduleID;
  m_fiberID = right.m_fiberID;
  m_hits.clear();
  for (unsigned int i = 0; i < right.m_hits.size(); i++) {
    m_hits.push_back(right.m_hits.at(i));
  }
  m_transverseMomentum = right.m_transverseMomentum;
  m_momentum = right.m_momentum;
  m_eventType = right.m_eventType;
  m_chi2 = right.m_chi2;
  m_dof = right.m_dof;
  return *this;
}

void RES_Event::AddHit(int module_ID, int fiber_ID, double x, double y, double z)
{
  m_moduleID.push_back(module_ID);
  m_fiberID.push_back(fiber_ID);
  m_hits.push_back(TVector3(x,y,z));
}
