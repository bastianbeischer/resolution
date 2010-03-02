// $Id: RES_Event.cc,v 1.10 2010/03/02 18:31:47 beischer Exp $

#include "RES_Event.hh"

ClassImp( RES_Event );

RES_Event::RES_Event() :
  m_ID(-1),
  m_moduleID(0),
  m_layerID(0),
  m_hits(0),
  m_smearedHits(0),
  m_transverseMomentum(0.),
  m_momentum(0.),
  m_eventType(generated),
  m_chi2(0.),
  m_dof(0),
  m_initialParameters(0),
  m_finalParameters(0)
{
}

RES_Event::RES_Event(const RES_Event& other)
{
  m_ID = other.m_ID;
  m_moduleID = other.m_moduleID;
  m_layerID = other.m_layerID;
  m_hits.clear();
  for (unsigned int i = 0; i < other.m_hits.size(); i++) {
    m_hits.push_back(other.m_hits.at(i));
  }
  m_smearedHits.clear();
  for (unsigned int i = 0; i < other.m_smearedHits.size(); i++) {
    m_smearedHits.push_back(other.m_smearedHits.at(i));
  }
  m_transverseMomentum = other.m_transverseMomentum;
  m_momentum = other.m_momentum;
  m_eventType = other.m_eventType;
  m_chi2 = other.m_chi2;
  m_dof = other.m_dof;

  m_initialParameters.clear();
  for (unsigned int i = 0; i < other.m_initialParameters.size(); i++)
    m_initialParameters.push_back(other.m_initialParameters.at(i));

  m_finalParameters.clear();
  for (unsigned int i = 0; i < other.m_finalParameters.size(); i++)
    m_finalParameters.push_back(other.m_finalParameters.at(i));
}

RES_Event::~RES_Event()
{
}

const RES_Event& RES_Event::operator=(const RES_Event& right)
{
  m_ID = right.m_ID;
  m_moduleID = right.m_moduleID;
  m_layerID = right.m_layerID;
  m_hits.clear();
  for (unsigned int i = 0; i < right.m_hits.size(); i++) {
    m_hits.push_back(right.m_hits.at(i));
  }
  m_smearedHits.clear();
  for (unsigned int i = 0; i < right.m_smearedHits.size(); i++) {
    m_smearedHits.push_back(right.m_smearedHits.at(i));
  }
  m_transverseMomentum = right.m_transverseMomentum;
  m_momentum = right.m_momentum;
  m_eventType = right.m_eventType;
  m_chi2 = right.m_chi2;
  m_dof = right.m_dof;

  m_initialParameters.clear();
  for (unsigned int i = 0; i < right.m_initialParameters.size(); i++)
    m_initialParameters.push_back(right.m_initialParameters.at(i));

  m_finalParameters.clear();
  for (unsigned int i = 0; i < right.m_finalParameters.size(); i++)
    m_finalParameters.push_back(right.m_finalParameters.at(i));

  return *this;
}

void RES_Event::SetInitialParameters(int nPar, double* par)
{
  m_initialParameters.clear();
  for (int i = 0; i < nPar; i++) {
    m_initialParameters.push_back(par[i]);
  }
}

void RES_Event::SetFinalParameters(int nPar, double* par)
{
  m_finalParameters.clear();
  for (int i = 0; i < nPar; i++) {
    m_finalParameters.push_back(par[i]);
  }
}

void RES_Event::AddHit(int module_ID, int layer_ID, double x, double y, double z)
{
  m_moduleID.push_back(module_ID);
  m_layerID.push_back(layer_ID);
  m_hits.push_back(TVector3(x,y,z));
}
void RES_Event::AddSmearedHit(double x, double y, double z)
{
  m_smearedHits.push_back(TVector3(x,y,z));
}
