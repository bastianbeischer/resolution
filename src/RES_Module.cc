#include "RES_Module.hh"

#include <cmath>

RES_Module::RES_Module()
{
  m_placement = G4ThreeVector(0,0,0);
  SetDefaultValuesForFiber();
}

RES_Module::RES_Module(G4ThreeVector placement)
{
  m_placement = placement;
  SetDefaultValuesForFiber();
}

RES_Module::~RES_Module()
{
}

void RES_Module::SetDefaultValuesForFiber()
{
  m_type          = fiber;
  m_angle         = 0.;
  m_internalAngle = 0.;
  m_length        = 30.*cm;
  m_width         = 6.912*cm;
  m_upperSigmaU   = m_length/sqrt(12);
  m_upperSigmaV   = 50.*um;
  m_upperSigmaZ   = 0.*cm;
  m_lowerSigmaU   = m_length/sqrt(12);
  m_lowerSigmaV   = 50.*um;
  m_lowerSigmaZ   = 0.*cm;
  m_efficiency    = 1.;
}

void RES_Module::SetDefaultValuesForSilicon()
{
  m_type          = silicon;
  m_angle         = 0.;
  m_internalAngle = M_PI/2.;
  m_length        = 30.*cm;
  m_width         = 10.*cm;
  m_upperSigmaU   = m_length/sqrt(12);
  m_upperSigmaV   = 10.*um;
  m_upperSigmaZ   = 0.*cm;
  m_lowerSigmaU   = m_length/sqrt(12);
  m_lowerSigmaV   = 30.*um;
  m_lowerSigmaZ   = 0.*cm;
  m_efficiency    = 1.;
}
