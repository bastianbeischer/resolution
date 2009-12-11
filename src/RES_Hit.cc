// $Id: RES_Hit.cc,v 1.1 2009/12/11 12:52:25 beischer Exp $

#include "RES_Hit.hh"

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"

RES_Hit::RES_Hit()
{
}

RES_Hit::~RES_Hit()
{
}

void RES_Hit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(m_position);
    circle.SetScreenSize(6.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void RES_Hit::Print()
{
  G4cout<< "HIT: at (" << m_position.x()/cm << ", " << m_position.y()/cm << ", " << m_position.z()/cm << ") cm" << G4endl;
}
