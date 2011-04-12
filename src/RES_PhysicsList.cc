// $Id: RES_PhysicsList.cc,v 1.5 2009/12/09 21:47:57 beischer Exp $

#include "RES_PhysicsList.hh"

#include "G4EmStandardPhysics.hh"

RES_PhysicsList::RES_PhysicsList()
{
  RegisterPhysics(new G4EmStandardPhysics);
}

RES_PhysicsList::~RES_PhysicsList()
{
}

void RES_PhysicsList::SetCuts()
{
  SetCutsWithDefault();
}

