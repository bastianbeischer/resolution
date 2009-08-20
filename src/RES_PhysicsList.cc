#include "RES_PhysicsList.hh"

#include "G4ParticleTypes.hh"

RES_PhysicsList::RES_PhysicsList()
{
}

RES_PhysicsList::~RES_PhysicsList()
{
}

void RES_PhysicsList::ConstructParticle()
{
  G4ChargedGeantino::ChargedGeantinoDefinition();
}

void RES_PhysicsList::ConstructProcess()
{
  AddTransportation();
}

void RES_PhysicsList::SetCuts()
{
  SetCutsWithDefault();
}

