#ifndef RES_PhysicsList_hh
#define RES_PhysicsList_hh

#include "G4VUserPhysicsList.hh"

class RES_PhysicsList : public G4VUserPhysicsList
{

public:
  RES_PhysicsList();
  ~RES_PhysicsList();

private:
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();

};

#endif /* RES_PhysicsList_hh */
