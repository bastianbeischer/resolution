// $Id: RES_PhysicsList.hh,v 1.3 2009/10/14 09:24:30 beischer Exp $

#ifndef RES_PhysicsList_hh
#define RES_PhysicsList_hh

#include "G4VModularPhysicsList.hh"

class RES_PhysicsList : public G4VModularPhysicsList
{

public:
  RES_PhysicsList();
  ~RES_PhysicsList();

private:
  void SetCuts();

};

#endif /* RES_PhysicsList_hh */
