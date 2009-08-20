#ifndef RES_FiberHit_hh
#define RES_FiberHit_hh

#include "G4VHit.hh"
#include "G4ThreeVector.hh"

class RES_FiberHit : public G4VHit
{

public:
  RES_FiberHit();
  ~RES_FiberHit();

public:
  void Draw();
  void Print();

  inline void SetPosition(G4ThreeVector position) {m_position = position;}

  inline G4ThreeVector GetPosition() {return m_position;}

private:
  G4ThreeVector m_position;

};

#endif /* RES_FiberHit_hh */
