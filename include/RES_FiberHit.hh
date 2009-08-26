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

  inline void          SetModuleID(G4int ID)               {m_moduleID = ID;}
  inline void          SetFiberID(G4int ID)                {m_fiberID = ID;}
  inline void          SetPosition(G4ThreeVector position) {m_position = position;}

  inline G4int         GetModuleID() {return m_moduleID;}
  inline G4int         GetFiberID()  {return m_fiberID;}
  inline G4ThreeVector GetPosition() {return m_position;}

private:
  G4int         m_moduleID;
  G4int         m_fiberID;
  G4ThreeVector m_position;

};

#endif /* RES_FiberHit_hh */
