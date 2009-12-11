// $Id: RES_Hit.hh,v 1.1 2009/12/11 12:52:23 beischer Exp $

#ifndef RES_Hit_hh
#define RES_Hit_hh

#include "G4VHit.hh"
#include "G4ThreeVector.hh"

class RES_Hit : public G4VHit
{

public:
  RES_Hit();
  ~RES_Hit();

public:
  void Draw();
  void Print();

  inline void          SetModuleID(G4int ID)               {m_moduleID = ID;}
  inline void          SetLayerID(G4int ID)                {m_layerID = ID;}
  inline void          SetPosition(G4ThreeVector position) {m_position = position;}

  inline G4int         GetModuleID() {return m_moduleID;}
  inline G4int         GetLayerID()  {return m_layerID;}
  inline G4ThreeVector GetPosition() {return m_position;}

private:
  G4int         m_moduleID;
  G4int         m_layerID;
  G4ThreeVector m_position;

};

#endif /* RES_Hit_hh */
