// $Id: RES_DummyField.hh,v 1.3 2010/01/04 09:47:00 beischer Exp $

#ifndef RES_DummyField_hh
#define RES_DummyField_hh

#include "G4MagneticField.hh"

#include "G4ThreeVector.hh"

class RES_DummyField : public G4MagneticField
{

public:
   RES_DummyField(G4ThreeVector fieldVector);
  ~RES_DummyField();

public:
  void GetFieldValue(const G4double* x, G4double* B) const;

public:
  G4double GetZ0() const {return m_z0;}
  G4double GetZ1() const {return m_z1;}

private:
  G4ThreeVector m_fieldVector;

  G4double m_z0;
  G4double m_z1;

};

#endif /* RES_DummyField_hh */
