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

private:
  G4ThreeVector m_fieldVector;

};

#endif /* RES_DummyField_hh */
