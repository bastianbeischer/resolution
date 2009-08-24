#include "RES_DummyField.hh"

#include "globals.hh"

RES_DummyField::RES_DummyField(G4ThreeVector fieldVector)
{
  m_fieldVector = fieldVector;
}

RES_DummyField::~RES_DummyField()
{
}

void RES_DummyField::GetFieldValue(const G4double* x, G4double* B) const
{
  B[0] = B[1] = B[2] = 0.;
  if ( (x[2] > -4.*cm) && (x[2] < 4.*cm) ) {
    B[0] = m_fieldVector.x();
    B[1] = m_fieldVector.y();
    B[2] = m_fieldVector.z();
  }
}
