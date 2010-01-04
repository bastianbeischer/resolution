// $Id: RES_DummyField.cc,v 1.3 2010/01/04 09:47:00 beischer Exp $

#include "RES_DummyField.hh"

#include "globals.hh"

RES_DummyField::RES_DummyField(G4ThreeVector fieldVector)
{
  m_fieldVector = fieldVector;
  m_z0 = -4*cm;
  m_z1 =  4*cm;
}

RES_DummyField::~RES_DummyField()
{
}

void RES_DummyField::GetFieldValue(const G4double* x, G4double* B) const
{
  B[0] = B[1] = B[2] = 0.;
  if ( (x[2] > m_z0) && (x[2] < m_z1) ) {
    B[0] = m_fieldVector.x();
    B[1] = m_fieldVector.y();
    B[2] = m_fieldVector.z();
  }
}
