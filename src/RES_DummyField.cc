// $Id: RES_DummyField.cc,v 1.4 2010/02/04 14:42:36 beischer Exp $

#include "RES_DummyField.hh"

#include "globals.hh"

RES_DummyField::RES_DummyField(G4ThreeVector fieldVector)
{
  m_fieldVector = fieldVector;
  m_displacement = G4ThreeVector(0,0,0);
  m_z0 = -4*cm;
  m_z1 =  4*cm;
}

RES_DummyField::~RES_DummyField()
{
}

void RES_DummyField::GetFieldValue(const G4double* x, G4double* B) const
{
  G4ThreeVector pos(x[0],x[1],x[2]);
  pos -= m_displacement;

  B[0] = B[1] = B[2] = 0.;
  if ( (pos.z() > m_z0) && (pos.z() < m_z1) ) {
    B[0] = m_fieldVector.x();
    B[1] = m_fieldVector.y();
    B[2] = m_fieldVector.z();
  }
}
