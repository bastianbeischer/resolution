// $Id: RES_UniformField.cc,v 1.2 2010/07/14 12:36:53 beischer Exp $

#include "RES_UniformField.hh"

#include "globals.hh"

RES_UniformField::RES_UniformField(G4ThreeVector fieldVector) :
  RES_MagneticField()
{
  m_z0 = -4.;
  m_z1 = 4.;
  m_fieldEstimate = fieldVector.x();
  m_fieldVector = fieldVector;
}

RES_UniformField::~RES_UniformField()
{
}

void RES_UniformField::GetFieldValue(const G4double* x, G4double* B) const
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
