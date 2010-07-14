// $Id: RES_UniformField.hh,v 1.1 2010/07/14 12:17:36 beischer Exp $

#ifndef RES_UniformField_hh
#define RES_UniformField_hh

#include "RES_MagneticField.hh"

#include "G4ThreeVector.hh"

class RES_UniformField : public RES_MagneticField
{

public:
  explicit RES_UniformField(G4ThreeVector fieldVector);
  ~RES_UniformField();

public:
  void GetFieldValue(const G4double* x, G4double* B) const;

private:
  G4ThreeVector m_fieldVector;

};

#endif /* RES_UniformField_hh */
