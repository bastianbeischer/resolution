// $Id: RES_InhomField.hh,v 1.4 2010/02/04 14:42:36 beischer Exp $

#ifndef RES_InhomField_hh
#define RES_InhomField_hh

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class RES_Axis;

class RES_InhomField : public G4MagneticField
{

public:
  explicit RES_InhomField(G4String dataFileName);
  ~RES_InhomField();

public:
  void GetFieldValue(const G4double* x, G4double* B) const;

  void SetDisplacement(G4ThreeVector displacement) {m_displacement = displacement;}

private:
  void ReadData();

public:
  G4String      m_dataFileName;
  G4ThreeVector m_displacement;

  RES_Axis*     m_axis_x;
  RES_Axis*     m_axis_y;
  RES_Axis*     m_axis_z;
  G4double***   m_field_x;
  G4double***   m_field_y;
  G4double***   m_field_z;
  
};

#endif /* RES_InhomField_hh */
