#ifndef RES_InhomField_hh
#define RES_InhomField_hh

#include "G4MagneticField.hh"
#include "globals.hh"

class TH3D;
class RES_InhomFieldMessenger;

class RES_InhomField : public G4MagneticField
{

public:
  explicit RES_InhomField(G4String dataFileName);
  ~RES_InhomField();

public:
  void GetFieldValue(const G4double* x, G4double* B) const;

private:
  void ReadData();

public:
  G4String  m_dataFileName;

  G4int nBins_x, nBins_y, nBins_z;           // bins for the histograms
  G4double x0, x1, y0, y1, z0, z1, R;        // dimensions

  TH3D*     m_field_x;
  TH3D*     m_field_y;
  TH3D*     m_field_z;

};

#endif /* RES_InhomField_hh */
