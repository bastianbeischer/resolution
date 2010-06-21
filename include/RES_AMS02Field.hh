// $Id: RES_AMS02Field.hh,v 1.2 2010/06/21 20:53:16 beischer Exp $

#ifndef RES_AMS02Field_hh
#define RES_AMS02Field_hh

#include "G4MagneticField.hh"
#include "globals.hh"

class MagField;

class RES_AMS02Field : public G4MagneticField
{

public:
  explicit RES_AMS02Field(G4String fileName);
  ~RES_AMS02Field();

public:
  void GetFieldValue(const G4double* x, G4double* B) const;

private:
  MagField* field;

};

#endif /* RES_AMS02Field_hh */
