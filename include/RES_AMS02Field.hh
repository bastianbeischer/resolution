// $Id: RES_AMS02Field.hh,v 1.3 2010/07/14 12:17:34 beischer Exp $

#ifndef RES_AMS02Field_hh
#define RES_AMS02Field_hh

#include "RES_MagneticField.hh"
#include "globals.hh"

class MagField;

class RES_AMS02Field : public RES_MagneticField
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
