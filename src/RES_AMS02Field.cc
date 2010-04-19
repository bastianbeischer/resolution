#include "RES_AMS02Field.hh"

#include "MagField.hh"

RES_AMS02Field::RES_AMS02Field(G4String fileName) : 
  G4MagneticField(),
  field(MagField::GetPtr())
{
  //set scale because fld07 map is at 460 A and we want 400 A
  field->SetScale(400./460.);
  field->Read(fileName);
  field->SetMagstat(1);
}

RES_AMS02Field::~RES_AMS02Field()
{
  delete field;
}

void RES_AMS02Field::GetFieldValue(const G4double* x, G4double* B) const
{
  if(field)
    field->GuFld((float*)x,(float*)B);
}
