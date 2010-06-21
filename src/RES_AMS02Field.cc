// $Id: RES_AMS02Field.cc,v 1.3 2010/06/21 20:53:15 beischer Exp $

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
  float position[3] = {x[0]/cm, x[1]/cm, x[2]/cm};
  float Bfield[3] = {0., 0., 0.};
  
  if(field)
    field->GuFld(position, Bfield);

  for(int i = 0; i < 3; i++)
    B[i] = Bfield[i];
}
