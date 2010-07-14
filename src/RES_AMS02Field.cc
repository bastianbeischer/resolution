// $Id: RES_AMS02Field.cc,v 1.4 2010/07/14 12:17:37 beischer Exp $

#include "RES_AMS02Field.hh"

#include "MagField.hh"

RES_AMS02Field::RES_AMS02Field(G4String fileName) : 
  RES_MagneticField(),
  field(MagField::GetPtr())
{
  m_z0 = 0;
  m_z1 = 0;
  m_fieldEstimate = 0;

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
