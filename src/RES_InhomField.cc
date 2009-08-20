#include "RES_InhomField.hh"

#include "TH3D.h"

#include <fstream>
#include "globals.hh"

RES_InhomField::RES_InhomField(G4String dataFileName) :
  m_dataFileName(dataFileName)
{
  ReadData();
}

RES_InhomField::~RES_InhomField()
{
  delete m_field_x;
  delete m_field_y;
  delete m_field_z;
}

void RES_InhomField::GetFieldValue(const G4double* x, G4double* B) const
{
  B[0] = 0.;
  B[1] = 0.;
  B[2] = 0.;

  // find the bin and store the data in B[0], B[1], B[2]
  int bin = m_field_x->FindBin(x[0]/cm,x[1]/cm,x[2]/cm);
  if ( (bin > 0) && (bin <= m_field_x->GetNbinsX() * m_field_x->GetNbinsY() * m_field_x->GetNbinsZ()) ) {
    B[0] = m_field_x->GetBinContent(bin);
    B[1] = m_field_y->GetBinContent(bin);
    B[2] = m_field_z->GetBinContent(bin);
  }
}

void RES_InhomField::ReadData()
{
  ifstream file((m_dataFileName).c_str());
  if (!file.is_open()) {
    G4cerr << "Error reading file : " << m_dataFileName.c_str() << G4endl;
    return;
  }
  
  double x0,x1,y0,y1,z0,z1,nBins_x,nBins_y,nBins_z;
  file >> x0 >> x1 >> y0 >> y1 >> z0 >> z1 >> nBins_x >> nBins_y >> nBins_z;
  this->x0 = x0;
  this->x1 = x1;
  this->y0 = y0;
  this->y1 = y1;
  this->z0 = z0;
  this->z1 = z1;
  this->nBins_x = nBins_x;
  this->nBins_y = nBins_y;
  this->nBins_z = nBins_z;

  m_field_x   = new TH3D((m_dataFileName + "_field_x").c_str(), "x component", nBins_x, x0, x1, nBins_y, y0, y1, nBins_z, z0, z1);
  m_field_y   = new TH3D((m_dataFileName + "_field_y").c_str(), "y component", nBins_x, x0, x1, nBins_y, y0, y1, nBins_z, z0, z1);
  m_field_z   = new TH3D((m_dataFileName + "_field_z").c_str(), "z component", nBins_x, x0, x1, nBins_y, y0, y1, nBins_z, z0, z1);

  // read the data (only contains all x, y >= 0, z >= 0)
  double x, y, z, f_x, f_y, f_z;
  while (!file.eof()) {
    file >> x >> y >> z >> f_x >> f_y >> f_z;
    m_field_x  ->Fill(x,y,z, f_x);
    m_field_y  ->Fill(x,y,z, f_y);
    m_field_z  ->Fill(x,y,z, f_z);
    // mirror to get the full world
    if (y != 0.0) {
      m_field_x  ->Fill(x,-y,z,  f_x);
      m_field_y  ->Fill(x,-y,z, -f_y);
      m_field_z  ->Fill(x,-y,z,  f_z);
    }
    if (z != 0.0) {
      m_field_x  ->Fill(x,y,-z,  f_x);
      m_field_y  ->Fill(x,y,-z,  f_y);
      m_field_z  ->Fill(x,y,-z, -f_z);
    }
    if ( (y != 0.0) && (z != 0) ) {
      m_field_x  ->Fill(x,-y,-z,  f_x);
      m_field_y  ->Fill(x,-y,-z, -f_y);
      m_field_z  ->Fill(x,-y,-z, -f_z);
    }
  }

}
