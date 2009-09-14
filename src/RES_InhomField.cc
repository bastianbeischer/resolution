#include "RES_InhomField.hh"

#include "TH3D.h"
#include "TAxis.h"

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

  G4int lowerBin[3];
  G4int upperBin[3];
  G4double d[3];

  G4int globalBin = m_field_x->FindBin(x[0]/cm,x[1]/cm,x[2]/cm);
  if ( (globalBin <= 0) && (globalBin > m_field_x->GetNbinsX() * m_field_x->GetNbinsY() * m_field_x->GetNbinsZ()) )
    return;

  m_field_x->GetBinXYZ(globalBin,lowerBin[0], lowerBin[1], lowerBin[2]);

  TAxis* axis[3] = {m_field_x->GetXaxis(), m_field_x->GetYaxis(), m_field_x->GetZaxis()};
  TH3D*  hist[3] = {m_field_x, m_field_y, m_field_z};

  for (int i = 0; i < 3; i++) {
    if (x[i]/cm < axis[i]->GetBinCenter(lowerBin[i])) {
      lowerBin[i]--;
    }
    if (lowerBin[i] > 0 && lowerBin[i] < axis[i]->GetNbins()) {
      upperBin[i] = lowerBin[i] + 1;
      d[i] = (x[i]/cm - axis[i]->GetBinCenter(lowerBin[i])) / (axis[i]->GetBinCenter(upperBin[i]) - axis[i]->GetBinCenter(lowerBin[i]));
    }
    else if (lowerBin[i] == 0) {
      lowerBin[i] = 1;
      upperBin[i] = lowerBin[i];
      d[i] = 0.;
    }
    else if (lowerBin[i] == axis[i]->GetNbins()) {
      lowerBin[i] = axis[i]->GetNbins();
      upperBin[i] = lowerBin[i];
      d[i] = 0.;
    }
  }

  for (int i = 0; i < 3;i ++) {
    G4double a1 = hist[i]->GetBinContent(lowerBin[0],lowerBin[1],lowerBin[2])*(1. - d[2]) + hist[i]->GetBinContent(lowerBin[0],lowerBin[1],upperBin[2])*d[2];
    G4double a2 = hist[i]->GetBinContent(lowerBin[0],upperBin[1],lowerBin[2])*(1. - d[2]) + hist[i]->GetBinContent(lowerBin[0],upperBin[1],upperBin[2])*d[2];
    G4double b1 = hist[i]->GetBinContent(upperBin[0],lowerBin[1],lowerBin[2])*(1. - d[2]) + hist[i]->GetBinContent(upperBin[0],lowerBin[1],upperBin[2])*d[2];
    G4double b2 = hist[i]->GetBinContent(upperBin[0],upperBin[1],lowerBin[2])*(1. - d[2]) + hist[i]->GetBinContent(upperBin[0],upperBin[1],upperBin[2])*d[2];
    G4double c1 = a1*(1. - d[1]) + a2*d[1]; 
    G4double c2 = b1*(1. - d[1]) + b2*d[1];
    B[i] = (c1*(1. - d[0]) + c2*d[0])*tesla;
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
