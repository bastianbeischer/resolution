// $Id: RES_InhomField.cc,v 1.16 2010/07/14 12:36:52 beischer Exp $

#include "RES_InhomField.hh"

#include "RES_Axis.hh"

#include <fstream>

RES_InhomField::RES_InhomField(G4String dataFileName) : 
  RES_MagneticField(),
  m_dataFileName(dataFileName),
  m_axis_x(0),
  m_axis_y(0),
  m_axis_z(0),
  m_field_x(0),
  m_field_y(0),
  m_field_z(0)
{
  ReadData();
}

RES_InhomField::~RES_InhomField()
{
  for (int i = 0; i < m_axis_x->GetNbins(); i++) {
    for (int j = 0; j < m_axis_y->GetNbins(); j++) {
      delete [] m_field_x[i][j];
      delete [] m_field_y[i][j];
      delete [] m_field_z[i][j];
    }
    delete [] m_field_x[i];
    delete [] m_field_y[i];
    delete [] m_field_z[i];
  }
  delete [] m_field_x;
  delete [] m_field_y;
  delete [] m_field_z;


  delete m_axis_x;
  delete m_axis_y;
  delete m_axis_z;
}

void RES_InhomField::GetFieldValue(const G4double* x, G4double* B) const
{
  G4double pos[3] =  {x[0],x[1],x[2]};
  pos[0] -= m_displacement.x();
  pos[1] -= m_displacement.y();
  pos[2] -= m_displacement.z();

  B[0] = 0.;
  B[1] = 0.;
  B[2] = 0.;

  G4int lowerBin[3] = {0,0,0};
  G4int upperBin[3] = {0,0,0};
  G4double d[3] = {0,0,0};
  RES_Axis* axis[3] = {m_axis_x, m_axis_y, m_axis_z};
  G4double*** values[3] = {m_field_x, m_field_y, m_field_z};

  // find the proper values for bin lower than, and above the positions on each axis
  for (int i = 0; i < 3; i++) {
    lowerBin[i] = axis[i]->GetBin(pos[i]/cm);

    // if the corresponding point is outside of the region where the magnetic field is defined, return
    if ( (lowerBin[i] < 0) || (lowerBin[i] >= axis[i]->GetNbins()) )
      return;

    // if we are in the first half of the bin, set the lowerBin to the next bin on the other side
    if (pos[i]/cm < axis[i]->GetBinCenter(lowerBin[i])) {
      lowerBin[i]--;
    }
    
    // if everything is normal, lowerBin[i] contains the number of the bin to the left of pos[i], while upperBin[i] will be set to the next bin, i.e. the one to the right
    if (lowerBin[i] >= 0 && lowerBin[i] < axis[i]->GetNbins()-1) {
      upperBin[i] = lowerBin[i] + 1;
      d[i] = (pos[i]/cm - axis[i]->GetBinCenter(lowerBin[i])) / (axis[i]->GetBinCenter(upperBin[i]) - axis[i]->GetBinCenter(lowerBin[i]));
    }

    // if lowerBin[i] is now "underflow" (this means that pos[i] is in the first half of the first bin), use the first bin for the calculation and don't interpolate
    else if (lowerBin[i] == -1) {
      lowerBin[i] = 0;
      upperBin[i] = lowerBin[i];
      d[i] = 0.;
    }

    // if lowerBin[i] is now "overflow" (this means that pos[i] is in the second half of the last bin), use the last bin for the calculation and don't interpolate
    else if (lowerBin[i] == axis[i]->GetNbins()-1) {
      lowerBin[i] = axis[i]->GetNbins() - 1;
      upperBin[i] = lowerBin[i];
      d[i] = 0.;
    }
  }

  // do the actual "trilinear interpolation"
  for (int i = 0; i < 3;i ++) {
    G4double a1 = values[i][lowerBin[0]][lowerBin[1]][lowerBin[2]]*(1. - d[2]) + values[i][lowerBin[0]][lowerBin[1]][upperBin[2]]*d[2];
    G4double a2 = values[i][lowerBin[0]][upperBin[1]][lowerBin[2]]*(1. - d[2]) + values[i][lowerBin[0]][upperBin[1]][upperBin[2]]*d[2];
    G4double b1 = values[i][upperBin[0]][lowerBin[1]][lowerBin[2]]*(1. - d[2]) + values[i][upperBin[0]][lowerBin[1]][upperBin[2]]*d[2];
    G4double b2 = values[i][upperBin[0]][upperBin[1]][lowerBin[2]]*(1. - d[2]) + values[i][upperBin[0]][upperBin[1]][upperBin[2]]*d[2];
    G4double c1 = a1*(1. - d[1]) + a2*d[1]; 
    G4double c2 = b1*(1. - d[1]) + b2*d[1];
    B[i] = (c1*(1. - d[0]) + c2*d[0])*tesla;
  }
}

void RES_InhomField::ReadData()
{
  std::ifstream file((m_dataFileName).c_str());
  if (!file.is_open()) {
    G4cerr << "Error reading file : " << m_dataFileName.c_str() << G4endl;
    return;
  }
  
  G4double disX, disY, disZ;
  file >> m_z0 >> m_z1 >> m_radius >> m_fieldEstimate >> disX >> disY >> disZ;
  
  m_z0 *= cm;
  m_z1 *= cm;
  m_radius *= cm;
  m_fieldEstimate *= tesla;
  m_displacement = G4ThreeVector(disX*cm, disY*cm, disZ*cm);

  G4double x0,x1,y0,y1,z0,z1;
  G4int nBins_x, nBins_y, nBins_z;
  G4int mirror;
  file >> x0 >> x1 >> y0 >> y1 >> z0 >> z1 >> nBins_x >> nBins_y >> nBins_z >> mirror;

  x0 /= 10.;
  x1 /= 10.;
  y0 /= 10.;
  y1 /= 10.;
  z0 /= 10.;
  z1 /= 10.;

  m_axis_x = new RES_Axis(x0,x1,nBins_x);
  m_axis_y = new RES_Axis(y0,y1,nBins_y);
  m_axis_z = new RES_Axis(z0,z1,nBins_z);

  m_field_x = new G4double**[nBins_x];
  m_field_y = new G4double**[nBins_x];
  m_field_z = new G4double**[nBins_x];
  for (int i = 0; i < nBins_x; i++) {
    m_field_x[i] = new G4double*[nBins_y];
    m_field_y[i] = new G4double*[nBins_y];
    m_field_z[i] = new G4double*[nBins_y];
    for (int j = 0; j < nBins_y; j++) {
      m_field_x[i][j] = new G4double[nBins_z];
      m_field_y[i][j] = new G4double[nBins_z];
      m_field_z[i][j] = new G4double[nBins_z];
      for (int k = 0; k < nBins_z; k++) {
        m_field_x[i][j][k] = 0.;
        m_field_y[i][j][k] = 0.;        
        m_field_z[i][j][k] = 0.;
      }
   } 
  }

  // read the data (only contains all x, y >= 0, z >= 0)
  G4double x, y, z, f_x, f_y, f_z, dummy;
  G4int nX, nY, nZ;
  while (true) {
    file >> x >> y >> z >> dummy >> f_x >> f_y >> f_z >> dummy;
    if (file.eof())
      break;
    
    x/=10.;
    y/=10.;
    z/=10.;

    //    file >> x >> y >> z >> f_x >> f_y >> f_z;

    nX = m_axis_x->GetBin(x);
    nY = m_axis_y->GetBin(y);
    nZ = m_axis_z->GetBin(z);
    m_field_x[nX][nY][nZ] = f_x;
    m_field_y[nX][nY][nZ] = f_y;
    m_field_z[nX][nY][nZ] = f_z;

    // mirror to get the full world
    if (mirror) {
      if (y != 0.0) {
        nX = m_axis_x->GetBin(x);
        nY = m_axis_y->GetBin(-y);
        nZ = m_axis_z->GetBin(z);
        m_field_x[nX][nY][nZ] = f_x;
        m_field_y[nX][nY][nZ] = -f_y;
        m_field_z[nX][nY][nZ] = f_z;
      }
      if (z != 0.0) {
        nX = m_axis_x->GetBin(x);
        nY = m_axis_y->GetBin(y);
        nZ = m_axis_z->GetBin(-z);
        m_field_x[nX][nY][nZ] = f_x;
        m_field_y[nX][nY][nZ] = f_y;
        m_field_z[nX][nY][nZ] = -f_z;
      }
      if ( (y != 0.0) && (z != 0) ) {
        nX = m_axis_x->GetBin(x);
        nY = m_axis_y->GetBin(-y);
        nZ = m_axis_z->GetBin(-z);
        m_field_x[nX][nY][nZ] = f_x;
        m_field_y[nX][nY][nZ] = -f_y;
        m_field_z[nX][nY][nZ] = -f_z;
      }
    }
  }
}
