// $Id: RES_Axis.cc,v 1.3 2009/10/14 09:24:26 beischer Exp $

#include "RES_Axis.hh"

RES_Axis::RES_Axis(G4double x0, G4double x1, G4int nBins)
{
  this->x0 = x0;
  this->x1 = x1;
  this->nBins = nBins;
}

RES_Axis::~RES_Axis()
{
}

G4double RES_Axis::GetBinCenter(G4int binX) 
{
  return x0 + (binX + 0.5) * (x1-x0)/((G4double)nBins);
}

G4int RES_Axis::GetBin(G4double x) 
{
  G4int value = (G4int) floor((x-x0)/(x1-x0) * nBins);
  if (value < 0) return -1;
  if (value > nBins) return nBins;
  return value;
}
