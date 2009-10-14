// $Id: RES_Axis.hh,v 1.2 2009/10/14 09:24:33 beischer Exp $

#ifndef RES_Axis_hh
#define RES_Axis_hh

#include "globals.hh"

class RES_Axis
{

public:
  explicit RES_Axis(G4double x0, G4double x1, G4int nBins);
  ~RES_Axis();

  G4double GetBinCenter(G4int binX);
  G4int GetBin(G4double x);
  G4int GetNbins() {return nBins;}

private:
  G4double x0,x1;
  G4int nBins;

};

#endif /* RES_Axis_hh */
