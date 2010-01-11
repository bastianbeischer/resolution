#ifndef RES_MagFieldInfo_hh
#define RES_MagFieldInfo_hh

#include "globals.hh"
#include "G4ThreeVector.hh"

class RES_MagFieldInfo
{

public:
  RES_MagFieldInfo();
  ~RES_MagFieldInfo();

public:
  G4double MeanFieldAlongTrack(G4ThreeVector startPoint, G4ThreeVector endPoint);

private:
  int m_nSteps;

};

#endif /* RES_MagFieldInfo_hh */
