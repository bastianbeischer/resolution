#ifndef RES_Layer_hh
#define RES_Layer_hh

#include "globals.hh"

#include <vector>
#include "RES_Module.hh"

class RES_Layer
{
  
public:
  RES_Layer();
  RES_Layer(unsigned int n, G4double z);
  ~RES_Layer();
  
  void SetType(RES_Module::ModuleType type);
  void SetLength(G4double length);
  void SetAngle(G4double angle);
  void SetInternalAngle(G4double internalAngle);
  void SetWidth(G4double width);
  void SetUpperSigmaU(G4double upperSigmaU);
  void SetUpperSigmaV(G4double upperSigmaV);
  void SetUpperSigmaZ(G4double upperSigmaZ);
  void SetLowerSigmaU(G4double lowerSigmaU);
  void SetLowerSigmaV(G4double lowerSigmaV);
  void SetLowerSigmaZ(G4double lowerSigmaZ);
  void SetUpperEfficiency(G4double upperEfficiency);
  void SetLowerEfficiency(G4double lowerEfficiency);
  void SetFoamThickness(G4double foamThickness);
  void SetCarbonFiberThickness(G4double carbonFiberThickness);
  void SetFiberThickness(G4double fiberThickness);
  void SetKaptonThickness(G4double kaptonThickness);
  void SetSiliconThickness(G4double siliconThickness);
  void SetGapSiliconThickness(G4double gapSiliconThickness);
  void SetSubtractHoles(G4bool subtractHoles);

private:
  std::vector<RES_Module*> m_modules;
  unsigned int m_firstModuleNumber;
  
};

#endif /* RES_Layer_hh */
