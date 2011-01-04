#include "RES_Layer.hh"

#include <G4RunManager.hh>

#include "RES_DetectorConstruction.hh"
#include "RES_Module.hh"

RES_Layer::RES_Layer() :
  m_firstModuleNumber(0)
{
}

RES_Layer::RES_Layer(unsigned int n, G4double z)
{
  RES_DetectorConstruction* det = (RES_DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  m_firstModuleNumber = det->GetNumberOfModules();
  for (unsigned int i = 0; i < n; i++) {
    G4double x = 0 * cm;
    G4double y = (-((int)n-1)/2. + (int)i) * 6.5*cm;
    G4ThreeVector where(x, y, z);
    det->AddModulePlacement(where);
    RES_Module* module = det->GetModule(m_firstModuleNumber + i);
    m_modules.push_back(module);
  }
}

RES_Layer::~RES_Layer()
{
}

void RES_Layer::SetType(RES_Module::ModuleType type)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetType(type);
  }
}

void RES_Layer::SetLength(G4double length)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetLength(length);
  }
}

void RES_Layer::SetAngle(G4double angle)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetAngle(angle);
  }
}

void RES_Layer::SetInternalAngle(G4double internalAngle)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetInternalAngle(internalAngle);
  }
}

void RES_Layer::SetWidth(G4double width)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetWidth(width);
  }
}

void RES_Layer::SetUpperSigmaU(G4double upperSigmaU)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetUpperSigmaU(upperSigmaU);
  }
}

void RES_Layer::SetUpperSigmaV(G4double upperSigmaV)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetUpperSigmaV(upperSigmaV);
  }
}

void RES_Layer::SetUpperSigmaZ(G4double upperSigmaZ)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetUpperSigmaZ(upperSigmaZ);
  }
}

void RES_Layer::SetLowerSigmaU(G4double lowerSigmaU)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetLowerSigmaU(lowerSigmaU);
  }
}

void RES_Layer::SetLowerSigmaV(G4double lowerSigmaV)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetLowerSigmaV(lowerSigmaV);
  }
}

void RES_Layer::SetLowerSigmaZ(G4double lowerSigmaZ)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetLowerSigmaZ(lowerSigmaZ);
  }
}

void RES_Layer::SetUpperEfficiency(G4double upperEfficiency)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetUpperEfficiency(upperEfficiency);
  }
}

void RES_Layer::SetLowerEfficiency(G4double lowerEfficiency)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetLowerEfficiency(lowerEfficiency);
  }
}

void RES_Layer::SetFoamThickness(G4double foamThickness)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetFoamThickness(foamThickness);
  }
}

void RES_Layer::SetCarbonFiberThickness(G4double carbonFiberThickness)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetCarbonFiberThickness(carbonFiberThickness);
  }
}

void RES_Layer::SetFiberThickness(G4double fiberThickness)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetFiberThickness(fiberThickness);
  }
}

void RES_Layer::SetKaptonThickness(G4double kaptonThickness)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetKaptonThickness(kaptonThickness);
  }
}

void RES_Layer::SetSiliconThickness(G4double siliconThickness)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetSiliconThickness(siliconThickness);
  }
}

void RES_Layer::SetGapSiliconThickness(G4double gapSiliconThickness)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetGapSiliconThickness(gapSiliconThickness);
  }
}

void RES_Layer::SetSubtractHoles(G4bool subtractHoles)
{
  for (std::vector<RES_Module*>::iterator it = m_modules.begin(); it != m_modules.end(); it++) {
    (*it)->SetSubtractHoles(subtractHoles);
  }
}
