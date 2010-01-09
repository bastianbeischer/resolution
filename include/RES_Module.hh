#ifndef RES_Module_hh
#define RES_Module_hh

#include "G4ThreeVector.hh"
#include "globals.hh"

class RES_Module
{

public:
  enum ModuleType {fiber, silicon};

public:
  RES_Module();
  RES_Module(G4ThreeVector placement);
  ~RES_Module();

  void SetType(ModuleType type)                 {type == fiber ? SetDefaultValuesForFiber() : SetDefaultValuesForSilicon();}
  void SetAngle(G4double angle)                 {m_angle = angle;}
  void SetInternalAngle(G4double internalAngle) {m_internalAngle = internalAngle;}
  void SetLength(G4double length)               {m_length = length;}
  void SetWidth(G4double width)                 {m_width = width;}
  void SetUpperSigmaU(G4double upperSigmaU)     {m_upperSigmaU = upperSigmaU;}
  void SetUpperSigmaV(G4double upperSigmaV)     {m_upperSigmaV = upperSigmaV;}
  void SetUpperSigmaZ(G4double upperSigmaZ)     {m_upperSigmaZ = upperSigmaZ;}
  void SetLowerSigmaU(G4double lowerSigmaU)     {m_lowerSigmaU = lowerSigmaU;}
  void SetLowerSigmaV(G4double lowerSigmaV)     {m_lowerSigmaV = lowerSigmaV;}
  void SetLowerSigmaZ(G4double lowerSigmaZ)     {m_lowerSigmaZ = lowerSigmaZ;}
  void SetEfficiency(G4double efficiency)       {m_efficiency = efficiency;}

  ModuleType    GetType()          {return m_type;}
  G4ThreeVector GetPlacement()     {return m_placement;}
  G4double      GetAngle()         {return m_angle;}
  G4double      GetInternalAngle() {return m_internalAngle;}
  G4double      GetLength()        {return m_length;}
  G4double      GetWidth()         {return m_width;}
  G4double      GetUpperSigmaU()   {return m_upperSigmaU;}
  G4double      GetUpperSigmaV()   {return m_upperSigmaV;}
  G4double      GetUpperSigmaZ()   {return m_upperSigmaZ;}
  G4double      GetLowerSigmaU()   {return m_lowerSigmaU;}
  G4double      GetLowerSigmaV()   {return m_lowerSigmaV;}
  G4double      GetLowerSigmaZ()   {return m_lowerSigmaZ;}
  G4double      GetEfficiency()    {return m_efficiency;}

private:
  void SetDefaultValuesForFiber();
  void SetDefaultValuesForSilicon();

private:
  ModuleType    m_type;
  G4ThreeVector m_placement;
  G4double      m_angle;
  G4double      m_internalAngle;
  G4double      m_length;
  G4double      m_width;
  G4double      m_upperSigmaU;
  G4double      m_upperSigmaV;
  G4double      m_upperSigmaZ;
  G4double      m_lowerSigmaU;
  G4double      m_lowerSigmaV;
  G4double      m_lowerSigmaZ;
  G4double      m_efficiency;

};

#endif /* RES_Module_hh */
