// $Id: RES_MagneticField.hh,v 1.2 2010/07/14 13:57:01 beischer Exp $

#ifndef RES_MagneticField_hh
#define RES_MagneticField_hh

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"

class RES_MagneticField : public G4MagneticField
{
  
public:
  RES_MagneticField();
  ~RES_MagneticField();
  
public:
  void SetDisplacement(G4ThreeVector displacement) {m_displacement = displacement;}

  G4double      GetZ0()            const {return m_z0;}
  G4double      GetZ1()            const {return m_z1;}
  G4double      GetFieldEstimate() const {return m_fieldEstimate;}
  G4ThreeVector GetDisplacement()  const {return m_displacement;}

  void SetZ0(G4double z0)                  {m_z0 = z0;}
  void SetZ1(G4double z1)                  {m_z1 = z1;}
  void SetFieldEstimate(G4double estimate) {m_fieldEstimate = estimate;}

  G4double MeanFieldAlongTrack(G4ThreeVector startPoint, G4ThreeVector endPoint);
  G4bool   CheckIfTrackIsInsideMagnet(G4ThreeVector position, G4ThreeVector direction);

protected:
  G4double      m_z0;
  G4double      m_z1;
  G4double      m_radius;
  G4double      m_fieldEstimate;
  G4ThreeVector m_displacement;
  G4int         m_nSteps;

};

#endif /* RES_MagneticField_hh */
