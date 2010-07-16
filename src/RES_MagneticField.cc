// $Id: RES_MagneticField.cc,v 1.2 2010/07/16 12:51:52 beischer Exp $

#include "RES_MagneticField.hh"

RES_MagneticField::RES_MagneticField() :
  G4MagneticField(),
  m_z0(0.),
  m_z1(0.),
  m_radius(0.),
  m_fieldEstimate(0.),
  m_displacement(G4ThreeVector(0.,0.,0.)),
  m_nSteps(1e4)
{
}

RES_MagneticField::~RES_MagneticField()
{
}

G4double RES_MagneticField::MeanFieldAlongTrack(G4ThreeVector startPoint, G4ThreeVector endPoint)
{
  G4ThreeVector direction = (endPoint - startPoint) / m_nSteps;

  G4double meanB[3] = {0,0,0};

  G4double maxLength = (endPoint - startPoint).mag();
  G4double currentLength = 0;
  G4ThreeVector point = startPoint;
  while (currentLength < maxLength) {
    G4double x[3] = {point.x(), point.y(), point.z()};
    G4double tempB[3];
    GetFieldValue(x, tempB);
    for(int i = 0; i < 3; i++) meanB[i] += tempB[i];

    point += direction;
    currentLength = (point - startPoint).mag();
  }

  for(int i = 0; i < 3; i++) meanB[i] /= m_nSteps;  

  return meanB[0];
}

G4bool RES_MagneticField::CheckIfTrackIsInsideMagnet(G4ThreeVector position, G4ThreeVector direction)
{
  G4double z[2] = {m_z0 + m_displacement.z(), m_z1 + m_displacement.z()};

  if (m_radius == 0.) return true;

  for (int i = 0; i < 2; i++) {
    G4double dz = z[i] - position.z();
    G4double l = dz / direction.z();
    G4ThreeVector currentPosition = position + l*direction;
    if (sqrt(pow(currentPosition.x(),2.) + pow(currentPosition.y(),2.)) > m_radius)
      return false;
  }

  return true;
}
