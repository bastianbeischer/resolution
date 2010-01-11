#include "RES_MagFieldInfo.hh"

#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"

RES_MagFieldInfo::RES_MagFieldInfo() :
  m_nSteps(1e4)
{
  
}

RES_MagFieldInfo::~RES_MagFieldInfo()
{

}

G4double RES_MagFieldInfo::MeanFieldAlongTrack(G4ThreeVector startPoint, G4ThreeVector endPoint)
{
  G4ThreeVector direction = (endPoint - startPoint) / m_nSteps;

  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  G4MagneticField* B = (G4MagneticField*) fieldMgr->GetDetectorField();

  G4double meanB[3] = {0,0,0};

  G4double maxLength = (endPoint - startPoint).mag();
  G4double currentLength = 0;
  G4ThreeVector point = startPoint;
  while (currentLength < maxLength) {
    G4double x[3] = {point.x(), point.y(), point.z()};
    G4double tempB[3];
    B->GetFieldValue(x, tempB);
    for(int i = 0; i < 3; i++) meanB[i] += tempB[i];

    point += direction;
    currentLength = (point - startPoint).mag();
  }

  for(int i = 0; i < 3; i++) meanB[i] /= m_nSteps;  

  return meanB[0]/tesla;
}
