#ifndef RES_PrimaryGeneratorAction_hh
#define RES_PrimaryGeneratorAction_hh

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class RES_PrimaryGeneratorMessenger;

class RES_PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{

public:
  RES_PrimaryGeneratorAction();
  ~RES_PrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* event);
  G4ParticleGun* GetParticleGun() const {return m_particleGun;}
  
public:
  inline G4bool GetRandomOrigin() {return m_randomOrigin;}
  inline G4bool GetRandomDirection() {return m_randomDirection;}
  inline void SetRandomOrigin(G4bool value) {m_randomOrigin = value;}
  inline void SetRandomDirection(G4bool value) {m_randomDirection = value;}

private:
  RES_PrimaryGeneratorMessenger* m_messenger;

  G4ParticleGun*                 m_particleGun;

  G4bool                         m_randomOrigin;
  G4bool                         m_randomDirection;
  G4double                       m_z_start;
};

#endif /* RES_PrimaryGeneratorAction_hh */
