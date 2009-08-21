#ifndef RES_RunManager_hh
#define RES_RunManager_hh

#include "G4RunManager.hh"

class RES_DataHandler;

class RES_RunManager : public G4RunManager
{

public:
  RES_RunManager();
  ~RES_RunManager();

public:
  void TestDataHandler();

private:
  RES_DataHandler* m_dataHandler;

};

#endif /* RES_RunManager_hh */
