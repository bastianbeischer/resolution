#ifndef RES_ApplicationManager_hh
#define RES_ApplicationManager_hh

#include "globals.hh"

class RES_ApplicationMessenger;
class RES_RunManager;

class RES_ApplicationManager
{

public:
  enum SessionType {Terminal=0,
#ifdef G4UI_USE_QT                    
                    Qt=1
#endif
  };

public:
  RES_ApplicationManager(int argc, char** argv);
  ~RES_ApplicationManager();

public:
  int RunBatchScript(G4String scriptName);
  void CreateSession(SessionType = Terminal, G4String = "");
  void SetSeedToSystemTime();

private:
  RES_ApplicationMessenger* m_messenger;
  RES_RunManager*           m_runManager;

  int                       m_argc;
  char**                    m_argv;

};

#endif /* RES_ApplicationManager_hh */
