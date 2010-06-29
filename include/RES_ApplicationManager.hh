#ifndef RES_ApplicationManager_hh
#define RES_ApplicationManager_hh

class G4String;
class RES_ApplicationMessenger;
class RES_RunManager;

class RES_ApplicationManager
{

public:
  RES_ApplicationManager(int argc, char** argv);
  ~RES_ApplicationManager();

public:
  int RunBatchScript(G4String scriptName);
  void CreateSession();

private:
  RES_ApplicationMessenger* m_messenger;
  RES_RunManager*           m_runManager;

  int                       m_argc;
  char**                    m_argv;

};

#endif /* RES_ApplicationManager_hh */
