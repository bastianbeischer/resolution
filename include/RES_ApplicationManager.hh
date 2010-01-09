#ifndef RES_ApplicationManager_hh
#define RES_ApplicationManager_hh

class G4String;
class RES_ApplicationMessenger;

class RES_ApplicationManager
{

public:
  RES_ApplicationManager();
  ~RES_ApplicationManager();

public:
  int RunBatchScript(G4String scriptName);
  void CreateSession();

private:
  RES_ApplicationMessenger* m_messenger;

};

#endif /* RES_ApplicationManager_hh */