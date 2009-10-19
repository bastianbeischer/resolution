// $Id: RES_AlignmentManager.hh,v 1.2 2009/10/19 12:27:45 beischer Exp $

#ifndef RES_AlignmentManager_hh
#define RES_AlignmentManager_hh

#include "globals.hh"

class RES_AlignmentMessenger;
class RES_DataHandler;

class RES_AlignmentManager
{

public:
  RES_AlignmentManager();
  ~RES_AlignmentManager();
  
public:
  void StartAlignment();
  void SetVerbose(G4int verbose) {m_verbose = verbose;}

private:

private:
  RES_AlignmentMessenger* m_messenger;
  RES_DataHandler*        m_dataHandler;

  G4float*                m_parameters;
  G4int                   m_verbose;

};

#endif /* RES_AlignmentManager_hh */
