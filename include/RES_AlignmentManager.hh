// $Id: RES_AlignmentManager.hh,v 1.1 2009/10/14 16:51:36 beischer Exp $

#ifndef RES_AlignmentManager_hh
#define RES_AlignmentManager_hh

class RES_AlignmentMessenger;
class RES_DataHandler;

class RES_AlignmentManager
{

public:
  RES_AlignmentManager();
  ~RES_AlignmentManager();
  
public:
  void StartAlignment();

private:

private:
  RES_AlignmentMessenger* m_messenger;
  RES_DataHandler*        m_dataHandler;

  float*                  m_parameters;

};

#endif /* RES_AlignmentManager_hh */
