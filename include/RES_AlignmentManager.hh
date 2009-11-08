// $Id: RES_AlignmentManager.hh,v 1.3 2009/11/08 17:09:11 beischer Exp $

#ifndef RES_AlignmentManager_hh
#define RES_AlignmentManager_hh

#include <map>
#include "globals.hh"

class RES_AlignmentMessenger;
class RES_DataHandler;

class RES_AlignmentManager
{

public:
  ~RES_AlignmentManager();
  
public:
  static RES_AlignmentManager* GetInstance();

  void StartAlignment();
  void SetXshift(unsigned int i, float shift) {m_xShifts[i] = shift;}
  void SetYshift(unsigned int i, float shift) {m_yShifts[i] = shift;}
  void SetVerbose(G4int verbose) {m_verbose = verbose;}

  G4float GetXshift(unsigned int i) {return m_xShifts[i];}
  G4float GetYshift(unsigned int i) {return m_yShifts[i];}
  G4float GetAngleShift(unsigned int i) {return m_angleShifts[i];}

private:
  RES_AlignmentManager();

private:
  static RES_AlignmentManager*  m_instance;

  RES_AlignmentMessenger*  m_messenger;
  RES_DataHandler*         m_dataHandler;

  std::map<G4int, G4float> m_xShifts;
  std::map<G4int, G4float> m_yShifts;
  std::map<G4int, G4float> m_angleShifts;
  G4float*                 m_parameters;
  G4int                    m_verbose;

};

#endif /* RES_AlignmentManager_hh */
