#ifndef Analysis_hh
#define Analysis_hh

class TFile;
class TTree;
class RES_Event;
class TH1D;

class Analysis
{
  
public:
  Analysis();
  ~Analysis();
  
public:
  void processFile(const char* filename);
  
private:
  double MS(double p, double m, double L, double X0);
  void   fitInChi2Range(TH1D* hist, double& sigma, double& error);

private:
  TFile*     m_file;
  TTree*     m_genTree;
  TTree*     m_recTree;
  RES_Event* m_genEvent;
  RES_Event* m_recEvent;

};

#endif /* Analysis_hh */
