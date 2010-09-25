#ifndef Analysis_hh
#define Analysis_hh

class TFile;
class TTree;
class RES_Event;
class TH1D;
class TH1I;

class Analysis
{
  
public:
  Analysis();
  ~Analysis();
  
public:
  void processFile(const char* filename);
  
private:
  void openFile(const char* filename);
  void setupHistograms();
  void fillData();
  void fit();
  void draw();

private:
  unsigned int m_nHits;
  double     m_genMom;

  TFile*     m_file;
  TTree*     m_genTree;
  TTree*     m_recTree;
  RES_Event* m_genEvent;
  RES_Event* m_recEvent;

  TH1D*      m_resHist;
  TH1D*      m_initialPHist;
  TH1D*      m_ptHist;
  TH1D**     m_xDeltaGenHist;
  TH1D**     m_yDeltaGenHist;
  TH1D**     m_xDeltaSmearedHist;
  TH1D**     m_yDeltaSmearedHist;
  TH1D*      m_xTotalHist;
  TH1D*      m_yTotalHist;
  TH1D*      m_chi2Hist;
  TH1D*      m_dofHist;
  TH1D*      m_angleHist;
  TH1D*      m_lHist;
  TH1D*      m_lOverAngleHist;
  TH1I*      m_innerOrOuterHist;

};

#endif /* Analysis_hh */
