#ifndef ResVsMom_hh
#define ResVsMom_hh

#include <vector>

class TGraphErrors;

class ResVsMom
{
  
public:
  ResVsMom();
  ~ResVsMom();
  
public:
  void processConfigFile();

private:
  void addGraph(const char* fileNameTemplate,double min,double max,double step,double guessA,double guessB);
  void draw();
  
private:
  static int m_markers[4];
  static int m_colors[4];

  std::vector<TGraphErrors*> m_graphs;

};

#endif /* ResVsMom_hh */
