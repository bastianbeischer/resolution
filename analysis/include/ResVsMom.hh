#ifndef ResVsMom_hh
#define ResVsMom_hh

#include <vector>

class TCanvas;
class TLatex;
class TLegend;
class TGraphErrors;

class ResVsMom
{
  
  enum ParticleType {electron = 0, proton = 1};

public:
  ResVsMom();
  ~ResVsMom();
  
public:
  void processConfigFile(const char* filename);

private:
  void addGraph(const char* fileNameTemplate, const char* title, ParticleType type, double min, double max, double step, double guessA, double guessB);
  void draw();
  
private:
  static int m_markers[4];
  static int m_colors[4];

  TCanvas* m_canvas;
  TLatex*  m_text;
  TLegend* m_legend;

  std::vector<TGraphErrors*> m_graphs;

};

#endif /* ResVsMom_hh */
