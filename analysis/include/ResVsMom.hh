#ifndef ResVsMom_hh
#define ResVsMom_hh

class TGraphErrors;

class ResVsMom
{
  
public:
  ResVsMom();
  ~ResVsMom();
  
public:
  void process();

private:
  void fillGraph(TGraphErrors& graph1,const char* fileNameTemplate,double min,double max,double step,double guessA,double guessB);
  
};

#endif /* ResVsMom_hh */
