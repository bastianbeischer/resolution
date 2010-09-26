// $Id: res_vs_mom.cc,v 1.22 2010/07/21 15:14:35 beischer Exp $

#include <iostream>
#include <cmath>

#include <MyROOTStyle.h>

#include <TApplication.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>

#include "fill_graph.hh"

double fitfunc(double* x, double* p)
{
  //  double m = 5.11e-4;
  double m = p[2];
  double beta = x[0] / sqrt(pow(x[0],2.) + pow(m,2.));
  return sqrt(pow(x[0]*p[0],2.) + pow(p[1]/beta, 2.));
}

int main(int argc, char** argv)
{
  TApplication* app = new TApplication("app", &argc, argv);

  MyROOTStyle* myStyle = new MyROOTStyle("myStyle");
  myStyle->cd();

  gStyle->SetOptFit(1111);


  app->Run();

  delete myStyle;
  delete app;

  return 1;

}
