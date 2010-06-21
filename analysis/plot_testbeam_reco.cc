#include "MyROOTStyle.h"
#include <TF1.h>
#include <TGraph.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>

#include <cmath>

// double fitfunc(double* x, double* p)
// {
//   //  double m = 5.11e-4;
//   //  double beta = x[0] / sqrt(pow(x[0],2.) + pow(m,2.));
//   return sqrt(pow(x[0]*0.051,2.) + pow(p[0], 2.));///beta, 2.) + pow(0.22, 2.));
// }

int main(int argc, char** argv)
{
  TApplication app("app", &argc, argv);

  MyROOTStyle style("style");
  style.cd();
  TGraph graph("testbeam_sigmaP_vs_p.txt");

  for (int i = 0; i < graph.GetN(); i++) {
    double x = graph.GetX()[i];
    double y = graph.GetY()[i];
    double y_new = sqrt(y*y - 0.2);
    graph.SetPoint(i, x, y_new);
  }

  graph.SetTitle("#sigma_p / p vs. reconstructed momentum for testbeam");
  graph.Draw("AP");
  graph.SetMarkerStyle(23);
  graph.SetMarkerSize(1.5);
  graph.SetMarkerColor(kRed);
  graph.GetXaxis()->SetTitle("p_{rec} / GeV");
  graph.GetYaxis()->SetTitle("#sigma_{p} / p");
  gStyle->SetOptFit(0);
  // TF1* fit = new TF1("fit", fitfunc, 0, 10, 1);
  // fit->FixParameter(0, 0.3);
  // fit->Draw("SAME");
  graph.Fit("pol1");

  TF1* fit =   graph.GetFunction("pol1");
  fit->SetLineColor(kRed);
  double a = fit->GetParameter(0);
  double b = fit->GetParameter(1);

  char text[128];
  sprintf(text, "slope #approx %.2f / GeV", b);
  TLatex latex(7.0, 0.7, text);
  latex.Draw("SAME");
  sprintf(text, "offset #approx %.2f", a);
  TLatex latex2(7.0, 0.65, text);
  latex2.Draw("SAME");


  app.Run();
  return 1;
}
