// $Id: res_vs_mom_pebs.cc,v 1.1 2010/04/25 19:25:27 beischer Exp $

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

#include "fill_graph.hh"

double fitfunc(double* x, double* p)
{
  return sqrt(pow(x[0]*p[0],2.) + pow(p[1], 2.));
}

int main(int argc, char** argv)
{
  TApplication* app = new TApplication("app", &argc, argv);

  MyROOTStyle* myStyle = new MyROOTStyle("myStyle");
  myStyle->cd();

  gStyle->SetOptFit(111);

  TGraphErrors graph1;
  graph1.SetMarkerStyle(23);
  graph1.SetMarkerColor(kRed);
  graph1.SetTitle("Momentum resolution for PEBS01 at 1^{o} stereo angle");

  int momMin = 10;
  int momMax = 1200;
  int momStep = 10;
  fillGraph(graph1, "../results/pebs01_%03d_GeV.root", momMin, momMax, momStep, 0.8e-3, 4e-3);

  TF1 fitMSC("fitMSC", fitfunc, momMin, momMax, 2);
  fitMSC.SetParNames("a", "b");

  TCanvas canvas("canvas", "canvas", 1024, 768);
  canvas.Draw();
  canvas.cd();
  canvas.SetGridx();
  canvas.SetGridy();
  graph1.Draw("AP");
  graph1.Fit("fitMSC", "E");
  graph1.SetMarkerSize(1.5);
  graph1.GetFunction("fitMSC")->SetLineColor(kRed);
  graph1.PaintStats(graph1.GetFunction("fitMSC"));
  graph1.GetXaxis()->SetTitle("p / GeV");
  graph1.GetYaxis()->SetTitle("#sigma_{p} / p");
  TPaveStats* pt1 = (TPaveStats*) graph1.GetListOfFunctions()->FindObject("stats");
  pt1->SetTextColor(kRed);
  pt1->SetX1NDC(0.6);
  pt1->SetX2NDC(0.95);
  pt1->SetY1NDC(0.3);
  pt1->SetY2NDC(0.45);

  TLatex text(120, 0.7, "#sigma_{p} / p = #sqrt{(ap)^{2} + b^{2}}");
  text.Draw("SAME");

  canvas.SaveAs("pebs01_1_deg.pdf");
  canvas.SaveAs("pebs01_1_deg.root");

  app->Run();

  delete myStyle;
  delete app;

  return 1;

}
