// $Id: res_vs_mom.cc,v 1.12 2010/04/24 19:17:53 beischer Exp $

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
  graph1.SetTitle("Momentum resolution for PERDaix at 1^{o} stereo angle");

  TGraphErrors graph2;
  graph2.SetMarkerStyle(22);
  graph2.SetMarkerColor(kBlue);

  TGraphErrors graph3;
  graph3.SetMarkerStyle(29);
  graph3.SetMarkerColor(kBlack);

  TGraphErrors graph4;
  graph4.SetMarkerStyle(27);
  graph4.SetMarkerColor(kViolet);

  double momMin = 0.25;
  double momMax = 9.0;
  double momStep = 0.25;
  fillGraph(graph1, "../results/perdaix_%.2f_GeV_1.00_deg_msc_inhom.root", momMin, momMax, momStep);
  fillGraph(graph2, "../results/perdaix_%.2f_GeV_1.00_deg_msc_hom.root", momMin, momMax, momStep);
  fillGraph(graph3, "../results/perdaix_%.2f_GeV_1.00_deg_nomsc_inhom.root", momMin, momMax, momStep, 0.08, 0.0);
  fillGraph(graph4, "../results/perdaix_%.2f_GeV_1.00_deg_nomsc_hom.root", momMin, momMax, momStep, 0.08, 0.0);

  TF1 fitMSC("fitMSC", fitfunc, 0., 10., 2);
  fitMSC.SetParNames("a", "b");
  TF1 fitNoMSC("fitNoMSC", "[0]*x + [1]", 0., 10.);
  fitNoMSC.SetParNames("a", "b");

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
  graph2.Draw("P");
  graph2.Fit("fitMSC", "E");
  graph2.SetMarkerSize(1.5);
  graph2.GetFunction("fitMSC")->SetLineColor(kBlue);
  graph2.PaintStats(graph2.GetFunction("fitMSC"));
  graph3.Draw("P");
  graph3.Fit("fitNoMSC", "E");
  graph3.SetMarkerSize(1.5);
  graph3.GetFunction("fitNoMSC")->SetLineColor(kBlack);
  graph3.PaintStats(graph3.GetFunction("fitNoMSC"));
  graph4.Draw("P");
  graph4.Fit("fitNoMSC", "E");
  graph4.SetMarkerSize(1.5);
  graph4.GetFunction("fitNoMSC")->SetLineColor(kViolet);
  graph4.PaintStats(graph4.GetFunction("fitNoMSC"));

  graph1.GetXaxis()->SetTitle("p / GeV");
  graph1.GetYaxis()->SetTitle("#sigma_{p} / p");
  graph1.GetYaxis()->SetRangeUser(0.0, 0.8);

  TPaveStats* pt1 = (TPaveStats*) graph1.GetListOfFunctions()->FindObject("stats");
  pt1->SetTextColor(kRed);
  TPaveStats* pt2 = (TPaveStats*) graph2.GetListOfFunctions()->FindObject("stats");
  pt2->SetTextColor(kBlue);
  TPaveStats* pt3 = (TPaveStats*) graph3.GetListOfFunctions()->FindObject("stats");
  pt3->SetTextColor(kBlack);
  TPaveStats* pt4 = (TPaveStats*) graph4.GetListOfFunctions()->FindObject("stats");
  pt4->SetTextColor(kViolet);
  TPaveStats* pt[4] = {pt1, pt2, pt3, pt4};
  for (int i = 0; i < 4; i++) {
    pt[i]->SetX1NDC(0.7);
    pt[i]->SetX2NDC(0.95);
    pt[i]->SetY1NDC(0.5 - i*0.12);
    pt[i]->SetY2NDC(0.6 - i*0.12);
  }

  TLatex text(3.1, 0.13, "#sigma_{p} / p = #sqrt{(ap)^{2} + b^{2}}");
  text.Draw("SAME");

  TLegend legend(0.12, 0.68, 0.6, 0.88);
  legend.AddEntry(&graph1, "Multiple Scattering activated + Measured inhom. field", "P");
  legend.AddEntry(&graph2, "Multiple Scattering activated + Hom. approximation in center", "P");
  legend.AddEntry(&graph3, "Multiple Scattering deactivated + Measured inhom. field", "P");
  legend.AddEntry(&graph4, "Multiple Scattering deactivated + Hom. approximation in center", "P");
  legend.Draw("SAME");

  canvas.SaveAs("perdaix_1_deg.pdf");
  canvas.SaveAs("perdaix_1_deg.root");

  app->Run();

  delete myStyle;
  delete app;

  return 1;

}
