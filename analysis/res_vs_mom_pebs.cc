// $Id: res_vs_mom_pebs.cc,v 1.3 2010/05/12 01:56:30 beischer Exp $

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

  TGraphErrors graph2;
  graph2.SetMarkerStyle(22);
  graph2.SetMarkerColor(kBlue);

  double momMin = 10;
  double momMax = 1200;
  double momStep = 10;
  fillGraph(graph1, "../results/pebs01_%03.0f_GeV.root", momMin, momMax, momStep, 0.8e-3, 4e-3);
  fillGraph(graph2, "../results/pebs01_%03.0f_GeV_protons.root", momMin, momMax, momStep, 0.8e-3, 4e-3);

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
  graph2.Draw("P");
  graph2.Fit("fitMSC", "E");
  graph2.SetMarkerSize(1.5);
  graph2.GetFunction("fitMSC")->SetLineColor(kBlue);
  graph2.PaintStats(graph2.GetFunction("fitMSC"));

  TLegend legend(0.22, 0.75, 0.4, 0.83);
  legend.AddEntry(&graph1, "electrons", "P");
  legend.AddEntry(&graph2, "protons", "P");
  legend.Draw("SAME");

  TPaveStats* pt1 = (TPaveStats*) graph1.GetListOfFunctions()->FindObject("stats");
  pt1->SetTextColor(kRed);
  TPaveStats* pt2 = (TPaveStats*) graph2.GetListOfFunctions()->FindObject("stats");
  pt2->SetTextColor(kBlue);

  TPaveText* pt[2] = {pt1, pt2};
  for(int i = 0; i < 2; i++) {
    pt[i]->SetX1NDC(0.55);
    pt[i]->SetX2NDC(0.80);
    pt[i]->SetY1NDC(0.35 - i*0.2);
    pt[i]->SetY2NDC(0.50 - i*0.2);
  }

  double MDR = graph1.GetFunction("fitMSC")->GetX(1.0);

  TLatex text(120, 0.7, "#sigma_{p} / p = #sqrt{(ap)^{2} + (b/#beta)^{2}}");
  text.Draw("SAME");

  char MDRtext[128];
  sprintf(MDRtext, "MDR #approx %.0f GV", ((int)MDR/10)*10.);
  TLatex text2(120, 0.55, MDRtext);
  text2.Draw("SAME");

  TLine yLine(0.0, 1.0, MDR, 1.0);
  yLine.Draw();
  yLine.SetLineColor(kBlack);
  yLine.SetLineStyle(2);
  yLine.SetLineWidth(3);


  TLine xLine(MDR, 0.0, MDR, 1.0);
  xLine.Draw();
  xLine.SetLineColor(kBlack);
  xLine.SetLineStyle(2);
  xLine.SetLineWidth(3);

  canvas.SaveAs("pebs01_1_deg.pdf");
  canvas.SaveAs("pebs01_1_deg.root");

  app->Run();

  delete myStyle;
  delete app;

  return 1;

}
