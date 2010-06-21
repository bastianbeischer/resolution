// $Id: res_vs_angle.cc,v 1.10 2010/06/21 12:20:40 beischer Exp $

#include <iostream>
#include <cmath>

#include <MyROOTStyle.h>

#include "RES_Event.hh"

#include <TApplication.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "fill_graph.hh"

int main(int argc, char** argv)
{
  TApplication* app = new TApplication("app", &argc, argv);

  MyROOTStyle* myStyle = new MyROOTStyle("myStyle");
  myStyle->cd();

  TGraphErrors graph1;
  graph1.SetMarkerStyle(23);
  graph1.SetMarkerColor(kRed);
  graph1.SetTitle("Momentum resolution for PERDaix for 5 GeV electrons");


  TGraphErrors graph2;
  graph2.SetMarkerStyle(23);
  graph2.SetMarkerColor(kRed);


  double angleMin = 0.1;
  double angleMax = 1;
  double angleStep = 0.1;
  fillGraph(graph1, "../results/perdaix_5.00_GeV_%.2f_deg_msc_inhom.root", angleMin, angleMax, angleStep);

  angleMin = 1.4;
  angleMax = 20;
  angleStep = 0.4;
  fillGraph(graph2, "../results/perdaix_5.00_GeV_%.2f_deg_msc_inhom.root", angleMin, angleMax, angleStep);

  TCanvas canvas("canvas", "canvas", 1024, 768);
  canvas.Draw();
  canvas.cd();
  canvas.SetGridx();
  canvas.SetGridy();
  graph1.Draw("AP");
  graph1.SetMarkerSize(1.5);
  graph1.GetXaxis()->SetTitle("stereo angle / deg");
  graph1.GetYaxis()->SetTitle("#sigma_{p} / p");
  graph2.Draw("P");
  graph2.SetMarkerSize(1.5);
  graph1.GetXaxis()->SetLimits(0., 16.);
  canvas.SaveAs("perdaix_angle_5GeV_zoomed.pdf");
  canvas.SaveAs("perdaix_angle_5GeV_zoomed.root");

  app->Run();

  delete myStyle;
  delete app;

  return 1;
}

