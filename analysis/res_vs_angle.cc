// $Id: res_vs_angle.cc,v 1.9 2010/04/25 19:25:27 beischer Exp $

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
  graph1.SetTitle("Momentum resolution for PERDaix for 9 GeV electrons");

  double angleMin = 1.0;
  double angleMax = 89.0;
  double angleStep = 4.0;
  fillGraph(graph1, "../results/perdaix_9.00_GeV_%.2f_deg_msc_inhom.root", angleMin, angleMax, angleStep);

  TCanvas canvas("canvas", "canvas", 1024, 768);
  canvas.Draw();
  canvas.cd();
  canvas.SetGridx();
  canvas.SetGridy();
  graph1.Draw("AP");
  graph1.SetMarkerSize(1.5);
  graph1.GetXaxis()->SetTitle("stereo angle / deg");
  graph1.GetYaxis()->SetTitle("#sigma_{p} / p");
  canvas.SaveAs("perdaix_angle_9GeV.pdf");
  canvas.SaveAs("perdaix_angle_9GeV.root");

  app->Run();

  delete myStyle;
  delete app;

  return 1;
}

