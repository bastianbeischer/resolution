// $Id: charge_id.cc,v 1.1 2010/05/02 22:59:30 beischer Exp $

#include <iostream>
#include <cmath>

#include <MyROOTStyle.h>

#include <TApplication.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMath.h>

#include "fill_graph.hh"

int main(int argc, char** argv)
{
  TApplication* app = new TApplication("app", &argc, argv);

  MyROOTStyle* myStyle = new MyROOTStyle("myStyle");
  myStyle->cd();

  gStyle->SetOptFit(1111);

  TGraphErrors graph1;
  TGraph chargeID;

  double momMin = 0.25;
  double momMax = 9.0;
  double momStep = 0.25;
  fillGraph(graph1, "../results/perdaix_%.2f_GeV_1.00_deg_msc_inhom.root", momMin, momMax, momStep);
  
  for (int i = 0; i < graph1.GetN(); i++) {
    double x, y;
    graph1.GetPoint(i,x,y);
    double miss = TMath::Erf(-1./(sqrt(2)*y))/2 + 0.5;
    chargeID.SetPoint(i, x, miss);
  }

  TCanvas canvas("canvas", "canvas");
  canvas.cd();
  canvas.SetGridx();
  canvas.SetGridy();
  chargeID.Draw("AP");
  chargeID.GetXaxis()->SetTitle("p / GeV");
  chargeID.GetYaxis()->SetTitle("probability");
  chargeID.SetTitle("PERDaix: probability for wrong reconstruction of charge sign");
  chargeID.SetMarkerStyle(23);
  chargeID.SetMarkerSize(1.5);

  canvas.SaveAs("perdaix_charge_id.pdf");
  canvas.SaveAs("perdaix_charge_id.root");

  app->Run();

  delete myStyle;
  delete app;

  return 1;

}
