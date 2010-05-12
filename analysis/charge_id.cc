// $Id: charge_id.cc,v 1.3 2010/05/12 01:56:30 beischer Exp $

#include <iostream>
#include <cmath>

#include <MyROOTStyle.h>

#include <TApplication.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
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

  int momMin = 10;
  int momMax = 1200;
  int momStep = 10;
  fillGraph(graph1, "../results/pebs01_%03.0f_GeV.root", momMin, momMax, momStep, 0.8e-3, 4e-3);
  
  for (int i = 0; i < graph1.GetN(); i++) {
    double x, y;
    graph1.GetPoint(i,x,y);
    double miss = TMath::Erf(-1./(sqrt(2)*y))/2 + 0.5;
    chargeID.SetPoint(i, x, miss);
  }

  double ymin = 1e-6;
  double ymax = 1.;

  TCanvas canvas("canvas", "canvas");
  canvas.cd();
  canvas.SetGridx();
  canvas.SetGridy();
  canvas.SetLogy();
  chargeID.Draw("AP");
  chargeID.GetXaxis()->SetTitle("p / GeV");
  chargeID.GetYaxis()->SetTitle("probability");
  chargeID.SetTitle("PEBS-01: probability for wrong reconstruction of charge sign");
  chargeID.SetMarkerStyle(23);
  chargeID.SetMarkerSize(1.5);
  chargeID.GetYaxis()->SetRangeUser(ymin, ymax);

  canvas.SaveAs("pebs_charge_id.pdf");
  canvas.SaveAs("pebs_charge_id.root");

  TGraph eContami;
  for(int i = 0; i < chargeID.GetN(); i++) {
    double x = chargeID.GetX()[i];
    double y = chargeID.GetY()[i];
    double frac = 0.5;
    eContami.SetPoint(i, x, y/frac/(1+2*y/frac));
  }


  ymin = 1e-5;
  ymax = 1.;
  TCanvas canvas2("canvas2", "canvas2");
  canvas2.cd();
  canvas2.SetGridx();
  canvas2.SetGridy();
  canvas2.SetLogy();
  eContami.Draw("AP");
  eContami.GetYaxis()->SetRangeUser(ymin, ymax);
  eContami.SetMarkerStyle(23);
  eContami.SetTitle("PEBS-01: electron contamination");
  eContami.GetXaxis()->SetTitle("p / GeV");
  eContami.GetYaxis()->SetTitle("electron contamination");

  double xmin = eContami.GetXaxis()->GetXmin();
  double xmax = eContami.GetXaxis()->GetXmax();  
  double xstep = (xmax - xmin) / 1000.;

  double cut = 1e-1;
  double x = xmin;
  double y = 0.;
  while (x < xmax) {
    y = eContami.Eval(x);
    if (y > cut)
      break;
    x += xstep;
  }
  
  TLine line(xmin, cut, x, cut);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.SetLineColor(kRed);
  line.Draw("SAME");

  TLine line2(x, ymin, x, y);
  line2.SetLineWidth(2);
  line2.SetLineStyle(2);
  line2.SetLineColor(kRed);
  line2.Draw("SAME");

  char pmax[256];
  sprintf(pmax, "p_{max} = %.1f GeV", x);
  TLatex text(5.5, 1e-2, pmax);
  text.SetTextColor(kRed);
  text.Draw("SAME");
      
  canvas2.SaveAs("pebs_econtami.pdf");
  canvas2.SaveAs("pebs_econtami.root");

  app->Run();

  delete myStyle;
  delete app;

  return 1;

}
