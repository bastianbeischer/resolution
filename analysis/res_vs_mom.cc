// $Id: res_vs_mom.cc,v 1.8 2010/02/26 21:08:23 beischer Exp $

#include <iostream>
#include <cmath>

#include <MyROOTStyle.h>

#include "RES_Event.hh"

#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>

double fitfunc(double*x, double* p)
{
  double mom = x[0];
  double a = p[0];
  double b = p[1];

  return sqrt(pow(a*mom,2.) + pow(b,2.));
}

double MS(double p, double m, double L, double X0) {
  double beta = p / sqrt(p*p + m*m);
  return 13.6e-3/(beta*p) * sqrt(X0) * L;
}

double analytical(double*x, double* p)
{
  double mom = x[0];
  double mass = p[0];
  double Lup = p[1];
  double Ldown = p[2];
  double Linner = p[3];
  double B = p[4];
  double X0 = p[5];
  double sigma = p[6];

  // double soUp   = sigma;
  // double siUp   = sqrt( pow(sigma, 2.) + pow(0.135e-3 / mom, 2.) );
  // double siDown = sqrt( pow(sigma, 2.) + pow(0.240e-3 / mom, 2.) );
  // double soDown = sqrt( pow(sigma, 2.) + pow(0.440e-3 / mom, 2.) );
  // if (X0 <= 0)  {
  //   soUp = sigma;
  //   siUp = sigma;
  //   siDown = sigma;
  //   soDown = sigma;
  // }
  double soUp = sigma;
  double siUp = sigma;
  double siDown = sqrt( pow(sigma, 2.0) + pow(MS(mom,mass,Linner,X0), 2.0) );
  double soDown = sqrt( pow(sigma, 2.0) + pow(MS(mom,mass,Ldown, X0), 2.0) );

  double part1 = mom*mom / (0.3*Linner*B);
  double part2 = 1.0/(Lup*Ldown) * sqrt(pow(soUp*Ldown, 2.0) +
                                        pow(siUp*Ldown, 2.0) +
                                        pow(siDown*Lup, 2.0) +
                                        pow(soDown*Lup, 2.0));

  double sigmaP = part1*part2;

  return sigmaP/mom;
}

double calculatePrediction(double* x, double* p)
{
  double sigma = 50e-6;
  double B = 0.3;
  double L = 91e-2;
  double N = 38;
  return sigma*x[0]/(0.3*B*L*L) * sqrt(720./(N+4));
}

int main(int argc, char** argv)
{
  TApplication* app = new TApplication("app", &argc, argv);

  MyROOTStyle* myStyle = new MyROOTStyle("myStyle");
  myStyle->cd();

  TTree* genTree;
  TTree* recTree;
  RES_Event* genEvent = new RES_Event;
  RES_Event* recEvent = new RES_Event;
  //  TF1 prediction("prediction", calculatePrediction, 0., 100., 0);
  TF1 analyticalFormula("analyticalFormula", analytical, 0., 100., 7);
  TF1 analyticalFormula2("analyticalFormula2", analytical, 0., 100., 7);
  double X0 = 0.0046875; // number of radiation lengths in fibers
  double Linner = 0.08;  // in m
  double Lup   = 0.14;   // in m
  double Ldown = 0.14;   // in m
  double magField = 0.27; // in T
  //double m = 0.511e-3; // electron mass in GeV
  double m = 0.0; // geantino mass in GeV
  double sigmaModule = 50e-6/sqrt(2);
  analyticalFormula.SetParameters(m, Lup, Ldown, Linner, magField, X0, sigmaModule);
  analyticalFormula2.SetParameters(m, Lup, Ldown, Linner, magField, 0., sigmaModule);

  TGraphErrors graph1;
  graph1.SetMarkerStyle(23);
  graph1.SetMarkerColor(kRed);
  graph1.GetXaxis()->SetTitle("p / GeV");
  graph1.GetYaxis()->SetTitle("#sigma_{p} / p");
  graph1.SetTitle("momentum resolution for perdaix");

  int momBins = 60;

  double sigmaLeft = 1.;
  double sigmaRight = 2.;

  int i = 0;
  double momMin = 0.25;
  double momMax = 9.;
  double momStep = 0.25;
  for (double mom = momMin; mom <= momMax; mom += momStep) {
    char filename[100];
    sprintf(filename, "../results/perdaix_%.2f_GeV_1.00_deg.root", mom);

    TFile file(filename);
    if (file.IsZombie())
      continue;

    std::cout << filename << std::endl;

    genTree = (TTree*) file.Get("resolution_gen_tree");
    recTree = (TTree*) file.Get("resolution_rec_tree");
    genTree->SetBranchAddress("event", &genEvent);
    recTree->SetBranchAddress("event", &recEvent);
    genTree->GetEntry(0);
    double genMom = genEvent->GetMomentum()/1000.;
    //    double momRes = analyticalFormula.Eval(genMom); //calculatePrediction(&genMom, 0);
    double momRes = sqrt(pow(genMom*0.12, 2.) + pow(0.25,2.));
    TH1D resHist("resHist", "resHist", momBins, 1-5*momRes, 1+5*momRes);
    //TH1D resHist("resHist", "resHist", 100, 0., 2.);
    for(int j = 0; j < genTree->GetEntries(); j++) {
      genTree->GetEntry(j);
      recTree->GetEntry(j);
      resHist.Fill(genEvent->GetMomentum()/recEvent->GetMomentum());
    }
    double rangeLower = 1-sigmaLeft*momRes;
    double rangeUpper = 1+sigmaRight*momRes;
    resHist.Fit("gaus", "EQR0", "", rangeLower, rangeUpper);
    double sigma = resHist.GetFunction("gaus")->GetParameter(2);
    double sigmaErr = resHist.GetFunction("gaus")->GetParError(2);

    graph1.SetPoint(i, mom, sigma);
    graph1.SetPointError(i, 0., sigmaErr);

    i++;

    file.Close();
  }

  TGraphErrors graph2;
  graph2.SetMarkerStyle(22);
  graph2.SetMarkerColor(kBlue);
  graph2.GetXaxis()->SetTitle("p / GeV");
  graph2.GetYaxis()->SetTitle("#sigma_{p} / p");
  graph2.SetTitle("momentum resolution for perdaix");

  i = 0;
  for (double mom = momMin; mom <= momMax; mom += momStep) {
    char filename[100];
    sprintf(filename, "../results/perdaix_%.2f_GeV_2.00_deg.root", mom);
    std::cout << filename << std::endl;
    TFile file(filename);

    if (file.IsZombie())
      continue;

    genTree = 0;
    recTree = 0;
    genEvent = 0;
    recEvent = 0;
    genTree = (TTree*) file.Get("resolution_gen_tree");
    recTree = (TTree*) file.Get("resolution_rec_tree");
    genTree->SetBranchAddress("event", &genEvent);
    recTree->SetBranchAddress("event", &recEvent);
    genTree->GetEntry(0);
    double genMom = genEvent->GetMomentum()/1000.;
    //    double momRes = analyticalFormula.Eval(genMom); //calculatePrediction(&genMom, 0);
    double momRes = sqrt(pow(genMom*0.12, 2.) + pow(0.25,2.));
    TH1D resHist("resHist", "resHist", momBins, 1-5*momRes, 1+5*momRes);
    //TH1D resHist("resHist", "resHist", 100, 0., 2.);
    for(int j = 0; j < genTree->GetEntries(); j++) {
      genTree->GetEntry(j);
      recTree->GetEntry(j);
      resHist.Fill(genEvent->GetMomentum()/recEvent->GetMomentum());
    }

    double rangeLower = 1-sigmaLeft*momRes;
    double rangeUpper = 1+sigmaRight*momRes;
    resHist.Fit("gaus", "EQR0", "", rangeLower, rangeUpper);
    double sigma = resHist.GetFunction("gaus")->GetParameter(2);
    double sigmaErr = resHist.GetFunction("gaus")->GetParError(2);

    graph2.SetPoint(i, mom, sigma);
    graph2.SetPointError(i, 0., sigmaErr);

    i++;

    file.Close();
  }

  //  prediction.SetLineWidth(2);

  TLegend legend(0.2, 0.6, 0.4, 0.8);
  legend.AddEntry(&graph1, "1.0deg stereo angle", "P");
  legend.AddEntry(&graph2, "2.0deg stereo angle", "P");

  TF1 fit("fit", fitfunc, 0., 10., 2);
  fit.SetParNames("a", "b");

  TCanvas canvas("canvas", "canvas", 1024, 768);
  canvas.Draw();
  canvas.cd();
  canvas.SetGridx();
  canvas.SetGridy();
  graph1.Draw("AP");
  graph1.Fit("fit", "E");
  graph1.SetMarkerSize(1.5);
  graph1.GetFunction("fit")->SetLineColor(kRed);
  graph1.GetFunction("fit")->SetParNames("a_{1 deg}", "b_{1 deg}");
  graph2.Draw("P");
  graph2.Fit("fit", "E");
  graph2.SetMarkerSize(1.5);
  graph2.GetFunction("fit")->SetLineColor(kBlue);
  graph2.GetFunction("fit")->SetParNames("a_{2 deg}", "b_{2 deg}");
  graph1.GetXaxis()->SetTitle("p / GeV");
  graph1.GetYaxis()->SetTitle("#sigma_{p} / p");

  legend.Draw("SAME");

  app->Run();

  delete genEvent;
  delete recEvent;
  delete myStyle;
  delete app;

  return 1;
}
