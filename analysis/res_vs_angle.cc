// $Id: res_vs_angle.cc,v 1.6 2009/10/14 16:51:37 beischer Exp $

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
  double X0 = 0.0046875; // number of radiation lengths in fibers
  double Linner = 0.08;  // in m
  double Lup   = 0.14;   // in m
  double Ldown = 0.14;   // in m
  double magField = 0.27; // in T
  //double m = 0.511e-3; // electron mass in GeV
  double m = 0.0; // geantino mass in GeV
  double sigmaModule = 50e-6/sqrt(2);
  analyticalFormula.SetParameters(m, Lup, Ldown, Linner, magField, 0., sigmaModule);
  
  TGraphErrors sigmaVsAngleGraph;
  sigmaVsAngleGraph.SetMarkerStyle(23);
  sigmaVsAngleGraph.SetMarkerColor(kRed);
  sigmaVsAngleGraph.GetXaxis()->SetTitle("stereo angle [deg]");
  sigmaVsAngleGraph.GetYaxis()->SetTitle("#sigma_{p} / p");
  sigmaVsAngleGraph.SetTitle("momentum resolution for perdaix");

  TGraphErrors muVsAngleGraph;
  muVsAngleGraph.SetMarkerStyle(23);
  muVsAngleGraph.SetMarkerColor(kRed);
  muVsAngleGraph.GetXaxis()->SetTitle("stereo angle [deg]");
  muVsAngleGraph.GetYaxis()->SetTitle("#mu_{p_{sim}/p_{rec}}");
  muVsAngleGraph.SetTitle("distribution of means");

  int i = 0;
  
  double angleMin = 0.1;
  double angleMax = 10.0;
  double angleStep = 0.1;
  for (double angle = angleMin; angle < angleMax; angle += angleStep) {
    char filename[100];
    sprintf(filename, "../results/perdaix_1.0_GeV_%.2f_deg_inhom.root", angle);
    TFile file(filename);

    std::cout << "Opening file: " << filename << std::endl;

    if (file.IsZombie())
      continue;

    genTree = (TTree*) file.Get("resolution_gen_tree");
    recTree = (TTree*) file.Get("resolution_rec_tree");
    genTree->SetBranchAddress("event", &genEvent);
    recTree->SetBranchAddress("event", &recEvent);
    genTree->GetEntry(0);
    double genMom = genEvent->GetMomentum()/1000.;
    double momRes = analyticalFormula.Eval(genMom);

    TH1D resHist("resHist", "resHist", 100, 1. - 5.*momRes, 1. + 5.*momRes);    
    for(int j = 0; j < genTree->GetEntries(); j++) {
      genTree->GetEntry(j);
      recTree->GetEntry(j);
      resHist.Fill(genEvent->GetMomentum()/recEvent->GetMomentum());
    }
    resHist.Fit("gaus", "Q0", "", 0.5, 2.5);

    double mu = resHist.GetFunction("gaus")->GetParameter(1);
    double muErr = resHist.GetFunction("gaus")->GetParError(1);
    double sigma = resHist.GetFunction("gaus")->GetParameter(2);
    double sigmaErr = resHist.GetFunction("gaus")->GetParError(2);

    sigmaVsAngleGraph.SetPoint(i, angle, sigma);
    sigmaVsAngleGraph.SetPointError(i, 0.0, sigmaErr);
    muVsAngleGraph.SetPoint(i, angle, mu);
    muVsAngleGraph.SetPointError(i, 0.0, muErr);
    i++;

    std::cout << "  --> " << genTree->GetEntries() << " events processed" << std::endl;

    file.Close();
  }

  TCanvas canvas("canvas", "canvas", 1024, 768);
  canvas.Draw();
  canvas.cd();
  canvas.SetGridx();
  canvas.SetGridy();
  sigmaVsAngleGraph.Draw("AP");

  TCanvas canvas2("canvas2", "canvas2", 1024, 768);
  canvas2.Draw();
  canvas2.cd();
  canvas2.SetGridx();
  canvas2.SetGridy();
  muVsAngleGraph.Draw("AP");

  app->Run();

  return 1;
}
