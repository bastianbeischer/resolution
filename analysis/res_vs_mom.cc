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

  TGraphErrors graph;
  graph.SetMarkerStyle(23);
  graph.SetMarkerColor(kRed);
  graph.GetXaxis()->SetTitle("p / GeV");
  graph.GetYaxis()->SetTitle("#sigma_{p} / p");
  graph.SetTitle("momentum resolution for perdaix");

  int i = 0;
  double momMin = 0.5;
  double momMax = 10.0;
  double momStep = 0.5;
  for (double mom = momMin; mom <= momMax; mom += momStep) {
    char filename[100];
    sprintf(filename, "../results/perdaix_%.1f_GeV_5.00_deg.root", mom);
    TFile file(filename);

    if (file.IsZombie())
      continue;

    genTree = (TTree*) file.Get("resolution_gen_tree");
    recTree = (TTree*) file.Get("resolution_rec_tree");
    genTree->SetBranchAddress("event", &genEvent);
    recTree->SetBranchAddress("event", &recEvent);
    genTree->GetEntry(0);
    double genMom = genEvent->GetMomentum()/1000.;
    double momRes = analyticalFormula.Eval(genMom); //calculatePrediction(&genMom, 0);

    TH1D resHist("resHist", "resHist", 100, 1. - 5.*momRes, 1. + 5.*momRes);    
    for(int j = 0; j < genTree->GetEntries(); j++) {
      genTree->GetEntry(j);
      recTree->GetEntry(j);
      resHist.Fill(genEvent->GetMomentum()/recEvent->GetMomentum());
    }
    resHist.Fit("gaus", "Q0", "", 0.5, 2.5);
    double sigma = resHist.GetFunction("gaus")->GetParameter(2);
    double sigmaErr = resHist.GetFunction("gaus")->GetParError(2);

    graph.SetPoint(i, mom, sigma);
    graph.SetPointError(i, 0., sigmaErr);

    i++;

    file.Close();
  }

  //  prediction.SetLineWidth(2);

  TLegend legend(0.2, 0.6, 0.4, 0.8);
  legend.AddEntry(&graph, "results of simulation", "P");
  //  legend.AddEntry(&prediction, "theoretical prediction", "L");
  legend.AddEntry(&analyticalFormula, "expectation", "L");
  legend.AddEntry(&analyticalFormula2, "expectation w/o multiple scattering", "L");

  analyticalFormula2.SetLineStyle(2);

  TF1 fit("fit", fitfunc, 0., 10., 2);

  TCanvas canvas("canvas", "canvas", 1024, 768);
  canvas.Draw();
  canvas.cd();
  canvas.SetGridx();
  canvas.SetGridy();
  graph.Draw("AP");
  graph.Fit("fit", "E");
  analyticalFormula.Draw("SAME");
  analyticalFormula2.Draw("SAME");
  legend.Draw("SAME");

  app->Run();

  delete genEvent;
  delete recEvent;
  delete myStyle;
  delete app;

  return 1;
}
