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

  int N = 10;

  TFile** file = new TFile*[N];
  double* momentum = new double[N];
  double* momentumErr = new double[N];
  double* results = new double[N];
  double* resultsErr = new double[N];

  TTree* genTree;
  TTree* recTree;
  RES_Event* genEvent = new RES_Event;
  RES_Event* recEvent = new RES_Event;
  TH1D resHist;
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


  for (int i = 0; i < N; i++) {
    char filename[100];
    sprintf(filename, "../results/perdaix_%d_GeV_5_deg.root", i+1);
    file[i] = new TFile(filename, "OPEN");
    genTree = (TTree*) file[i]->Get("resolution_gen_tree");
    recTree = (TTree*) file[i]->Get("resolution_rec_tree");
    genTree->SetBranchAddress("event", &genEvent);
    recTree->SetBranchAddress("event", &recEvent);
    genTree->GetEntry(0);
    double genMom = genEvent->GetMomentum()/1000.;
    double momRes = analyticalFormula2.Eval(genMom); //calculatePrediction(&genMom, 0);
    momentum[i] = genMom;
    momentumErr[i] = 0.;

    resHist = TH1D("resHist", "resHist", 100, 1. - 5.*momRes, 1. + 5.*momRes);    
    //    resHist = TH1D("resHist", "resHist", 100, -0.5, 2.5);
    for(int j = 0; j < genTree->GetEntries(); j++) {
      genTree->GetEntry(j);
      recTree->GetEntry(j);
      resHist.Fill(genEvent->GetMomentum()/recEvent->GetMomentum());
    }
    resHist.Fit("gaus", "Q0", "", 0.5, 2.5);
    results[i] = resHist.GetFunction("gaus")->GetParameter(2);
    resultsErr[i] = resHist.GetFunction("gaus")->GetParError(2);
    file[i]->Close();
  }

  TGraphErrors graph(N, momentum, results, momentumErr, resultsErr);
  graph.SetMarkerStyle(23);
  graph.SetMarkerColor(kRed);
  graph.GetXaxis()->SetTitle("p / GeV");
  graph.GetYaxis()->SetTitle("#sigma_{p} / p");
  graph.SetTitle("momentum resolution for perdaix - homogeneous field");
  //  prediction.SetLineWidth(2);

  TLegend legend(0.2, 0.6, 0.4, 0.8);
  legend.AddEntry(&graph, "results of simulation", "P");
  //  legend.AddEntry(&prediction, "theoretical prediction", "L");
  legend.AddEntry(&analyticalFormula, "expectation", "L");
  legend.AddEntry(&analyticalFormula2, "expectation w/o multiple scattering", "L");

  analyticalFormula2.SetLineStyle(2);

  TCanvas canvas("canvas", "canvas", 1024, 768);
  canvas.Draw();
  canvas.cd();
  canvas.SetGridx();
  canvas.SetGridy();
  graph.Draw("AP");
  analyticalFormula.Draw("SAME");
  analyticalFormula2.Draw("SAME");
  legend.Draw("SAME");

  app->Run();

  for (int i = 0; i < N; i++)
    delete file[i];
  delete[] file;
  delete[] momentum;
  delete[] results;
  delete[] momentumErr;
  delete[] resultsErr;

  return 1;
}
