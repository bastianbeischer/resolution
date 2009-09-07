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

  int N = 20;

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
  TF1 prediction("prediction", calculatePrediction, 0., 100., 0);

  for (int i = 0; i < N; i++) {
    char filename[100];
    sprintf(filename, "../results/res_%d.root", i+1000);
    file[i] = new TFile(filename, "OPEN");
    genTree = (TTree*) file[i]->Get("resolution_gen_tree");
    recTree = (TTree*) file[i]->Get("resolution_rec_tree");
    genTree->SetBranchAddress("event", &genEvent);
    recTree->SetBranchAddress("event", &recEvent);
    genTree->GetEntry(0);
    double genMom = genEvent->GetMomentum()/1000.;
    double momRes = calculatePrediction(&genMom, 0);
    momentum[i] = genMom;
    momentumErr[i] = 0.;

    resHist = TH1D("resHist", "resHist", 30, 1. - 5.*momRes, 1. + 5.*momRes);    
    for(int j = 0; j < genTree->GetEntries(); j++) {
      genTree->GetEntry(j);
      recTree->GetEntry(j);
      resHist.Fill(genEvent->GetMomentum()/recEvent->GetMomentum());
    }
    resHist.Fit("gaus", "Q0");
    results[i] = resHist.GetFunction("gaus")->GetParameter(2);
    resultsErr[i] = resHist.GetFunction("gaus")->GetParError(2);
    file[i]->Close();
  }

  TGraphErrors graph(N, momentum, results, momentumErr, resultsErr);
  graph.SetMarkerStyle(23);
  graph.SetMarkerColor(kRed);
  graph.GetXaxis()->SetTitle("p / GeV");
  graph.GetYaxis()->SetTitle("#sigma_{p} / p");
  graph.SetTitle("momentum resolution for toy detector");
  prediction.SetLineWidth(2);

  TLegend legend(0.2, 0.6, 0.4, 0.8);
  legend.AddEntry(&graph, "results of simulation", "P");
  legend.AddEntry(&prediction, "theoretical prediction", "L");

  TCanvas canvas("canvas", "canvas", 1024, 768);
  canvas.Draw();
  canvas.cd();
  graph.Draw("AP");
  prediction.Draw("SAME");
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
