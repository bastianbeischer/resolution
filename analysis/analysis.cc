#include <iostream>
#include <cmath>

#include <MyROOTStyle.h>

#include "RES_Event.hh"

#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>

double calculatePrediction(double* x, double* p)
{
  double sigma = 50e-6; // m
  double B = 0.3;       // T
  double L = 91e-2;     // m
  double N = 38;
  return sigma*x[0]/(0.3*B*L*L) * sqrt(720./(N+4));
}

int main(int argc, char** argv)
{
  if (!argc == 2)
    std::cout << "please provide root file for analysis" << std::endl;
  const char* filename = argv[1];

  TApplication app("app", &argc, argv);
  
  MyROOTStyle myStyle("myStyle");
  myStyle.cd();

  TFile file(filename, "READ");
  TTree* genTree = (TTree*) file.Get("resolution_gen_tree");
  TTree* recTree = (TTree*) file.Get("resolution_rec_tree");

  RES_Event* genEvent = new RES_Event;
  RES_Event* recEvent = new RES_Event;
  genTree->SetBranchAddress("event", &genEvent);
  recTree->SetBranchAddress("event", &recEvent);

  genTree->GetEntry(0);
  double genMom = genEvent->GetMomentum()/1000.;
  // double momRes = calculatePrediction(&genMom, 0);
  // TH1D resHist("resHist", "resHist", 100, 1. - 5.*momRes, 1. + 5.*momRes);    
  TH1D resHist("resHist", "resHist", 100, 1. - 5.*0.1, 1. + 5.*0.1);    
  TH1D chi2Hist("chi2Hist", "chi2Hist", 100, 0.0, 10.0);
  for(int i = 0; i < genTree->GetEntries(); i++) {
    genTree->GetEntry(i);
    recTree->GetEntry(i);
    std::cout << "rec mom: " << recEvent->GetMomentum() << "  --> frac: " << genEvent->GetMomentum()/recEvent->GetMomentum() << std::endl;

    resHist.Fill(genEvent->GetMomentum()/recEvent->GetMomentum());
    chi2Hist.Fill(recEvent->GetChi2OverDof());
  }

  TCanvas canvas("canvas", "canvas", 1024, 768);
  canvas.Draw();
  resHist.Draw();
  resHist.GetXaxis()->SetTitle("p_{sim}/p_{rec}");
  resHist.GetYaxis()->SetTitle("N");

  TCanvas canvas2("canvas2", "canvas2", 1024, 768);
  canvas2.Draw();
  chi2Hist.Draw();
  chi2Hist.GetXaxis()->SetTitle("#chi^{2}/dof");
  chi2Hist.GetYaxis()->SetTitle("N");

  app.Run();

  delete genEvent;
  delete recEvent;
  file.Close();


  return 1;
}
