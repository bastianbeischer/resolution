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
  recTree->GetEntry(0);
  double genMom = genEvent->GetMomentum()/1000.;
  // double momRes = calculatePrediction(&genMom, 0);
  // TH1D resHist("resHist", "resHist", 100, 1. - 5.*momRes, 1. + 5.*momRes);    
  TH1D resHist("resHist", "resHist", 100, 0.5, 1.5);    
  TH1D xHist("xHist", "xHist", 100, -5, 5);
  TH1D yHist("yHist", "yHist", 100, -0.2, 0.2);
  TH1D chi2Hist("chi2Hist", "chi2Hist", 100, 0.0, 10.0);

  char title[128];
  sprintf(title, "#chi^{2} Distribution (dof = %d)", recEvent->GetDof());
  chi2Hist.SetTitle(title);

  for(int i = 0; i < 3500; i++) {
    genTree->GetEntry(i);
    recTree->GetEntry(i);
    //    std::cout << "rec mom: " << recEvent->GetMomentum() << "  --> frac: " << genEvent->GetMomentum()/recEvent->GetMomentum() << std::endl;
    std::cout << i << " --> rec pos: " << recEvent->GetHitPosition(0).x() << "  <--> sim pos: " << genEvent->GetHitPosition(0).x() << std::endl;
    resHist.Fill(genEvent->GetMomentum()/recEvent->GetMomentum());
    xHist.Fill(genEvent->GetHitPosition(0).x() - recEvent->GetHitPosition(0).x());
    yHist.Fill(genEvent->GetHitPosition(0).y() - recEvent->GetHitPosition(0).y());
    chi2Hist.Fill(recEvent->GetChi2());
  }


  TCanvas canvas("canvas", "canvas", 1024, 768);
  canvas.Draw();
  resHist.Draw();
  resHist.GetXaxis()->SetTitle("p_{sim}/p_{rec}");
  resHist.GetYaxis()->SetTitle("N");

  TCanvas canvas2("canvas2", "canvas2", 1024, 768);
  canvas2.Divide(1,2);
  canvas2.Draw();

  canvas2.cd(1);
  xHist.Draw();
  xHist.GetXaxis()->SetTitle("x_{0,sim} - x_{0,rec}");
  xHist.GetYaxis()->SetTitle("N");

  canvas2.cd(2);
  yHist.Draw();
  yHist.GetYaxis()->SetTitle("y_{0,sim} - y_{0,rec}");
  yHist.GetYaxis()->SetTitle("N");

  TCanvas canvas3("canvas3", "canvas3", 1024, 768);
  canvas3.Draw();
  chi2Hist.Draw();
  chi2Hist.GetXaxis()->SetTitle("#chi^{2}/dof");
  chi2Hist.GetYaxis()->SetTitle("N");

  app.Run();

  delete genEvent;
  delete recEvent;
  file.Close();


  return 1;
}
