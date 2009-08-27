#include <iostream>

#include "RES_Event.hh"

#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>

int main(int argc, char** argv)
{
  TApplication* app = new TApplication("app", &argc, argv);
  
  TFile* file = new TFile("../results/equi.root", "READ");
  TTree* genTree = (TTree*) file->Get("resolution_gen_tree");
  TTree* recTree = (TTree*) file->Get("resolution_rec_tree");

  RES_Event* genEvent = new RES_Event;
  RES_Event* recEvent = new RES_Event;
  genTree->SetBranchAddress("event", &genEvent);
  recTree->SetBranchAddress("event", &recEvent);

  TH1D* resHist = new TH1D("resHist", "resHist", 200, 0.90, 1.10);
  TH1D* chi2Hist = new TH1D("chi2Hist", "chi2Hist", 100, 0.0, 10.0);
  for(int i = 0; i < genTree->GetEntries(); i++) {
    genTree->GetEntry(i);
    recTree->GetEntry(i);
    resHist->Fill(genEvent->GetMomentum()/recEvent->GetMomentum());
    chi2Hist->Fill(recEvent->GetChi2OverDof());
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas", 1024, 768);
  canvas->Draw();
  resHist->Draw();

  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 1024, 768);
  canvas2->Draw();
  chi2Hist->Draw();

  app->Run();

  return 1;
}
