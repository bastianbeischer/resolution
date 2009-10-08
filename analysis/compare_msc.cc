#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <RES_Event.hh>
#include <MyROOTStyle.h>
#include <iostream>

int main(int argc, char** argv)
{
  if (argc != 3) {
    std::cout << "please provide two root files to compare" << std::endl;
    return -1;
  }

  const char* filename1 = argv[1];
  const char* filename2 = argv[2];
  
  TApplication app("app", &argc, argv);
  MyROOTStyle myStyle("myStyle");
  myStyle.cd();

  TFile file1(filename1);
  TFile file2(filename2);

  TTree* genTree1 = (TTree*) file1.Get("resolution_gen_tree");
  TTree* genTree2 = (TTree*) file2.Get("resolution_gen_tree");

  RES_Event* genEvent1 = new RES_Event;
  RES_Event* genEvent2 = new RES_Event;
  genTree1->SetBranchAddress("event", &genEvent1);
  genTree2->SetBranchAddress("event", &genEvent2);

  genTree1->GetEntry(0);
  int nHits = genEvent1->GetNbOfHits();

  TH1D** xHist = new TH1D*[nHits];
  for (int i = 0;i < nHits; i++) {
    char title[256];
    sprintf(title, "xHist%d", i);
    xHist[i] = new TH1D(title,title, 100, -1.0, 1.0);
  }
  TH1D** yHist = new TH1D*[nHits];
  for (int i = 0;i < nHits; i++) {
    char title[256];
    sprintf(title, "yHist%d", i);
    yHist[i] = new TH1D(title,title, 100, -1.0, 1.0);
  }

  for (int i = 0; i < genTree1->GetEntries(); i++) {
    genTree1->GetEntry(i);
    genTree2->GetEntry(i);
    if ( (genEvent1->GetNbOfHits() < nHits) || (genEvent2->GetNbOfHits() < nHits) )
      continue;
    for (int j = 0; j < nHits; j++) {
      xHist[j]->Fill(genEvent1->GetHitPosition(j).x() - genEvent2->GetHitPosition(j).x());
      yHist[j]->Fill(genEvent1->GetHitPosition(j).y() - genEvent2->GetHitPosition(j).y());
    }
  }

  TCanvas canvas2("canvas2", "canvas2", 1024, 768);
  canvas2.Divide(nHits/2,2);
  canvas2.Draw();
  for (int i = 0; i < nHits; i++) {
    char xtitle[256];
    sprintf(xtitle, "(x_{%d,sim} - x_{%d,rec}) / mm", i, i);
    canvas2.cd(i+1);
    xHist[i]->Draw();
    xHist[i]->GetXaxis()->SetTitle(xtitle);
    xHist[i]->GetYaxis()->SetTitle("N");
    xHist[i]->Fit("gaus", "Q");
    TF1* fitFunc = xHist[i]->GetFunction("gaus");
    std::cout << "x" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }

  TCanvas canvas3("canvas3", "canvas3", 1024, 768);
  canvas3.Divide(nHits/2,2);
  canvas3.Draw();
  for (int i = 0; i < nHits; i++) {
    char ytitle[256];
    sprintf(ytitle, "(y_{%d,sim} - y_{%d,rec}) / mm", i, i);
    canvas3.cd(i+1);
    yHist[i]->Draw();
    yHist[i]->GetXaxis()->SetTitle(ytitle);
    yHist[i]->GetYaxis()->SetTitle("N");
    yHist[i]->Fit("gaus", "Q");
    TF1* fitFunc = yHist[i]->GetFunction("gaus");
    std::cout << "y" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }

  app.Run();

  file1.Close();
  file2.Close();

  delete genEvent1;
  delete genEvent2;

  for (int i = 0; i < nHits; i++)
    delete xHist[i];
  delete xHist;

  for (int i = 0; i < nHits; i++)
    delete yHist[i];
  delete yHist;

  return 0;
}
