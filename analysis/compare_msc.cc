// $Id: compare_msc.cc,v 1.8 2010/01/04 09:47:00 beischer Exp $

#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <RES_Event.hh>
#include <MyROOTStyle.h>
#include <iostream>
#include <cstdlib>

int main(int argc, char** argv)
{
  if (argc != 4) {
    std::cout << "please provide two root files to compare" << std::endl;
    return -1;
  }

  const char* filename1 = argv[1];
  const char* filename2 = argv[2];
  int printLayer = atoi(argv[3]);
  
  TApplication app("app", &argc, argv);
  MyROOTStyle myStyle("myStyle");
  myStyle.cd();
  gStyle->SetOptFit(111);

  TFile file1(filename1);
  TFile file2(filename2);

  TTree* genTree1 = (TTree*) file1.Get("resolution_gen_tree");
  TTree* genTree2 = (TTree*) file2.Get("resolution_gen_tree");
  TTree* recTree1 = (TTree*) file1.Get("resolution_rec_tree");
  TTree* recTree2 = (TTree*) file2.Get("resolution_rec_tree");

  RES_Event* genEvent1 = new RES_Event;
  RES_Event* genEvent2 = new RES_Event;
  RES_Event* recEvent1 = new RES_Event;
  RES_Event* recEvent2 = new RES_Event;
  genTree1->SetBranchAddress("event", &genEvent1);
  genTree2->SetBranchAddress("event", &genEvent2);
  recTree1->SetBranchAddress("event", &recEvent1);
  recTree2->SetBranchAddress("event", &recEvent2);

  genTree1->GetEntry(0);
  int nHits = genEvent1->GetNbOfHits();

  TH1D** xHist1 = new TH1D*[nHits];
  TH1D** yHist1 = new TH1D*[nHits];
  TH1D** xHist2 = new TH1D*[nHits];
  TH1D** yHist2 = new TH1D*[nHits];
  TH1D** xHist3 = new TH1D*[nHits];
  TH1D** yHist3 = new TH1D*[nHits];
  TH1D** xHist4 = new TH1D*[nHits];
  TH1D** yHist4 = new TH1D*[nHits];
  int nBins = 100;
  double lowerLimit = -0.3;
  double upperLimit =  0.3;
  for (int i = 0;i < nHits; i++) {
    char title[256];
    sprintf(title, "xHist1%d", i);
    xHist1[i] = new TH1D(title,title, nBins, lowerLimit, upperLimit);
    sprintf(title, "yHist1%d", i);
    yHist1[i] = new TH1D(title,title, nBins, lowerLimit, upperLimit);
    sprintf(title, "xHist2%d", i);
    xHist2[i] = new TH1D(title,title, nBins, lowerLimit, upperLimit);
    sprintf(title, "yHist2%d", i);
    yHist2[i] = new TH1D(title,title, nBins, lowerLimit, upperLimit);
    sprintf(title, "xHist3%d", i);
    xHist3[i] = new TH1D(title,title, nBins, lowerLimit, upperLimit);
    sprintf(title, "yHist3%d", i);
    yHist3[i] = new TH1D(title,title, nBins, lowerLimit, upperLimit);
    sprintf(title, "xHist4%d", i);
    xHist4[i] = new TH1D(title,title, nBins, lowerLimit, upperLimit);
    sprintf(title, "yHist4%d", i);
    yHist4[i] = new TH1D(title,title, nBins, lowerLimit, upperLimit);
  }

  for (int i = 0; i < genTree1->GetEntries(); i++) {
    genTree1->GetEntry(i);
    genTree2->GetEntry(i);
    recTree1->GetEntry(i);
    recTree2->GetEntry(i);
    if ( (genEvent1->GetNbOfHits() < nHits) || (genEvent2->GetNbOfHits() < nHits) || (recEvent1->GetNbOfHits() < nHits) || (recEvent1->GetNbOfHits() < nHits) )
      continue;
    for (int j = 0; j < nHits; j++) {
      xHist1[j]->Fill(genEvent2->GetHitPosition(j).x() - genEvent1->GetHitPosition(j).x());
      yHist1[j]->Fill(genEvent2->GetHitPosition(j).y() - genEvent1->GetHitPosition(j).y());
      xHist2[j]->Fill(recEvent2->GetHitPosition(j).x() - genEvent1->GetHitPosition(j).x());
      yHist2[j]->Fill(recEvent2->GetHitPosition(j).y() - genEvent1->GetHitPosition(j).y());
      xHist3[j]->Fill(recEvent2->GetHitPosition(j).x() - genEvent2->GetSmearedHitPosition(j).x());
      yHist3[j]->Fill(recEvent2->GetHitPosition(j).y() - genEvent2->GetSmearedHitPosition(j).y());
      xHist4[j]->Fill(recEvent2->GetHitPosition(j).x() - genEvent2->GetHitPosition(j).x());
      yHist4[j]->Fill(recEvent2->GetHitPosition(j).y() - genEvent2->GetHitPosition(j).y());
    }
  }

  // TCanvas canvas2("canvas2", "canvas2", 1024, 768);
  // canvas2.Divide(nHits/2,2);
  // canvas2.Draw();
  // for (int i = 0; i < nHits; i++) {
  //   char xtitle[256];
  //   sprintf(xtitle, "(x_{%d,sim} - x_{%d,rec}) / mm", i, i);
  //   canvas2.cd(i+1);
  //   xHist[i]->Draw();
  //   xHist[i]->GetXaxis()->SetTitle(xtitle);
  //   xHist[i]->GetYaxis()->SetTitle("N");
  //   xHist[i]->Fit("gaus", "Q");
  //   TF1* fitFunc = xHist[i]->GetFunction("gaus");
  //   std::cout << "x" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  // }

  TCanvas canvas("gengen_can", "msc_gen - no_msc_gen", 1024, 768);
  canvas.Divide(nHits/2,2);
  //  canvas.Divide(2,2);
  canvas.Draw();

  for (int i = 0; i < nHits; i++) {
    char ytitle[256];
    sprintf(ytitle, "(y_{%d,sim} - y_{%d,rec}) / mm", i, i);
    // switch(i) {
    // case 0:
    //   canvas.cd(1);
    //   break;
    // case 3:
    //   canvas.cd(2);
    //   break;
    // case 6:
    //   canvas.cd(3);
    //   break;
    // case 10:
    //   canvas.cd(4);
    //   break;
    // default:
    //   continue;
    // }
    canvas.cd(i+1);
    yHist1[i]->Draw();
    yHist1[i]->GetXaxis()->SetTitle(ytitle);
    yHist1[i]->GetYaxis()->SetTitle("N");
    yHist1[i]->Fit("gaus", "Q");
    TF1* fitFunc = yHist1[i]->GetFunction("gaus");
    //    std::cout << "y" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }
  TCanvas canvas2("recgen_can", "rec - no_msc_gen", 1024, 768);
  canvas2.Divide(nHits/2,2);
  canvas2.Draw();
  for (int i = 0; i < nHits; i++) {
    char ytitle[256];
    sprintf(ytitle, "(y_{%d,sim} - y_{%d,rec}) / mm", i, i);
    canvas2.cd(i+1);
    yHist2[i]->Draw();
    yHist2[i]->GetXaxis()->SetTitle(ytitle);
    yHist2[i]->GetYaxis()->SetTitle("N");
    yHist2[i]->Fit("gaus", "Q");
    TF1* fitFunc = yHist2[i]->GetFunction("gaus");
    //std::cout << "y" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }
  TCanvas canvas3("genrec_can", "rec - msc_gen_smeared (RESIDUALS)", 1024, 768);
  canvas3.Divide(nHits/2,2);
  canvas3.Draw();
  for (int i = 0; i < nHits; i++) {
    char ytitle[256];
    sprintf(ytitle, "(y_{%d,sim} - y_{%d,rec}) / mm", i, i);
    canvas3.cd(i+1);
    yHist3[i]->Draw();
    yHist3[i]->GetXaxis()->SetTitle(ytitle);
    yHist3[i]->GetYaxis()->SetTitle("N");
    yHist3[i]->Fit("gaus", "Q");
    TF1* fitFunc = yHist3[i]->GetFunction("gaus");
    //std::cout << "y" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }
  TCanvas canvas4("recrec_can", "rec - msc_gen (REAL ERROR OF TRACKFIT)", 1024, 768);
  canvas4.Divide(nHits/2,2);
  canvas4.Draw();
  for (int i = 0; i < nHits; i++) {
    char ytitle[256];
    sprintf(ytitle, "(y_{%d,sim} - y_{%d,rec}) / mm", i, i);
    canvas4.cd(i+1);
    yHist4[i]->Draw();
    yHist4[i]->GetXaxis()->SetTitle(ytitle);
    yHist4[i]->GetYaxis()->SetTitle("N");
    yHist4[i]->Fit("gaus", "Q");
    TF1* fitFunc = yHist4[i]->GetFunction("gaus");
    if (i==printLayer) std::cout << fitFunc->GetParameter(2) << std::endl;
    //std::cout << "y" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }

  app.Run();

  file1.Close();
  file2.Close();

  delete genEvent1;
  delete genEvent2;
  delete recEvent1;
  delete recEvent2;

  for (int i = 0; i < nHits; i++) {
    delete xHist1[i];
    delete yHist1[i];
    delete xHist2[i];
    delete yHist2[i];
    delete xHist3[i];
    delete yHist3[i];
    delete xHist4[i];
    delete yHist4[i];
  }
  delete [] xHist1;
  delete [] yHist1;
  delete [] xHist2;
  delete [] yHist2;
  delete [] xHist3;
  delete [] yHist3;
  delete [] xHist4;
  delete [] yHist4;

  return 0;
}
