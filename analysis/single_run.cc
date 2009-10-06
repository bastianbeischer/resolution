#include <iostream>
#include <cmath>

#include <MyROOTStyle.h>

#include "RES_Event.hh"

#include <TApplication.h>
#include <TFile.h>
#include <TF1.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>

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

double chi2dist(double*x, double*p)
{
  double amplitude = p[0];
  double ndf = p[1];
  double nom = pow(x[0],ndf/2. - 1.) * exp(-x[0]/2.);
  double denom = pow(2., ndf/2.) * TMath::Gamma(ndf/2.);
  return amplitude*nom/denom;
}

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
  TF1 analyticalFormula("analyticalFormula", analytical, 0., 100., 7);
  double X0 = 0.0046875; // number of radiation lengths in fibers
  double Linner = 0.08;  // in m
  double Lup   = 0.14;   // in m
  double Ldown = 0.14;   // in m
  double magField = 0.25; // in T
  //double m = 0.511e-3; // electron mass in GeV
  double m = 0.0; // geantino mass in GeV
  double sigmaModule = 50e-6/sqrt(2);
  analyticalFormula.SetParameters(m, Lup, Ldown, Linner, magField, 0., sigmaModule);
  double momRes = analyticalFormula.Eval(genMom);
  TH1D resHist("resHist", "resHist", 100, 1. - 5.*momRes, 1. + 5.*momRes);
  //  TH1D resHist("resHist", "resHist", 100, 0.5, 1.5);
  TH1D ptHist("ptHist", "ptHist", 100, 1. - 5.*momRes, 1. + 5.*momRes);
  int nHits = genEvent->GetNbOfHits();
  TH1D** xHist = new TH1D*[nHits];
  for (int i = 0;i < nHits; i++) {
    char title[256];
    sprintf(title, "xHist%d", i);
    xHist[i] = new TH1D(title,title, 100, -20, 20);
  }
  TH1D** yHist = new TH1D*[nHits];
  for (int i = 0;i < nHits; i++) {
    char title[256];
    sprintf(title, "yHist%d", i);
    yHist[i] = new TH1D(title,title, 100, -1.0, 1.0);
  }

  TH1D totalXhist("totalXhist", "totalXhist", 100, -20, 20);
  TH1D totalYhist("totalYhist", "totalYhist", 100, -1.0, 1.0);

  TH1D chi2Hist("chi2Hist", "chi2Hist", 500, 0.0, 100.0);

  char title[128];
  sprintf(title, "#chi^{2} Distribution (dof = %d)", recEvent->GetDof());
  chi2Hist.SetTitle(title);

  for(int i = 0; i < genTree->GetEntries(); i++) {
    //  for(int i = 0; i < 90; i++) {
    genTree->GetEntry(i);
    recTree->GetEntry(i);
    int nHitsGen = genEvent->GetNbOfHits();
    int nHitsRec = recEvent->GetNbOfHits();
    if (nHits == 0 || nHitsGen != nHitsRec) continue;
    resHist.Fill(genEvent->GetMomentum()/recEvent->GetTransverseMomentum());
    ptHist.Fill(genEvent->GetTransverseMomentum()/recEvent->GetTransverseMomentum());
    for (int i = 0; i < nHitsRec; i++) {
      xHist[i]->Fill(genEvent->GetHitPosition(i).x() - recEvent->GetHitPosition(i).x());
      yHist[i]->Fill(genEvent->GetHitPosition(i).y() - recEvent->GetHitPosition(i).y());
      totalXhist.Fill(genEvent->GetHitPosition(i).x() - recEvent->GetHitPosition(i).x());
      totalYhist.Fill(genEvent->GetHitPosition(i).y() - recEvent->GetHitPosition(i).y());
    }
    chi2Hist.Fill(recEvent->GetChi2());
    char title[128];
    sprintf(title, "#chi^{2} Distribution (dof = %d)", recEvent->GetDof());
    chi2Hist.SetTitle(title);
  }


  TCanvas canvas("canvas", "canvas", 1024, 768);
  canvas.Draw();
  canvas.Divide(1,2);
  canvas.cd(1);
  resHist.Draw();
  resHist.Fit("gaus", "EQR", "", 0.1, 1.+5*momRes);
  resHist.GetXaxis()->SetTitle("p_{sim}/p_{rec}");
  resHist.GetYaxis()->SetTitle("N");
  canvas.cd(2);
  ptHist.Draw();
  ptHist.Fit("gaus", "EQR", "", 0.1, 1.+5*momRes);
  ptHist.GetXaxis()->SetTitle("pt_{sim}/pt_{rec}");
  ptHist.GetYaxis()->SetTitle("N");

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

  TCanvas canvas4("canvas4", "canvas4", 1024, 768);
  canvas4.Draw();
  canvas4.Divide(1,2);
  canvas4.cd(1);
  totalXhist.Draw();
  totalXhist.GetXaxis()->SetTitle("(x_{sim,total} - x_{rec,total}) / mm");
  totalXhist.GetYaxis()->SetTitle("N");
  canvas4.cd(2);
  totalYhist.Draw();
  totalYhist.GetXaxis()->SetTitle("(y_{sim,total} - y_{rec,total}) / mm");
  totalYhist.GetYaxis()->SetTitle("N");

  TF1 chi2Dist("chi2Dist", chi2dist, 0.0, 100.0, 2);
  //  chi2Dist.SetParameter(0,recEvent->GetDof());

  chi2Dist.SetParameters(1000, 3);

  TCanvas canvas5("canvas5", "canvas5", 1024, 768);
  canvas5.Draw();
  chi2Hist.Draw();
  chi2Hist.GetXaxis()->SetTitle("#chi^{2}");
  chi2Hist.GetYaxis()->SetTitle("N");
  chi2Hist.Fit(&chi2Dist, "E");
  chi2Dist.Draw("SAME");
  chi2Dist.SetLineColor(kRed);

  app.Run();

  delete genEvent;
  delete recEvent;
  file.Close();


  return 1;
}
