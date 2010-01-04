// $Id: single_run.cc,v 1.11 2010/01/04 09:47:00 beischer Exp $

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
#include <TStyle.h>

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
  if (argc != 2) {
    std::cout << "please provide root file for analysis" << std::endl;
    return -1;
  }
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
  double m = 0.511e-3; // electron mass in GeV
  //double m = 0.0; // geantino mass in GeV
  double sigmaModule = 50e-6/sqrt(2);
  analyticalFormula.SetParameters(m, Lup, Ldown, Linner, magField, X0, sigmaModule);
  double momRes = analyticalFormula.Eval(genMom);
  TH1D resHist("resHist", "resHist", 100, 1. - 5.*momRes, 1. + 5.*momRes);
  //  TH1D resHist("resHist", "resHist", 100, 0.5, 1.5);
  TH1D ptHist("ptHist", "ptHist", 100, 1. - 5.*momRes, 1. + 5.*momRes);
  int nHits = genEvent->GetNbOfHits();
  int nBins = 100;
  TH1D** xDeltaGenHist = new TH1D*[nHits];
  for (int i = 0;i < nHits; i++) {
    char title[256];
    sprintf(title, "xDeltaGenHist%d", i);
    xDeltaGenHist[i] = new TH1D(title,title, nBins, -20, 20);
  }
  TH1D** yDeltaGenHist = new TH1D*[nHits];
  for (int i = 0;i < nHits; i++) {
    char title[256];
    sprintf(title, "yDeltaGenHist%d", i);
    if (i == 0 || i == nHits - 2) yDeltaGenHist[i] = new TH1D(title,title, 500, -0.3, 0.3);
    else yDeltaGenHist[i] = new TH1D(title,title, nBins, -0.3, 0.3);
  }
  TH1D** xDeltaSmearedHist = new TH1D*[nHits];
  for (int i = 0;i < nHits; i++) {
    char title[256];
    sprintf(title, "xDeltaSmearedHist%d", i);
    xDeltaSmearedHist[i] = new TH1D(title,title, nBins, -20, 20);
  }
  TH1D** yDeltaSmearedHist = new TH1D*[nHits];
  for (int i = 0;i < nHits; i++) {
    char title[256];
    sprintf(title, "yDeltaSmearedHist%d", i);
    if (i == 0 || i == nHits - 2) yDeltaSmearedHist[i] = new TH1D(title,title, 500, -0.3, 0.3);
    else yDeltaSmearedHist[i] = new TH1D(title,title, nBins, -0.3, 0.3);
  }

  TH1D totalXhist("totalXhist", "totalXhist", 500, -20, 20);
  TH1D totalYhist("totalYhist", "totalYhist", 500, -0.3, 0.3);

  TH1D chi2Hist("chi2Hist", "chi2Hist", 500, 0.0, 100.0);

  TH1D angleHist("angleHist", "angleHist", 500, -5e-3, 5e-3);

  char title[128];
  sprintf(title, "#chi^{2} Distribution (dof = %d)", recEvent->GetDof());
  chi2Hist.SetTitle(title);

  for(int i = 0; i < genTree->GetEntries(); i++) {
    //  for(int i = 0; i < 90; i++) {
    genTree->GetEntry(i);
    recTree->GetEntry(i);
    int nHitsGen = genEvent->GetNbOfHits();
    int nHitsRec = recEvent->GetNbOfHits();
    if (nHits == 0 || (nHitsGen != nHitsRec)) continue;
    if (recEvent->GetChi2() / recEvent->GetDof() < 5)
      resHist.Fill(genEvent->GetMomentum()/recEvent->GetMomentum());
    ptHist.Fill(genEvent->GetTransverseMomentum()/recEvent->GetTransverseMomentum());
    for (int i = 0; i < nHitsRec; i++) {
      xDeltaGenHist[i]->Fill(genEvent->GetHitPosition(i).x() - recEvent->GetHitPosition(i).x());
      yDeltaGenHist[i]->Fill(genEvent->GetHitPosition(i).y() - recEvent->GetHitPosition(i).y());
      // xDeltaSmearedHist[i]->Fill(genEvent->GetSmearedHitPosition(i).x() - recEvent->GetHitPosition(i).x());
      // yDeltaSmearedHist[i]->Fill(genEvent->GetSmearedHitPosition(i).y() - recEvent->GetHitPosition(i).y());
      totalXhist.Fill(genEvent->GetHitPosition(i).x() - recEvent->GetHitPosition(i).x());
      totalYhist.Fill(genEvent->GetHitPosition(i).y() - recEvent->GetHitPosition(i).y());
    }

    double angle1 = (genEvent->GetHitPosition(1).y() - genEvent->GetHitPosition(0).y())/(genEvent->GetHitPosition(1).z() - genEvent->GetHitPosition(0).z());
    double angle2 = (genEvent->GetHitPosition(nHitsGen-1).y() - genEvent->GetHitPosition(nHitsGen-2).y())/(genEvent->GetHitPosition(nHitsGen-1).z() - genEvent->GetHitPosition(nHitsGen-2).z());
    angleHist.Fill(angle2-angle1);

    chi2Hist.Fill(recEvent->GetChi2());
    char title[128];
    sprintf(title, "#chi^{2} Distribution (dof = %d)", recEvent->GetDof());
    chi2Hist.SetTitle(title);
  }


  gStyle->SetOptFit(11);
  TCanvas canvas("canvas", "Momentum resolution", 1024, 768);
  canvas.Draw();
  canvas.Divide(1,2);
  canvas.cd(1);
  resHist.Draw();
  resHist.Fit("gaus", "EQR", "", 0.1, 1.+5*momRes);
  resHist.GetXaxis()->SetTitle("p_{gen}/p_{rec}");
  resHist.GetYaxis()->SetTitle("N");
  canvas.cd(2);
  ptHist.Draw();
  ptHist.Fit("gaus", "EQR", "", 0.1, 1.+5*momRes);
  ptHist.GetXaxis()->SetTitle("pt_{gen}/pt_{rec}");
  ptHist.GetYaxis()->SetTitle("N");

  TCanvas canvas2("canvas2", "x: Reconstructed vs generated position", 1024, 768);
  canvas2.Divide(nHits/2,2);
  canvas2.Draw();

  for (int i = 0; i < nHits; i++) {
    char xtitle[256];
    sprintf(xtitle, "(x_{%d,gen} - x_{%d,rec}) / mm", i, i);
    canvas2.cd(i+1);
    xDeltaGenHist[i]->Draw();
    xDeltaGenHist[i]->GetXaxis()->SetTitle(xtitle);
    xDeltaGenHist[i]->GetYaxis()->SetTitle("N");
    xDeltaGenHist[i]->Fit("gaus", "Q");
    TF1* fitFunc = xDeltaGenHist[i]->GetFunction("gaus");
    std::cout << "x" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }

  TCanvas canvas3("canvas3", "y: Reconstructed vs generated position", 1024, 768);
  canvas3.Divide(nHits/2,2);
  canvas3.Draw();
  for (int i = 0; i < nHits; i++) {
    char ytitle[256];
    sprintf(ytitle, "(y_{%d,gen} - y_{%d,rec}) / mm", i, i);
    canvas3.cd(i+1);
    yDeltaGenHist[i]->Draw();
    yDeltaGenHist[i]->GetXaxis()->SetTitle(ytitle);
    yDeltaGenHist[i]->GetYaxis()->SetTitle("N");
    yDeltaGenHist[i]->Fit("gaus", "Q");
    TF1* fitFunc = yDeltaGenHist[i]->GetFunction("gaus");
    std::cout << "y" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }

  TCanvas canvas4("canvas4", "x: Reconstructed vs measured position", 1024, 768);
  canvas4.Divide(nHits/2,2);
  canvas4.Draw();

  for (int i = 0; i < nHits; i++) {
    char xtitle[256];
    sprintf(xtitle, "(x_{%d,meas} - x_{%d,rec}) / mm", i, i);
    canvas4.cd(i+1);
    xDeltaSmearedHist[i]->Draw();
    xDeltaSmearedHist[i]->GetXaxis()->SetTitle(xtitle);
    xDeltaSmearedHist[i]->GetYaxis()->SetTitle("N");
    //    xDeltaSmearedHist[i]->Fit("gaus", "Q");
    // TF1* fitFunc = xDeltaSmearedHist[i]->GetFunction("gaus");
    // std::cout << "x" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }

  TCanvas canvas5("canvas5", "y: Reconstructed vs measured Position", 1024, 768);
  canvas5.Divide(nHits/2,2);
  canvas5.Draw();
  for (int i = 0; i < nHits; i++) {
    char ytitle[256];
    sprintf(ytitle, "(y_{%d,meas} - y_{%d,rec}) / mm", i, i);
    canvas5.cd(i+1);
    yDeltaSmearedHist[i]->Draw();
    yDeltaSmearedHist[i]->GetXaxis()->SetTitle(ytitle);
    yDeltaSmearedHist[i]->GetYaxis()->SetTitle("N");
    // yDeltaSmearedHist[i]->Fit("gaus", "Q");
    // TF1* fitFunc = yDeltaSmearedHist[i]->GetFunction("gaus");
    // std::cout << "y" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }

  TCanvas canvas6("canvas6", "Sum of reconstructed vs generated position histograms", 1024, 768);
  canvas6.Draw();
  canvas6.Divide(1,2);
  canvas6.cd(1);
  totalXhist.Draw();
  totalXhist.GetXaxis()->SetTitle("(x_{sim,total} - x_{rec,total}) / mm");
  totalXhist.GetYaxis()->SetTitle("N");
  canvas6.cd(2);
  totalYhist.Draw();
  totalYhist.GetXaxis()->SetTitle("(y_{sim,total} - y_{rec,total}) / mm");
  totalYhist.GetYaxis()->SetTitle("N");

  TF1 chi2Dist("chi2Dist", chi2dist, 0.0, 100.0, 2);
  chi2Dist.SetNpx(1000);
  
  chi2Dist.FixParameter(0, 1);
  chi2Dist.FixParameter(1, recEvent->GetDof());
  TCanvas canvas7("canvas7", "Chi2 distribution", 1024, 768);
  canvas7.Draw();
  chi2Hist.Draw();
  chi2Hist.Scale(1./(chi2Hist.Integral("WIDTH")));
  chi2Hist.GetXaxis()->SetTitle("#chi^{2}");
  chi2Hist.GetYaxis()->SetTitle("N");
  //  chi2Hist.Fit(&chi2Dist, "E");
  chi2Dist.Draw("SAME");
  chi2Dist.SetLineColor(kRed);

  TCanvas canvas8("canvas8", "Deflection angle distribution", 1024, 768);
  canvas8.Draw();
  angleHist.Draw();
  angleHist.GetXaxis()->SetTitle("#Delta #theta [rad]");
  angleHist.GetYaxis()->SetTitle("N");

  app.Run();

  file.Close();
  delete genEvent;
  delete recEvent;

  for (int i = 0; i < nHits; i++)
    delete xDeltaGenHist[i];
  delete xDeltaGenHist;

  for (int i = 0; i < nHits; i++)
    delete yDeltaGenHist[i];
  delete yDeltaGenHist;


  return 1;
}
