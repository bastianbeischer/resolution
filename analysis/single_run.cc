// $Id: single_run.cc,v 1.36 2010/06/29 15:09:21 beischer Exp $

#include <iostream>
#include <cmath>

#include <MyROOTStyle.h>

#include "RES_Event.hh"

#include <TApplication.h>
#include <TFile.h>
#include <TF1.h>
#include <TTree.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <THistPainter.h>
#include <TPaveStats.h>

double MS(double p, double m, double L, double X0) {
  double beta = p / sqrt(p*p + m*m);
  return 13.6e-3/(beta*p) * sqrt(X0) * L;
}


void fitInChi2Range(TH1D* hist, double& sigma, double& error)
{
  if ( hist ) {
    double chi2 = DBL_MAX;
    double ndf = 1;
    int iBin = 1;
    int nBins = hist->GetNbinsX();
    TF1* function;
    while(chi2/ndf > 1 && iBin < nBins/2) {
      double center = hist->GetBinCenter(nBins/2); 
      double range = center - hist->GetBinLowEdge(iBin);
      hist->Fit("gaus", "QR", "", center-range, center+range);
      function = hist->GetFunction("gaus");
      if (function) {
        chi2 = function->GetChisquare();
        ndf = function->GetNDF();
      }
      iBin++;
    }
    if (function) {
      sigma = function->GetParameter(2);
      error = function->GetParError(2);
    }
    return;
  }
  return;
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
  // TF1 analyticalFormula("analyticalFormula", analytical, 0., 100., 7);
  // double X0 = 0.0046875; // number of radiation lengths in fibers
  // double Linner = 0.08;  // in m
  // double Lup   = 0.14;   // in m
  // double Ldown = 0.14;   // in m
  // double magField = 0.25; // in T
  // double m = 0.511e-3; // electron mass in GeV
  // //double m = 0.0; // geantino mass in GeV
  // double sigmaModule = 50e-6/sqrt(2);
  // analyticalFormula.SetParameters(m, Lup, Ldown, Linner, magField, X0, sigmaModule);
  // double momRes = analyticalFormula.Eval(genMom);
  //  TH1D resHist("resHist", "resHist", 100, 1. - 5.*momRes, 1. + 5.*momRes);
  double momRes = sqrt(pow(genMom*0.08, 2.) + pow(0.21,2.));
  //double momRes = 0.8;
  //  double momRes = sqrt(pow(genMom*.8e-3, 2.) + pow(0.04,2.));
  char title[256];
  sprintf(title, "Data with nominal energy of %.2f GeV", genMom);
  TH1D resHist("Mom. Resolution", title, 100, 1-10*momRes, 1+10*momRes);
  //TH1D resHist("Mom. Resolution", title, 100, -2.5, 3.0);
  sprintf(title, "Initial values for %.2f GeV", genMom);
  TH1D initialP("Inital values", title, 100, 1-10*momRes, 1+10*momRes);
  TH1D ptHist("ptHist", "ptHist", 50, 1-5*momRes, 1+5*momRes);
  int nHits = recEvent->GetNbOfHits();
  //int nHits = 12;
  int nBins = 300;
  TH1D** xDeltaGenHist = new TH1D*[nHits];
  for (int i = 0;i < nHits; i++) {
    char title[256];
    sprintf(title, "xLayer %d", i+1);
    int nBins = 100;
    xDeltaGenHist[i] = new TH1D(title,title, nBins, -20, 20);
  }
  TH1D** yDeltaGenHist = new TH1D*[nHits];
  for (int i = 0;i < nHits; i++) {
    char title[256];
    sprintf(title, "Layer %d", i+1);
    if (i == 0 || i == nHits - 2) yDeltaGenHist[i] = new TH1D(title,title, 100, -0.2, 0.2);
    else yDeltaGenHist[i] = new TH1D(title,title, 100, -0.2, 0.2);
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
    if (i == 0 || i == nHits - 2) yDeltaSmearedHist[i] = new TH1D(title,title, 500, -1.0, 1.0);
    else yDeltaSmearedHist[i] = new TH1D(title,title, nBins, -1.0, 1.0);
  }

  TH1D totalXhist("totalXhist", "totalXhist", 500, -20, 20);
  TH1D totalYhist("totalYhist", "totalYhist", 500, -1.0, 1.0);
  sprintf(title, "#chi^{2} Distribution (dof = %d)", recEvent->GetDof());
  //sprintf(title, "#chi^{2} Distribution (dof = %d)", 11);
  TH1D chi2Hist("#chi^{2} Distribution", title, 200, 0.0, 200.0);
  TH1D dofHist("n_{dof} Disribution", title, 8, 0, 8);
  TH1D angleHist("angleHist", "angleHist", 500, -100e-3, 100e-3);
  TH1D lHist("lHist", "lHist", 100, 0.07, 0.11);
  TH1D lOverAngleHist("lOverAngleHist", "lOverAngleHist", 100, -5., -1.);

  int total = 0;
  int chi2Passed = 0;
  for(int iEvent = 0; iEvent < genTree->GetEntries(); iEvent++) {
    //  for(int i = 0; i < 90; i++) {
    genTree->GetEntry(iEvent);
    recTree->GetEntry(iEvent);
    int nHitsGen = genEvent->GetNbOfHits();
    int nHitsRec = recEvent->GetNbOfHits();

    double chi2Cut = 200;
    double chi2 = recEvent->GetChi2();
    total++;
    // if (chi2 > chi2Cut)
    //   continue;
    chi2Passed++;
    
    //if (nHitsGen < 8 || nHitsRec == 0 || nHitsGen > nHitsRec || chi2 > chi2Cut) continue;
    if (nHitsGen < 8 || nHitsRec == 0 || nHitsGen > nHitsRec) continue;

    double angle1 = (genEvent->GetHitPosition(3).y() - genEvent->GetHitPosition(0).y())/(genEvent->GetHitPosition(3).z() - genEvent->GetHitPosition(0).z());
    double angle2 = (genEvent->GetHitPosition(nHitsGen-1).y() - genEvent->GetHitPosition(4).y())/(genEvent->GetHitPosition(nHitsGen-1).z() - genEvent->GetHitPosition(4).z());

    double y0 = genEvent->GetHitPosition(4).y();;
    double y1 = genEvent->GetHitPosition(3).y();;
    double z0 = genEvent->GetHitPosition(4).z();;
    double z1 = genEvent->GetHitPosition(3).z();;

    double L = sqrt(pow(y1 - y0, 2.) + pow (z1 - z0, 2.)) / 1000.; //mm -> m
    lHist.Fill(L);

    // if (angle2 - angle1  -3e-3) continue;
    // angleHist.Fill(angle2-angle1);

    lOverAngleHist.Fill(L/(angle2-angle1));

    resHist.Fill(genEvent->GetMomentum()/recEvent->GetMomentum());
    ptHist.Fill(genEvent->GetTransverseMomentum()/recEvent->GetTransverseMomentum());

    for (int i = 0; i < nHitsGen; i++) {
      unsigned int genUniqueLayer = 2*genEvent->GetModuleID(i) + genEvent->GetLayerID(i);
      unsigned int iRec = genUniqueLayer;
      xDeltaGenHist[genUniqueLayer]->Fill(genEvent->GetHitPosition(i).x() - recEvent->GetHitPosition(iRec).x());
      yDeltaGenHist[genUniqueLayer]->Fill(genEvent->GetHitPosition(i).y() - recEvent->GetHitPosition(iRec).y());
      xDeltaGenHist[genUniqueLayer]->Fill(genEvent->GetHitPosition(i).x());
      yDeltaGenHist[genUniqueLayer]->Fill(genEvent->GetHitPosition(i).y());
      xDeltaSmearedHist[genUniqueLayer]->Fill(genEvent->GetSmearedHitPosition(i).x() - recEvent->GetHitPosition(iRec).x());
      yDeltaSmearedHist[genUniqueLayer]->Fill(genEvent->GetSmearedHitPosition(i).y() - recEvent->GetHitPosition(iRec).y());
      totalXhist.Fill(genEvent->GetHitPosition(i).x() - recEvent->GetHitPosition(iRec).x());
      totalYhist.Fill(genEvent->GetHitPosition(i).y() - recEvent->GetHitPosition(iRec).y());
    }


    initialP.Fill(recEvent->GetInitialParameter(0) * genEvent->GetMomentum());

    double dof = genEvent->GetNbOfHits() - 5;
     dofHist.Fill(dof);
    // if (dof == 5)
    chi2Hist.Fill(chi2);
    chi2Hist.SetTitle(title);
  }

  for (int i = 0; i < nHits; i++) {
    yDeltaGenHist[i]->Fit("gaus", "Q0");
    TF1* fitFunc = yDeltaGenHist[i]->GetFunction("gaus");
    if (fitFunc)
      std::cout << "y" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }

  gStyle->SetOptFit(11111);
  TCanvas canvas("canvas", "Momentum resolution", 1024, 768);
  canvas.Draw();
  //  canvas.Divide(1,2);
  //  canvas.cd(1);
  resHist.Draw();
  
  // double rangeLower = 1-2*momRes;
  // double rangeUpper = 1+2*momRes;
  // double rangeLower = 0.25;
  // double rangeUpper = 1.25;
  double sigma, error;
  //fitInChi2Range(&resHist, sigma, error);
  // resHist.Fit("gaus", "EQR", "", rangeLower, rangeUpper);
  resHist.Fit("gaus", "EQ");
  TF1* gausFit = resHist.GetFunction("gaus");
  //  gausFit->SetParName(2, "#sigma_{p}/p");
  double muFit = gausFit->GetParameter(1);
  double muFitErr= gausFit->GetParError(1);
  double sigmaFit = gausFit->GetParameter(2);
  double sigmaFitErr= gausFit->GetParError(2);
  double p_rec = genMom / muFit;
  double p_rec_err = genMom * muFitErr / (muFit*muFit);
  double sigmaP_overP = sigmaFit / muFit;
  double sigmaP_overP_err = sqrt( pow(sigmaFitErr/muFit, 2.) + pow(muFitErr*sigmaFit/(muFit*muFit), 2.) );
  char text1[128];

  double max = resHist.GetMaximum();
  sprintf(text1, "#sigma_{p}/p = %.2f #pm %.2f", sigmaP_overP, sigmaP_overP_err);
  TLatex latex1(-2, 0.7*max, text1);
  latex1.SetTextColor(kRed);
  latex1.SetTextFont(42);
  latex1.Draw("SAME");
  char text2[128];
  sprintf(text2, "p_{rec} = (%.2f #pm %.2f) GeV", p_rec, p_rec_err);
  TLatex latex2(-2, 0.85*max, text2);
  latex2.SetTextColor(kRed);
  latex2.SetTextFont(42);
  latex2.Draw("SAME");

  resHist.GetXaxis()->SetTitle("p_{nom}/p_{rec}");
  resHist.GetYaxis()->SetTitle("N");
  // canvas.cd(2);
  // ptHist.Draw();
  // ptHist.Fit("gaus", "EQR", "", 0.1, 1.+5*momRes);
  // ptHist.GetXaxis()->SetTitle("pt_{gen}/pt_{rec}");
  // ptHist.GetYaxis()->SetTitle("N");
  char stem[256];
  sprintf(stem,"neg_%.0fGeV", genMom);
  char saveName[256];
  sprintf(saveName, "%s_resolution.%s", stem, "pdf");
  canvas.SaveAs(saveName);
  sprintf(saveName, "%s_resolution.%s", stem, "root");
  canvas.SaveAs(saveName);


  TCanvas canvas2("canvas2", "x: Reconstructed vs generated position", 1024, 768);
  canvas2.Divide(2,nHits/2);
  canvas2.Draw();

  for (int i = 0; i < nHits; i++) {
    char xtitle[256];
    sprintf(xtitle, "(x_{%d,gen} - x_{%d,rec}) / mm", i+1, i+1);
    canvas2.cd(i+1);
    xDeltaGenHist[i]->Draw();
    xDeltaGenHist[i]->GetXaxis()->SetTitle(xtitle);
    xDeltaGenHist[i]->GetYaxis()->SetTitle("N");
    xDeltaGenHist[i]->Fit("gaus", "Q");
    TF1* fitFunc = xDeltaGenHist[i]->GetFunction("gaus");
    THistPainter* painter = (THistPainter*) xDeltaGenHist[i]->GetPainter();
    painter->PaintStat(1, fitFunc);
    TPaveStats* pt = (TPaveStats*) xDeltaGenHist[i]->GetListOfFunctions()->FindObject("stats");
    pt->SetY1NDC(0.45);
    if (fitFunc)
      std::cout << "x" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }

  TCanvas canvas3("canvas3", "y: Reconstructed vs generated position", 1024, 768);
  canvas3.Divide(2, nHits/2);
  canvas3.Draw();
  for (int i = 0; i < nHits; i++) {
    char ytitle[256];
    sprintf(ytitle, "(y_{%d,gen} - y_{%d,rec}) / mm", i+1, i+1);
    canvas3.cd(i+1);
    yDeltaGenHist[i]->Draw();
    yDeltaGenHist[i]->GetXaxis()->SetTitle(ytitle);
    yDeltaGenHist[i]->GetYaxis()->SetTitle("N");
    yDeltaGenHist[i]->Fit("gaus", "Q");
    TF1* fitFunc = yDeltaGenHist[i]->GetFunction("gaus");
    THistPainter* painter = (THistPainter*) yDeltaGenHist[i]->GetPainter();
    painter->PaintStat(1, fitFunc);
    TPaveStats* pt = (TPaveStats*) yDeltaGenHist[i]->GetListOfFunctions()->FindObject("stats");
    pt->SetY1NDC(0.45);
    if (fitFunc)
      std::cout << "y" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }

  TCanvas canvas4("canvas4", "x: Reconstructed vs measured position", 1024, 768);
  canvas4.Divide(nHits/2,2);
  canvas4.Draw();

  for (int i = 0; i < nHits; i++) {
    char xtitle[256];
    sprintf(xtitle, "(x_{%d,meas} - x_{%d,rec}) / mm", i+1, i+1);
    canvas4.cd(i+1);
    xDeltaSmearedHist[i]->Draw();
    xDeltaSmearedHist[i]->GetXaxis()->SetTitle(xtitle);
    xDeltaSmearedHist[i]->GetYaxis()->SetTitle("N");
    xDeltaSmearedHist[i]->Fit("gaus", "Q");
    TF1* fitFunc = xDeltaSmearedHist[i]->GetFunction("gaus");
    if (fitFunc)
      std::cout << "x" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
   }

  TCanvas canvas5("canvas5", "y: Reconstructed vs measured Position", 1024, 768);
  canvas5.Divide(nHits/2,2);
  canvas5.Draw();
  for (int i = 0; i < nHits; i++) {
    char ytitle[256];
    sprintf(ytitle, "(y_{%d,meas} - y_{%d,rec}) / mm", i+1, i+1);
    canvas5.cd(i+1);
    yDeltaSmearedHist[i]->Draw();
    yDeltaSmearedHist[i]->GetXaxis()->SetTitle(ytitle);
    yDeltaSmearedHist[i]->GetYaxis()->SetTitle("N");
    yDeltaSmearedHist[i]->Fit("gaus", "Q");
    TF1* fitFunc = yDeltaSmearedHist[i]->GetFunction("gaus");
    if (fitFunc)
      std::cout << "y" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
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
  chi2Dist.FixParameter(1, 5);
  //chi2Dist.FixParameter(1, recEvent->GetDof());
  //chi2Dist.FixParameter(1, 11);
  TCanvas canvas7("canvas7", "Chi2 distribution", 1024, 768);
  canvas7.Draw();
  chi2Hist.Draw();
  chi2Hist.Scale(1./(chi2Hist.Integral("WIDTH")));
  chi2Hist.GetXaxis()->SetTitle("#chi^{2}");
  chi2Hist.GetYaxis()->SetTitle("N");
  //  chi2Hist.Fit(&chi2Dist, "E");
  chi2Dist.Draw("SAME");
  chi2Dist.SetLineColor(kRed);
  sprintf(saveName, "%s_chi2.%s", stem, "pdf");
  //  canvas7.SaveAs(saveName);
  sprintf(saveName, "%s_chi2.%s", stem, "root");
  //  canvas7.SaveAs(saveName);

  TCanvas canvasDof("canvasDof", "Degrees of freedom");
  canvasDof.cd();
  dofHist.Draw();

  TCanvas canvas8("canvas8", "Deflection angle distribution", 1024, 768);
  canvas8.Draw();
  canvas8.Divide(3,1);
  canvas8.cd(1);
  angleHist.Draw();
  angleHist.GetXaxis()->SetTitle("#Delta #theta [rad]");
  angleHist.GetYaxis()->SetTitle("N");
  canvas8.cd(2);  
  lHist.Draw();
  canvas8.cd(3);  
  lOverAngleHist.Draw();

  TCanvas canvas9("canvas9", "initialP", 1024, 768);
  canvas9.Draw();
  initialP.Draw();
  initialP.Fit("gaus", "EQ");
  initialP.GetFunction("gaus")->SetParName(2, "#sigma_{p}/p");
  sprintf(saveName, "%s_initial.%s", stem, "pdf");
  //  canvas9.SaveAs(saveName);
  sprintf(saveName, "%s_initial.%s", stem, "root");
  //  canvas9.SaveAs(saveName);

  std::cout << (double)(chi2Passed) / (double)total << std::endl;
  std::cout << genTree->GetEntries() << std::endl;

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
