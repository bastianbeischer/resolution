#include <TGraphErrors.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>

#include "RES_Event.hh"

void fillGraph(TGraphErrors& graph1, const char* fileNameTemplate, double min, double max, double step, double guessA = 0.08, double guessB = 0.21)
{
  TTree* genTree;
  TTree* recTree;
  RES_Event* genEvent = new RES_Event;
  RES_Event* recEvent = new RES_Event;

  int momBins = 120;
  int i = 0;
  for (double value = min; value <= max; value += step) {
    char filename[100];
    sprintf(filename, fileNameTemplate, value);

    TFile file(filename);
    if (file.IsZombie())
      continue;

    std::cout << filename << std::endl;

    genTree = (TTree*) file.Get("resolution_gen_tree");
    recTree = (TTree*) file.Get("resolution_rec_tree");
    genTree->SetBranchAddress("event", &genEvent);
    recTree->SetBranchAddress("event", &recEvent);
    genTree->GetEntry(0);
    double genMom = genEvent->GetMomentum()/1000.;
    double momRes = sqrt(pow(genMom*guessA, 2.) + pow(guessB,2.));
    TH1D resHist("resHist", "resHist", momBins, 1-10*momRes, 1+10*momRes);
    for(int j = 0; j < genTree->GetEntries(); j++) {
      genTree->GetEntry(j);
      recTree->GetEntry(j);
      resHist.Fill(genEvent->GetMomentum()/recEvent->GetMomentum());
    }
    // double sigmaLeft = 2.5;
    // double sigmaRight = 2.5;
    // double rangeLower = 1-sigmaLeft*momRes;
    // double rangeUpper = 1+sigmaRight*momRes;
    //    resHist.Fit("gaus", "EQR0", "", rangeLower, rangeUpper);
    resHist.Fit("gaus", "EQ0");
    double sigma = resHist.GetFunction("gaus")->GetParameter(2);
    double sigmaErr = resHist.GetFunction("gaus")->GetParError(2);

    graph1.SetPoint(i, value, sigma);
    graph1.SetPointError(i, 0., sigmaErr);

    i++;

    file.Close();
  }
}

double chi2dist(double*x, double*p)
{
  double amplitude = p[0];
  double ndf = p[1];
  double nom = pow(x[0],ndf/2. - 1.) * exp(-x[0]/2.);
  double denom = pow(2., ndf/2.) * TMath::Gamma(ndf/2.);
  return amplitude*nom/denom;
}

double MS(double p, double m, double L, double X0) 
{
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

