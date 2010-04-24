#include <TGraphErrors.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>

#include "RES_Event.hh"

void fillGraph(TGraphErrors& graph1, const char* fileNameTemplate, double min, double max, double step, double guessA = 0.08, double guessB = 0.21)
{
  TTree* genTree;
  TTree* recTree;
  RES_Event* genEvent = new RES_Event;
  RES_Event* recEvent = new RES_Event;

  int momBins = 60;
  double sigmaLeft = 2.5;
  double sigmaRight = 2.5;
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
    TH1D resHist("resHist", "resHist", momBins, 1-5*momRes, 1+5*momRes);
    for(int j = 0; j < genTree->GetEntries(); j++) {
      genTree->GetEntry(j);
      recTree->GetEntry(j);
      resHist.Fill(genEvent->GetMomentum()/recEvent->GetMomentum());
    }
    double rangeLower = 1-sigmaLeft*momRes;
    double rangeUpper = 1+sigmaRight*momRes;
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
