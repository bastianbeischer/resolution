#include "ResVsMom.hh"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>

#include "RES_Event.hh"

#include <cmath>
#include <iostream>
#include <fstream>

int ResVsMom::m_markers[4] = {23, 22, 29, 27};
int ResVsMom::m_colors[4] = {kRed, kBlue, kBlack, kViolet};

double fitfunc(double* x, double* p)
{
  double m = p[2];
  double beta = x[0] / sqrt(pow(x[0],2.) + pow(m,2.));
  return sqrt(pow(x[0]*p[0],2.) + pow(p[1]/beta, 2.));
}

ResVsMom::ResVsMom() :
  m_canvas(0),
  m_text(0),
  m_legend(0)
{
}

ResVsMom::~ResVsMom()
{
  delete m_canvas;
  delete m_text;
  delete m_legend;

  for (unsigned int i = 0; i < m_graphs.size(); i++) {
    delete m_graphs.at(i);
  }
  m_graphs.clear();
}

void ResVsMom::addGraph(const char* fileNameTemplate, const char* title, ParticleType type, double min, double max, double step, double guessA = 0.08, double guessB = 0.21)
{
  TGraphErrors* graph = new TGraphErrors;

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

    resHist.Fit("gaus", "EQ0");
    double sigma = resHist.GetFunction("gaus")->GetParameter(2);
    double sigmaErr = resHist.GetFunction("gaus")->GetParError(2);

    graph->SetPoint(i, value, sigma);
    graph->SetPointError(i, 0., sigmaErr);

    i++;

    file.Close();
  }

  delete genEvent;
  delete recEvent;

  graph->SetMarkerStyle(m_markers[m_graphs.size()%4]);
  graph->SetMarkerColor(m_colors[m_graphs.size()%4]);
  graph->SetMarkerSize(1.5);
  graph->SetTitle(title);

  m_graphs.push_back(graph);

  // fit
  double m_electron = 5.11e-4;
  double m_proton   = 9.38e-1;

  TF1 fitMSC("fitMSC", fitfunc, min, max, 3);
  fitMSC.SetParNames("a", "b");
  fitMSC.SetParNames("a", "b", "m");
  // TF1 fitNoMSC("fitNoMSC", "[0]*x + [1]", min, max);
  // fitNoMSC.SetParNames("a", "b");

  if (type == electron) {
    fitMSC.FixParameter(2, m_electron);
  }
  else if (type == proton) {
    fitMSC.FixParameter(2, m_proton);
  }
  graph->Fit("fitMSC", "E");
  graph->GetFunction("fitMSC")->SetLineColor(graph->GetMarkerColor());
  
  graph->GetXaxis()->SetTitle("p / GeV");
  graph->GetYaxis()->SetTitle("#sigma_{p} / p");
  graph->GetYaxis()->SetRangeUser(0.0, 1.2);

}

void ResVsMom::processConfigFile(const char* filename)
{
  ifstream infile(filename);
  if (!infile.is_open()) {
    std::cout << "Error opening config file: " << filename << std::endl;
    return;
  }

  while (!infile.eof()) {
    
    char fileTemplate[128];
    char title[128];
    double momMin, momMax, momStep, guessA, guessB;
    int typeInt;
    infile >> fileTemplate >> title >> typeInt >> momMin >> momMax >> momStep >> guessA >> guessB;
    ParticleType type = (ParticleType)typeInt;
    addGraph(fileTemplate, title, type, momMin, momMax, momStep, guessA, guessB);

  }

  draw();
}

void ResVsMom::draw()
{
  if (m_canvas) delete m_canvas;
  m_canvas = new TCanvas("canvas", "canvas", 1024, 768);
  m_canvas->Draw();
  m_canvas->cd();
  m_canvas->SetGridx();
  m_canvas->SetGridy();

  for (unsigned int i = 0; i < m_graphs.size(); i++) {
    if (i == 0) m_graphs.at(i)->Draw("AP");
    else        m_graphs.at(i)->Draw("P");
    m_graphs.at(i)->PaintStats(m_graphs.at(i)->GetFunction("fitMSC"));
  }

  TPaveStats* pt[m_graphs.size()];

  if(m_text) delete m_text;
  m_text = new TLatex(3.1, 0.13, "#sigma_{p} / p = #sqrt{(ap)^{2} + (b/#beta)^{2}}");

  if(m_legend) delete m_legend;
  m_legend = new TLegend(0.25, 0.75, 0.45, 0.88);

  for (unsigned int i = 0; i < m_graphs.size(); i++) {
    TPaveStats* stats = (TPaveStats*) m_graphs.at(i)->GetListOfFunctions()->FindObject("stats");
    stats->SetTextColor(m_graphs.at(i)->GetMarkerColor());
    pt[i] = stats;
    pt[i]->SetX1NDC(0.7);
    pt[i]->SetX2NDC(0.95);
    pt[i]->SetY1NDC(0.38 - i*0.12);
    pt[i]->SetY2NDC(0.48 - i*0.12);

    m_legend->AddEntry(m_graphs.at(i), m_graphs.at(i)->GetTitle(), "P");
  }

  m_text->Draw("SAME");
  m_legend->Draw("SAME");

  // canvas->SaveAs("perdaix_1_deg.pdf");
  // canvas->SaveAs("perdaix_1_deg.root");
}
