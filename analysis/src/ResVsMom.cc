#include "ResVsMom.hh"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>

#include "RES_Event.hh"

#include <cmath>
#include <iostream>

int ResVsMom::m_markers[4] = {23, 22, 29, 27};
int ResVsMom::m_colors[4] = {kRed, kBlue, kBlack, kViolet};

ResVsMom::ResVsMom()
{
}

ResVsMom::~ResVsMom()
{
}

void ResVsMom::addGraph(const char* fileNameTemplate, double min, double max, double step, double guessA = 0.08, double guessB = 0.21)
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
    // double sigmaLeft = 2.5;
    // double sigmaRight = 2.5;
    // double rangeLower = 1-sigmaLeft*momRes;
    // double rangeUpper = 1+sigmaRight*momRes;
    //    resHist.Fit("gaus", "EQR0", "", rangeLower, rangeUpper);
    resHist.Fit("gaus", "EQ0");
    double sigma = resHist.GetFunction("gaus")->GetParameter(2);
    double sigmaErr = resHist.GetFunction("gaus")->GetParError(2);

    graph->SetPoint(i, value, sigma);
    graph->SetPointError(i, 0., sigmaErr);

    i++;

    file.Close();
  }

  graph->SetMarkerStyle(m_markers[m_graphs.size()%4]);
  graph->SetMarkerColor(m_colors[m_graphs.size()%4]);
  graph->SetTitle("Momentum resolution for PERDaix at 1^{o} stereo angle");

  m_graphs.push_back(graph);
}

void ResVsMom::processConfigFile()
{
  double momMin = 0.25;
  double momMax = 10.0;
  double momStep = 0.25;
  addGraph("../results/test.root", 0, 10, 1, 0.077, 0.18);
  //  addGraph("../results/test.root", 0, 10, 1, 0.077, 0.18);
  // addGraph("../results/perdaix_%.2f_GeV_1.00_deg_msc_inhom_electrons_windows.root", momMin, momMax, momStep, 0.077, 0.18);
  // addGraph("../results/perdaix_%.2f_GeV_1.00_deg_msc_inhom_protons_windows.root", momMin, momMax, momStep, 0.077, 0.23);
  //  addGraph("../results/perdaix_%.2f_GeV_1.00_deg_msc_inhom_protons_nowindows.root", momMin, momMax, momStep, 0.08, 0.30);
  // addGraph("../results/perdaix_%.2f_GeV_1.00_deg_nomsc_hom.root", momMin, momMax, momStep, 0.08, 0.0);

  draw();
}

void ResVsMom::draw()
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 1024, 768);
  canvas->Draw();
  canvas->cd();
  canvas->SetGridx();
  canvas->SetGridy();

  for (int i = 0; i < m_graphs.size(); i++) {
    if (i == 0) m_graphs.at(i)->Draw("AP");
    else        m_graphs.at(i)->Draw("P");
  }

  // double m_electron = 5.11e-4;
  // double m_proton   = 9.38e-1;

  // TF1 fitMSC("fitMSC", fitfunc, momMin, momMax, 3);
  // fitMSC.SetParNames("a", "b");
  // fitMSC.SetParNames("a", "b", "m");
  // TF1 fitNoMSC("fitNoMSC", "[0]*x + [1]", momMin, momMax);
  // fitNoMSC.SetParNames("a", "b");

  // graph1.Draw("AP");

  // //  fitMSC.FixParameter(2, .938);

  // fitMSC.FixParameter(2, m_electron);
  // graph1.Fit("fitMSC", "E");
  // graph1.SetMarkerSize(1.5);
  // graph1.GetFunction("fitMSC")->SetLineColor(kRed);
  // graph1.PaintStats(graph1.GetFunction("fitMSC"));
  // graph2.Draw("P");
  // fitMSC.FixParameter(2, m_proton);
  // graph2.Fit("fitMSC", "E");
  // graph2.SetMarkerSize(1.5);
  // graph2.GetFunction("fitMSC")->SetLineColor(kBlue);
  // graph2.PaintStats(graph2.GetFunction("fitMSC"));
  // // graph3.Draw("P");
  // // graph3.Fit("fitMSC", "E");
  // // graph3.SetMarkerSize(1.5);
  // // graph3.GetFunction("fitMSC")->SetLineColor(kBlack);
  // // graph3.PaintStats(graph3.GetFunction("fitMSC"));
  // // graph4.Draw("P");
  // // graph4.Fit("fitNoMSC", "E");
  // // graph4.SetMarkerSize(1.5);
  // // graph4.GetFunction("fitNoMSC")->SetLineColor(kViolet);
  // // graph4.PaintStats(graph4.GetFunction("fitNoMSC"));

  // graph1.GetXaxis()->SetTitle("p / GeV");
  // graph1.GetYaxis()->SetTitle("#sigma_{p} / p");
  // graph1.GetYaxis()->SetRangeUser(0.0, 1.2);
  // double upperRange = 10.0;
  // graph1.GetXaxis()->SetLimits(0., upperRange);
  // graph1.GetFunction("fitMSC")->SetRange(0., upperRange);
  // graph2.GetFunction("fitMSC")->SetRange(0., upperRange);
  // //  graph3.GetFunction("fitMSC")->SetRange(0., upperRange);
  // // graph4.GetFunction("fitNoMSC")->SetRange(0., upperRange);

  // TPaveStats* pt1 = (TPaveStats*) graph1.GetListOfFunctions()->FindObject("stats");
  // pt1->SetTextColor(kRed);
  // TPaveStats* pt2 = (TPaveStats*) graph2.GetListOfFunctions()->FindObject("stats");
  // pt2->SetTextColor(kBlue);
  // // TPaveStats* pt3 = (TPaveStats*) graph3.GetListOfFunctions()->FindObject("stats");
  // // pt3->SetTextColor(kBlack);
  // // TPaveStats* pt4 = (TPaveStats*) graph4.GetListOfFunctions()->FindObject("stats");
  // // pt4->SetTextColor(kViolet);
  // TPaveStats* pt[2] = {pt1, pt2};//, pt3, pt4};
  // for (int i = 0; i < 2; i++) {
  //   pt[i]->SetX1NDC(0.7);
  //   pt[i]->SetX2NDC(0.95);
  //   pt[i]->SetY1NDC(0.38 - i*0.12);
  //   pt[i]->SetY2NDC(0.48 - i*0.12);
  // }

  // TLatex text(3.1, 0.13, "#sigma_{p} / p = #sqrt{(ap)^{2} + (b/#beta)^{2}}");
  // text.Draw("SAME");

  //  TLegend legend(0.25, 0.75, 0.45, 0.88);
  //  legend.AddEntry(&graph1, "electrons", "P");
  //  legend.AddEntry(&graph2, "protons", "P");
  // // legend.AddEntry(&graph3, "module without windows", "P");
  // // legend.AddEntry(&graph4, "Multiple Scattering deactivated + Hom. approximation in center", "P");
  //  legend.Draw("SAME");

  // canvas.SaveAs("perdaix_1_deg.pdf");
  // canvas.SaveAs("perdaix_1_deg.root");
}
