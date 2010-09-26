#include "SingleFile.hh"

#include <cmath>

#include <TCanvas.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <THistPainter.h>
#include <TPaveStats.h>
#include <TMath.h>

#include <MyROOTStyle.h>
#include "RES_Event.hh"

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

SingleFile::SingleFile() :
  m_nHits(0),
  m_genMom(0),
  m_file(0),
  m_genTree(0),
  m_recTree(0),
  m_genEvent(0),
  m_recEvent(0),
  m_resHist(0),
  m_initialPHist(0),
  m_ptHist(0),
  m_xDeltaGenHist(0),
  m_yDeltaGenHist(0),
  m_xDeltaSmearedHist(0),
  m_yDeltaSmearedHist(0),
  m_xTotalHist(0),
  m_yTotalHist(0),
  m_chi2Hist(0),
  m_dofHist(0),
  m_angleHist(0),
  m_lHist(0),
  m_lOverAngleHist(0),
  m_innerOrOuterHist(0)
{
  MyROOTStyle myStyle("myStyle");
  myStyle.cd();
}

SingleFile::~SingleFile()
{
  delete m_file;
  delete m_genTree;
  delete m_recTree;
  delete m_genEvent;
  delete m_recEvent;

  deleteHistograms();
}

void SingleFile::processFile(const char* filename)
{
  openFile(filename);
  setupHistograms();
  fillData();
  draw();
}

void SingleFile::openFile(const char* filename)
{
  m_file = new TFile(filename, "READ");
  m_genTree = (TTree*) m_file->Get("resolution_gen_tree");
  m_recTree = (TTree*) m_file->Get("resolution_rec_tree");

  m_genEvent = new RES_Event;
  m_recEvent = new RES_Event;
  m_genTree->SetBranchAddress("event", &m_genEvent);
  m_recTree->SetBranchAddress("event", &m_recEvent);

  m_genTree->GetEntry(0);
  m_recTree->GetEntry(0);
}

void SingleFile::deleteHistograms()
{
  if (m_resHist) delete m_resHist;
  if (m_initialPHist) delete m_initialPHist;
  if (m_ptHist) delete m_ptHist;
  for (int i = 0; i < m_nHits; i++) {
    if (m_xDeltaGenHist && m_xDeltaGenHist[i]) delete m_xDeltaGenHist[i];
    if (m_yDeltaGenHist && m_yDeltaGenHist[i]) delete m_yDeltaGenHist[i];
    if (m_xDeltaSmearedHist && m_xDeltaSmearedHist[i]) delete m_xDeltaSmearedHist[i];
    if (m_yDeltaSmearedHist && m_yDeltaSmearedHist[i]) delete m_yDeltaSmearedHist[i];
  }
  delete [] m_xDeltaGenHist;
  delete [] m_yDeltaGenHist;
  delete [] m_xDeltaSmearedHist;
  delete [] m_yDeltaSmearedHist;
  if (m_xTotalHist) delete m_xTotalHist;
  if (m_yTotalHist) delete m_yTotalHist;
  if (m_chi2Hist) delete m_chi2Hist;
  if (m_dofHist) delete m_dofHist;
  if (m_angleHist) delete m_angleHist;
  if (m_lHist) delete m_lHist;
  if (m_lOverAngleHist) delete m_lOverAngleHist;
  if (m_innerOrOuterHist) delete m_innerOrOuterHist;
}

void SingleFile::setupHistograms()
{
  deleteHistograms();

  m_genMom = m_genEvent->GetMomentum()/1000.;
  double momRes = sqrt(pow(m_genMom*0.08, 2.) + pow(0.21,2.));
  //double momRes = 0.8;
  //double momRes = sqrt(pow(genMom*.8e-3, 2.) + pow(0.04,2.));

  char title[256];
  sprintf(title, "Data with nominal energy of %.2f GeV", m_genMom);
  m_resHist = new TH1D("Mom. Resolution", title, 100, 1-10*momRes, 1+10*momRes);
  m_ptHist = new TH1D("ptHist", "ptHist", 50, 1-5*momRes, 1+5*momRes);
  
  sprintf(title, "Initial values for %.2f GeV", m_genMom);
  m_initialPHist = new TH1D("Inital values", title, 100, 1-10*momRes, 1+10*momRes);

  m_nHits = m_genEvent->GetNbOfHits();
  //int m_nHits = 12;
  int nBins = 300;
  m_xDeltaGenHist = new TH1D*[m_nHits];
  for (int i = 0;i < m_nHits; i++) {
    char title[256];
    sprintf(title, "xLayer %d", i+1);
    int nBins = 100;
    m_xDeltaGenHist[i] = new TH1D(title,title, nBins, -20, 20);
  }
  m_yDeltaGenHist = new TH1D*[m_nHits];
  for (int i = 0;i < m_nHits; i++) {
    char title[256];
    sprintf(title, "Layer %d", i+1);
    //   if (i == 0 || i == m_nHits - 2) m_yDeltaGenHist[i] = new TH1D(title,title, 100, -0.5, 0.5);
    m_yDeltaGenHist[i] = new TH1D(title,title, 100, -0.5, 0.5);
  }
  m_xDeltaSmearedHist = new TH1D*[m_nHits];
  for (int i = 0;i < m_nHits; i++) {
    char title[256];
    sprintf(title, "xDeltaSmearedHist%d", i);
    m_xDeltaSmearedHist[i] = new TH1D(title,title, nBins, -20, 20);
  }
  m_yDeltaSmearedHist = new TH1D*[m_nHits];
  for (int i = 0;i < m_nHits; i++) {
    char title[256];
    sprintf(title, "yDeltaSmearedHist%d", i);
    //    if (i == 0 || i == m_nHits - 2) m_yDeltaSmearedHist[i] = new TH1D(title,title, 500, -1.0, 1.0);
    m_yDeltaSmearedHist[i] = new TH1D(title,title, nBins, -1.0, 1.0);
  }

  m_xTotalHist = new TH1D("totalXhist", "totalXhist", 500, -20, 20);
  m_yTotalHist = new TH1D("totalYhist", "totalYhist", 500, -1.0, 1.0);
  sprintf(title, "#chi^{2} Distribution (dof = %d)", m_recEvent->GetDof());
  m_chi2Hist = new TH1D("#chi^{2} Distribution", title, 200, 0.0, 200.0);
  m_dofHist = new TH1D("n_{dof} Disribution", title, 8, 0, 8);
  m_angleHist = new TH1D("angleHist", "angleHist", 500, -100e-3, 100e-3);
  m_lHist = new TH1D("lHist", "lHist", 100, 0.07, 0.11);
  m_lOverAngleHist = new TH1D("lOverAngleHist", "lOverAngleHist", 100, -5., -1.);
  m_innerOrOuterHist = new TH1I("innerOrOuterHist", "innerOrOuterHist", 2, 0, 2);
}

void SingleFile::fillData()
{
  int total = 0;
  int chi2Passed = 0;
  for(int iEvent = 0; iEvent < m_genTree->GetEntries(); iEvent++) {
    m_genTree->GetEntry(iEvent);
    m_recTree->GetEntry(iEvent);
    int nHitsGen = m_genEvent->GetNbOfHits();
    int nHitsRec = m_recEvent->GetNbOfHits();

    double chi2 = m_recEvent->GetChi2();
    total++;
    // double chi2Cut = 200;
    // if (chi2 > chi2Cut)
    //   continue;
    chi2Passed++;
    
    if (nHitsGen < 8 || nHitsRec == 0 || nHitsGen > nHitsRec) continue;

    // double angle1 = (m_genEvent->GetHitPosition(3).y() - m_genEvent->GetHitPosition(0).y())/(m_genEvent->GetHitPosition(3).z() - m_genEvent->GetHitPosition(0).z());
    // double angle2 = (m_genEvent->GetHitPosition(nHitsM_Gen-1).y() - m_genEvent->GetHitPosition(4).y())/(m_genEvent->GetHitPosition(nHitsM_Gen-1).z() - m_genEvent->GetHitPosition(4).z());

    // double y0 = m_genEvent->GetHitPosition(4).y();;
    // double y1 = m_genEvent->GetHitPosition(3).y();;
    // double z0 = m_genEvent->GetHitPosition(4).z();;
    // double z1 = m_genEvent->GetHitPosition(3).z();;

    // double L = sqrt(pow(y1 - y0, 2.) + pow (z1 - z0, 2.)) / 1000.; //mm -> m
    // lHist->Fill(L);

    // // if (angle2 - angle1  -3e-3) continue;
    // // angleHist->Fill(angle2-angle1);

    // lOverAngleHist->Fill(L/(angle2-angle1));

    m_resHist->Fill(m_genEvent->GetMomentum()/m_recEvent->GetMomentum());
    m_ptHist->Fill(m_genEvent->GetTransverseMomentum()/m_recEvent->GetTransverseMomentum());

    for (int i = 0; i < nHitsGen; i++) {

      int iModule = m_genEvent->GetModuleID(i);
      unsigned int genUniqueLayer;
      if (iModule >= 0 && iModule <= 2)
        genUniqueLayer  = 0;
      else if(iModule >= 3 && iModule <= 4)
        genUniqueLayer  = 1;
      else if(iModule >= 5 && iModule <= 6)
        genUniqueLayer  = 2;
      else if(iModule >= 7 && iModule <= 9)
        genUniqueLayer  = 3;

      genUniqueLayer = 2*genUniqueLayer + m_genEvent->GetLayerID(i);
      unsigned int iRec = genUniqueLayer;

      m_xDeltaGenHist[genUniqueLayer]->Fill(m_genEvent->GetHitPosition(i).x() - m_recEvent->GetHitPosition(iRec).x());
      m_yDeltaGenHist[genUniqueLayer]->Fill(m_genEvent->GetHitPosition(i).y() - m_recEvent->GetHitPosition(iRec).y());
      // xDeltaGenHist[genUniqueLayer]->Fill(m_genEvent->GetHitPosition(i).x());
      // yDeltaGenHist[genUniqueLayer]->Fill(m_genEvent->GetHitPosition(i).y());
      m_xDeltaSmearedHist[genUniqueLayer]->Fill(m_genEvent->GetSmearedHitPosition(i).x() - m_recEvent->GetHitPosition(iRec).x());
      m_yDeltaSmearedHist[genUniqueLayer]->Fill(m_genEvent->GetSmearedHitPosition(i).y() - m_recEvent->GetHitPosition(iRec).y());
      m_xTotalHist->Fill(m_genEvent->GetHitPosition(i).x() - m_recEvent->GetHitPosition(iRec).x());
      m_yTotalHist->Fill(m_genEvent->GetHitPosition(i).y() - m_recEvent->GetHitPosition(iRec).y());
      //  std::cout << genUniqueLayer << "  --> " <<m_genEvent->GetHitPosition(i).z() << std::endl;
    }

    TVector3 direction = m_genEvent->GetHitPosition(nHitsGen-1) - m_genEvent->GetHitPosition(0);
    direction = 1./direction.Mag() * direction;
    TVector3 start = m_genEvent->GetHitPosition(0);
    
    double l = 35 - m_genEvent->GetHitPosition(0).z();
    TVector3 currentPos = start + l*direction;
    double r1 = sqrt(pow(currentPos.x(), 2.) + pow(currentPos.y(), 2.))/10.;

    l = -35 - m_genEvent->GetHitPosition(0).z();
    currentPos = start + l*direction;
    double r2 = sqrt(pow(currentPos.x(), 2.) + pow(currentPos.y(), 2.))/10.;

    //    std::cout << "r1: " << r1 << "  ---- r2: " << r2 << std::endl;
    if (r1 < 7.5 && r2 < 7.5)
      m_innerOrOuterHist->Fill(1);
    else
      m_innerOrOuterHist->Fill(0);

    m_initialPHist->Fill(m_recEvent->GetInitialParameter(0) * m_genEvent->GetMomentum());

    double dof = m_genEvent->GetNbOfHits() - 5;
    m_dofHist->Fill(dof);
    // if (dof == 5)
    m_chi2Hist->Fill(chi2);
    char title[256];
    sprintf(title, "#chi^{2} Distribution (dof = %d)", m_recEvent->GetDof());
    m_chi2Hist->SetTitle(title);
  }
}

void SingleFile::draw()
{
  for (int i = 0; i < m_nHits; i++) {
    m_yDeltaGenHist[i]->Fit("gaus", "Q0");
    //TF1* fitFunc = m_yDeltaGenHist[i]->GetFunction("gaus");
    //if (fitFunc)
      //      std::cout << "y" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }

  gStyle->SetOptFit(11111);
  TCanvas* canvas = new TCanvas("canvas", "Momentum resolution", 1024, 768);
  canvas->Draw();
  //  canvas.Divide(1,2);
  //  canvas.cd(1);
  m_resHist->Draw();
  
  // double rangeLower = 1-2*momRes;
  // double rangeUpper = 1+2*momRes;
  // double rangeLower = 0.25;
  // double rangeUpper = 1.25;
  //double sigma, error;
  //fitInChi2Range(&resHist, sigma, error);
  // resHist.Fit("gaus", "EQR", "", rangeLower, rangeUpper);
  m_resHist->Fit("gaus", "EQ");
  TF1* gausFit = m_resHist->GetFunction("gaus");
  //  gausFit->SetParName(2, "#sigma_{p}/p");
  double muFit = gausFit->GetParameter(1);
  double muFitErr= gausFit->GetParError(1);
  double sigmaFit = gausFit->GetParameter(2);
  double sigmaFitErr= gausFit->GetParError(2);
  double p_rec = m_genMom / muFit;
  double p_rec_err = m_genMom * muFitErr / (muFit*muFit);
  double sigmaP_overP = sigmaFit / muFit;
  double sigmaP_overP_err = sqrt( pow(sigmaFitErr/muFit, 2.) + pow(muFitErr*sigmaFit/(muFit*muFit), 2.) );
  char text1[128];

  double max = m_resHist->GetMaximum();
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

  m_resHist->GetXaxis()->SetTitle("p_{nom}/p_{rec}");
  m_resHist->GetYaxis()->SetTitle("N");
  // canvas.cd(2);
  // ptHist.Draw();
  // ptHist.Fit("gaus", "EQR", "", 0.1, 1.+5*momRes);
  // ptHist.GetXaxis()->SetTitle("pt_{gen}/pt_{rec}");
  // ptHist.GetYaxis()->SetTitle("N");
  // char stem[256];
  // sprintf(stem,"neg_%.0fGeV", m_genMom);
  // char saveName[256];
  // sprintf(saveName, "%s_resolution.%s", stem, "pdf");
  // canvas->SaveAs(saveName);
  // sprintf(saveName, "%s_resolution.%s", stem, "root");
  // canvas->SaveAs(saveName);


  // TCanvas* canvas2 = new TCanvas("canvas2", "x: Reconstructed vs generated position", 1024, 768);
  // canvas2->Divide(2,m_nHits/2);
  // canvas2->Draw();

  // for (int i = 0; i < m_nHits; i++) {
  //   char xtitle[256];
  //   sprintf(xtitle, "(x_{%d,gen} - x_{%d,rec}) / mm", i+1, i+1);
  //   canvas2->cd(i+1);
  //   xDeltaGenHist[i]->Draw();
  //   xDeltaGenHist[i]->GetXaxis()->SetTitle(xtitle);
  //   xDeltaGenHist[i]->GetYaxis()->SetTitle("N");
  //   xDeltaGenHist[i]->Fit("gaus", "Q");
  //   TF1* fitFunc = xDeltaGenHist[i]->GetFunction("gaus");
  //   THistPainter* painter = (THistPainter*) xDeltaGenHist[i]->GetPainter();
  //   painter->PaintStat(1, fitFunc);
  //   TPaveStats* pt = (TPaveStats*) xDeltaGenHist[i]->GetListOfFunctions()->FindObject("stats");
  //   pt->SetY1NDC(0.45);
  //   if (fitFunc)
  //     std::cout << "x" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  // }

  TCanvas* canvas3 = new TCanvas("canvas3", "y: Reconstructed vs generated position", 1024, 768);
  canvas3->Divide(2, m_nHits/2);
  canvas3->Draw();
  for (int i = 0; i < m_nHits; i++) {
    char ytitle[256];
    sprintf(ytitle, "(y_{%d,gen} - y_{%d,rec}) / mm", i+1, i+1);
    canvas3->cd(i+1);
    m_yDeltaGenHist[i]->Draw();
    m_yDeltaGenHist[i]->GetXaxis()->SetTitle(ytitle);
    m_yDeltaGenHist[i]->GetYaxis()->SetTitle("N");
    m_yDeltaGenHist[i]->Fit("gaus", "Q");
    TF1* fitFunc = m_yDeltaGenHist[i]->GetFunction("gaus");
    THistPainter* painter = (THistPainter*) m_yDeltaGenHist[i]->GetPainter();
    painter->PaintStat(1, fitFunc);
    TPaveStats* pt = (TPaveStats*) m_yDeltaGenHist[i]->GetListOfFunctions()->FindObject("stats");
    pt->SetY1NDC(0.45);
    //  if (fitFunc)
      //      std::cout << "y" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }

  TCanvas* canvas4 = new TCanvas("canvas4", "x: Reconstructed vs measured position", 1024, 768);
  canvas4->Divide(m_nHits/2,2);
  canvas4->Draw();

  for (int i = 0; i < m_nHits; i++) {
    char xtitle[256];
    sprintf(xtitle, "(x_{%d,meas} - x_{%d,rec}) / mm", i+1, i+1);
    canvas4->cd(i+1);
    m_xDeltaSmearedHist[i]->Draw();
    m_xDeltaSmearedHist[i]->GetXaxis()->SetTitle(xtitle);
    m_xDeltaSmearedHist[i]->GetYaxis()->SetTitle("N");
    m_xDeltaSmearedHist[i]->Fit("gaus", "Q");
    // TF1* fitFunc = m_xDeltaSmearedHist[i]->GetFunction("gaus");
    // if (fitFunc)
      //      std::cout << "x" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
   }

  TCanvas* canvas5 = new TCanvas("canvas5", "y: Reconstructed vs measured Position", 1024, 768);
  canvas5->Divide(m_nHits/2,2);
  canvas5->Draw();
  for (int i = 0; i < m_nHits; i++) {
    char ytitle[256];
    sprintf(ytitle, "(y_{%d,meas} - y_{%d,rec}) / mm", i+1, i+1);
    canvas5->cd(i+1);
    m_yDeltaSmearedHist[i]->Draw();
    m_yDeltaSmearedHist[i]->GetXaxis()->SetTitle(ytitle);
    m_yDeltaSmearedHist[i]->GetYaxis()->SetTitle("N");
    m_yDeltaSmearedHist[i]->Fit("gaus", "Q");
    //TF1* fitFunc = m_yDeltaSmearedHist[i]->GetFunction("gaus");
    //if (fitFunc)
      //      std::cout << "y" << i  << " --> mu = " << fitFunc->GetParameter(1) << ", rms = " << fitFunc->GetParameter(2) << std::endl;
  }

  TCanvas* canvas6 = new TCanvas("canvas6", "Sum of reconstructed vs generated position histograms", 1024, 768);
  canvas6->Draw();
  canvas6->Divide(1,2);
  canvas6->cd(1);
  m_xTotalHist->Draw();
  m_xTotalHist->GetXaxis()->SetTitle("(x_{sim,total} - x_{rec,total}) / mm");
  m_xTotalHist->GetYaxis()->SetTitle("N");
  canvas6->cd(2);
  m_yTotalHist->Draw();
  m_yTotalHist->GetXaxis()->SetTitle("(y_{sim,total} - y_{rec,total}) / mm");
  m_yTotalHist->GetYaxis()->SetTitle("N");

  TF1 chi2Dist("chi2Dist", chi2dist, 0.0, 100.0, 2);
  chi2Dist.SetNpx(1000);
  
  chi2Dist.FixParameter(0, 1);
  chi2Dist.FixParameter(1, 5);
  //chi2Dist.FixParameter(1, recEvent->GetDof());
  //chi2Dist.FixParameter(1, 11);
  TCanvas* canvas7 = new TCanvas("canvas7", "Chi2 distribution", 1024, 768);
  canvas7->Draw();
  m_chi2Hist->Draw();
  m_chi2Hist->Scale(1./(m_chi2Hist->Integral("WIDTH")));
  m_chi2Hist->GetXaxis()->SetTitle("#chi^{2}");
  m_chi2Hist->GetYaxis()->SetTitle("N");
  //  chi2Hist.Fit(&chi2Dist, "E");
  chi2Dist.Draw("SAME");
  chi2Dist.SetLineColor(kRed);
  //sprintf(saveName, "%s_chi2.%s", stem, "pdf");
  //  canvas7->SaveAs(saveName);
  //sprintf(saveName, "%s_chi2.%s", stem, "root");
  //  canvas7->SaveAs(saveName);

  TCanvas* canvasDof = new TCanvas("canvasDof", "Degrees of freedom");
  canvasDof->cd();
  m_dofHist->Draw();

  TCanvas* canvas8 = new TCanvas("canvas8", "Deflection angle distribution", 1024, 768);
  canvas8->Draw();
  canvas8->Divide(3,1);
  canvas8->cd(1);
  m_angleHist->Draw();
  m_angleHist->GetXaxis()->SetTitle("#Delta #theta [rad]");
  m_angleHist->GetYaxis()->SetTitle("N");
  canvas8->cd(2);  
  m_lHist->Draw();
  canvas8->cd(3);  
  m_lOverAngleHist->Draw();

  TCanvas* canvas9 = new TCanvas("canvas9", "initialP", 1024, 768);
  canvas9->Draw();
  m_initialPHist->Draw();
  m_initialPHist->Fit("gaus", "EQ");
  m_initialPHist->GetFunction("gaus")->SetParName(2, "#sigma_{p}/p");
  //  sprintf(saveName, "%s_initial.%s", stem, "pdf");
  //  canvas9->SaveAs(saveName);
  //sprintf(saveName, "%s_initial.%s", stem, "root");
  //  canvas9->uSaveAs(saveName);

  TCanvas* canvas10 = new TCanvas("canvas10", "innerOuter", 1024, 768);
  canvas10->Draw();
  m_innerOrOuterHist->Draw();
  // innerOrOuterHist.Fit("gaus", "EQ");
  // innerOrOuterHist.GetFunction("gaus")->SetParName(2, "#sigma_{p}/p");
  //  sprintf(saveName, "%s_initial.%s", stem, "pdf");
  //  canvas10.SaveAs(saveName);
  // sprintf(saveName, "%s_initial.%s", stem, "root");
  //  canvas10.SaveAs(saveName);


  // std::cout << (double)(chi2Passed) / (double)total << std::endl;
  // std::cout << genTree->GetEntries() << std::endl;
}

