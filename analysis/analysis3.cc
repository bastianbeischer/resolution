#include <iostream>
#include <fstream>

#include <MyROOTStyle.h>

#include <TApplication.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH2D.h>

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Need exactly one argument" << std::endl;
    return 0;
  }

  std::ifstream infile(argv[1]);
  if (!infile.is_open()) {
    std::cerr << "Error opening file" << std::endl;
    return 0;
  }

  TApplication app("app", &argc, argv);

  MyROOTStyle myStyle("myStyle");
  myStyle.SetOptStat(0);
  myStyle.cd();


  int n = 0;
  double x0,x1,theta0,theta1;
  infile >> n >> x0 >> x1 >> theta0 >> theta1;

  TH2D hist("hist","#chi^{2} distribution after fitting in bending plane", n+1, x0 - (x1-x0)/(2.*n), x1 + (x1-x0)/(2.*n), n+1, theta0 - (theta1-theta0)/(2.*n),theta1 + (theta1-theta0)/(2.*n));
  for (int i = 0; i <= n; i++) {
    for (int j = 0; j <= n; j++) {
      double x,theta, chi2;
      infile >> x >> theta >> chi2;
      //      std::cout << x << " " << theta << " " << chi2 << " " << std::endl;
      int bin = hist.FindBin(x,theta);
      hist.SetBinContent(bin,chi2);
    }
  }

  TCanvas can("can", "chi2 scan", 1024, 768);
  can.Draw();
  can.SetLogz(1);
  can.SetRightMargin(0.15);
  hist.Draw("cont4 z");
  hist.SetContour(99);
  hist.GetXaxis()->SetTitle("x [cm]");
  hist.GetYaxis()->SetTitle("#theta [rad]");
  hist.GetZaxis()->SetTitle("#chi^{2}");


  app.Run();
  return 1;
}
