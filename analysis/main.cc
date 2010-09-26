#include <TApplication.h>

#include <iostream>

#include "ResVsMom.hh"
#include "SingleFile.hh"
#include <MyROOTStyle.h>

int main(int argc, char* argv[])
{
  // if (argc != 2) {
  //   std::cerr << "please provide root file for analysis" << std::endl;
  //   return -1;
  // }

  const int nArgs = argc;
  const char* filename = argv[1];

  TApplication app("app", &argc, argv);
  MyROOTStyle myStyle("myStyle");
  myStyle.cd();

  SingleFile singleFile;
  ResVsMom resVsMom; 
    
  // if (nArgs == 2) {
  //   singleFile.processFile(filename);
  // }
  // else {
  resVsMom.processConfigFile(filename);
  // }

  app.Run();

  return 0;
}
