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

  TApplication app("app", &argc, argv);
  MyROOTStyle myStyle("myStyle");
  myStyle.cd();

  if (argc == 2) {
    const char* filename = argv[1];
    SingleFile singleFile;
    singleFile.processFile(filename);
  }
  else {
    ResVsMom resVsMom;
    resVsMom.processConfigFile();
  }

  app.Run();

  return 0;
}
