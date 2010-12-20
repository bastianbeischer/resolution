#include <TApplication.h>

#include <iostream>
#include <string>

#include "ResVsMom.hh"
#include "SingleFile.hh"
#include <MyROOTStyle.h>

int main(int argc, char** argv)
{
  if (argc != 2) {
    std::cerr << "please provide ROOT file or conf file for analysis" << std::endl;
    return -1;
  }

  std::string filename = argv[1];

  int bargc = 1;
  char* bargv[1] = {(char*)"not_analysis"};
  TApplication app("app", &bargc, bargv);

  MyROOTStyle myStyle("myStyle");
  myStyle.cd();

  if (filename.compare(filename.size()-4, 4, "root") == 0) {
    SingleFile singleFile;
    singleFile.processFile(filename.c_str());
    app.Run();
  }
  else if (filename.compare(filename.size()-4, 4, "conf") == 0) {
    ResVsMom resVsMom; 
    resVsMom.processConfigFile(filename.c_str());
    app.Run();
  }

  return 0;
}
