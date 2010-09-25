#include <TApplication.h>

#include <iostream>

#include "Analysis.hh"


int main(int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "please provide root file for analysis" << std::endl;
    return -1;
  }
  const char* filename = argv[1];

  TApplication app("app", &argc, argv);

  Analysis ana;
  ana.processFile(filename);

  app.Run();

  return 0;
}
