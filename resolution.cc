//test
// $Id: resolution.cc,v 1.14 2010/06/29 13:26:14 beischer Exp $

#include "RES_ApplicationManager.hh"

#include <iostream>
#include <globals.hh>

int main(int argc, char** argv)
{
  if (argc >= 3) {
    std::cout << "Please run without command line arguments to create a session, or give the name of a macro file." << std::endl;
    return -1;
  } 

  RES_ApplicationManager appManager(argc, argv);

  if (argc == 1) {
    appManager.CreateSession();
  }
  else if (argc == 2) {
    G4String string(argv[1]);
    appManager.RunBatchScript(argv[1]);
  }

  return 0;
}
