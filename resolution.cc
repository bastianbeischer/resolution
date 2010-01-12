// $Id: resolution.cc,v 1.13 2010/01/12 14:32:23 beischer Exp $

#include "RES_ApplicationManager.hh"

#include <iostream>
#include <globals.hh>

int main(int argc, char** argv)
{
  if (argc >= 3) {
    std::cout << "Please run without command line arguments to create a session, or give the name of a macro file." << std::endl;
    return -1;
  } 

  RES_ApplicationManager appManager;

  if (argc == 1) {
    appManager.CreateSession();
  }
  else if (argc == 2) {
    G4String string(argv[1]);
    appManager.RunBatchScript(argv[1]);
  }

  return 0;
}
