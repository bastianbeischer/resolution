// $Id: resolution.cc,v 1.11 2009/12/14 08:52:50 beischer Exp $

#include "RES_ApplicationManager.hh"

#include <iostream>
#include <globals.hh>

int main(int argc, char** argv)
{
  if (!argc == 2) {
    std::cout << "This program is meant for batch mode only. Please provide a script to run" << std::endl;
  } 

  G4String string(argv[1]);

  RES_ApplicationManager wrapper;
  wrapper.RunBatchScript(argv[1]);

  return 0;
}
