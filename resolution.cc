// $Id: resolution.cc,v 1.12 2010/01/11 14:47:40 beischer Exp $

#include "RES_ApplicationManager.hh"

#include <iostream>
#include <globals.hh>

int main(int argc, char** argv)
{
  if (argc != 2) {
    std::cout << "This program is meant for batch mode only. Please provide a script to run" << std::endl;
    return -1;
  } 

  G4String string(argv[1]);

  RES_ApplicationManager wrapper;
  wrapper.RunBatchScript(argv[1]);

  return 0;
}
