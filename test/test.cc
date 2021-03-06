// $Id: test.cc,v 1.15 2010/05/12 01:56:31 beischer Exp $

#include "RES_Event.hh"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

int main(int argc, char** argv)
{
  if (!argc == 2)
    std::cout << "please provide root file for analysis" << std::endl;
  const char* filename = argv[1];

  TFile* file = new TFile(filename, "READ");

  TTree* genTree = (TTree*) file->Get("resolution_gen_tree");
  if (genTree) {
    RES_Event* genEvent = new RES_Event();
      genTree->SetBranchAddress("event", &genEvent);
    std::cout << "Generated events:" << std::endl;
    std::cout << genTree->GetEntries() << std::endl;
    for (unsigned int i = 0; i < genTree->GetEntries(); i++) {
      //for (unsigned int i = 0; i < 1; i++) {
      genTree->GetEntry(i);
      std::cout << "--------------------------------------------------" << std::endl;
      std::cout << "ID: " << genEvent->GetID() << " --> "
                << " p: " << genEvent->GetMomentum()
                << " -----> chi2: " <<  genEvent->GetChi2() << std::endl;
      for (unsigned int j = 0; j < genEvent->GetNbOfHits(); j++) {
        std::cout << " (i,j): (" << genEvent->GetModuleID(j) << ", " << genEvent->GetLayerID(j) << ")"
                  << " x: " << genEvent->GetSmearedHitPosition(j).x() 
                  << " y: " << genEvent->GetSmearedHitPosition(j).y()
                  << " z: " << genEvent->GetSmearedHitPosition(j).z() << std::endl;
      }
    }
  }

  TTree* recTree = (TTree*) file->Get("resolution_rec_tree");
  if (recTree) {
    RES_Event* recEvent = new RES_Event();
    recTree->SetBranchAddress("event", &recEvent);
    std::cout << "Reconstructed events:" << std::endl;
    for (unsigned int i = 0; i < recTree->GetEntries(); i++) {
      //for (unsigned int i = 0; i < 1; i++) {
      recTree->GetEntry(i);
      std::cout << "--------------------------------------------------" << std::endl;
      std::cout << "ID: " << recEvent->GetID() << " --> "
                << " p: " << recEvent->GetMomentum()
                << " -----> chi2: " <<  recEvent->GetChi2() << std::endl;
      for (unsigned int j = 0; j < recEvent->GetNbOfHits(); j++) {
        std::cout << " nHits: " << recEvent->GetNbOfHits();
        std::cout << " (i,j): (" << recEvent->GetModuleID(j) << ", " << recEvent->GetLayerID(j) << ")"
                  << " x: " << recEvent->GetHitPosition(j).x() 
                  << " y: " << recEvent->GetHitPosition(j).y()
                  << " z: " << recEvent->GetHitPosition(j).z() << std::endl;
      }
    }
  }

  return 0;
}
