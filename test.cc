#include "RES_Event.hh"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

int main(int argc, char** argv)
{
  TFile* file = new TFile("results/test.root", "READ");
  TTree* genTree = (TTree*) file->Get("resolution_gen_tree");
  TTree* recTree = (TTree*) file->Get("resolution_rec_tree");
  RES_Event* genEvent = new RES_Event();
  RES_Event* recEvent = new RES_Event();
  genTree->SetBranchAddress("event", &genEvent);
  recTree->SetBranchAddress("event", &recEvent);

  std::cout << "Generated events:" << std::endl;
  for (unsigned int i = 0; i < genTree->GetEntries(); i++) {
    genTree->GetEntry(i);
    for (unsigned int j = 0; j < genEvent->GetNbOfHits(); j++) {
      std::cout << "ID: " << genEvent->GetID() << " --> "
                << " p: " << genEvent->GetMomentum()
                << " x: " << genEvent->GetHit(j).x() 
                << " y: " << genEvent->GetHit(j).y()
                << " z: " << genEvent->GetHit(j).z() << std::endl;

    }
  }

  std::cout << "Reconstructed events:" << std::endl;
  for (unsigned int i = 0; i < recTree->GetEntries(); i++) {
    recTree->GetEntry(i);
    for (unsigned int j = 0; j < recEvent->GetNbOfHits(); j++) {
      std::cout << "ID: " << recEvent->GetID() << " --> "
                << " p: " << recEvent->GetMomentum()
                << " x: " << recEvent->GetHit(j).x() 
                << " y: " << recEvent->GetHit(j).y()
                << " z: " << recEvent->GetHit(j).z() << std::endl;
    }
  }

  return 0;
}
