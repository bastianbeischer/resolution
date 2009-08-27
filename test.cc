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
                << " (i,j): (" << genEvent->GetModuleID(j) << ", " << genEvent->GetFiberID(j) << ")"
                << " x: " << genEvent->GetHitPosition(j).x() 
                << " y: " << genEvent->GetHitPosition(j).y()
                << " z: " << genEvent->GetHitPosition(j).z() << std::endl;

    }
  }

  std::cout << "Reconstructed events:" << std::endl;
  for (unsigned int i = 0; i < recTree->GetEntries(); i++) {
    recTree->GetEntry(i);
    for (unsigned int j = 0; j < recEvent->GetNbOfHits(); j++) {
      std::cout << "ID: " << recEvent->GetID() << " --> "
                << " p: " << recEvent->GetMomentum()
                << " (i,j): (" << recEvent->GetModuleID(j) << ", " << recEvent->GetFiberID(j) << ")"
                << " x: " << recEvent->GetHitPosition(j).x() 
                << " y: " << recEvent->GetHitPosition(j).y()
                << " z: " << recEvent->GetHitPosition(j).z()
                << " -----> chi2/dof: " <<  recEvent->GetChi2OverDof() << std::endl;
    }
  }

  return 0;
}
