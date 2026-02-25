#include "TClonesArray.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <cstdio> 
#include <ctime>
#include <iostream> 
#include <time.h> 
#include <valarray>
#include <chrono>

float getradialdistance(float phi1, float eta1, float phi2, float eta2)
{
  float dphi = std::abs(phi1 - phi2);
  if (dphi > M_PI)
    dphi = 2 * M_PI - dphi;
  float deta = eta1 - eta2;
  return std::sqrt(dphi * dphi + deta * deta);
}
float getpairmass(TLorentzVector part1,
                  TLorentzVector part2)
{
  TLorentzVector Sum = part1 + part2;
  return Sum.M();
}

struct Candidate
{
  TLorentzVector tvec;
  float eta;
  float phi;
  bool isprimary;
  bool ispair;
  int positioninvector; 
  vector<int> indexpair;
};

void buildHistos(int nEvents = 1000)
{

  std::vector<std::vector<Candidate>> lambdas;
  std::vector<std::vector<Candidate>> antiLambdas;

  // mass, Sigma pt, proton pt, mother flag, kstar, sigma charge, proton charge
  Int_t bins[7] = {700, 250, 250, 2, 300, 3, 3};
  Double_t xmin[7] = {1.0, 0.0, 0.0, 0.0, 0.0, -1.5, -1.5};
  Double_t xmax[7] = {1.4, 10, 10, 2.0, 3.0, 1.5, 1.5};
  double value[7] = {-10.0, -10.0, -10.0, -5.0, -10.0, 0.0, 0.0};
  THnSparseD hsparse_same("hsparse_same", "hsparse_same", 7, bins, xmin, xmax);
  //THnSparseD hsparse_mix("hsparse_mix", "hsparse_mix", 7, bins, xmin, xmax);
  
  // --- open tree file ---
  TString Sinputfile = Form("DoubleLambdatree_%d.root", nEvents);
  TFile *f = new TFile(Sinputfile, "");
  TTree *t = (TTree *)f->Get("fTreeEvent");

  // --- branches (main tree) ---
  TBranch *b_TreemultLambda = nullptr;
  TBranch *b_TreemultAntiLambda = nullptr;
  TBranch *b_TreemultTotalLambda = nullptr;
  TBranch *b_TreechargeLambda = nullptr;
  TBranch *b_TreePrimaryLambda = nullptr;
  TBranch *b_TreeIsLambdaFromPair = nullptr;
  TBranch *b_TreeLambdaMother1 = nullptr;
  TBranch *b_TreeLambdaMother2 = nullptr;
  TBranch *b_TreeLambdaMother1ID = nullptr;
  TBranch *b_TreeLambdaMother2ID = nullptr;
  TBranch *b_TreeLambdaPx = nullptr;
  TBranch *b_TreeLambdaPy = nullptr;
  TBranch *b_TreeLambdaPz = nullptr;
  TBranch *b_TreeLambdaMass = nullptr;
  TBranch *b_TreeLambdaEta = nullptr;
  TBranch *b_TreeLambdaPhi = nullptr;
  TBranch *b_TreechargeAntiLambda = nullptr;
  TBranch *b_TreePrimaryAntiLambda = nullptr;
  TBranch *b_TreeIsAntiLambdaFromPair = nullptr;
  TBranch *b_TreeAntiLambdaMother1 = nullptr;
  TBranch *b_TreeAntiLambdaMother2 = nullptr;
  TBranch *b_TreeAntiLambdaMother1ID = nullptr;
  TBranch *b_TreeAntiLambdaMother2ID = nullptr;
  TBranch *b_TreeAntiLambdaPx = nullptr;
  TBranch *b_TreeAntiLambdaPy = nullptr;
  TBranch *b_TreeAntiLambdaPz = nullptr;
  TBranch *b_TreeAntiLambdaMass = nullptr;
  TBranch *b_TreeAntiLambdaEta = nullptr;
  TBranch *b_TreeAntiLambdaPhi = nullptr;

  // leaf types
  Int_t TreemultLambda;
  Int_t TreechargeLambda[1000];
  Int_t TreeAncestorIndexLambda[1000];
  Float_t TreeLambdaPx[1000];
  Float_t TreeLambdaPy[1000];
  Float_t TreeLambdaPz[1000];
  Float_t TreeLambdaMass[1000];
  Float_t TreeLambdaEta[1000];
  Float_t TreeLambdaPhi[1000];
  Int_t TreeLambdaMother1[10000];
  Int_t TreeLambdaMother2[10000];
  Int_t TreeLambdaMother1ID[1000];
  Int_t TreeLambdaMother2ID[1000];
  Bool_t TreePrimaryLambda[1000];
  Bool_t TreeIsLambdaFromPair[1000];
  vector<int> TreeIndexPair[1000];

  Int_t TreemultAntiLambda;
  Float_t TreechargeAntiLambda[1000];
  Int_t TreeAncestorIndexAntiLambda[1000];
  Float_t TreeAntiLambdaPx[1000];
  Float_t TreeAntiLambdaPy[1000];
  Float_t TreeAntiLambdaPz[1000];
  Float_t TreeAntiLambdaMass[1000];
  Float_t TreeAntiLambdaEta[1000];
  Float_t TreeAntiLambdaPhi[1000];
  Int_t TreeAntiLambdaMother1[10000];
  Int_t TreeAntiLambdaMother2[10000];
  Int_t TreeAntiLambdaMother1ID[1000];
  Int_t TreeAntiLambdaMother2ID[1000];
  Bool_t TreePrimaryAntiLambda[1000];
  Bool_t TreeIsAntiLambdaFromPair[1000];
  vector<int> TreeIndexPairAnti[1000];

  Int_t TreemultTotalLambda;

  // set branch addresses (main tree)
  t->SetBranchAddress("TreemultLambda", &TreemultLambda, &b_TreemultLambda);
  t->SetBranchAddress("TreemultAntiLambda", &TreemultAntiLambda, &b_TreemultAntiLambda);
  t->SetBranchAddress("TreemultTotalLambda", &TreemultTotalLambda, &b_TreemultTotalLambda);
  t->SetBranchAddress("TreechargeLambda", &TreechargeLambda, &b_TreechargeLambda);
  t->SetBranchAddress("TreePrimaryLambda", &TreePrimaryLambda, &b_TreePrimaryLambda);
  t->SetBranchAddress("TreeIsLambdaFromPair", &TreeIsLambdaFromPair, &b_TreeIsLambdaFromPair);
  t->SetBranchAddress("TreeLambdaMother1", &TreeLambdaMother1, &b_TreeLambdaMother1);
  t->SetBranchAddress("TreeLambdaMother2", &TreeLambdaMother2, &b_TreeLambdaMother2);
  t->SetBranchAddress("TreeLambdaMother1ID", &TreeLambdaMother1ID, &b_TreeLambdaMother1ID);
  t->SetBranchAddress("TreeLambdaMother2ID", &TreeLambdaMother2ID, &b_TreeLambdaMother2ID);
  t->SetBranchAddress("TreeLambdaPx", &TreeLambdaPx, &b_TreeLambdaPx);
  t->SetBranchAddress("TreeLambdaPy", &TreeLambdaPy, &b_TreeLambdaPy);
  t->SetBranchAddress("TreeLambdaPz", &TreeLambdaPz, &b_TreeLambdaPz);
  t->SetBranchAddress("TreeLambdaMass", &TreeLambdaMass, &b_TreeLambdaMass);
  t->SetBranchAddress("TreeLambdaEta", &TreeLambdaEta, &b_TreeLambdaEta);
  t->SetBranchAddress("TreeLambdaPhi", &TreeLambdaPhi, &b_TreeLambdaPhi);
  t->SetBranchAddress("TreechargeAntiLambda", &TreechargeAntiLambda, &b_TreechargeAntiLambda);
  t->SetBranchAddress("TreePrimaryAntiLambda", &TreePrimaryAntiLambda, &b_TreePrimaryAntiLambda);
  t->SetBranchAddress("TreeIsAntiLambdaFromPair", &TreeIsAntiLambdaFromPair, &b_TreeIsAntiLambdaFromPair);
  t->SetBranchAddress("TreeAntiLambdaMother1", &TreeAntiLambdaMother1, &b_TreeAntiLambdaMother1);
  t->SetBranchAddress("TreeAntiLambdaMother2", &TreeAntiLambdaMother2, &b_TreeAntiLambdaMother2);
  t->SetBranchAddress("TreeAntiLambdaMother1ID", &TreeAntiLambdaMother1ID, &b_TreeAntiLambdaMother1ID);
  t->SetBranchAddress("TreeAntiLambdaMother2ID", &TreeAntiLambdaMother2ID, &b_TreeAntiLambdaMother2ID);
  t->SetBranchAddress("TreeAntiLambdaPx", &TreeAntiLambdaPx, &b_TreeAntiLambdaPx);
  t->SetBranchAddress("TreeAntiLambdaPy", &TreeAntiLambdaPy, &b_TreeAntiLambdaPy);
  t->SetBranchAddress("TreeAntiLambdaPz", &TreeAntiLambdaPz, &b_TreeAntiLambdaPz);
  t->SetBranchAddress("TreeAntiLambdaMass", &TreeAntiLambdaMass, &b_TreeAntiLambdaMass);
  t->SetBranchAddress("TreeAntiLambdaEta", &TreeAntiLambdaEta, &b_TreeAntiLambdaEta);
  t->SetBranchAddress("TreeAntiLambdaPhi", &TreeAntiLambdaPhi, &b_TreeAntiLambdaPhi);

  int nevent = t->GetEntries();

  for (Long64_t jentry = 0; jentry < nevent; ++jentry)
  {

    t->GetEntry(jentry);
    int lambdamult = TreemultLambda;
    int antilambdamult = TreemultAntiLambda;

    Candidate lambdaCand;
    std::vector<Candidate> lambdaVec;
    for (int i = 0; i < lambdamult; ++i)
    {
      lambdaCand.tvec.SetXYZM(TreeLambdaPx[i], TreeLambdaPy[i], TreeLambdaPz[i], TreeLambdaMass[i]);
      lambdaCand.eta = TreeLambdaEta[i];
      lambdaCand.phi = TreeLambdaPhi[i];
      lambdaCand.isprimary = TreePrimaryLambda[i];
      lambdaCand.ispair = TreeIsLambdaFromPair[i];
      lambdaCand.indexpair = TreeIndexPair[i];
      lambdaCand.positioninvector = i;
      lambdaVec.push_back(lambdaCand);
    }
    lambdas.push_back(lambdaVec);

    Candidate antiLambdaCand;
    std::vector<Candidate> antiLambdaVec;
    for (int j = 0; j < antilambdamult; ++j)
    {
      antiLambdaCand.tvec.SetXYZM(TreeAntiLambdaPx[j], TreeAntiLambdaPy[j], TreeAntiLambdaPz[j], TreeAntiLambdaMass[j]);
      antiLambdaCand.eta = TreeAntiLambdaEta[j];
      antiLambdaCand.phi = TreeAntiLambdaPhi[j];
      antiLambdaCand.isprimary = TreePrimaryAntiLambda[j];
      antiLambdaCand.ispair = TreeIsAntiLambdaFromPair[j];
      antiLambdaCand.indexpair = TreeIndexPairAnti[j];
      antiLambdaCand.positioninvector = j;
      antiLambdaVec.push_back(antiLambdaCand);
    }
    antiLambdas.push_back(antiLambdaVec);
  }
  std::cout << "Finished loading trees, now processing pairs, SE" << std::endl;
  std::cout << "Lambda sizes: " << lambdas.size() << ", AntiLambda sizes: " << antiLambdas.size() << std::endl;

  for (Long64_t jentry = 0; jentry < nevent; ++jentry)
  {
    auto lambdaVec = lambdas[jentry];
    auto antiLambdaVec = antiLambdas[jentry];
    for (auto &lambda : lambdaVec)
    {
      for (auto &antiLambda : antiLambdaVec)
      {
        // discard if mother or daughter out of eta acceptance
        // if (std::abs(pidau1->Eta()) > 0.8 || std::abs(lambda.tvec.Eta()) > 0.8 || std::abs(antiLambda.tvec.Eta()) > 0.8)
        //  continue;
        float relative_distance = getradialdistance(lambda.phi, lambda.eta, antiLambda.phi, antiLambda.eta);
        float pair_invmass = getpairmass(lambda.tvec, antiLambda.tvec);
        bool pair_isprimary = lambda.isprimary && antiLambda.isprimary;
        bool pair_ispair = (lambda.ispair && antiLambda.ispair && lambda.indexpair == antiLambda.positioninvector && antiLambda.indexpair == lambda.positioninvector);
        value[0] = relative_distance;
        value[1] = pair_invmass;
        value[2] = pair_isprimary ? 1.0 : 0.0;
        value[3] = pair_ispair ? 1.0 : 0.0;
        hsparse_same.Fill(value);
      }
    }
  }

  /*
  std::cout << "Finished same-event pairs, now processing pairs, ME" << std::endl;

  const int kMaxMixDepth = 4;
  for (Long64_t jentry = 0; jentry < nevent; ++jentry)
  {
    auto &lambdaVec = lambdas[jentry];
    int nMixPartners = std::min<Long64_t>(kMaxMixDepth, nevent - 1);
    for (int idepth = 1; idepth <= nMixPartners; ++idepth)
    {
      Long64_t mixEntry = (jentry + idepth) % nevent;
      auto &antiLambdaVec = antiLambdas[mixEntry];
      for (auto &lambda : lambdaVec)
      {
        for (auto &antiLambda : antiLambdaVec)
        {
          // create acceptance effect here if needed
          phasespaceevent.SetDecay(lambda.tvec, 2, masses);
          phasespaceevent.Generate();
          auto pidau1 = phasespaceevent.GetDecay(0); // pion

          if (std::abs(pidau1->Eta()) > 0.8 ||
              std::abs(lambda.tvec.Eta()) > 0.8 ||
              std::abs(antiLambda.tvec.Eta()) > 0.8)
            continue;

          float relative_momentum = getkstar(lambda.tvec, antiLambda.tvec);
          // mass, Lambda pt, AntiLambda pt, mother flag, kstar, Lambda charge, AntiLambda charge
          value[0] = lambda.tvec.M();
          value[1] = lambda.tvec.Pt();
          value[2] = antiLambda.tvec.Pt();
          value[3] = 1.5; // mixed event: uncorrelated
          value[4] = relative_momentum;
          value[5] = lambda.charge > 0 ? 1.0 : -1.0;
          value[6] = antiLambda.charge > 0 ? 1.0 : -1.0;

          if (std::abs(antiLambda.mother1Pdg) == 3122 || std::abs(antiLambda.mother1Pdg) == 3222)
          {
            continue;
          }

          if (reweightPt)
          {
            double wLambda = lambdaPtWeight->GetBinContent(lambdaPtWeight->FindBin(lambda.tvec.Pt()));
            double wAntiLambda = antiLambdaPtWeight->GetBinContent(antiLambdaPtWeight->FindBin(antiLambda.tvec.Pt()));
            hsparse_mix.Fill(value, wLambda * wAntiLambda);
          }
          else
          {
            hsparse_mix.Fill(value);
          }
        }
      }
    }
  }
*/
  // dump to output file
  TString Sfileout = Form("DoubleLambdaHistos_%d.root", nEvents);
  TFile *fout = new TFile(Sfileout, "recreate");
  fout->cd();
  hsparse_same.Write();
  //hsparse_mix.Write();
  fout->Close();
}