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

float getradialdistance(float phi1, float y1, float phi2, float y2)
{
  float dphi = std::abs(phi1 - phi2);
  if (dphi > M_PI)
    dphi = 2 * M_PI - dphi;
  float dy = y1 - y2;
  return std::sqrt(dphi * dphi + dy * dy);
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
  float y;
  float pt;
  bool isprimary;
  bool ispair;
  int positioninvector;
  vector<int> indexpair;
};

float maxY = 0.5;
void buildHistos(int cmEnergyGeV = 13600, int nEvents = 1000, float ptmin = 0.6, float ptmax = 4.0, float ymax = 0.5)
{

  if (ymax > maxY){
    cout <<"You set a looser y cut then what was set to produce the input file" << endl;
    return;
  }
  std::vector<std::vector<Candidate>> lambdas;
  std::vector<std::vector<Candidate>> antiLambdas;

  // relatiuve mass, relative distance, isPrimaryPair, hasCommonMother
  Int_t bins[4] = {200, 200, 2, 2};
  Double_t xmin[4] = {0.0, 0.0, 0.0, 0.0};
  Double_t xmax[4] = {6, 10, 2, 2};
  double value[4] = {-10.0, -10.0, -10.0, -10.0};
  THnSparseD hsparse_same("hsparse_same", "hsparse_same", 4, bins, xmin, xmax);

  // --- open tree file ---
  TString Sinputfile = Form("../DoubleLambdatree_%d_cmEnergy%d_maxY%.2f.root", nEvents, cmEnergyGeV, maxY);
  TFile *f = new TFile(Sinputfile, "");
  if (!f || f->IsZombie())
  {
    std::cerr << "Error opening file: " << Sinputfile << std::endl;
    return;
  }
  cout << "Opened file: " << Sinputfile << std::endl;
  TTree *t = (TTree *)f->Get("fTreeEvent");

  // --- branches (main tree) ---
  TBranch *b_TreemultLambda = nullptr;
  TBranch *b_TreemultAntiLambda = nullptr;
  TBranch *b_TreemultTotalLambda = nullptr;
  TBranch *b_TreePrimaryLambda = nullptr;
  TBranch *b_TreeIsLambdaFromPair = nullptr;
  TBranch *b_TreeLambdaPx = nullptr;
  TBranch *b_TreeLambdaPy = nullptr;
  TBranch *b_TreeLambdaPz = nullptr;
  TBranch *b_TreeLambdaMass = nullptr;
  TBranch *b_TreeLambdaEta = nullptr;
  TBranch *b_TreeLambdaPhi = nullptr;
  TBranch *b_TreePrimaryAntiLambda = nullptr;
  TBranch *b_TreeIsAntiLambdaFromPair = nullptr;
  TBranch *b_TreeAntiLambdaPx = nullptr;
  TBranch *b_TreeAntiLambdaPy = nullptr;
  TBranch *b_TreeAntiLambdaPz = nullptr;
  TBranch *b_TreeAntiLambdaMass = nullptr;
  TBranch *b_TreeAntiLambdaEta = nullptr;
  TBranch *b_TreeAntiLambdaPhi = nullptr;

  // leaf types
  Int_t TreemultLambda;
  Float_t TreeLambdaPx[1000];
  Float_t TreeLambdaPy[1000];
  Float_t TreeLambdaPz[1000];
  Float_t TreeLambdaMass[1000];
  Float_t TreeLambdaEta[1000];
  Float_t TreeLambdaPhi[1000];
  Bool_t TreePrimaryLambda[1000];
  Bool_t TreeIsLambdaFromPair[1000];
  vector<int> TreeIndexPair[1000];

  Int_t TreemultAntiLambda;
  Float_t TreeAntiLambdaPx[1000];
  Float_t TreeAntiLambdaPy[1000];
  Float_t TreeAntiLambdaPz[1000];
  Float_t TreeAntiLambdaMass[1000];
  Float_t TreeAntiLambdaEta[1000];
  Float_t TreeAntiLambdaPhi[1000];
  Bool_t TreePrimaryAntiLambda[1000];
  Bool_t TreeIsAntiLambdaFromPair[1000];
  vector<int> TreeIndexPairAnti[1000];

  Int_t TreemultTotalLambda;

  // set branch addresses (main tree)
  t->SetBranchAddress("TreemultLambda", &TreemultLambda, &b_TreemultLambda);
  t->SetBranchAddress("TreemultAntiLambda", &TreemultAntiLambda, &b_TreemultAntiLambda);
  t->SetBranchAddress("TreemultTotalLambda", &TreemultTotalLambda, &b_TreemultTotalLambda);
  t->SetBranchAddress("TreePrimaryLambda", &TreePrimaryLambda, &b_TreePrimaryLambda);
  t->SetBranchAddress("TreeIsLambdaFromPair", &TreeIsLambdaFromPair, &b_TreeIsLambdaFromPair);
  t->SetBranchAddress("TreeLambdaPx", &TreeLambdaPx, &b_TreeLambdaPx);
  t->SetBranchAddress("TreeLambdaPy", &TreeLambdaPy, &b_TreeLambdaPy);
  t->SetBranchAddress("TreeLambdaPz", &TreeLambdaPz, &b_TreeLambdaPz);
  t->SetBranchAddress("TreeLambdaMass", &TreeLambdaMass, &b_TreeLambdaMass);
  t->SetBranchAddress("TreeLambdaEta", &TreeLambdaEta, &b_TreeLambdaEta);
  t->SetBranchAddress("TreeLambdaPhi", &TreeLambdaPhi, &b_TreeLambdaPhi);
  t->SetBranchAddress("TreePrimaryAntiLambda", &TreePrimaryAntiLambda, &b_TreePrimaryAntiLambda);
  t->SetBranchAddress("TreeIsAntiLambdaFromPair", &TreeIsAntiLambdaFromPair, &b_TreeIsAntiLambdaFromPair);
  t->SetBranchAddress("TreeAntiLambdaPx", &TreeAntiLambdaPx, &b_TreeAntiLambdaPx);
  t->SetBranchAddress("TreeAntiLambdaPy", &TreeAntiLambdaPy, &b_TreeAntiLambdaPy);
  t->SetBranchAddress("TreeAntiLambdaPz", &TreeAntiLambdaPz, &b_TreeAntiLambdaPz);
  t->SetBranchAddress("TreeAntiLambdaMass", &TreeAntiLambdaMass, &b_TreeAntiLambdaMass);
  t->SetBranchAddress("TreeAntiLambdaEta", &TreeAntiLambdaEta, &b_TreeAntiLambdaEta);
  t->SetBranchAddress("TreeAntiLambdaPhi", &TreeAntiLambdaPhi, &b_TreeAntiLambdaPhi);

  int nevent = t->GetEntries();
  cout << "Number of events in tree: " << nevent << endl;
  cout << "This is the number of events where at least two (anti-)Lambdas are found" << endl;

  for (Long64_t jentry = 0; jentry < nevent; ++jentry)
  {

    t->GetEntry(jentry);
    int lambdamult = TreemultLambda;
    int antilambdamult = TreemultAntiLambda;

    Candidate lambdaCand;
    std::vector<Candidate> lambdaVec;
    int numberLambdas = 0;
    for (int i = 0; i < lambdamult; ++i)
    {
      lambdaCand.tvec.SetXYZM(TreeLambdaPx[i], TreeLambdaPy[i], TreeLambdaPz[i], TreeLambdaMass[i]);
      lambdaCand.eta = TreeLambdaEta[i];
      lambdaCand.phi = TreeLambdaPhi[i];
      lambdaCand.y = lambdaCand.tvec.Rapidity();
      lambdaCand.isprimary = TreePrimaryLambda[i];
      lambdaCand.ispair = TreeIsLambdaFromPair[i];
      lambdaCand.indexpair = TreeIndexPair[i];
      //cout << "TreeIndexPair " << TreeIndexPair[i].size() << endl;
      //I do not understand why it is always zero
      lambdaCand.positioninvector = i;
      lambdaCand.pt = lambdaCand.tvec.Pt();
      if (std::abs(lambdaCand.y) > ymax)
        continue; // discard if Lambda out of acceptance
      if (lambdaCand.pt < ptmin || lambdaCand.pt > ptmax)
        continue;
      numberLambdas++;
      lambdaVec.push_back(lambdaCand);
    }
    if (numberLambdas > 0)
      lambdas.push_back(lambdaVec);
    else
      lambdas.push_back(std::vector<Candidate>()); // push empty vector if no Lambdas after selections

    Candidate antiLambdaCand;
    std::vector<Candidate> antiLambdaVec;
    int numberAntiLambdas = 0;
    for (int j = 0; j < antilambdamult; ++j)
    {
      antiLambdaCand.tvec.SetXYZM(TreeAntiLambdaPx[j], TreeAntiLambdaPy[j], TreeAntiLambdaPz[j], TreeAntiLambdaMass[j]);
      antiLambdaCand.eta = TreeAntiLambdaEta[j];
      antiLambdaCand.phi = TreeAntiLambdaPhi[j];
      antiLambdaCand.y = antiLambdaCand.tvec.Rapidity();
      antiLambdaCand.isprimary = TreePrimaryAntiLambda[j];
      antiLambdaCand.ispair = TreeIsAntiLambdaFromPair[j];
      antiLambdaCand.indexpair = TreeIndexPairAnti[j];
      antiLambdaCand.positioninvector = j;
      antiLambdaCand.pt = antiLambdaCand.tvec.Pt();
      //cout << "TreeIndexPairAnti " << TreeIndexPairAnti[j].size() << endl;
      //I do not understand why it is always zero
      if (std::abs(antiLambdaCand.y) > ymax)
        continue; // discard if AntiLambda out of acceptance
      if (antiLambdaCand.pt < ptmin || antiLambdaCand.pt > ptmax)
        continue;
      numberAntiLambdas++;
      antiLambdaVec.push_back(antiLambdaCand);
    }
    if (numberAntiLambdas > 0)
      antiLambdas.push_back(antiLambdaVec);
    else
      antiLambdas.push_back(std::vector<Candidate>()); // push empty vector if no AntiLambdas after selections
  }

  std::cout << "Finished loading trees, now processing pairs" << std::endl;

  for (Long64_t jentry = 0; jentry < nevent; ++jentry)
  {
    auto lambdaVec = lambdas[jentry];
    auto antiLambdaVec = antiLambdas[jentry];
    for (auto &lambda : lambdaVec)
    {
      for (auto &antiLambda : antiLambdaVec)
      {
        float relative_distance = getradialdistance(lambda.phi, lambda.y, antiLambda.phi, antiLambda.y);
        float pair_invmass = getpairmass(lambda.tvec, antiLambda.tvec);
        bool pair_isprimary = lambda.isprimary && antiLambda.isprimary;
        bool pair_ispair = 0;
        //if (lambda.indexpair.size()!=0 && antiLambda.indexpair.size()!=0){
        //  pair_ispair = (pair_isprimary && lambda.indexpair[0] == antiLambda.positioninvector && antiLambda.indexpair[0] == lambda.positioninvector);
        //}
        pair_ispair = pair_isprimary && lambda.ispair && antiLambda.ispair;
        //BE CAREFUL: this guarantees that the lambda shares a common mother with an anti-lambda and the anti-lambda shares a common mother with a lambda. 
        //However, this pair might not come from the same mother. 
        //I need to check the indexpair to be sure, but it is always zero for some reason. 
        value[0] = relative_distance;
        value[1] = pair_invmass;
        value[2] = pair_isprimary ? 1.0 : 0.0;
        value[3] = pair_ispair ? 1.0 : 0.0;
        hsparse_same.Fill(value);
      }
    }
  }
  // dump to output file
  TString Sfileout = Form("../DoubleLambdaHistos_%d_cmEnergy%d_ymax%.2f.root", nEvents, cmEnergyGeV, ymax);
  TFile *fout = new TFile(Sfileout, "recreate");
  fout->cd();
  hsparse_same.Write();
  fout->Close();
  cout << "Histograms saved to file: " << Sfileout << endl;
}