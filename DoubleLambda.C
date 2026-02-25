
// R__LOAD_LIBRARY(EvtGen)
// R__ADD_INCLUDE_PATH($EVTGEN_ROOT/include)

#include "Pythia8/Pythia.h"
// #include "Pythia8Plugins/ColourReconnectionHooks.h"
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
// #include "Pythia8Plugins/EvtGen.h"
// #include "EvtGen/EvtGen.hh"
// #include "fastjet/PseudoJet.hh"
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
// #include "fastjet/ClusterSequence.hh"
// #include "fastjet/ClusterSequenceArea.hh"
#include <cstdio> // needed for io
#include <ctime>
#include <iostream> // needed for io
#include <time.h>   /* time */
#include <valarray>
#include <chrono>
// #include <yaml.h>
// #include <stdio.h>
// #include <glib.h>
// #include <yaml-cpp/yaml.h>

using namespace Pythia8;

int DoubleLambda(int nEvents = 1000)
{

    auto start = std::chrono::high_resolution_clock::now();
    // Read JOBID from environment to set unique random seeds per job
    int jobid = 0;
    if (const char *env_p = std::getenv("JOBID"))
    {
        jobid = std::atoi(env_p);
    }

    // Derive a unique seed for Pythia
    int uniqueSeed = 101 + jobid; // 101 is your old seed
    std::cout << "Using jobid = " << jobid << ", seed = " << uniqueSeed << std::endl;

    Pythia8::Pythia pythia;
    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("Tune:pp = 14");
    // pythia.readString("HardQCD:all = on");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 13600");
    pythia.readString("PartonLevel:FSR = on");
    pythia.readString("PartonLevel:ISR = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(uniqueSeed));
    pythia.init();
    /*EvtGenDecays *evtgen = 0;
    evtgen = new EvtGenDecays(&pythia, "./DECAY_2010.DEC", "./evt.pdl");
    evtgen->readDecayFile("./DstarDecay.dec");*/

    const int kMaxParticles = 20;
    TTree *fTreeEvent;
    Int_t TreemultLambda{};
    Float_t TreeLambdaPx[kMaxParticles]{};
    Float_t TreeLambdaPy[kMaxParticles]{};
    Float_t TreeLambdaPz[kMaxParticles]{};
    Float_t TreeLambdaEta[kMaxParticles]{};
    Float_t TreeLambdaPhi[kMaxParticles]{};
    Bool_t TreePrimaryLambda[kMaxParticles]{};
    Bool_t TreeIsLambdaFromPair[kMaxParticles]{};
    vector<int> motherListLambda[kMaxParticles];
    vector<int> TreeIndexPair[kMaxParticles]; // to store the index of the AntiLambda that forms a pair with each Lambda, initialized to -1
    Int_t TreeNumberOfPairsLambda[kMaxParticles]{};

    Int_t TreemultAntiLambda{};
    Float_t TreeAntiLambdaPx[kMaxParticles]{};
    Float_t TreeAntiLambdaPy[kMaxParticles]{};
    Float_t TreeAntiLambdaPz[kMaxParticles]{};
    Float_t TreeAntiLambdaEta[kMaxParticles]{};
    Float_t TreeAntiLambdaPhi[kMaxParticles]{};
    Bool_t TreePrimaryAntiLambda[kMaxParticles]{};
    Bool_t TreeIsAntiLambdaFromPair[kMaxParticles]{};
    vector<int> motherListAntiLambda[kMaxParticles];
    vector<int> TreeIndexPairAnti[kMaxParticles]; // to store the index of the Lambda that forms a pair with each AntiLambda, initialized to -1
    Int_t TreeNumberOfPairsAntiLambda[kMaxParticles]{};

    Int_t TreemultTotalLambda{};

    fTreeEvent = new TTree("fTreeEvent", "Event");
    fTreeEvent->Branch("TreemultLambda", &TreemultLambda, "TreemultLambda/I");
    fTreeEvent->Branch("TreemultAntiLambda", &TreemultAntiLambda, "TreemultAntiLambda/I");
    fTreeEvent->Branch("TreemultTotalLambda", &TreemultTotalLambda, "TreemultTotalLambda/I");
    fTreeEvent->Branch("TreePrimaryLambda", &TreePrimaryLambda, "TreePrimaryLambda[TreemultLambda]/O");
    fTreeEvent->Branch("TreeIsLambdaFromPair", &TreeIsLambdaFromPair, "TreeIsLambdaFromPair[TreemultLambda]/O");
    fTreeEvent->Branch("TreeNumberOfPairsLambda", &TreeNumberOfPairsLambda, "TreeNumberOfPairsLambda[TreemultLambda]/I");
    fTreeEvent->Branch("TreeIndexPair", &TreeIndexPair, "TreeIndexPair[TreemultLambda]/I");
    fTreeEvent->Branch("TreeLambdaPx", &TreeLambdaPx, "TreeLambdaPx[TreemultLambda]/F");
    fTreeEvent->Branch("TreeLambdaPy", &TreeLambdaPy, "TreeLambdaPy[TreemultLambda]/F");
    fTreeEvent->Branch("TreeLambdaPz", &TreeLambdaPz, "TreeLambdaPz[TreemultLambda]/F");
    fTreeEvent->Branch("TreeLambdaEta", &TreeLambdaEta, "TreeLambdaEta[TreemultLambda]/F");
    fTreeEvent->Branch("TreeLambdaPhi", &TreeLambdaPhi, "TreeLambdaPhi[TreemultLambda]/F");
    fTreeEvent->Branch("TreePrimaryAntiLambda", &TreePrimaryAntiLambda, "TreePrimaryAntiLambda[TreemultAntiLambda]/O");
    fTreeEvent->Branch("TreeIsAntiLambdaFromPair", &TreeIsAntiLambdaFromPair, "TreeIsAntiLambdaFromPair[TreemultAntiLambda]/O");
    fTreeEvent->Branch("TreeIndexPairAnti", &TreeIndexPairAnti, "TreeIndexPairAnti[TreemultAntiLambda]/I");
    fTreeEvent->Branch("TreeAntiLambdaPx", &TreeAntiLambdaPx, "TreeAntiLambdaPx[TreemultAntiLambda]/F");
    fTreeEvent->Branch("TreeAntiLambdaPy", &TreeAntiLambdaPy, "TreeAntiLambdaPy[TreemultAntiLambda]/F");
    fTreeEvent->Branch("TreeAntiLambdaPz", &TreeAntiLambdaPz, "TreeAntiLambdaPz[TreemultAntiLambda]/F");
    fTreeEvent->Branch("TreeAntiLambdaEta", &TreeAntiLambdaEta, "TreeAntiLambdaEta[TreemultAntiLambda]/F");
    fTreeEvent->Branch("TreeAntiLambdaPhi", &TreeAntiLambdaPhi, "TreeAntiLambdaPhi[TreemultAntiLambda]/F");

    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        int lambdaNumber = 0;
        int antiLambdaNumber = 0;
        int totalLambdaNumber = 0;
        if (!pythia.next())
        {
            cout << "check wrong" << "\n";
            continue;
        }
        // evtgen->decay();
        for (int i = 0; i < pythia.event.size(); ++i)
        {
            Particle &part = pythia.event[i];
            Int_t ist = part.status();
            Int_t pid = part.id();
            if (std::fabs(part.eta()) > 0.8)
                continue;

            auto mother1 = part.mother1();
            auto mother2 = part.mother2();

            int motherId1 = -999;
            int motherId2 = -999;
            bool isPrimaryLambda = false;

            if (mother1 > 0)
                motherId1 = pythia.event[mother1].id();
            if (mother2 >= 0)
                motherId2 = pythia.event[mother2].id();

            if (std::abs(pid) != 3122)
                continue; // only keep Lambdas and AntiLambdas

            if (mother1 > 0 && mother2 == 0)
            {
                // cout << "Lambda from decay" << endl;
                isPrimaryLambda = false;
            }
            else if (mother1 > 0 && mother2 > 0 && mother1 < mother2 && std::abs(ist) >= 81 && std::abs(ist) <= 89)
            {
                //    cout << "Lambda from hard process" << endl;
                isPrimaryLambda = true;
            }
            else
            {
                // cout << "Lambda from other source" << endl;
                isPrimaryLambda = false;
            }

            /*
            Particle &partmother1 = pythia.event[mother1];
            Particle &partmother2 = pythia.event[mother2];
            auto gmother1 = partmother1.mother1();
            auto gmother2 = partmother1.mother2(); // these two look the same
            auto gmother3 = partmother2.mother1();
            auto gmother4 = partmother2.mother2(); // these two look the same
            Particle &partgmother1 = pythia.event[gmother1];
            Particle &partgmother2 = pythia.event[gmother3];

            if (std::abs(motherId1) <= 8 || std::abs(motherId2) <= 8)
            {
                // cout << "This Lambda is the daughter of light quarks.\n";
                auto ggmother1 = partgmother1.mother1();
                auto ggmother2 = partgmother2.mother1();
                // cout << "The two mothers are " << mother1 << " (PID = " << motherId1 << ") and " << mother2 << " (PID = " << motherId2 << "). " << endl;
                // cout << "The gmothers of mother 1 are: index " << gmother1 << " with PID " << pythia.event[gmother1].id() << " and " << gmother2 << " with PID " << pythia.event[gmother2].id() << "\n";
                // cout << "The gmothers of mother 2 are: index " << gmother3 << " with PID " << pythia.event[gmother3].id() << " and " << gmother4 << " with PID " << pythia.event[gmother4].id() << "\n";
            }
            else if (std::abs(motherId1) == 9 || std::abs(motherId2) == 9)
            {
                // cout << "This Lambda is the daughter of a gluon.\n";
                // cout << "Particle " << i << " with PID " << pid << " has mother1 index " << mother1 << " with PID " << motherId1 << "\n";
                // cout << "Particle " << i << " with PID " << pid << " has mother2 index " << mother2 << " with PID " << motherId2 << "\n";
            }
            else
            {
                // cout << "Particle " << i << " with PID " << pid << " has mother1 index " << mother1 << " with PID " << motherId1 << "\n";
                // cout << "Particle " << i << " with PID " << pid << " has mother2 index " << mother2 << " with PID " << motherId2 << "\n";
            }
                */

            if (pid == 3122) // Lambda
            {
                TreePrimaryLambda[lambdaNumber] = isPrimaryLambda;
                TreeLambdaPx[lambdaNumber] = part.px();
                TreeLambdaPy[lambdaNumber] = part.py();
                TreeLambdaPz[lambdaNumber] = part.pz();
                TreeLambdaEta[lambdaNumber] = part.eta();
                TreeLambdaPhi[lambdaNumber] = part.phi();
                motherListLambda[lambdaNumber] = part.motherList();
                lambdaNumber++;
                totalLambdaNumber++;
            }
            else if (pid == -3122) // AntiLambda
            {
                TreePrimaryAntiLambda[antiLambdaNumber] = isPrimaryLambda;
                TreeAntiLambdaPx[antiLambdaNumber] = part.px();
                TreeAntiLambdaPy[antiLambdaNumber] = part.py();
                TreeAntiLambdaPz[antiLambdaNumber] = part.pz();
                TreeAntiLambdaEta[antiLambdaNumber] = part.eta();
                TreeAntiLambdaPhi[antiLambdaNumber] = part.phi();
                motherListAntiLambda[antiLambdaNumber] = part.motherList();
                antiLambdaNumber++;
                totalLambdaNumber++;
            }
        }
        // check which Lambda - AntiLambda pairs have common mothers

        if (antiLambdaNumber >= 1 && lambdaNumber >= 1)
        {
            cout << "\nIn this event, there are " << lambdaNumber << " Lambdas and " << antiLambdaNumber << " AntiLambdas. " << endl;
            for (int i = 0; i < lambdaNumber; ++i)
            {
                int numberofPairs = 0;
                cout << "The Lambda n. " << i << endl;
                TreeIsLambdaFromPair[i] = false;
                cout << "AAA number of pairs " << numberofPairs << endl;
                if (TreePrimaryAntiLambda[i] != 1)
                {
                    TreeNumberOfPairsLambda[i] = numberofPairs;
                    continue; // only consider primary Lambdas for pairing
                }
                for (int j = 0; j < antiLambdaNumber; ++j)
                {
                    cout << "AA number of pairs " << numberofPairs << endl;
                    bool isPair = false;
                    if (TreePrimaryAntiLambda[j] != 1)
                        continue; // only consider primary AntiLambdas for pairing
                    for (const auto &a : motherListLambda[i])
                    {
                        for (const auto &b : motherListAntiLambda[j])
                        {
                            if (a == b && std::abs(pythia.event[a].status()) > 19) // check they have a common mother with status code indicating it's not a beam particle
                            {
                                isPair = true;
                                break;
                            }
                        }
                        if (isPair)
                            break;
                    }
                    cout << "A number of pairs " << numberofPairs << endl;
                    if (isPair)
                    {
                        cout << "isPair " << isPair << endl;
                        numberofPairs++;
                        cout << "number of pairs " << numberofPairs << endl;
                        TreeIsLambdaFromPair[i] = isPair;
                        TreeIsAntiLambdaFromPair[j] = isPair;
                        TreeIndexPair[i].push_back(j);     // store the index of the AntiLambda that forms a pair with this Lambda
                        TreeIndexPairAnti[j].push_back(i); // store the index of the Lambda that forms a pair with this AntiLambda
                    }
                }
                TreeNumberOfPairsLambda[i] = numberofPairs; // store the number of pairs that this Lambda forms with AntiLambdas in the event
            }
            for (int i = 0; i < lambdaNumber; ++i)
            {
                cout << "The Lambda n. " << i << " forms " << TreeIndexPair[i].size() << " pairs." << endl;
                cout << "isPair : " << TreeIsLambdaFromPair[i] << endl;
                cout << TreeNumberOfPairsLambda[i] << " pairs." << endl;
                if (TreeIsLambdaFromPair[i] && TreeNumberOfPairsLambda[i] == 0)
                { // non succede
                    cout << "Error: Lambda n. " << i << " is marked as from pair but has no pairs in TreeIndexPair." << endl;
                }
                if (!TreeIsLambdaFromPair[i] && TreeNumberOfPairsLambda[i] > 0)
                {
                    cout << "Error: Lambda n. " << i << " is marked as not from pair but has " << TreeNumberOfPairsLambda[i] << " pairs in TreeIndexPair." << endl;
                }
                cout << "with the antilambdas: " << endl;
                for (const auto &indexAnti : TreeIndexPair[i])
                {
                    cout << "AntiLambda n. " << indexAnti << endl;
                }
            }
            for (int j = 0; j < antiLambdaNumber; ++j)
            {
                cout << "The AntiLambda n. " << j << " forms " << TreeIndexPairAnti[j].size() << " pairs." << endl;
                cout << "with the lambdas: " << endl;
                for (const auto &indexLambda : TreeIndexPairAnti[j])
                {
                    cout << "Lambda n. " << indexLambda << endl;
                }
            }
        }
        if (totalLambdaNumber >= 2)
        {
            TreemultLambda = lambdaNumber;
            TreemultAntiLambda = antiLambdaNumber;
            TreemultTotalLambda = totalLambdaNumber;
            fTreeEvent->Fill();
        }
        for (int i = 0; i < lambdaNumber; ++i)
            TreeIndexPair[i].clear();
        for (int j = 0; j < antiLambdaNumber; ++j)
            TreeIndexPairAnti[j].clear();
        for (int i = 0; i < kMaxParticles; i++)
        {
            TreeIsLambdaFromPair[i] = false;
            TreeIsAntiLambdaFromPair[i] = false;
        }
    }

    TString Sfileout = Form("DoubleLambdatree_%d.root", nEvents);
    TFile *fout = new TFile(Sfileout, "recreate");
    fout->cd();
    fTreeEvent->Write();
    fout->Close();
    pythia.stat();

    auto end = std::chrono::high_resolution_clock::now();
    // Compute the elapsed time
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";

    cout << "I produced the file " << Sfileout << " with " << nEvents << " events." << endl;

    return 0;
}