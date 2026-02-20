
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
// #include <yaml.h>
// #include <stdio.h>
// #include <glib.h>
// #include <yaml-cpp/yaml.h>

using namespace Pythia8;

int DoubleLambda(int nEvents = 1000)
{

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

    TTree *fTreeEvent;
    Int_t TreemultLambda;
    Float_t TreeLambdaPx[1000];
    Float_t TreeLambdaPy[1000];
    Float_t TreeLambdaPz[1000];
    Float_t TreeLambdaMass[1000];
    Int_t TreeLambdaMother1[10000];
    Int_t TreeLambdaMother2[10000];
    Int_t TreeLambdaMother1ID[1000];
    Int_t TreeLambdaMother2ID[1000];
    Int_t TreechargeLambda[1000];
    Int_t TreemultAntiLambda;
    Float_t TreechargeAntiLambda[1000];
    Float_t TreeAntiLambdaPx[1000];
    Float_t TreeAntiLambdaPy[1000];
    Float_t TreeAntiLambdaPz[1000];
    Float_t TreeAntiLambdaMass[1000];
    Int_t TreeAntiLambdaMother1[10000];
    Int_t TreeAntiLambdaMother2[10000];
    Int_t TreeAntiLambdaMother1ID[1000];
    Int_t TreeAntiLambdaMother2ID[1000];
    Int_t TreemultTotalLambda;
    fTreeEvent = new TTree("fTreeEvent", "Event");
    fTreeEvent->Branch("TreemultLambda", &TreemultLambda, "TreemultLambda/I");
    fTreeEvent->Branch("TreemultAntiLambda", &TreemultAntiLambda, "TreemultAntiLambda/I");
    fTreeEvent->Branch("TreemultTotalLambda", &TreemultTotalLambda, "TreemultTotalLambda/I");
    fTreeEvent->Branch("TreeLambdaMother1", &TreeLambdaMother1, "TreeLambdaMother1[TreemultLambda]/I");
    fTreeEvent->Branch("TreeLambdaMother2", &TreeLambdaMother2, "TreeLambdaMother2[TreemultLambda]/I");
    fTreeEvent->Branch("TreeLambdaMother1ID", &TreeLambdaMother1ID, "TreeLambdaMother1ID[TreemultLambda]/I");
    fTreeEvent->Branch("TreeLambdaMother2ID", &TreeLambdaMother2ID, "TreeLambdaMother2ID[TreemultLambda]/I");
    fTreeEvent->Branch("TreeLambdaPx", &TreeLambdaPx, "TreeLambdaPx[TreemultLambda]/F");
    fTreeEvent->Branch("TreeLambdaPy", &TreeLambdaPy, "TreeLambdaPy[TreemultLambda]/F");
    fTreeEvent->Branch("TreeLambdaPz", &TreeLambdaPz, "TreeLambdaPz[TreemultLambda]/F");
    fTreeEvent->Branch("TreeLambdaMass", &TreeLambdaMass, "TreeLambdaMass[TreemultLambda]/F");
    fTreeEvent->Branch("TreechargeLambda", &TreechargeLambda, "TreechargeLambda[TreemultLambda]/I");
    fTreeEvent->Branch("TreechargeAntiLambda", &TreechargeAntiLambda, "TreechargeAntiLambda[TreemultAntiLambda]/F");
    fTreeEvent->Branch("TreeAntiLambdaPx", &TreeAntiLambdaPx, "TreeAntiLambdaPx[TreemultAntiLambda]/F");
    fTreeEvent->Branch("TreeAntiLambdaPy", &TreeAntiLambdaPy, "TreeAntiLambdaPy[TreemultAntiLambda]/F");
    fTreeEvent->Branch("TreeAntiLambdaPz", &TreeAntiLambdaPz, "TreeAntiLambdaPz[TreemultAntiLambda]/F");
    fTreeEvent->Branch("TreeAntiLambdaMother1", &TreeAntiLambdaMother1, "TreeAntiLambdaMother1[TreemultAntiLambda]/I");
    fTreeEvent->Branch("TreeAntiLambdaMother2", &TreeAntiLambdaMother2, "TreeAntiLambdaMother2[TreemultAntiLambda]/I");
    fTreeEvent->Branch("TreeAntiLambdaMother1ID", &TreeAntiLambdaMother1ID, "TreeAntiLambdaMother1ID[TreemultAntiLambda]/I");
    fTreeEvent->Branch("TreeAntiLambdaMother2ID", &TreeAntiLambdaMother2ID, "TreeAntiLambdaMother2ID[TreemultAntiLambda]/I");

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
            if (mother1 > 0)
            {
                motherId1 = pythia.event[mother1].id();
            }
            if (mother2 >= 0)
            {
                motherId2 = pythia.event[mother2].id();
            }

            if (pid == 3122) // Lambda
            {
                TreechargeLambda[lambdaNumber] = part.id();
                TreeLambdaPx[lambdaNumber] = part.px();
                TreeLambdaPy[lambdaNumber] = part.py();
                TreeLambdaPz[lambdaNumber] = part.pz();
                TreeLambdaMass[lambdaNumber] = part.m();
                TreeLambdaMother1[lambdaNumber] = part.mother1();
                TreeLambdaMother2[lambdaNumber] = part.mother2();
                TreeLambdaMother1ID[lambdaNumber] = motherId1;
                TreeLambdaMother2ID[lambdaNumber] = motherId2;
                lambdaNumber++;
                totalLambdaNumber++;
            }
            else if (pid == -3122) // AntiLambda
            {
                TreechargeAntiLambda[antiLambdaNumber] = part.id();
                TreeAntiLambdaPx[antiLambdaNumber] = part.px();
                TreeAntiLambdaPy[antiLambdaNumber] = part.py();
                TreeAntiLambdaPz[antiLambdaNumber] = part.pz();
                TreeAntiLambdaMass[antiLambdaNumber] = part.m();
                TreeAntiLambdaMother1[antiLambdaNumber] = part.mother1();
                TreeAntiLambdaMother2[antiLambdaNumber] = part.mother2();
                TreeAntiLambdaMother1ID[antiLambdaNumber] = motherId1;
                TreeAntiLambdaMother2ID[antiLambdaNumber] = motherId2;
                antiLambdaNumber++;
                totalLambdaNumber++;
            }
        }
        if (totalLambdaNumber >= 2)
        {
            TreemultLambda = lambdaNumber;
            TreemultAntiLambda = antiLambdaNumber;
            TreemultTotalLambda = totalLambdaNumber;
            fTreeEvent->Fill();
        }
    }

    TFile *fout = new TFile("DoubleLambdatree_test.root", "recreate");
    fout->cd();
    fTreeEvent->Write();
    fout->Close();
    pythia.stat();
    return 0;
}