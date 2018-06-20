#include<iostream>
#include<vector>
#include"TH1.h"
#include"TH2.h"
#include"TF1.h"
#include"TMath.h"
#include"TCanvas.h"
#include"TFile.h"
#include"TLatex.h"
#include"TStyle.h"
#include"TPad.h"
#include"TLegend.h"
#include"TPaveText.h"
#include"TAxis.h"
#include"TTree.h"
#include"TFitResultPtr.h"
#include"TFitResult.h"
#include"TString.h"

using namespace std;

void Resave_Ntuple_for_TMVA_signal()
{

    TFile *soubor = new TFile("../myOutput/2018-04-20_02-38_new_HFT_pT_bins_final/merge/output.root", "READ"); //original production (signal + background)
    TFile *out_file = new TFile("./output/Dpm_TMVA_signal.root", "RECREATE"); //input for TMVA from data (siglnal)

    //original TTree
    TTree *tree = (TTree*)soubor->Get("nt"); //TNtuple from PYTHIA+fast-sim 

    //define new TTrees for TMVA
    //separate tree for each pT and centrality bin
    int nCentBins = 4;
    const int nPtBins = 12;

    double pT_bins[nPtBins+1] = { 1., 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10. }; //pT binning

    TTree *TMVA_tree[nCentBins][nPtBins];

    for(unsigned int i = 0; i<nCentBins; i++)
    {
      for(unsigned int j = 0; j<nPtBins; j++)
      {
        TMVA_tree[i][j] = new TTree(Form("Dpm_TMVA_signal_cent_%i_pT_%i", i, j), Form("Dpm_TMVA_signal_cent_%i_pT_%i", i, j));
      }
    }


    //-----------------ORIGINAL (fast-sim)-------------------------------------
    //Pion1
	Float_t p1RPhi, p1REta, p1RPt, p1RDca; //, pi1_dedx_old, pi1_nSigma_old;

  Float_t p1Hft, p1Tpc;
	

	//Pion2
  Float_t p2RPhi, p2REta, p2RPt, p2RDca;

  Float_t p2Hft, p2Tpc;
	

	//Kaon
	Float_t kRPhi, kREta, kRPt, kRDca;

  Float_t kHft, kTpc;

	//dca, flag, prim. vertex
	Float_t dcaDaughters;
	//Int_t flag_old;
	//Float_t primVz_old;

	//D meson
	Float_t cosTheta, decayLength, rPhi, rEta, rY, rPt, rM, mdV0Max;

	//centrality
	Float_t cent;
  Float_t w;

	//---------------------FOR NEW (for TMVA input)-------------------------------

	//Pion1
	//Float_t pi1_runId, pi1_eventId;
	Float_t pi1_phi, pi1_eta, pi1_pt, pi1_dca; //, pi1_dedx, pi1_nSigma;
	//Float_t pi1_nHitFit, pi1_nHitsMax, pi1_nHitdedx;
	//Float_t pi1_TOFinvbeta, pi1_betaBase;

	//Pion2
	//Float_t pi2_runId, pi2_eventId;
	Float_t pi2_phi, pi2_eta, pi2_pt, pi2_dca; //, pi2_dedx, pi2_nSigma;
	//Float_t pi2_nHitFit, pi2_nHitsMax, pi2_nHitdedx;
	//Float_t pi2_TOFinvbeta, pi2_betaBase;

	//Kaon
	//Float_t k_runId, k_eventId;
	Float_t k_phi, k_eta, k_pt, k_dca; //, k_dedx, k_nSigma;
	//Float_t k_nHitFit, k_nHitsMax, k_nHitdedx;
	Float_t k_TOFinvbeta, k_betaBase;

	//dca, flag, prim. vertex
	Float_t mdcaMax;
	//Float_t flag;
	Float_t v0z;

	//D meson
	Float_t D_cos_theta, D_decayL, D_phi, D_eta, D_pt, D_mass, D_dV0Max;

	//centrality, refMult
	//Float_t mcentrality;
  Float_t mreweight;
  



    //---------------LOAD OLD BRANCHES------------------------
    //change old branches to fast-sim 
  //Pion1
	tree->SetBranchAddress("p1RPhi", &p1RPhi);			 			//Float_t pi1_phi
	tree->SetBranchAddress("p1REta", &p1REta);			 			//Float_t pi1_eta
	tree->SetBranchAddress("p1RPt", &p1RPt);				 			//Float_t pi1_pt
	tree->SetBranchAddress("p1RDca", &p1RDca);			 			//Float_t pi1_dca

  tree->SetBranchAddress("p1Hft", &p1Hft);
  tree->SetBranchAddress("p1Tpc", &p1Tpc);
	
	//Pion2
	tree->SetBranchAddress("p2RPhi", &p2RPhi);
	tree->SetBranchAddress("p2REta", &p2REta);
	tree->SetBranchAddress("p2RPt", &p2RPt);
	tree->SetBranchAddress("p2RDca", &p2RDca);

  tree->SetBranchAddress("p2Hft", &p2Hft);
  tree->SetBranchAddress("p2Tpc", &p2Tpc);

	//Kaon
	tree->SetBranchAddress("kRPhi", &kRPhi);
	tree->SetBranchAddress("kREta", &kREta);
	tree->SetBranchAddress("kRPt", &kRPt);
	tree->SetBranchAddress("kRDca", &kRDca);

  tree->SetBranchAddress("kHft", &kHft);
  tree->SetBranchAddress("kTpc", &kTpc);

	//dca, flag, prim. vertex
	tree->SetBranchAddress("dcaDaughters", &dcaDaughters); //Float_t mdcaMax
  tree->SetBranchAddress("v0z", &v0z);

	//D meson cosTheta, decayLength, rPhi, rEta, rPt, rM, mdV0Max
	tree->SetBranchAddress("cosTheta", &cosTheta); 	//Float_t D_theta
	tree->SetBranchAddress("decayLength", &decayLength);	//Float_t D_decayL
	tree->SetBranchAddress("rPhi", &rPhi);			//Float_t D_phi
	tree->SetBranchAddress("rEta", &rEta);			//Float_t pi2_eta
  tree->SetBranchAddress("rY", &rY);			//Float_t pi2_eta
	tree->SetBranchAddress("rPt", &rPt);				//Float_t pi2_pt
	tree->SetBranchAddress("rM", &rM);		//Float_t D_mass
	tree->SetBranchAddress("mdV0Max", &mdV0Max);	//Float_t D_dV0Max

	//centrality, refmult
	tree->SetBranchAddress("cent", &cent);		//Float_t cmentrality
  tree->SetBranchAddress("w", &w); //weight from pp 200 GeV



   //--------------------FOR NEW----------------------------------------------------------------
   for(unsigned int i = 0; i < nCentBins; i++)
   {
     for(unsigned int j=0; j < nPtBins; j++)
     {
        //Pion1
        TMVA_tree[i][j]->Branch("pi1_phi", &pi1_phi, "pi1_phi/F");			 			//Float_t pi1_phi
        TMVA_tree[i][j]->Branch("pi1_eta", &pi1_eta, "pi1_eta/F");			 			//Float_t pi1_eta
        TMVA_tree[i][j]->Branch("pi1_pt", &pi1_pt, "pi1_pt/F");				 			//Float_t pi1_pt
        TMVA_tree[i][j]->Branch("pi1_dca", &pi1_dca, "pi1_dca/F");			 			//Float_t pi1_dca

        //Pion2
        TMVA_tree[i][j]->Branch("pi2_phi", &pi2_phi, "pi2_phi/F");			 			//Float_t pi2_phi
        TMVA_tree[i][j]->Branch("pi2_eta", &pi2_eta, "pi2_eta/F");			 			//Float_t pi2_eta
        TMVA_tree[i][j]->Branch("pi2_pt", &pi2_pt, "pi2_pt/F");				 			//Float_t pi2_pt
        TMVA_tree[i][j]->Branch("pi2_dca", &pi2_dca, "pi2_dca/F");			 			//Float_t pi2_dca
        
        //Kaon
        
        TMVA_tree[i][j]->Branch("k_phi", &k_phi, "k_phi/F");			 			//Float_t k_phi
        TMVA_tree[i][j]->Branch("k_eta", &k_eta, "k_eta/F");			 			//Float_t k_eta
        TMVA_tree[i][j]->Branch("k_pt", &k_pt, "k_pt/F");				 			//Float_t k_pt
        TMVA_tree[i][j]->Branch("k_dca", &k_dca, "k_dca/F");			 			//Float_t k_dca
 
        //dca, flag, prim. vertex
        TMVA_tree[i][j]->Branch("mdcaMax", &mdcaMax, "mdcaMax/F"); //Float_t mdcaMax
        

        //D meson
        TMVA_tree[i][j]->Branch("D_cos_theta", &D_cos_theta, "D_cos_theta/F"); 	//Float_t D_theta
        TMVA_tree[i][j]->Branch("D_decayL", &D_decayL, "D_decayL/F");	//Float_t D_decayL
        TMVA_tree[i][j]->Branch("D_phi", &D_phi, "D_phi/F");			//Float_t D_phi
        TMVA_tree[i][j]->Branch("D_eta", &D_eta, "D_eta/F");			//Float_t pi2_eta
        TMVA_tree[i][j]->Branch("D_pt", &D_pt, "D_pt/F");				//Float_t pi2_pt
        TMVA_tree[i][j]->Branch("D_mass", &D_mass, "D_mass/F");		//Float_t D_mass
        TMVA_tree[i][j]->Branch("D_dV0Max", &D_dV0Max, "D_dV0Max/F");	//Float_t D_dV0Max

        TMVA_tree[i][j]->Branch("mreweight", &mreweight, "mreweight/F");	//Float_t D_dV0Max

     }
   }




  //event, trackquality and PID cuts
  //Dpm
/*    Float_t D_y_cut = 1.; //check, if this is needed

    //TPC track quality cuts
    Int_t pi1_nHitFit_cut, pi2_nHitFit_cut, k_nHitFit_cut;
    pi1_nHitFit_cut = 20;
    pi2_nHitFit_cut = 20;
    k_nHitFit_cut = 20;

    Float_t pi1_nHitsMax_cut, pi2_nHitsMax_cut, k_nHitsMax_cut;
    pi1_nHitsMax_cut = 0.52;
    pi2_nHitsMax_cut = 0.52;
    k_nHitsMax_cut = 0.52;

    //PID cuts
    Float_t k_nSigma_cut; //pi_nSigma set in production code
    k_nSigma_cut = 2;

    Float_t pi1_TOFinvbeta_cut, pi2_TOFinvbeta_cut, k_TOFinvbeta_cut;
    pi1_TOFinvbeta_cut = 0.03;
    pi2_TOFinvbeta_cut = 0.03;
    k_TOFinvbeta_cut = 0.03;

    //double D_restMass = 1.86959; //rest mass in GeV
*/
    Long64_t nEntries = tree->GetEntries();

    cout<<nEntries<<endl;

    for(Long64_t i = 0; i < nEntries; i++)
    {
        if(i%100000 == 0)
        {
            cout<<i<<endl;
        }
        tree->GetEntry(i);

        //change PID, track and event quality cuts to fast-sim

        //Dmp rapidity cut
        if( !( TMath::Abs(rY) < 1. ) ) continue;

        if (TMath::Abs(v0z) > 60000.) continue;
    		
    		// |TPC Vz - VPD Vz| missing
    		//
    		// cent
    		//if (cent != con_cent) continue;
    		if (cent > 9) continue;
    		if (cent < -1) continue;
    		
    		if (fabs(kREta) > 1 || fabs(p1REta) > 1 || fabs(p2REta) > 1) continue;

/*      //probably do not need HFT and TPC matching here     		
        // HFT
    		if (kHft != 1 || p1Hft != 1 || p2Hft != 1) continue;    		
    		// TPC
    		if (kTpc != 1 || p1Tpc != 1 || p2Tpc != 1) continue;
    		//if (kTof != 1 || p1Tof != 1 || p2Tof != 1) continue;
*/    		
    		// kRPt
    		if (kRPt < 0.5) continue;
    		
    		// p1RPt
    		if (p1RPt < 0.5) continue;
    		
    		// p2RPt
    		if (p2RPt < 0.5) continue;


        //topological pre-cuts - same as in analysis production
        //to match variables ranges in simulation and background sample
        if (cosTheta < 0.997) continue;
    		
    		// dcaDaughters
    		if (dcaDaughters > 90) continue; //kuba
    		
    		// decayLength
    		if (decayLength < 30) continue;
    		
    		// kDca
    		if (kRDca < 70) continue; //orig. kDca
    		
    		// p1Dca
		    if (p1RDca < 90) continue; //orig. pi1Dca
		    
		    // p2Dca
		    if (p2RDca < 90) continue; //orig. pi2Dca

        if (mdV0Max > 220) continue;
    		

        
        pi1_phi = p1RPhi;
        pi1_eta = p1REta;
        pi1_pt  = p1RPt;
        pi1_dca = p1RDca/10000; // !!!!In simulation in mum, in data in cm!!!!

        pi2_phi = p2RPhi;
        pi2_eta = p2REta;
        pi2_pt  = p2RPt;
        pi2_dca = p2RDca/10000; // !!!!In simulation in mum, in data in cm!!!!

        k_phi = kRPhi;
        k_eta = kREta;
        k_pt  = kRPt;
        k_dca = kRDca/10000; // !!!!In simulation in mum, in data in cm!!!!

        mdcaMax = dcaDaughters/10000; // !!!!In simulation in mum, in data in cm!!!!
      //flag already saved on lines 554 - 564
        //primVz = primVz_old;

        D_cos_theta = cosTheta; //save as cos theta - in data (base production) is just theta!!!
        D_decayL = decayLength/10000; // !!!!In simulation in mum, in data in cm!!!!
        D_phi = rPhi;
        D_eta = rEta;
        D_pt = rPt;
        D_mass = rM;
        D_dV0Max = mdV0Max/10000; // !!!!In simulation in mum, in data in cm!!!!
        //mcentrality = mcentrality_old;

        mreweight = w;
        

        int centrality = -1; //default value
        int pT_bin = -1; //default value

        for(unsigned int l = 0; l < nPtBins; l++) //determine pT bin
        {
            if(rPt > pT_bins[l] && rPt <= pT_bins[l+1])
            {
                pT_bin = l;
            }
        }

        if( pT_bin < 0 ) continue;

        if(cent == 7 || cent == 8)  //centrality 0-10%
        {
          centrality = 0;
          TMVA_tree[centrality][pT_bin]->Fill();
        }

        if(cent == 4 || cent == 5 ||  cent == 6)  //centrality 10-40%
        {
          centrality = 1;
          TMVA_tree[centrality][pT_bin]->Fill();
        }

        if(cent == 0 || cent == 1 || cent == 2 || cent == 3)  //centrality 40-80%
        {
          centrality = 2;
          TMVA_tree[centrality][pT_bin]->Fill();
        }

        if(cent != -1 ) //centralitu 0-80%
        {
          centrality = 3;
          TMVA_tree[centrality][pT_bin]->Fill();
        }

    }//end of loop over TTree events


    for(unsigned int i = 0; i < nCentBins; i++ )
    {
      for(unsigned int j = 0; j< nPtBins; j++)
      {
        TMVA_tree[i][j]->Write();
      }
    }

//------------------------------------------------------------------------------------------------------------

    soubor->Close();
}
