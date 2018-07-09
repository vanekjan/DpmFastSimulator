#define nt_TMVA_Cuts_cxx
#include "nt_TMVA_Cuts.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>

//#include "./TMVA_cuts/analysis/tmvaCuts.h"
#include "./TMVA_cuts/_iteration_2/tmvaCuts.h"

using namespace tmvaCuts;

void nt_TMVA_Cuts::Loop()
{
	//   In a ROOT session, you can do:
	//      Root > .L nt.C
	//      Root > nt t
	//      Root > t.GetEntry(12); // Fill t data members with entry number 12
	//      Root > t.Show();       // Show values of entry 12
	//      Root > t.Show(16);     // Read and show values of entry 16
	//      Root > t.Loop();       // Loop on all entries
	//

	//     This is the loop skeleton where:
	//    jentry is the global entry number in the chain
	//    ientry is the entry number in the current Tree
	//  Note that the argument to GetEntry must be:
	//    jentry for TChain::GetEntry
	//    ientry for TTree::GetEntry and TBranch::GetEntry
	//
	//       To read only selected branches, Insert statements like:
	// METHOD1:
	//    fChain->SetBranchStatus("*",0);  // disable all branches
	//    fChain->SetBranchStatus("branchname",1);  // activate branchname
	// METHOD2: replace line
	//    fChain->GetEntry(jentry);       //read all branches
	//by  b_branchname->GetEntry(ientry); //read only this branch
	if (fChain == 0) return;

  //double pT_bins[13] = { 1., 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10. }; //defined in tmvaCuts.h

	// generated Dpm histograms
  TH1F *h_dpm_pt[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_dpm_pt[i] = new TH1F(Form("h_dpm_pt%i", i),Form("Generated_Dpm_cent_%i", i),nPtBins_TMVA, pT_bins_TMVA);
    h_dpm_pt[i]->Sumw2();

  }
	TH1F *h_dpm_phi = new TH1F("h_dpm_phi","Generated_Dpm_phi",64, -3.2, 3.2);
	TH2F *h_dpm_eta_phi = new TH2F("h_dpm_eta_phi","Generated_Dpm_eta_phi", 140, -1.2, 1.2, 64, -3.2, 3.2);

	h_dpm_phi->Sumw2();
	h_dpm_eta_phi->Sumw2();
  //---------------------------------------------------------------------------------------------
	// reconstructed Dpm histograms without TOF matching
  TH1F *h_r_dpm_pt[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_r_dpm_pt[i] = new TH1F(Form("h_r_dpm_pt%i", i),Form("Reconstructed_Dpm_no_TOF_cent_%i", i),nPtBins_TMVA, pT_bins_TMVA);
    h_r_dpm_pt[i]->Sumw2();

  }
	TH1F *h_r_dpm_phi = new TH1F("h_r_dpm_phi","Reconstructed_Dpm_no_TOF_phi",64, -3.2, 3.2);
	TH2F *h_r_dpm_eta_phi = new TH2F("h_r_dpm_eta_phi","Reconstructed_Dpm_no_TOF_eta_phi", 140, -1.2, 1.2, 64, -3.2, 3.2);

	h_r_dpm_phi->Sumw2();
	h_r_dpm_eta_phi->Sumw2();
  //----------------------------------------------------------------------------------------------------
  

	// now cuts
	TH1D *n_cuts = new TH1D("n_cuts","n_cuts",100,0,100);

	Long64_t nentries = fChain->GetEntries();
	cout << "#events: " << nentries << endl;

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
    if(jentry%(nentries/10)==0) std::cout<<((jentry+10)*100/nentries)<<" % done..."<<std::endl;		
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
		// if (Cut(ientry) < 0) continue;

		// generated Dpm
		if(cent > -0.5 and cent < 0.5) h_dpm_pt[0]->Fill(pt);
		if(cent > 0.5 and cent < 1.5) h_dpm_pt[1]->Fill(pt);
		if(cent > 1.5 and cent < 2.5) h_dpm_pt[2]->Fill(pt);
		if(cent > 2.5 and cent < 3.5) h_dpm_pt[3]->Fill(pt);
		if(cent > 3.5 and cent < 4.5) h_dpm_pt[4]->Fill(pt);
		if(cent > 4.5 and cent < 5.5) h_dpm_pt[5]->Fill(pt);
		if(cent > 5.5 and cent < 6.5) h_dpm_pt[6]->Fill(pt);
		if(cent > 6.5 and cent < 7.5) h_dpm_pt[7]->Fill(pt);
		if(cent > 7.5 and cent < 8.5) h_dpm_pt[8]->Fill(pt);
		h_dpm_phi->Fill(phi);
		h_dpm_eta_phi->Fill(eta, phi);

		//if (TMath::Abs(rV0z) > con_rV0z) continue;
		// |TPC Vz|
		if (TMath::Abs(v0z) > con_v0z) continue;
		n_cuts->Fill(1);
		// |TPC Vz - VPD Vz| missing
		//
		// cent
		//if (cent != con_cent) continue;
		if (cent > con_cent_up) continue;
		if (cent < con_cent_down) continue;
		n_cuts->Fill(2);
		if (fabs(kREta) > 1 || fabs(p1REta) > 1 || fabs(p2REta) > 1) continue;
		// HFT
		if (kHft != 1 || p1Hft != 1 || p2Hft != 1) continue;
		n_cuts->Fill(3);
		// TPC
		if (kTpc != 1 || p1Tpc != 1 || p2Tpc != 1) continue;
		//if (kTof != 1 || p1Tof != 1 || p2Tof != 1) continue;
		n_cuts->Fill(4);

    // kRPt
		if (kRPt < con_kRPt) continue;
		n_cuts->Fill(11);
		// p1RPt
		if (p1RPt < con_pRPt) continue;
		n_cuts->Fill(12);
		// p2RPt
		if (p2RPt < con_pRPt) continue;
		n_cuts->Fill(13);
		

    
    int ptBin = h_dpm_pt[0]->FindBin(rPt)-1; //find pT bin
    int centBin_9 = (int)cent; //set centrality bin

    int centBin = -1;

    if( cent == 7 || cent == 8 ) centBin = 0; //0-10%
    if( cent == 4 || cent == 5 || cent == 6  ) centBin = 1; //10-40%
    if( cent == 0 || cent == 1 || cent ==  2 || cent ==  3 ) centBin = 2; //40-80%
    if( cent != -1 ) centBin = 3; //0-80%


    if(ptBin < 0 || ptBin > nPtBins_TMVA-1) continue; //reject Dpm with pT out of range

    //centra and pT dependent topological cuts from TMVA		
    if (mdV0Max > D_dV0Max_cut[centBin][ptBin]*10000) continue; //simulation is in mum, data in cm - has to convert units to use the same file for cuts
		n_cuts->Fill(14);
    // cosTheta
		if (cosTheta < D_cos_theta_cut[centBin][ptBin]) continue;
		n_cuts->Fill(5);
		// dcaDaughters
		if (dcaDaughters > mdcaMax_cut[centBin][ptBin]*10000) continue; //kuba
		n_cuts->Fill(6);
		// decayLength
		if (decayLength < D_decayL_min_cut[centBin][ptBin]*10000) continue;
		n_cuts->Fill(7);
		// kDca
		if (kRDca < k_dca_cut[centBin][ptBin]*10000) continue; //orig. kDca
		n_cuts->Fill(8);
		// p1Dca
		if (p1RDca < pi1_dca_cut[centBin][ptBin]*10000) continue; //orig. pi1Dca
		n_cuts->Fill(9);
		// p2Dca
		if (p2RDca < pi2_dca_cut[centBin][ptBin]*10000) continue; //orig. pi2Dca
		n_cuts->Fill(10);
		
		
    // reconstructed Dpm without TOF matching
  	h_r_dpm_pt[centBin_9]->Fill(pt);		

		h_r_dpm_phi->Fill(phi);
		h_r_dpm_eta_phi->Fill(eta, phi);

   
		
	} //end entries loop
  //_____________________________________________________________________________

  //write histograms to a new file
	TFile *f = new TFile(out_file_name, "recreate");

  //generated Dpm
  for(unsigned int j = 0; j<9; j++ )
  {
    h_dpm_pt[j]->Write();
  }
	h_dpm_phi->Write();
	h_dpm_eta_phi->Write();

  //reconstructed Dpm, no TOF matching
  for(unsigned int j = 0; j<9; j++ )
  {
    h_r_dpm_pt[j]->Write();
  }
	h_r_dpm_phi->Write();
	h_r_dpm_eta_phi->Write();


  n_cuts->Write();


	f->Close();
}
