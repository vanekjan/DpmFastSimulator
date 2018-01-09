#define nt_cxx
#include "nt.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>

void nt::Loop()
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

	// generated Dpm histograms
    TH1F *h_dpm_pt[9];
	h_dpm_pt[0] = new TH1F("h_dpm_pt0","h_dpm_pt0",150, 0, 15);
	h_dpm_pt[1] = new TH1F("h_dpm_pt1","h_dpm_pt1",150, 0, 15);
	h_dpm_pt[2] = new TH1F("h_dpm_pt2","h_dpm_pt2",150, 0, 15);
	h_dpm_pt[3] = new TH1F("h_dpm_pt3","h_dpm_pt3",150, 0, 15);
	h_dpm_pt[4] = new TH1F("h_dpm_pt4","h_dpm_pt4",150, 0, 15);
	h_dpm_pt[5] = new TH1F("h_dpm_pt5","h_dpm_pt5",150, 0, 15);
	h_dpm_pt[6] = new TH1F("h_dpm_pt6","h_dpm_pt6",150, 0, 15);
	h_dpm_pt[7] = new TH1F("h_dpm_pt7","h_dpm_pt7",150, 0, 15);
	h_dpm_pt[8] = new TH1F("h_dpm_pt8","h_dpm_pt8",150, 0, 15);
	TH1F *h_dpm_phi = new TH1F("h_dpm_phi","h_dpm_phi",64, -3.2, 3.2);
	TH2F *h_dpm_eta_phi = new TH2F("h_dpm_eta_phi","h_dpm_eta_phi", 140, -1.2, 1.2, 64, -3.2, 3.2);
	h_dpm_pt[0]->Sumw2();
	h_dpm_pt[1]->Sumw2();
	h_dpm_pt[2]->Sumw2();
	h_dpm_pt[3]->Sumw2();
	h_dpm_pt[4]->Sumw2();
	h_dpm_pt[5]->Sumw2();
	h_dpm_pt[6]->Sumw2();
	h_dpm_pt[7]->Sumw2();
	h_dpm_pt[8]->Sumw2();
	h_dpm_phi->Sumw2();
	h_dpm_eta_phi->Sumw2();

	// reconstructed Dpm histograms
    TH1F *h_r_dpm_pt[9];
	h_r_dpm_pt[0] = new TH1F("h_r_dpm_pt0","h_r_dpm_pt0",150, 0, 15);
	h_r_dpm_pt[1] = new TH1F("h_r_dpm_pt1","h_r_dpm_pt1",150, 0, 15);
	h_r_dpm_pt[2] = new TH1F("h_r_dpm_pt2","h_r_dpm_pt2",150, 0, 15);
	h_r_dpm_pt[3] = new TH1F("h_r_dpm_pt3","h_r_dpm_pt3",150, 0, 15);
	h_r_dpm_pt[4] = new TH1F("h_r_dpm_pt4","h_r_dpm_pt4",150, 0, 15);
	h_r_dpm_pt[5] = new TH1F("h_r_dpm_pt5","h_r_dpm_pt5",150, 0, 15);
	h_r_dpm_pt[6] = new TH1F("h_r_dpm_pt6","h_r_dpm_pt6",150, 0, 15);
	h_r_dpm_pt[7] = new TH1F("h_r_dpm_pt7","h_r_dpm_pt7",150, 0, 15);
	h_r_dpm_pt[8] = new TH1F("h_r_dpm_pt8","h_r_dpm_pt8",150, 0, 15);
	TH1F *h_r_dpm_phi = new TH1F("h_r_dpm_phi","h_r_dpm_phi",64, -3.2, 3.2);
	TH2F *h_r_dpm_eta_phi = new TH2F("h_r_dpm_eta_phi","h_r_dpm_eta_phi", 140, -1.2, 1.2, 64, -3.2, 3.2);
	h_r_dpm_pt[0]->Sumw2();
	h_r_dpm_pt[1]->Sumw2();
	h_r_dpm_pt[2]->Sumw2();
	h_r_dpm_pt[3]->Sumw2();
	h_r_dpm_pt[4]->Sumw2();
	h_r_dpm_pt[5]->Sumw2();
	h_r_dpm_pt[6]->Sumw2();
	h_r_dpm_pt[7]->Sumw2();
	h_r_dpm_pt[8]->Sumw2();
	h_r_dpm_phi->Sumw2();
	h_r_dpm_eta_phi->Sumw2();

	// now cuts
	TH1D *n_cuts = new TH1D("n_cuts","n_cuts",100,0,100);

	Long64_t nentries = fChain->GetEntries();
	cout << "#events: " << nentries << endl;

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
if(jentry%(nentries/10)==0) std::cout<<((jentry+10)*100/nentries)<<" % done..."<<std::endl;		
nb = fChain->GetEntry(jentry);   nbytes += nb;
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
		// cosTheta
		if (cosTheta < con_cosTheta) continue;
		n_cuts->Fill(5);
		// dcaDaughters
		if (dcaDaughters > con_dcaDaughters) continue; //kuba
		n_cuts->Fill(6);
		// decayLength
		if (decayLength < con_decayLength) continue;
		n_cuts->Fill(7);
		// kDca
		if (kDca < con_kDca) continue;
		n_cuts->Fill(8);
		// p1Dca
		if (p1Dca < con_pDca) continue;
		n_cuts->Fill(9);
		// p2Dca
		if (p2Dca < con_pDca) continue;
		n_cuts->Fill(10);
		// V0max missing
		//
		// dcaDpmToPv not in our cuts, maybe incude
		//
		// kRPt
		if (kRPt < con_kRPt) continue;
		n_cuts->Fill(11);
		// p1RPt
		if (p1RPt < con_pRPt) continue;
		n_cuts->Fill(12);
		// p2RPt
		if (p2RPt < con_pRPt) continue;
		n_cuts->Fill(13);
		if (mdV0Max > con_dV0Max) continue;
		n_cuts->Fill(14);
		// dca daughters
		/*
		if (dcakp1 < con_dca_daughters) continue;
		n_cuts->Fill(14);
		if (dcap1p2 < con_dca_daughters) continue;
		n_cuts->Fill(15);
		if (dcap2k < con_dca_daughters) continue;
		n_cuts->Fill(16);
		*/

		// reconstructed Dpm 
		if(cent > -0.5 and cent < 0.5) h_r_dpm_pt[0]->Fill(pt);
		if(cent > 0.5 and cent < 1.5) h_r_dpm_pt[1]->Fill(pt);
		if(cent > 1.5 and cent < 2.5) h_r_dpm_pt[2]->Fill(pt);
		if(cent > 2.5 and cent < 3.5) h_r_dpm_pt[3]->Fill(pt);
		if(cent > 3.5 and cent < 4.5) h_r_dpm_pt[4]->Fill(pt);
		if(cent > 4.5 and cent < 5.5) h_r_dpm_pt[5]->Fill(pt);
		if(cent > 5.5 and cent < 6.5) h_r_dpm_pt[6]->Fill(pt);
		if(cent > 6.5 and cent < 7.5) h_r_dpm_pt[7]->Fill(pt);
		if(cent > 7.5 and cent < 8.5) h_r_dpm_pt[8]->Fill(pt);


		h_r_dpm_phi->Fill(phi);
		h_r_dpm_eta_phi->Fill(eta, phi);
		
	}
	TFile *f = new TFile(out_file_name, "recreate");
	h_dpm_pt[0]->Write();
	h_dpm_pt[1]->Write();
	h_dpm_pt[2]->Write();
	h_dpm_pt[3]->Write();
	h_dpm_pt[4]->Write();
	h_dpm_pt[5]->Write();
	h_dpm_pt[6]->Write();
	h_dpm_pt[7]->Write();
	h_dpm_pt[8]->Write();
	h_dpm_phi->Write();
	h_dpm_eta_phi->Write();
	h_r_dpm_pt[0]->Write();
	h_r_dpm_pt[1]->Write();
	h_r_dpm_pt[2]->Write();
	h_r_dpm_pt[3]->Write();
	h_r_dpm_pt[4]->Write();
	h_r_dpm_pt[5]->Write();
	h_r_dpm_pt[6]->Write();
	h_r_dpm_pt[7]->Write();
	h_r_dpm_pt[8]->Write();
	h_r_dpm_phi->Write();
	h_r_dpm_eta_phi->Write();
	n_cuts->Write();
	f->Close();
}
