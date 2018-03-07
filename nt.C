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
  for(unsigned int i = 0; i<9; i++)
  {
    h_dpm_pt[i] = new TH1F(Form("h_dpm_pt%i", i),Form("Generated_Dpm_cent_%i", i),150, 0, 15);
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
    h_r_dpm_pt[i] = new TH1F(Form("h_r_dpm_pt%i", i),Form("Reconstructed_Dpm_no_TOF_cent_%i", i),150, 0, 15);
    h_r_dpm_pt[i]->Sumw2();

  }
	TH1F *h_r_dpm_phi = new TH1F("h_r_dpm_phi","Reconstructed_Dpm_no_TOF_phi",64, -3.2, 3.2);
	TH2F *h_r_dpm_eta_phi = new TH2F("h_r_dpm_eta_phi","Reconstructed_Dpm_no_TOF_eta_phi", 140, -1.2, 1.2, 64, -3.2, 3.2);

	h_r_dpm_phi->Sumw2();
	h_r_dpm_eta_phi->Sumw2();
  //----------------------------------------------------------------------------------------------------

  // reconstructed Dpm histograms with TOF matching:
  //at least 1 daughter has to have TOF
  TH1F *h_r_dpm_pt_TOF_match[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_r_dpm_pt_TOF_match[i] = new TH1F(Form("h_r_dpm_pt_TOF_match%i", i),Form("Reconstructed_Dpm_hybrid_TOF_cent_%i", i),150, 0, 15);
    h_r_dpm_pt_TOF_match[i]->Sumw2();

  }
	TH1F *h_r_dpm_phi_TOF_match = new TH1F("h_r_dpm_phi_TOF_match","Reconstructed_Dpm_hybrid_TOF_phi",64, -3.2, 3.2);
	TH2F *h_r_dpm_eta_phi_TOF_match = new TH2F("h_r_dpm_eta_phi_TOF_match","Reconstructed_Dpm_hybrid_TOF_eta_phi", 140, -1.2, 1.2, 64, -3.2, 3.2);

	h_r_dpm_phi_TOF_match->Sumw2();
	h_r_dpm_eta_phi_TOF_match->Sumw2();
  //----------------------------------------------------------------------------------------------------------------------------------

  //individual possible combiunations seperately:
  //only one daughter has TOF
  TH1F *h_r_dpm_pt_TOF_match_one[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_r_dpm_pt_TOF_match_one[i] = new TH1F(Form("h_r_dpm_pt_TOF_match_one%i", i),Form("Reconstructed_Dpm_one_TOF_cent_%i", i),150, 0, 15);
    h_r_dpm_pt_TOF_match_one[i]->Sumw2();

  }
	TH1F *h_r_dpm_phi_TOF_match_one = new TH1F("h_r_dpm_phi_TOF_match_one","Reconstructed_Dpm_one_TOF_phi",64, -3.2, 3.2);
	TH2F *h_r_dpm_eta_phi_TOF_match_one = new TH2F("h_r_dpm_eta_phi_TOF_match_one","Reconstructed_Dpm_one_TOF_eta_phi", 140, -1.2, 1.2, 64, -3.2, 3.2);

	h_r_dpm_phi_TOF_match_one->Sumw2();
	h_r_dpm_eta_phi_TOF_match_one->Sumw2();
  //______________________________________________________________________________________________________________________________________________
  //two of three daughter have TOF
  TH1F *h_r_dpm_pt_TOF_match_two[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_r_dpm_pt_TOF_match_two[i] = new TH1F(Form("h_r_dpm_pt_TOF_match_two%i", i),Form("Reconstructed_Dpm_two_TOF_cent_%i", i),150, 0, 15);
    h_r_dpm_pt_TOF_match_two[i]->Sumw2();

  }
	TH1F *h_r_dpm_phi_TOF_match_two = new TH1F("h_r_dpm_phi_TOF_match_two","Reconstructed_Dpm_two_TOF_phi",64, -3.2, 3.2);
	TH2F *h_r_dpm_eta_phi_TOF_match_two = new TH2F("h_r_dpm_eta_phi_TOF_match_two","Reconstructed_Dpm_two_TOF_eta_phi", 140, -1.2, 1.2, 64, -3.2, 3.2);

	h_r_dpm_phi_TOF_match_two->Sumw2();
	h_r_dpm_eta_phi_TOF_match_two->Sumw2();
  //______________________________________________________________________________________________________________________________________________
  //all three daughter have TOF
  TH1F *h_r_dpm_pt_TOF_match_three[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_r_dpm_pt_TOF_match_three[i] = new TH1F(Form("h_r_dpm_pt_TOF_match_three%i", i),Form("Reconstructed_Dpm_three_TOF_cent_%i", i),150, 0, 15);
    h_r_dpm_pt_TOF_match_three[i]->Sumw2();

  }
	TH1F *h_r_dpm_phi_TOF_match_three = new TH1F("h_r_dpm_phi_TOF_match_three","Reconstructed_Dpm_three_TOF_phi",64, -3.2, 3.2);
	TH2F *h_r_dpm_eta_phi_TOF_match_three = new TH2F("h_r_dpm_eta_phi_TOF_match_three","Reconstructed_Dpm_three_TOF_eta_phi", 140, -1.2, 1.2, 64, -3.2, 3.2);

	h_r_dpm_phi_TOF_match_three->Sumw2();
	h_r_dpm_eta_phi_TOF_match_three->Sumw2();
  //______________________________________________________________________________________________________________________________________________
  //none of daughters has TOF - reject all that have TOF
  TH1F *h_r_dpm_pt_TOF_match_zero[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_r_dpm_pt_TOF_match_zero[i] = new TH1F(Form("h_r_dpm_pt_TOF_match_zero%i", i),Form("Reconstructed_Dpm_zero_TOF_cent_%i", i),150, 0, 15);
    h_r_dpm_pt_TOF_match_zero[i]->Sumw2();

  }
	TH1F *h_r_dpm_phi_TOF_match_zero = new TH1F("h_r_dpm_phi_TOF_match_zero","Reconstructed_Dpm_zero_TOF_phi",64, -3.2, 3.2);
	TH2F *h_r_dpm_eta_phi_TOF_match_zero = new TH2F("h_r_dpm_eta_phi_TOF_match_zero","Reconstructed_Dpm_zero_TOF_eta_phi", 140, -1.2, 1.2, 64, -3.2, 3.2);

	h_r_dpm_phi_TOF_match_zero->Sumw2();
	h_r_dpm_eta_phi_TOF_match_zero->Sumw2();
  //------------------------------------------------------------------------------------------------------------------------------------------------

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

		// reconstructed Dpm without TOF matching
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

    //TOF matching, "full" TOF - all particles have to have TOF
    //if (kTof != 1 || p1Tof != 1 || p2Tof != 1) continue;
    //n_cuts->Fill(15);

    //TOF matching, "soft" TOF - one or more particles has to have TOF information
    if (kTof == 1 || p1Tof == 1 || p2Tof == 1)
    {
    n_cuts->Fill(15);

    // reconstructed Dpm with "soft" TOF matching
		if(cent > -0.5 and cent < 0.5) h_r_dpm_pt_TOF_match[0]->Fill(pt);
		if(cent > 0.5 and cent < 1.5) h_r_dpm_pt_TOF_match[1]->Fill(pt);
		if(cent > 1.5 and cent < 2.5) h_r_dpm_pt_TOF_match[2]->Fill(pt);
		if(cent > 2.5 and cent < 3.5) h_r_dpm_pt_TOF_match[3]->Fill(pt);
		if(cent > 3.5 and cent < 4.5) h_r_dpm_pt_TOF_match[4]->Fill(pt);
		if(cent > 4.5 and cent < 5.5) h_r_dpm_pt_TOF_match[5]->Fill(pt);
		if(cent > 5.5 and cent < 6.5) h_r_dpm_pt_TOF_match[6]->Fill(pt);
		if(cent > 6.5 and cent < 7.5) h_r_dpm_pt_TOF_match[7]->Fill(pt);
		if(cent > 7.5 and cent < 8.5) h_r_dpm_pt_TOF_match[8]->Fill(pt);

    h_r_dpm_phi_TOF_match->Fill(phi);
		h_r_dpm_eta_phi_TOF_match->Fill(eta, phi);
    }

    //TOF matching, fill only if only one of daughters has TOF
    if ( (kTof == 1 && p1Tof == 0 && p2Tof == 0) || (kTof == 0 && p1Tof == 1 && p2Tof == 0) || (kTof == 0 && p1Tof == 0 && p2Tof == 1) )
    {
    n_cuts->Fill(16);

    // reconstructed Dpm with TOF matching
		if(cent > -0.5 and cent < 0.5) h_r_dpm_pt_TOF_match_one[0]->Fill(pt);
		if(cent > 0.5 and cent < 1.5) h_r_dpm_pt_TOF_match_one[1]->Fill(pt);
		if(cent > 1.5 and cent < 2.5) h_r_dpm_pt_TOF_match_one[2]->Fill(pt);
		if(cent > 2.5 and cent < 3.5) h_r_dpm_pt_TOF_match_one[3]->Fill(pt);
		if(cent > 3.5 and cent < 4.5) h_r_dpm_pt_TOF_match_one[4]->Fill(pt);
		if(cent > 4.5 and cent < 5.5) h_r_dpm_pt_TOF_match_one[5]->Fill(pt);
		if(cent > 5.5 and cent < 6.5) h_r_dpm_pt_TOF_match_one[6]->Fill(pt);
		if(cent > 6.5 and cent < 7.5) h_r_dpm_pt_TOF_match_one[7]->Fill(pt);
		if(cent > 7.5 and cent < 8.5) h_r_dpm_pt_TOF_match_one[8]->Fill(pt);

    h_r_dpm_phi_TOF_match_one->Fill(phi);
		h_r_dpm_eta_phi_TOF_match_one->Fill(eta, phi);
    }

    //TOF matching, fill only if two of daughters have TOF
    if ( (kTof == 1 && p1Tof == 1 && p2Tof == 0) || (kTof == 1 && p1Tof == 0 && p2Tof == 1) || (kTof == 0 && p1Tof == 1 && p2Tof == 1) )
    {
    n_cuts->Fill(17);

    // reconstructed Dpm with TOF matching
		if(cent > -0.5 and cent < 0.5) h_r_dpm_pt_TOF_match_two[0]->Fill(pt);
		if(cent > 0.5 and cent < 1.5) h_r_dpm_pt_TOF_match_two[1]->Fill(pt);
		if(cent > 1.5 and cent < 2.5) h_r_dpm_pt_TOF_match_two[2]->Fill(pt);
		if(cent > 2.5 and cent < 3.5) h_r_dpm_pt_TOF_match_two[3]->Fill(pt);
		if(cent > 3.5 and cent < 4.5) h_r_dpm_pt_TOF_match_two[4]->Fill(pt);
		if(cent > 4.5 and cent < 5.5) h_r_dpm_pt_TOF_match_two[5]->Fill(pt);
		if(cent > 5.5 and cent < 6.5) h_r_dpm_pt_TOF_match_two[6]->Fill(pt);
		if(cent > 6.5 and cent < 7.5) h_r_dpm_pt_TOF_match_two[7]->Fill(pt);
		if(cent > 7.5 and cent < 8.5) h_r_dpm_pt_TOF_match_two[8]->Fill(pt);

    h_r_dpm_phi_TOF_match_two->Fill(phi);
		h_r_dpm_eta_phi_TOF_match_two->Fill(eta, phi);
    }

    //TOF matching, fill only if all three daughters have TOF
    if ( kTof == 1 && p1Tof == 1 && p2Tof == 1 )
    {
    n_cuts->Fill(18);

    // reconstructed Dpm with TOF matching
		if(cent > -0.5 and cent < 0.5) h_r_dpm_pt_TOF_match_three[0]->Fill(pt);
		if(cent > 0.5 and cent < 1.5) h_r_dpm_pt_TOF_match_three[1]->Fill(pt);
		if(cent > 1.5 and cent < 2.5) h_r_dpm_pt_TOF_match_three[2]->Fill(pt);
		if(cent > 2.5 and cent < 3.5) h_r_dpm_pt_TOF_match_three[3]->Fill(pt);
		if(cent > 3.5 and cent < 4.5) h_r_dpm_pt_TOF_match_three[4]->Fill(pt);
		if(cent > 4.5 and cent < 5.5) h_r_dpm_pt_TOF_match_three[5]->Fill(pt);
		if(cent > 5.5 and cent < 6.5) h_r_dpm_pt_TOF_match_three[6]->Fill(pt);
		if(cent > 6.5 and cent < 7.5) h_r_dpm_pt_TOF_match_three[7]->Fill(pt);
		if(cent > 7.5 and cent < 8.5) h_r_dpm_pt_TOF_match_three[8]->Fill(pt);

    h_r_dpm_phi_TOF_match_three->Fill(phi);
		h_r_dpm_eta_phi_TOF_match_three->Fill(eta, phi);
    }

    //TOF matching, fill only if all three daughters don't have TOF
    if ( kTof == 0 && p1Tof == 0 && p2Tof == 0 )
    {
    n_cuts->Fill(19);

    // reconstructed Dpm with TOF matching
		if(cent > -0.5 and cent < 0.5) h_r_dpm_pt_TOF_match_zero[0]->Fill(pt);
		if(cent > 0.5 and cent < 1.5) h_r_dpm_pt_TOF_match_zero[1]->Fill(pt);
		if(cent > 1.5 and cent < 2.5) h_r_dpm_pt_TOF_match_zero[2]->Fill(pt);
		if(cent > 2.5 and cent < 3.5) h_r_dpm_pt_TOF_match_zero[3]->Fill(pt);
		if(cent > 3.5 and cent < 4.5) h_r_dpm_pt_TOF_match_zero[4]->Fill(pt);
		if(cent > 4.5 and cent < 5.5) h_r_dpm_pt_TOF_match_zero[5]->Fill(pt);
		if(cent > 5.5 and cent < 6.5) h_r_dpm_pt_TOF_match_zero[6]->Fill(pt);
		if(cent > 6.5 and cent < 7.5) h_r_dpm_pt_TOF_match_zero[7]->Fill(pt);
		if(cent > 7.5 and cent < 8.5) h_r_dpm_pt_TOF_match_zero[8]->Fill(pt);

    h_r_dpm_phi_TOF_match_zero->Fill(phi);
		h_r_dpm_eta_phi_TOF_match_zero->Fill(eta, phi);
    }
    
		
	} //end entries loop


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

  //reconstructed Dpm, "soft" TOF matching
  for(unsigned int j = 0; j<9; j++ )
  {	
    h_r_dpm_pt_TOF_match[j]->Write();
	}
	h_r_dpm_phi_TOF_match->Write();
	h_r_dpm_eta_phi_TOF_match->Write();

  //reconstructed Dpm, one daughter has TOF
  for(unsigned int j = 0; j<9; j++ )
  {	
    h_r_dpm_pt_TOF_match_one[j]->Write();
	}
	h_r_dpm_phi_TOF_match_one->Write();
	h_r_dpm_eta_phi_TOF_match_one->Write();

  //reconstructed Dpm, two daughters have TOF
  for(unsigned int j = 0; j<9; j++ )
  {	
    h_r_dpm_pt_TOF_match_two[j]->Write();
	}
	h_r_dpm_phi_TOF_match_two->Write();
	h_r_dpm_eta_phi_TOF_match_two->Write();

  //reconstructed Dpm, three daughters have TOF
  for(unsigned int j = 0; j<9; j++ )
  {	
    h_r_dpm_pt_TOF_match_three[j]->Write();
	}
	h_r_dpm_phi_TOF_match_three->Write();
	h_r_dpm_eta_phi_TOF_match_three->Write();

  //reconstructed Dpm, no daughters have TOF
  for(unsigned int j = 0; j<9; j++ )
  {	
    h_r_dpm_pt_TOF_match_zero[j]->Write();
	}
	h_r_dpm_phi_TOF_match_zero->Write();
	h_r_dpm_eta_phi_TOF_match_zero->Write();

  n_cuts->Write();


	f->Close();
}
