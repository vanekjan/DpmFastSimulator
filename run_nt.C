#include <TString.h>
//#include <TFile.h>
//#include <TTree.h>
//#include "nt.h"
//#include <iostream>
//#include"TROOT.h"


void run_nt(TString infile = "./myOutput/2018-04-04_05-33_1Sigma_TOF_match/merge/output.root", TString outfile = "./myOutput/Histo_output/Dpm.out_1Sigma_TOF_match_full_stat.toyMc.root") { //for output from submit
//void run_nt(TString infile = "Dpm.toyMc.root", TString outfile = "Dpm.out_ana_cuts.toyMc.root") { //for output from local test
  std::cout << "start " << std::endl;
	gROOT->ProcessLine(".L nt.C+");
	//gROOT->ProcessLine(".L nt.C");
	TFile *f_input = new TFile(infile, "open");
	TTree *tree_nt = (TTree*)f_input->Get("nt");
	nt n(tree_nt);
	n.Set_out_file_name(outfile);
	//n.Set_con_cent(8);
	n.Set_con_cent_up(9);
	n.Set_con_cent_down(-1);
	n.Set_con_v0z(60000.);
	n.Set_con_cosTheta(0.998);
	n.Set_con_dcaDaughters(80);
	n.Set_con_decayLength(30);
	n.Set_con_kDca(80);
	n.Set_con_pDca(100);
	/*
	n.Set_con_kRPt(0.6);
	n.Set_con_pRPt(0.8);
	*/
	n.Set_con_kRPt(0.5);
	n.Set_con_pRPt(0.5);
	n.Set_con_dV0Max(200.);
	n.Loop();
	std::cout << "end " << std::endl;
	f_input->Close();
}
