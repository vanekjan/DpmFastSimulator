#include <TString.h>
//#include <TFile.h>
//#include <TTree.h>
//#include "nt.h"
//#include <iostream>
//#include"TROOT.h"


void run_nt_TMVA_Cuts(TString infile = "./myOutput/2018-04-20_02-38_new_HFT_pT_bins_final/merge/output.root", TString outfile = "./myOutput/Histo_output/TMVA/Dpm.out_eff_TMVA_6_pT_bins_iter_2.toyMc.root") { //for output from submit
//void run_nt(TString infile = "Dpm.toyMc.root", TString outfile = "Dpm.out_ana_cuts.toyMc.root") { //for output from local test
  std::cout << "start " << std::endl;
	gROOT->ProcessLine(".L nt_TMVA_Cuts.C+");
	//gROOT->ProcessLine(".L nt.C");
	TFile *f_input = new TFile(infile, "open");
	TTree *tree_nt = (TTree*)f_input->Get("nt");
	nt_TMVA_Cuts n(tree_nt);
	n.Set_out_file_name(outfile);
	//n.Set_con_cent(8);
	n.Set_con_cent_up(9);
	n.Set_con_cent_down(-1);
	n.Set_con_v0z(60000.);
	
	n.Set_con_kRPt(0.3);
	n.Set_con_pRPt(0.3);

  //set topological cuts in tmvaCuts.h
	//_____________________________________________________________________________
	
	n.Loop();
	std::cout << "end " << std::endl;
	f_input->Close();
}
