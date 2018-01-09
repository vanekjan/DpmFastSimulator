#include <TString.h>
//#include <TFile.h>
//#include <TTree.h>
//#include "nt.h"
//#include <iostream>


//void run_nt(TString infile = "Dpm.B47A2BF7B0B04BE7838C9638213FD0E8_99.toyMc.root", TString outfile = "test.root") {
void run_nt(TString infile = "Dpm.toyMc.root", TString outfile = "test_output.root") {
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
