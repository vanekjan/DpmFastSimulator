/* *********************************************************************
 *  ROOT macro - Toy Monte Carlo Simulation for D0 decay
 *  Includes Momentum Resolution, DCA, hft ration, TPC efficiency ...
 *  Example for D0 --> Kpi
 *
 *  Authors:
 *            Guannan Xie (guannanxie@lbl.gov)
 *            **Mustafa Mustafa (mmustafa@lbl.gov)
 *            Hao Qiu (hqiu@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
 */

#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TClonesArray.h"
#include "TPythia6.h"
#include "TPythia6Decayer.h"
#include "TRandom3.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMath.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "TStopwatch.h"
// #include "TSystem.h"
// #include "TMemStat.h"

using namespace std;

void setDecayChannels(int const mdme);
void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters);
void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& p1Mom, TLorentzVector const& p2Mom, TVector3 v00);
void getKinematics(TLorentzVector& b, double const mass);
TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution);
TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos);
TVector3 smearPosData(int iParticleIndex, double vz, int cent, TLorentzVector const& rMom, TVector3 const& pos);
float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0);
TVector3 getVertex(int centrality);
bool matchHft(int iParticleIndex, double vz, int cent, TLorentzVector const& mom);
bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom);
bool matchTOF(int const iParticleIndex, TLorentzVector const& mom);
void bookObjects();
void write();
int getPtIndexDca(double);
int getEtaIndexDca(double);
int getVzIndexDca(double);
int getPhiIndexDca(double);

int getPtIndexHftRatio(double);
int getEtaIndexHftRatio(double);
int getVzIndexHftRatio(double);
int getPhiIndexHftRatio(double);

TPythia6Decayer* pydecay;
TNtuple* nt;
TFile* result;

TF1* fKaonMomResolution = NULL;
TF1* fPionMomResolution = NULL;
TF1* fWeightFunction = NULL;
const Int_t nParticles = 2;
const Int_t nCentHftRatio = 9;

// HFT ratio binning
const Int_t nEtasHftRatio = 10;
const Int_t nVzsHftRatio = 6;
const Int_t nPtBinsHftRatio = 36;
const Double_t EtaEdgeHftRatio[nEtasHftRatio + 1] =
{
	-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0
};
const Double_t VzEdgeHftRatio[nVzsHftRatio + 1] = { -6., -4., -2.,  0., 2.,  4.,  6.};
/*{
	-6.0e4, -4.0e4, -2.0e4, 0.0, 2.0e4, 4.0e4, 6.0e4
};
*/
const Double_t ptEdgeHftRatio[nPtBinsHftRatio + 1] =
   {
      0.3, 0.4, 0.5, 0.6 , 0.7 , 0.8 , 0.9 ,
      1. , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 ,
      2. , 2.2 , 2.4 , 2.6 , 2.8 , 3.0 ,
      3.4 , 3.8 , 4.2 , 4.6 , 5.0 ,  5.5 ,
      6. , 6.5 , 7.0 , 8.0 , 9.0 , 10. , 11,  12.0
   };


/*{
	0.2, 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 ,
	1. , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 ,
	2. , 2.2 , 2.4 , 2.6 , 2.8 , 3.0 , 3.2 , 3.4 , 3.6 , 3.8 ,
	4. , 4.2 , 4.4 , 4.6 , 4.8 , 5.0 , 5.2 , 5.4 , 5.6 , 5.8 ,
	6. , 6.5 , 7.0 , 7.5 , 8.0 , 8.5 , 9.0 , 9.5 , 10. , 10.5,
	11.0 , 11.5 , 12.0
};
*/
const Int_t nPhisHftRatio = 11;
const Double_t PhiEdgeHftRatio[nPhisHftRatio + 1] =
{
	-3.14159 , -2.80359 , -2.17527 , -1.54696 , -0.918637 , -0.290319 , 0.338 , 0.966319 , 1.59464 , 2.22296 , 2.85127 , 3.14159 //Sector by Sector  // sector number 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3
};

const Int_t nPhisDca = 11;
const Double_t PhiEdgeDca[nPhisDca + 1] =
{
	-3.14159 , -2.80359 , -2.17527 , -1.54696 , -0.918637 , -0.290319 , 0.338 , 0.966319 , 1.59464 , 2.22296 , 2.85127 , 3.14159 //Sector by Sector  // sector number 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3
};

// DCA binning
int const nVzsDca = 4;
float const VzEdgeDca[nVzsDca + 1] = {   -6.e4, -3.e4, 0, 3.e4, 6.e4};

int const nEtasDca = 5;
float const EtaEdgeDca[nEtasDca + 1] = { -1.0, -0.6, -0.2, 0.2, 0.6, 1.0};

const Int_t nPtBinsDca = 19;
const Double_t ptEdgeDca[nPtBinsDca + 1] =
   {
      0.3, 0.4, 0.5,
      0.6,  0.7 , 0.8 , 0.9 ,
      1. , 1.25 , 1.5 , 1.75 , 2.  , 2.25 , 2.5 , 2.75 , 3.0 , 3.5,
      // 3. , 3.5 , 4.  , 4.5 , 5. , 6. , 8.0 , 10. , 12.0
      4.  , 6. , 12.0
   };

/*{
	0.2, 0.3, 0.4, 0.5, 0.6,  0.7 , 0.8 , 0.9 ,
	1. , 1.25 , 1.5 , 1.75 , 2.  , 2.25 , 2.5 , 2.75 , 3.0 , 3.5,
	// 3. , 3.5 , 4.  , 4.5 , 5. , 6. , 8.0 , 10. , 12.0
	4.  , 4.5 , 5. , 12.0
};
*/

TH1D* h1Vz[nCentHftRatio];

TH1D* hHftRatio1[nParticles][nEtasHftRatio][nVzsHftRatio][nPhisHftRatio][nCentHftRatio];
int const nCentDca = 9;
TH2D* h2Dca[nParticles][nEtasDca][nVzsDca][nCentDca][nPtBinsDca];

TH1D* hTpcPiPlus[nCentHftRatio];
TH1D* hTpcPiMinus[nCentHftRatio];
TH1D* hTpcKPlus[nCentHftRatio];
TH1D* hTpcKMinus[nCentHftRatio];

string outFileName = "Dpm.toyMc.root";
std::pair<int, int> const decayChannels(673, 735);
std::pair<float, float> const momentumRange(0, 15);

float const gVzCut = 6.0e4;
float const acceptanceRapidity = 1.0;
float const sigmaPos0 = 15.2;
float const pxlLayer1Thickness = 0.00486;
float const sigmaVertexCent[nCentHftRatio] = {31., 18.1, 12.8, 9.3, 7.2, 5.9, 5., 4.6, 4.};

TF1* fPionTofEff= NULL;
TF1* fKaonTofEff= NULL;
TH1D* h_pi_tof_eff;
TH1D* h_k_tof_eff;

//============== main  program ==================
void toyMcDpm(int npart = 10000)
{
	// TMemStat mem;
	// mem.Enable();
	TStopwatch*   stopWatch = new TStopwatch();
	stopWatch->Start();
	gRandom->SetSeed();
	bookObjects();

	pydecay = TPythia6Decayer::Instance();
	pydecay->Init();

	setDecayChannels(719); // D+ --> K- pi+ pi+
	// setDecayChannels(696); // D+ --> K*bar0 + pi+
	// setDecayChannels(617); // D+ --> K*bar0 + pi+ // --> K*bar0 -->K- pi+


	TLorentzVector* b_d = new TLorentzVector;
	TClonesArray ptl("TParticle", 10);
	for (int ipart = 0; ipart < npart; ipart++)
	{
		if (!(ipart % 1000))
    {
			cout << "____________ ipart = " << ipart / static_cast<float>(npart) << " ________________" << endl;
    }

		getKinematics(*b_d, M_D_PLUS);

		decayAndFill(411, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
		decayAndFill(-411, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
		if (ipart % 1000 == 1) nt->AutoSave("SaveSelf");
	}

	write();
	// mem.Show();
	stopWatch->Stop();
	stopWatch->Print();
}

void setDecayChannels(int const mdme)
{
	for (int idc = decayChannels.first; idc < decayChannels.second + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0);
	TPythia6::Instance()->SetMDME(mdme, 1, 1);
}

void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters)
{
	pydecay->Decay(kf, b);
	pydecay->ImportParticles(&daughters);

	TLorentzVector kMom;
	TLorentzVector p1Mom;
	TLorentzVector p2Mom;
	TVector3 v00;

	int nTrk = daughters.GetEntriesFast();
	for (int iTrk = 0; iTrk < nTrk; ++iTrk)
	{
		TParticle* ptl0 = (TParticle*)daughters.At(iTrk);

		switch (abs(ptl0->GetPdgCode()))
		{
			case 321:
				ptl0->Momentum(kMom);
				// v00.SetXYZ(0,0,0);
				v00.SetXYZ(ptl0->Vx() * 1000., ptl0->Vy() * 1000., ptl0->Vz() * 1000.); // converted to μm
				//if (ptl0->Momentum(kMom).mag() < 0.15) return;
				break;
			case 211:
				//if (ptl0->Momentum(kMom).mag() < 0.15) return;
				if (!p1Mom.P()) ptl0->Momentum(p1Mom);
				else ptl0->Momentum(p2Mom);
				break;
			default:
				break;
		}
	}
	daughters.Clear();

	fill(kf, b, weight, kMom, p1Mom, p2Mom, v00);
}

void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& p1Mom, TLorentzVector const& p2Mom, TVector3 v00)
{
	int const centrality = floor(nCentHftRatio * gRandom->Rndm());

	TVector3 const vertex = getVertex(centrality);
	// smear primary vertex
	// float const sigmaVertex = sigmaVertexCent[cent];
	// TVector3 const vertex(gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex));

	v00 += vertex;

	// smear momentum
	TLorentzVector const kRMom  = smearMom(kMom,  fKaonMomResolution);
	TLorentzVector const p1RMom = smearMom(p1Mom, fPionMomResolution);
	TLorentzVector const p2RMom = smearMom(p2Mom, fPionMomResolution);

	// smear position
	TVector3 const kRPos  = smearPosData(1, vertex.z(), centrality, kRMom,  v00);
	TVector3 const p1RPos = smearPosData(0, vertex.z(), centrality, p1RMom, v00);
	TVector3 const p2RPos = smearPosData(0, vertex.z(), centrality, p2RMom, v00);
	// TVector3 const kRPos = smearPos(kMom, kRMom, v00);
	// TVector3 const p1RPos = smearPos(p1Mom, p1RMom, v00);
	// TVector3 const p2RPos = smearPos(p2Mom, p1RMom, v00);

	// reconstruct
	TLorentzVector const rMom = kRMom + p1RMom + p2RMom;
	float const kDca = dca(kMom.Vect(), v00, vertex);
	float const p1Dca = dca(p1Mom.Vect(), v00, vertex);
	float const p2Dca = dca(p2Mom.Vect(), v00, vertex);
	float const kRDca = dca(kRMom.Vect(), kRPos, vertex);
	float const kRSDca = dcaSigned(kRMom.Vect(), kRPos, vertex);
	float const kRDcaXY = dcaXY(kRMom.Vect(), kRPos, vertex);
	float const kRDcaZ = dcaZ(kRMom.Vect(), kRPos, vertex);
	float const p1RDca = dca(p1RMom.Vect(), p1RPos, vertex);
	float const p1RSDca = dcaSigned(p1RMom.Vect(), p1RPos, vertex);
	float const p1RDcaXY = dcaXY(p1RMom.Vect(), p1RPos, vertex);
	float const p1RDcaZ = dcaZ(p1RMom.Vect(), p1RPos, vertex);
	float const p2RDca = dca(p2RMom.Vect(), p2RPos, vertex);
	float const p2RSDca = dcaSigned(p2RMom.Vect(), p2RPos, vertex);
	float const p2RDcaXY = dcaXY(p2RMom.Vect(), p2RPos, vertex);
	float const p2RDcaZ = dcaZ(p2RMom.Vect(), p2RPos, vertex);

	TVector3 v12, v23, v31, v0;
	float const dca12 = dca1To2(kRMom.Vect(), kRPos, p1RMom.Vect(), p1RPos, v12);
	float const dca23 = dca1To2(p1RMom.Vect(), p1RPos, p2RMom.Vect(), p2RPos, v23);
	float const dca31 = dca1To2(p2RMom.Vect(), p2RPos, kRMom.Vect(), kRPos, v31);
	float dcaDaughters = dca12 > dca31 ? dca12 : dca31;
	dcaDaughters = dcaDaughters > dca23 ? dcaDaughters : dca23;
	v0 = (v12 + v23 + v31) * 0.333333333;
	float const decayLength = (v0 - vertex).Mag();
	float const dcaDpmToPv = dca(rMom.Vect(), v0, vertex);
	float const cosTheta = (v0 - vertex).Unit().Dot(rMom.Vect().Unit());
	float const angle12 = kRMom.Vect().Angle(p1RMom.Vect());

	// Distance between v12 and v23
	//float const v12 = (p1AtDcaToP2 + p2AtDcaToP1 - p2AtDcaToP3 - p3AtDcaToP2).mag()/2.0;
	float const v12_23 = (v12 - v23).Mag()/2.0;
	// Distance between v23 and v31
	//float const v23 = (p2AtDcaToP3 + p3AtDcaToP2 - p3AtDcaToP1 - p1AtDcaToP3).mag()/2.0;
	float const v23_31 = (v23 - v31).Mag()/2.0;
	// Distance between v31 and v12
	//float const v31 = (p3AtDcaToP1 + p1AtDcaToP3 - p1AtDcaToP2 - p2AtDcaToP1).mag()/2.0;
	float const v31_12 = (v31 - v12).Mag()/2.0;
	//maximum dist between reo v0's to be averaging
	float const max12 =  v12_23 > v23_31 ? v12_23 : v23_31 ;
	float const mdV0Max = max12 > v31_12 ? max12 : v31_12;

	TLorentzVector kRMomRest = kRMom;
	TVector3 beta;
	beta.SetMagThetaPhi(rMom.Beta(), rMom.Theta(), rMom.Phi());
	kRMomRest.Boost(-beta);
	float const cosThetaStar = rMom.Vect().Unit().Dot(kRMomRest.Vect().Unit());

	int const charge = kf > 0 ? 1 : -1;

	/*

	cout << "kaon " << endl;
	//cout << hk_tof_eff->GetEntries() << endl;
	//int const bin = hk_tof_eff->FindBin(kRMom.Perp());
	int const bin = hk_tof_eff->FindBin(5.);
	cout << "bin " << bin << endl;
	if (gRandom->Rndm() < hk_tof_eff->GetBinContent(bin)) {
		cout << "true " << endl;
	} else {
		cout << "false" << endl;
	}
	*/

	// save
	float arr[110];
	int iArr = 0;
	arr[iArr++] = centrality;
	arr[iArr++] = vertex.X();
	arr[iArr++] = vertex.Y();
	arr[iArr++] = vertex.Z();
	arr[iArr++] = getVzIndexDca(vertex.Z());

	arr[iArr++] = kf;
	arr[iArr++] = weight;
	arr[iArr++] = b->M();
	arr[iArr++] = b->Perp();
	arr[iArr++] = b->PseudoRapidity();
	arr[iArr++] = b->Rapidity();
	arr[iArr++] = b->Phi();
	arr[iArr++] = v00.X();
	arr[iArr++] = v00.Y();
	arr[iArr++] = v00.Z();

	arr[iArr++] = rMom.M();
	arr[iArr++] = rMom.Perp();
	arr[iArr++] = rMom.PseudoRapidity();
	arr[iArr++] = rMom.Rapidity();
	arr[iArr++] = rMom.Phi();
	arr[iArr++] = v0.X();
	arr[iArr++] = v0.Y();
	arr[iArr++] = v0.Z();

	arr[iArr++] = dcaDaughters;
	arr[iArr++] = decayLength;
	arr[iArr++] = dcaDpmToPv;
	arr[iArr++] = cosTheta;
	arr[iArr++] = angle12;
	arr[iArr++] = cosThetaStar;

	arr[iArr++] = kMom.M();
	arr[iArr++] = kMom.Perp();
	arr[iArr++] = kMom.PseudoRapidity();
	arr[iArr++] = kMom.Rapidity();
	arr[iArr++] = kMom.Phi();
	arr[iArr++] = kDca;

	arr[iArr++] = kRMom.M();
	arr[iArr++] = kRMom.Perp();
	arr[iArr++] = kRMom.PseudoRapidity();
	arr[iArr++] = kRMom.Rapidity();
	arr[iArr++] = kRMom.Phi();
	arr[iArr++] = kRPos.X();
	arr[iArr++] = kRPos.Y();
	arr[iArr++] = kRPos.Z();
	arr[iArr++] = kRDca;
	arr[iArr++] = kRSDca;
	arr[iArr++] = kRDcaXY;
	arr[iArr++] = kRDcaZ;
	arr[iArr++] = getEtaIndexDca(kRMom.PseudoRapidity());
	arr[iArr++] = getPtIndexDca(kRMom.Perp());
	arr[iArr++] = tpcReconstructed(1, -1 * charge, centrality, kRMom);

	arr[iArr++] = p1Mom.M();
	arr[iArr++] = p1Mom.Perp();
	arr[iArr++] = p1Mom.PseudoRapidity();
	arr[iArr++] = p1Mom.Rapidity();
	arr[iArr++] = p1Mom.Phi();
	arr[iArr++] = p1Dca;

	arr[iArr++] = p1RMom.M();
	arr[iArr++] = p1RMom.Perp();
	arr[iArr++] = p1RMom.PseudoRapidity();
	arr[iArr++] = p1RMom.Rapidity();
	arr[iArr++] = p1RMom.Phi();
	arr[iArr++] = p1RPos.X();
	arr[iArr++] = p1RPos.Y();
	arr[iArr++] = p1RPos.Z();
	arr[iArr++] = p1RDca;
	arr[iArr++] = p1RSDca;
	arr[iArr++] = p1RDcaXY;
	arr[iArr++] = p1RDcaZ;
	arr[iArr++] = getEtaIndexDca(p1RMom.PseudoRapidity());
	arr[iArr++] = getPtIndexDca(p1RMom.Perp());
	arr[iArr++] = tpcReconstructed(0, charge, centrality, p1RMom);

	arr[iArr++] = p2Mom.M();
	arr[iArr++] = p2Mom.Perp();
	arr[iArr++] = p2Mom.PseudoRapidity();
	arr[iArr++] = p2Mom.Rapidity();
	arr[iArr++] = p2Mom.Phi();
	arr[iArr++] = p2Dca;

	arr[iArr++] = p2RMom.M();
	arr[iArr++] = p2RMom.Perp();
	arr[iArr++] = p2RMom.PseudoRapidity();
	arr[iArr++] = p2RMom.Rapidity();
	arr[iArr++] = p2RMom.Phi();
	arr[iArr++] = p2RPos.X();
	arr[iArr++] = p2RPos.Y();
	arr[iArr++] = p2RPos.Z();
	arr[iArr++] = p2RDca;
	arr[iArr++] = p2RSDca;
	arr[iArr++] = p2RDcaXY;
	arr[iArr++] = p2RDcaZ;
	arr[iArr++] = getEtaIndexDca(p2RMom.PseudoRapidity());
	arr[iArr++] = getPtIndexDca(p2RMom.Perp());
	arr[iArr++] = tpcReconstructed(0, charge, centrality, p2RMom);

	arr[iArr++] = matchHft(1, vertex.z(), centrality, kRMom);
	arr[iArr++] = matchHft(0, vertex.z(), centrality, p1RMom);
	arr[iArr++] = matchHft(0, vertex.z(), centrality, p2RMom);
	arr[iArr++] = matchTOF(1, kRMom);
	arr[iArr++] = matchTOF(0, p1RMom);
	arr[iArr++] = matchTOF(0, p2RMom);
	arr[iArr++] = pow((kRMom + p1RMom).M(),2);
	arr[iArr++] = pow((kRMom + p2RMom).M(),2);
	arr[iArr++] = dca12; //DCA kp1
	arr[iArr++] = dca23; //DCA p1p2 
	arr[iArr++] = dca31; //DCA p2k
	arr[iArr++] = mdV0Max;
	//arr[iArr++] = ;
	//cout << "iArr: " << iArr << endl;

	nt->Fill(arr);
}

void getKinematics(TLorentzVector& b, double const mass)
{
	float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
	float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
	float const phi = TMath::TwoPi() * gRandom->Rndm();

	float const mT = sqrt(mass * mass + pt * pt);
	float const pz = mT * sinh(y);
	float const E = mT * cosh(y);

	b.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
	TVector3 posDiff = pos - vertex;
	return fabs(p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff));
}

float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
	TVector3 posDiff = pos - vertex;
	float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;

	return sign * p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff);
}

float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
	TVector3 newPos(pos);
	newPos.SetZ(0);

	TVector3 newP(p);
	newP.SetZ(0);

	TVector3 newVertex(vertex);
	newVertex.SetZ(0);

	TVector3 posDiff = newPos - newVertex;
	float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
	return sign * newP.Cross(posDiff.Cross(newP)).Unit().Dot(posDiff);
}

float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
	TVector3 posDiff = pos - vertex;
	if (sin(p.Theta()) == 0) return 0;
	else return (-posDiff.x() * cos(p.Phi()) - posDiff.y() * sin(p.Phi())) * cos(p.Theta()) / sin(p.Theta()) + posDiff.z();
}

float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0)
{
	TVector3 posDiff = pos2 - pos1;
	TVector3 pu1 = p1.Unit();
	TVector3 pu2 = p2.Unit();
	double pu1Pu2 = pu1.Dot(pu2);
	double g = posDiff.Dot(pu1);
	double k = posDiff.Dot(pu2);
	double s2 = (k - pu1Pu2 * g) / (pu1Pu2 * pu1Pu2 - 1.);
	double s1 = g + s2 * pu1Pu2;
	TVector3 posDca1 = pos1 + pu1 * s1;
	TVector3 posDca2 = pos2 + pu2 * s2;
	v0 = 0.5 * (posDca1 + posDca2);
	return (posDca1 - posDca2).Mag();
}

TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution)
{
	float const pt = b.Perp();
	float const sPt = gRandom->Gaus(pt, pt * fMomResolution->Eval(pt));

	TLorentzVector sMom;
	sMom.SetXYZM(sPt * cos(b.Phi()), sPt * sin(b.Phi()), sPt * sinh(b.PseudoRapidity()), b.M());
	return sMom;
}

TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos)
{
	float thetaMCS = 13.6 / mom.Beta() / rMom.P() / 1000 * sqrt(pxlLayer1Thickness / fabs(sin(mom.Theta())));
	float sigmaMCS = thetaMCS * 28000 / fabs(sin(mom.Theta()));
	float sigmaPos = sqrt(pow(sigmaMCS, 2) + pow(sigmaPos0, 2));

	return TVector3(gRandom->Gaus(pos.X(), sigmaPos), gRandom->Gaus(pos.Y(), sigmaPos), gRandom->Gaus(pos.Z(), sigmaPos));
}

int getPtIndexDca(double pT)
{
	for (int i = 0; i < nPtBinsDca; i++)
	{
		if ((pT >= ptEdgeDca[i]) && (pT < ptEdgeDca[i + 1]))
			return i;
	}
	return nPtBinsDca - 1 ;
}

int getEtaIndexDca(double Eta)
{
	for (int i = 0; i < nEtasDca; i++)
	{
		if ((Eta >= EtaEdgeDca[i]) && (Eta < EtaEdgeDca[i + 1]))
			return i;
	}
	return nEtasDca - 1 ;
}

int getVzIndexDca(double Vz)
{
	for (int i = 0; i < nVzsDca; i++)
	{
		if ((Vz >= VzEdgeDca[i]) && (Vz < VzEdgeDca[i + 1]))
			return i;
	}
	return nVzsDca - 1 ;
}

int getPhiIndexDca(double Phi)
{
	for (int i = 0; i < nPhisDca; i++)
	{
		if ((Phi >= PhiEdgeDca[i]) && (Phi < PhiEdgeDca[i + 1]))
			return i;
	}
	return nPhisDca - 1 ;
}

int getPtIndexHftRatio(double pT)
{
	for (int i = 0; i < nPtBinsHftRatio; i++)
	{
		if ((pT >= ptEdgeHftRatio[i]) && (pT < ptEdgeHftRatio[i + 1]))
			return i;
	}
	return nPtBinsHftRatio - 1 ;
}

int getEtaIndexHftRatio(double Eta)
{
	for (int i = 0; i < nEtasHftRatio; i++)
	{
		if ((Eta >= EtaEdgeHftRatio[i]) && (Eta < EtaEdgeHftRatio[i + 1]))
			return i;
	}
	return nEtasHftRatio - 1 ;
}

int getVzIndexHftRatio(double Vz)
{
	for (int i = 0; i < nVzsHftRatio; i++)
	{
		if ((Vz >= VzEdgeHftRatio[i]) && (Vz < VzEdgeHftRatio[i + 1]))
			return i;
	}
	return nVzsHftRatio - 1 ;
}

int getPhiIndexHftRatio(double Phi)
{
	for (int i = 0; i < nPhisHftRatio; i++)
	{
		if ((Phi >= PhiEdgeHftRatio[i]) && (Phi < PhiEdgeHftRatio[i + 1]))
			return i;
	}
	return nPhisHftRatio - 1 ;
}

TVector3 smearPosData(int const iParticleIndex, double const vz, int cent, TLorentzVector const& rMom, TVector3 const& pos)
{
	int const iEtaIndex = getEtaIndexDca(rMom.PseudoRapidity());
	int const iVzIndex = getVzIndexDca(vz);
	// int const iPhiIndex = getPhiIndexDca(rMom.Phi());
	int const iPtIndex = getPtIndexDca(rMom.Perp());

	double sigmaPosZ = 0;
	double sigmaPosXY = 0;

	// if (cent == 8) cent = 7;
	// All the centrality position smear was based on 0-10% centrality input
	// changed to 0-80%

	h2Dca[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom2(sigmaPosXY, sigmaPosZ);
  if(sigmaPosXY==0 || sigmaPosZ==0)
  {
    cout<<iParticleIndex<<" "<<iEtaIndex<<" "<<iVzIndex<<" "<<cent<<" "<<iPtIndex<<endl;
  }
	// h2Dca[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][iPtIndex]->GetRandom2(sigmaPosXY, sigmaPosZ);
	sigmaPosZ *= 1.e4;
	sigmaPosXY *= 1.e4;
	/*if (h1DcaZ1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->ComputeIntegral())
	  {
	  do sigmaPosZ = h1DcaZ1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom() * 1e4;
	  while (fabs(sigmaPosZ) > 1.e3);
	  }

	  if (h1DcaXY1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->ComputeIntegral())
	  {
	  do sigmaPosXY = h1DcaXY1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom() * 1e4;
	  while (fabs(sigmaPosXY) > 1.e3);
	  }
	  */

	TVector3 newPos(pos);
	newPos.SetZ(0);
	TVector3 momPerp(-rMom.Vect().Y(), rMom.Vect().X(), 0.0);
	newPos -= momPerp.Unit() * sigmaPosXY;

	return TVector3(newPos.X(), newPos.Y(), pos.Z() + sigmaPosZ);
}

TVector3 getVertex(int const centrality)
{
	double rdmVz;

	if (h1Vz[centrality]->GetEntries() == 0) rdmVz = 0.;
	else
	{
		do rdmVz = h1Vz[centrality]->GetRandom() * 1e4;
		while (fabs(rdmVz) > gVzCut);
	}

	return TVector3(0., 0., rdmVz);
}

bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom)
{
	TH1D* h = NULL;

	if (iParticleIndex == 0)
	{
		if (charge > 0) h = hTpcPiPlus[cent];
		else h = hTpcPiMinus[cent];
	}
	else
	{
		if (charge > 0) h = hTpcKPlus[cent];
		else h = hTpcKMinus[cent];
	}

	int const bin = h->FindBin(mom.Perp());

	return gRandom->Rndm() < h->GetBinContent(bin);
}

bool matchHft(int const iParticleIndex, double const vz, int const cent, TLorentzVector const& mom)
{
	int const iEtaIndex = getEtaIndexHftRatio(mom.PseudoRapidity());
	int const iVzIndex = getVzIndexHftRatio(vz);
	int const iPhiIndex = getPhiIndexHftRatio(mom.Phi());

	int const bin = hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][cent]->FindBin(mom.Perp());
	return gRandom->Rndm() < hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][cent]->GetBinContent(bin);
}

bool matchTOF(int const iParticleIndex, TLorentzVector const& mom)
{

	//cout << "tof eff ... " << endl;
	if (iParticleIndex == 0) { // pion
		//return gRandom->Rndm() < fPionTofEff->Eval(mom.Perp()); //from function
    
    return gRandom->Rndm() < h_pi_tof_eff->GetBinContent(h_pi_tof_eff->FindBin(mom.Perp())); //from histogram
    
		/*
		if (gRandom->Rndm() < fKaonTofEff->Eval(5.)) {
		cout << "pion " << endl;
		if (mom.Perp() > 10.) {
			int const bin = h_pi_tof_eff->FindBin(9.9);
			cout << "bin " << bin << endl;
			return gRandom->Rndm() < h_pi_tof_eff->GetBinContent(bin);
		} else {
			int const bin = h_pi_tof_eff->FindBin(mom.Perp());
			cout << "bin " << bin << endl;
			return gRandom->Rndm() < h_pi_tof_eff->GetBinContent(bin);
		}
		*/
	} 
  else if (iParticleIndex == 1) { // kaon
		//return gRandom->Rndm() < fKaonTofEff->Eval(mom.Perp()); //from function

    return gRandom->Rndm() < h_k_tof_eff->GetBinContent(h_k_tof_eff->FindBin(mom.Perp())); //from histogram
		/*
		cout << "kaon " << endl;
		if (mom.Perp() > 10.) {
			//Double_t p = 9.9;
			//int const bin = hk_tof_eff->FindBin(p);
			cout << "<" << endl;
			int const bin = hk_tof_eff->GetNbinsX() - 1;
			cout << "bin " << bin << endl;
			return gRandom->Rndm() < hk_tof_eff->GetBinContent(bin);
		} else {
			cout << ">" << endl;
			cout << hk_tof_eff->GetEntries() << endl;
			int const bin = hk_tof_eff->FindBin(mom.Perp());
			cout << "bin " << bin << endl;
			return gRandom->Rndm() < hk_tof_eff->GetBinContent(bin);
		}
	*/
	} else {
		return false;
	}
}

//___________
void bookObjects()
{
	cout << "Loading input momentum resolution ..." << endl;
	TFile f("./input/momentum_resolution.root");
	fPionMomResolution = (TF1*)f.Get("fPion")->Clone("fPion");
	fKaonMomResolution = (TF1*)f.Get("fKaon")->Clone("fKaon");
	f.Close();

	cout << "Loading TOF eff ..." << endl;
	TFile f_tof("./input/tof_eff_Dmp_run16_HFT_1sig_20_hist.root"); //tof_eff_Dmp_run16_HFT_1sig_20_hist.root contains histograms only!!!
	//fPionTofEff = (TF1*)f_tof.Get("fPion")->Clone("fPion"); //commented for tof_eff_Dmp_run16_HFT_1sig_20_hist.root input
	//fKaonTofEff = (TF1*)f_tof.Get("fKaon")->Clone("fKaon");
  h_pi_tof_eff = (TH1D*)f_tof.Get("h_pi_tof_eff");
  h_pi_tof_eff->SetDirectory(0);
  h_k_tof_eff = (TH1D*)f_tof.Get("h_k_tof_eff");
  h_k_tof_eff->SetDirectory(0);
	/*
	float pp = 5.;
	cout << "at 5 " << fKaonTofEff->Eval(pp) << endl;
	*/
	f_tof.Close();
	/*
	cout << "at 5 " << fKaonTofEff->Eval(5.) << endl;
	if (gRandom->Rndm() < fKaonTofEff->Eval(5.)) {
		cout << "true " << endl;
	} else {
		cout << "false" << endl;
	}
	*/


  //all files moved into "input" folder
	cout << "Loading input spectra ..." << endl;
	TFile fPP("./input/pp200_spectra.root");
	fWeightFunction = (TF1*)fPP.Get("run12/f1Levy")->Clone("f1Levy");
	fPP.Close();

	TFile fVertex("./input/Vz_Cent_Run16.root");

	for (int ii = 0; ii < nCentHftRatio; ++ii)
	{
		h1Vz[ii]      = (TH1D*)(fVertex.Get(Form("mh1Vz_%i", ii)));
		h1Vz[ii]->SetDirectory(0);
	}

	fVertex.Close();

	cout << "Loading input HFT ratios and DCA ..." << endl;
	
	TFile *fHftRatio1 = new TFile("./input/HFT_Ratio_VsPt_Centrality_Eta_Phi_Vz_Zdcx_new_cuts_new.root", "read"); //strict nSigma cuts (1Sigma)
	TFile fDca1("./input/2DProjection_simCent_NoBinWidth_3D_Dca_VsPt_Centrality_Eta_Phi_Vz_Zdcx_new_cuts_new.root");
	
	//TFile* fHftRatio1 = new TFile("./input/HFT_Ratio_VsPt_Centrality_Eta_Phi_Vz_Zdcx_ana_cuts_new.root", "read"); //nSigma cuts as in data production
        //TFile fDca1("./input/2DProjection_simCent_NoBinWidth_3D_Dca_VsPt_Centrality_Eta_Phi_Vz_Zdcx_ana_cuts_new.root");

//cout<<"test"<<endl;
	
	for (int iParticle = 0; iParticle < nParticles; ++iParticle)
	{
		for (int iCent = 0; iCent < nCentHftRatio; ++iCent)
		{
			// HFT ratio
			for (int iEta = 0; iEta < nEtasHftRatio; ++iEta)
			{
				for (int iVz = 0; iVz < nVzsHftRatio; ++iVz)
				{
					for (int iPhi = 0; iPhi < nPhisHftRatio; ++iPhi)
					{
            //cout<<iParticle<<" "<<iEta<<" "<<iVz<<" "<<iPhi<<" "<<iCent<<endl;
						hHftRatio1[iParticle][iEta][iVz][iPhi][iCent] = (TH1D*)fHftRatio1->Get(Form("mh1HFT1PtCentPartEtaVzPhi_%i_%i_%i_%i_%i", iParticle, iEta, iVz, iPhi, iCent));
						hHftRatio1[iParticle][iEta][iVz][iPhi][iCent]->SetDirectory(0);
					}
				}
			}
		}
		cout << "Finished loading HFT Ratio: " <<  endl;

		for (int iCent = 0; iCent < nCentDca; ++iCent)
		{
			// DCA
			for (int iEta = 0; iEta < nEtasDca; ++iEta)
			{
				for (int iVz = 0; iVz < nVzsDca; ++iVz)
				{
					// for (int iPhi = 0; iPhi < nPhisDca; ++iPhi)
					for (int iPt = 0; iPt < nPtBinsDca; ++iPt)
					{
						h2Dca[iParticle][iEta][iVz][iCent][iPt] = (TH2D*)((fDca1.Get(Form("mh2DcaPtCentPartEtaVzPhi_%i_%i_%i_%i_%i", iParticle, iEta, iVz, iCent, iPt))));
						h2Dca[iParticle][iEta][iVz][iCent][iPt]->SetDirectory(0);
					}
				}
			}
		}
		// cout << "Finished loading centrality: " << iCent << endl;
	}
	cout << "Finished loading Dca: " <<  endl;

	fHftRatio1->Close();
	fDca1.Close();

	cout << " Loading TPC tracking efficiencies " << endl;

	TFile fTpcPiPlus("./input/Eff_PionPlus_embedding_v2.root");
	TFile fTpcPiMinus("./input/Eff_PionMinus_embedding_v2.root");
	TFile fTpcKPlus("./input/Eff_KaonPlus_embedding_v2.root");
	TFile fTpcKMinus("./input/Eff_KaonMinus_embedding_v2.root");

	for (int iCent = 0; iCent < nCentHftRatio; ++iCent)
	{
		hTpcPiPlus[iCent] = (TH1D*)fTpcPiPlus.Get(Form("h1Ratiocent_%i", iCent));
		hTpcPiPlus[iCent]->SetDirectory(0);
		hTpcPiMinus[iCent] = (TH1D*)fTpcPiMinus.Get(Form("h1Ratiocent_%i", iCent));
		hTpcPiMinus[iCent] ->SetDirectory(0);
		hTpcKPlus[iCent] = (TH1D*)fTpcKPlus.Get(Form("h1Ratiocent_%i", iCent));
		hTpcKPlus[iCent]->SetDirectory(0);
		hTpcKMinus[iCent] = (TH1D*)fTpcKMinus.Get(Form("h1Ratiocent_%i", iCent));
		hTpcKMinus[iCent]->SetDirectory(0);
	}

	fTpcPiPlus.Close();
	fTpcPiMinus.Close();
	fTpcKPlus.Close();
	fTpcKMinus.Close();

	cout << "Done with loading all files ..." << endl;

	result = new TFile(outFileName.c_str(), "recreate");
	result->SetCompressionLevel(1);
	result->cd();

	int BufSize = (int)pow(2., 16.);
	// int Split = 1;
	nt = new TNtuple("nt", "", "cent:vx:vy:vz:vzIdx:"
			"pid:w:m:pt:eta:y:phi:v0x:v0y:v0z:" // MC Dpm
			"rM:rPt:rEta:rY:rPhi:rV0x:rV0y:rV0z:" // Rc Dpm
			"dcaDaughters:decayLength:dcaDpmToPv:cosTheta:angle12:cosThetaStar:" // Rc pair
			"kM:kPt:kEta:kY:kPhi:kDca:" // MC Kaon
			"kRM:kRPt:kREta:kRY:kRPhi:kRVx:kRVy:kRVz:kRDca:kRSDca:kRDcaXY:kRDcaZ:kEtaIdx:kPtIdx:kTpc:" // Rc Kaon
			"p1M:p1Pt:p1Eta:p1Y:p1Phi:p1Dca:" // MC Pion1
			"p1RM:p1RPt:p1REta:p1RY:p1RPhi:p1RVx:p1RVy:p1RVz:p1RDca:p1RSDca:p1RDcaXY:p1RDcaZ:p1EtaIdx:p1PtIdx:p1Tpc:" // Rc Pion1
			"p2M:p2Pt:p2Eta:p2Y:p2Phi:p2Dca:" // MC Pion2
			"p2RM:p2RPt:p2REta:p2RY:p2RPhi:p2RVx:p2RVy:p2RVz:p2RDca:p2RSDca:p2RDcaXY:p2RDcaZ:p2EtaIdx:p2PtIdx:p2Tpc:" // Rc Pion2
			"kHft:p1Hft:p2Hft:kTof:p1Tof:p2Tof:sA:sB:dcakp1:dcap1p2:dcap2k:mdV0Max", BufSize);
	// nt->SetAutoSave(-500000); // autosave every 1 Mbytes
}
	
//___________
void write()
{
	result->cd();
	nt->Write();
	result->Close();
}
