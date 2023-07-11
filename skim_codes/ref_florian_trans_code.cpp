#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include <iostream>
#include <fstream>
#include <cmath>

// (https://twiki.cern.ch/twiki/bin/viewauth/CMS/UpsilonPolarizationInPbPb5TeV)

void UpsilonRefFrameReco2() { //version2 (In this version, I used TLorentzVector Boost class in stead of the explicit maxtrix)

	// ******** Start measuring time ******** //
	clock_t start, end, cpu_time;
	start = clock();

	// ******** Open OniaTree file ******** //
	// (To get the file, type the command below on the CERN server)
	// (xrdcp root://cms-xrd-global.cern.ch//store/user/fdamas//UpsilonPolarizationPbPb/MC/UpsilonEmbeddedMC_2018PbPb_oniatree_10_3_2/Upsilon1S_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/crab_UpsilonEmbeddedMC_2018PbPb_oniatree_10_3_2/220912_133418/0001/Oniatree_MC_numEvent1000_1342.root .)
	// TFile *infile = TFile::Open("Oniatree_MC_numEvent1000_1342.root"); //(a MC file from Florian (applied cuts are unknown))
	TFile* infile = TFile::Open("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/MiniAOD/OniaTree_MiniAOD_2018DoubleMuonPD_GlbAndTrkMuon_MuonJSON_merged.root"); //(a data file made by JaeBeom)
	TDirectoryFile* hionia = (TDirectoryFile*)gDirectory->Get("hionia");
	TTree* OniaTree = (TTree*)hionia->Get("myTree");

	// ******** Select Upsilon mass region bits ******** //
	// 2018
	// Bit1: HLT_HIL1DoubleMuOpen_v1       (Double muon inclusive)
	// Bit13: HLT_HIL3MuONHitQ10_L2MuO_MAXdR3p5_M1to5_v1  (J/psi region)
	// Bit14: HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v1 (Upsilon + high masses)
	const Int_t NTriggers = 3;
	const Int_t Bits[NTriggers] = {1, 13, 14};
	Int_t SelectedBit = 2; //(This will be used in the loop for HLTrigger and Reco_QQ_Trig)

	// ******** Define variables in the tree ******** //
	ULong64_t HLTriggers;
	ULong64_t Reco_QQ_trig[66];
	Int_t Centrality;
	TClonesArray* CloneArr_QQ;
	TClonesArray* CloneArr_mu;
	Short_t Reco_QQ_size;
	Short_t Reco_QQ_sign[1000];
	Short_t Reco_QQ_mupl_idx[1000];
	Short_t Reco_QQ_mumi_idx[1000];

	Int_t Reco_mu_SelectionType[1000];
	//(parameters for quality cuts)
	float_t Reco_QQ_VtxProb[66];
	Int_t Reco_mu_nPixWMea[66];
	Int_t Reco_mu_nTrkWMea[66];
	float_t Reco_mu_dxy[1000];
	float_t Reco_mu_dz[1000];

	CloneArr_QQ = 0;
	CloneArr_mu = 0;

	OniaTree->SetBranchAddress("HLTriggers", &HLTriggers);
	OniaTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
	OniaTree->SetBranchAddress("Centrality", &Centrality);
	OniaTree->SetBranchAddress("Reco_QQ_4mom", &CloneArr_QQ);
	OniaTree->SetBranchAddress("Reco_mu_4mom", &CloneArr_mu);
	OniaTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
	OniaTree->SetBranchAddress("Reco_QQ_sign", &Reco_QQ_sign);
	OniaTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
	OniaTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
	OniaTree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);

	OniaTree->SetBranchAddress("Reco_QQ_VtxProb", &Reco_QQ_VtxProb);
	OniaTree->SetBranchAddress("Reco_mu_nPixWMea", &Reco_mu_nPixWMea);
	OniaTree->SetBranchAddress("Reco_mu_nTrkWMea", &Reco_mu_nTrkWMea);
	OniaTree->SetBranchAddress("Reco_mu_dxy", &Reco_mu_dxy);
	OniaTree->SetBranchAddress("Reco_mu_dz", &Reco_mu_dz);

	Double_t Reco_QQ_phi, Reco_QQ_costheta, Reco_QQ_pt, Reco_QQ_y, Reco_QQ_eta, Reco_QQ_px, Reco_QQ_py, Reco_QQ_pz, Reco_QQ_E, Reco_QQ_m;
	Double_t Reco_mupl_phi, Reco_mupl_costheta, Reco_mupl_pt, Reco_mupl_y, Reco_mupl_eta, Reco_mupl_px, Reco_mupl_py, Reco_mupl_pz, Reco_mupl_E;
	Double_t Reco_mumi_phi, Reco_mumi_costheta, Reco_mumi_pt, Reco_mumi_y, Reco_mumi_eta, Reco_mumi_px, Reco_mumi_py, Reco_mumi_pz, Reco_mumi_E;

	// ******** Create a Ntuple to store kinematics of Upsilon and daughter muons ******** //
	gROOT->cd();
	TString varlist = "upsM:upsRap:upsPt:upsPz:upsEta:upsPhi:upsCosTheta:muplPt:muplPz:muplEta:muplPhi:muplCosTheta:mumiPt:mumiPz:mumiEta:mumiPhi:muplM:muplCosThetaPrimeHX:muplPhiPrimeHX:muplCosThetaPrimeCS:muplPhiPrimeCS";
	TNtuple* UpsMuNTuple = new TNtuple("UpsMuKinematics", "Upsilon in the lab frame ntuple", varlist);

	// ******** Set beam energy for the Collins-Soper reference frame ******** //
	double sqrt_S_NN = 5.02;                 //(Center of mass Energy per nucleon pair in TeV)
	double beam1_p = sqrt_S_NN * 1000. / 2.; //(in GeV) (Note. sqrt_S_NN = sqrt(2*E1*E2+2*p1*p2) = 2E1 when two beams have the same E)
	double beam1_E = beam1_p;
	double beam2_p = -beam1_p;
	double beam2_E = beam1_E;
	double delta = 0; //(Angle between ZHX(Z-axis in the Helicity frame) and ZCS(Z-axis in the Collins-Soper frame))

	// ******** Print out the total number of entries ******** //
	double totEntries = OniaTree->GetEntries();
	cout << "Total Entries: " << totEntries << endl;

	// ******** Start the event loop - read onia tree and save values for bit 14 in Ntuples ******** //
	for (int EveNum = 0; EveNum < (totEntries); EveNum++) {
		// cout << "*********************************************" << endl;

		// ******** Show how much % of the process has been completed ******** //
		if (EveNum % 100000 == 0) {
			cout << "********* " << (double)(100. * EveNum / (totEntries)) << "% completed"
			     << " ********" << endl;
		}

		// // (when wanting to test with a small number of events)
		// if(EveNum>5){
		// 	break;
		// }

		// ******** Load the values ******** //
		OniaTree->GetEntry(EveNum);

		if (!((HLTriggers & (ULong64_t)(1 << (Bits[SelectedBit] - 1))) == (ULong64_t)(1 << (Bits[SelectedBit] - 1)))) continue;
		// cout<< "Cen:" << Centrality << endl;

		for (int QQEveNum = 0; QQEveNum < Reco_QQ_size; QQEveNum++) {
			TLorentzVector* Reco_QQ_4mom = (TLorentzVector*)CloneArr_QQ->At(QQEveNum);
			TLorentzVector* Reco_mupl_4mom = (TLorentzVector*)CloneArr_mu->At(Reco_QQ_mupl_idx[QQEveNum]);
			TLorentzVector* Reco_mumi_4mom = (TLorentzVector*)CloneArr_mu->At(Reco_QQ_mumi_idx[QQEveNum]);

			// ******** Apply cuts ******** //

			bool passMuonTypePl = true;
			passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[QQEveNum]] & ((int)pow(2, 1)));
			passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[QQEveNum]] & ((int)pow(2, 3)));

			bool passMuonTypeMi = true;
			passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[QQEveNum]] & ((int)pow(2, 1)));
			passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[QQEveNum]] & ((int)pow(2, 3)));

			//(2018 Hybrid soft id cut)
			bool muplSoft = (                                          //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
			  (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[QQEveNum]] > 5) &&    // (at least 6 hits in the silicon strip layers)
			  (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[QQEveNum]] > 0) &&    // (at least 1 hit in the pixel detectors)
			  (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[QQEveNum]]) < 0.3) && // (distance btw the track and the event vertex_xy <0.3cm)
			  (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[QQEveNum]]) < 20.) &&  // (distance btw the track and the event vertex_z <20cm)
			  passMuonTypePl                                           //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
			);

			bool mumiSoft = (                                          //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
			  (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[QQEveNum]] > 5) &&    // (at least 6 hits in the silicon strip layers)
			  (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[QQEveNum]] > 0) &&    // (at least 1 hit in the pixel detectors)
			  (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[QQEveNum]]) < 0.3) && // (distance btw the track and the event vertex_xy <0.3cm)
			  (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[QQEveNum]]) < 20.) &&  // (distance btw the track and the event vertex_z <20cm)
			  passMuonTypeMi                                           //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
			);

			if (!(muplSoft && mumiSoft)) continue;

			if (Reco_QQ_VtxProb[QQEveNum] < 0.01) continue;

			if (((Reco_QQ_trig[QQEveNum] & (ULong64_t)(1 << (Bits[SelectedBit] - 1))) == (ULong64_t)(1 << (Bits[SelectedBit] - 1))) && (Reco_QQ_sign[QQEveNum]) == 0
			    // && (Reco_QQ_type[QQEveNum]==1)
			    // && ((Reco_QQ_4mom->Pt())<50)
			    // && (abs(Reco_QQ_4mom->Rapidity())<2.40)
			    && ((Reco_mupl_4mom->Pt()) > 3.5) && ((Reco_mumi_4mom->Pt()) > 3.5) && (abs(Reco_mupl_4mom->Eta()) < 2.4) && (abs(Reco_mumi_4mom->Eta()) < 2.4)
			    // && (Centrality/2. >= 10 && Centrality/2. < 90)
			) {
				// ******** Store kinematics of upsilon and muons (Lab Frame) into variables ******** //
				Reco_QQ_phi = Reco_QQ_4mom->Phi();
				Reco_QQ_costheta = Reco_QQ_4mom->CosTheta();
				Reco_QQ_pt = Reco_QQ_4mom->Pt();
				Reco_QQ_y = Reco_QQ_4mom->Rapidity();
				Reco_QQ_eta = Reco_QQ_4mom->Eta();
				Reco_QQ_px = Reco_QQ_4mom->Px();
				Reco_QQ_py = Reco_QQ_4mom->Py();
				Reco_QQ_pz = Reco_QQ_4mom->Pz();
				Reco_QQ_E = Reco_QQ_4mom->Energy();
				Reco_QQ_m = Reco_QQ_4mom->M();

				Reco_mupl_phi = Reco_mupl_4mom->Phi();
				Reco_mupl_costheta = Reco_mupl_4mom->CosTheta();
				Reco_mupl_pt = Reco_mupl_4mom->Pt();
				Reco_mupl_y = Reco_mupl_4mom->Rapidity();
				Reco_mupl_eta = Reco_mupl_4mom->Eta();
				Reco_mupl_px = Reco_mupl_4mom->Px();
				Reco_mupl_py = Reco_mupl_4mom->Py();
				Reco_mupl_pz = Reco_mupl_4mom->Pz();
				Reco_mupl_E = Reco_mupl_4mom->Energy();

				Reco_mumi_phi = Reco_mumi_4mom->Phi();
				Reco_mumi_costheta = Reco_mumi_4mom->CosTheta();
				Reco_mumi_pt = Reco_mumi_4mom->Pt();
				Reco_mumi_y = Reco_mumi_4mom->Rapidity();
				Reco_mumi_eta = Reco_mumi_4mom->Eta();
				Reco_mumi_px = Reco_mumi_4mom->Px();
				Reco_mumi_py = Reco_mumi_4mom->Py();
				Reco_mumi_pz = Reco_mumi_4mom->Pz();
				Reco_mumi_E = Reco_mumi_4mom->Energy();

				// ******** Construct 4-momentum vector of upsilon and muons (Lab Frame) ******** //
				// (documetation of TVector3 and TLorentzVector: https://root.cern.ch/root/html534/guides/users-guide/PhysicsVectors.html#lorentz-boost)
				TVector3 upsPvecLab(Reco_QQ_px, Reco_QQ_py, Reco_QQ_pz);
				TLorentzVector ups4MomLab(upsPvecLab, Reco_QQ_E);

				TVector3 muplPvecLab(Reco_mupl_px, Reco_mupl_py, Reco_mupl_pz);
				TLorentzVector mupl4MomLab(muplPvecLab, Reco_mupl_E);

				TVector3 mumiPvecLab(Reco_mumi_px, Reco_mumi_py, Reco_mumi_pz);
				TLorentzVector mumi4MomLab(mumiPvecLab, Reco_mumi_E);

				TVector3 beam1PvecLab(0, 0, beam1_p);
				TLorentzVector beam14MomLab(beam1PvecLab, beam1_E);

				TVector3 beam2PvecLab(0, 0, beam2_p);
				TLorentzVector beam24MomLab(beam2PvecLab, beam2_E);

				// cout << "<<In the lab frame>>" << endl;
				// cout << "ups: p = (" << upsPvecLab.Px() << ", " << upsPvecLab.Py()  << ", " << upsPvecLab.Pz() << ")" << endl;
				// cout << "mu+: p = (" << muplPvecLab.Px() << ", " << muplPvecLab.Py()  << ", " << muplPvecLab.Pz() << ")" << endl;
				// cout << "mu-: p = (" << mumiPvecLab.Px() << ", " << mumiPvecLab.Py()  << ", " << mumiPvecLab.Pz() << ")" << endl;

				// ******** Transform variables of muons from the lab frame to the upsilon's rest frame ******** //
				TLorentzVector ups4MomBoosted(upsPvecLab, Reco_QQ_E);
				TLorentzVector mupl4MomBoosted(muplPvecLab, Reco_mupl_E);
				TLorentzVector mumi4MomBoosted(mumiPvecLab, Reco_mumi_E);

				//(Note. TLorentzVector.BoostVector() gives beta(=px/E,py/E,pz/E) of the parents)
				//(TLorentzVector.Boost() boosts from the rod frame to the lab frame, so plug in -beta to get lab to rod)
				ups4MomBoosted.Boost(-ups4MomLab.BoostVector());
				mupl4MomBoosted.Boost(-ups4MomLab.BoostVector());
				mumi4MomBoosted.Boost(-ups4MomLab.BoostVector());

				// ******** Print out momentums of upsilon and daughter muons in the upsilon's rest frame ******** //
				// cout << endl;
				// cout << "<<Boosted to the quarkonium rest frame>>" << endl;
				// cout << "ups: p = (" << ups4MomBoosted.Px() << ", " << ups4MomBoosted.Py()  << ", " << ups4MomBoosted.Pz() << ")" << endl;
				// cout << "mu+: p = (" << mupl4MomBoosted.Px() << ", " << mupl4MomBoosted.Py()  << ", " << mupl4MomBoosted.Pz() << ")" << endl;
				// cout << "mu-: p = (" << mumi4MomBoosted.Px() << ", " << mumi4MomBoosted.Py()  << ", " << mumi4MomBoosted.Pz() << ")" << endl;

				// ******** Rotate the coordinate ******** //
				TVector3 muplPvecBoosted(mupl4MomBoosted.Px(), mupl4MomBoosted.Py(), mupl4MomBoosted.Pz());
				TVector3 mumiPvecBoosted(mumi4MomBoosted.Px(), mumi4MomBoosted.Py(), mumi4MomBoosted.Pz());

				//(Note. TVector3.Ratate() rotates the vectors not the coordinates, so should rotate -phi and -theta)
				muplPvecBoosted.RotateZ(-upsPvecLab.Phi());
				muplPvecBoosted.RotateY(-upsPvecLab.Theta());
				mumiPvecBoosted.RotateZ(-upsPvecLab.Phi());
				mumiPvecBoosted.RotateY(-upsPvecLab.Theta());

				// ******** Print out momentums of daughter muons in the upsilon's rest frame after coordinate rotation ******** //
				// cout << endl;
				// cout << "<<Rotated the quarkonium rest frame>>" << endl;
				// cout << "mu+: p = (" << muplPvecBoosted.Px() << ", " << muplPvecBoosted.Py()  << ", " << muplPvecBoosted.Pz() << ")" << endl;
				// cout << "mu-: p = (" << mumiPvecBoosted.Px() << ", " << mumiPvecBoosted.Py()  << ", " << mumiPvecBoosted.Pz() << ")" << endl;

				TLorentzVector mupl4MomBoostedRot(muplPvecBoosted, mupl4MomBoosted.E());

				// ******** HX to CS (rotation from HX frame to CS frame) ******** //
				// (1. Boost two beams to upsilon's rest frame)
				// (2. Rotate the coordinate)
				// (3. Get angle between two beams(b1 and -b2), and between b1 and ZHX in the upsilon's rest frame)
				// (4. Calculate delta (angle btw ZHX and ZCS))

				// ******** Transform variables of beams from the lab frame to the upsilon's rest frame ******** //
				TLorentzVector beam14MomBoosted(beam1PvecLab, beam1_E);
				TLorentzVector beam24MomBoosted(beam2PvecLab, beam2_E);

				// ups4MomLab.SetX(-1);
				// ups4MomLab.SetY(0);
				// ups4MomLab.SetZ(0);

				beam14MomBoosted.Boost(-ups4MomLab.BoostVector());
				beam24MomBoosted.Boost(-ups4MomLab.BoostVector());

				// ******** Print out momentums of two beams in the upsilon's rest frame ******** //
				// cout << endl;
				// cout << "<<Boosted to the quarkonium rest frame>>" << endl;
				// cout << "ups: p = (" << ups4MomBoosted.Px() << ", " << ups4MomBoosted.Py()  << ", " << ups4MomBoosted.Pz() << ")" << endl;
				// cout << "beam1: p = (" << beam14MomBoosted.Px() << ", " << beam14MomBoosted.Py()  << ", " << beam14MomBoosted.Pz() << ")" << endl;
				// cout << "beam2: p = (" << beam24MomBoosted.Px() << ", " << beam24MomBoosted.Py()  << ", " << beam24MomBoosted.Pz() << ")" << endl;

				// ******** Rotate the coordinate ******** //
				TVector3 beam1PvecBoosted(beam14MomBoosted.Px(), beam14MomBoosted.Py(), beam14MomBoosted.Pz());
				TVector3 beam2PvecBoosted(beam24MomBoosted.Px(), beam24MomBoosted.Py(), beam24MomBoosted.Pz());

				// upsPvecLab.SetX(-1);
				// upsPvecLab.SetY(0);
				// upsPvecLab.SetZ(0);

				beam1PvecBoosted.RotateZ(-upsPvecLab.Phi());
				beam1PvecBoosted.RotateY(-upsPvecLab.Theta());
				beam2PvecBoosted.RotateZ(-upsPvecLab.Phi());
				beam2PvecBoosted.RotateY(-upsPvecLab.Theta());

				// ******** Print out momentums of daughter muons in the upsilon's rest frame after coordinate rotation ******** //
				// cout << endl;
				// cout << "<<Rotated the quarkonium rest frame>>" << endl;
				// cout << "beam1: p = (" << beam1PvecBoosted.Px() << ", " << beam1PvecBoosted.Py()  << ", " << beam1PvecBoosted.Pz() << ")" << endl;
				// cout << "beam2: p = (" << beam2PvecBoosted.Px() << ", " << beam2PvecBoosted.Py()  << ", " << beam2PvecBoosted.Pz() << ")" << endl;

				// ******** Calculate the angle between z_HX and z_CS ******** //
				TVector3 ZHXunitVec(0, 0, 1);                                    //(define z_HX unit vector)
				double Angle_B1ZHX = beam1PvecBoosted.Angle(ZHXunitVec);         //(angle between beam1 and z_HX)
				double Angle_B2ZHX = beam2PvecBoosted.Angle(-ZHXunitVec);        //(angle between beam2 and -z_HX =(-beam2 and z_HX) )
				double Angle_B1miB2 = beam1PvecBoosted.Angle(-beam2PvecBoosted); //(angle between beam1 and -beam2)

				double delta = 0; //(define and initialize the angle between z_HX and z_CS)

				// // (The math for caculating the angle between z_HX and z_CS is different depending on the sign of the beam1's z-coordinate)
				// if(beam1PvecBoosted.Pz()>0) delta = Angle_B1ZHX + Angle_B1miB2/2.;
				// else if(beam1PvecBoosted.Pz()<0) delta = Angle_B1ZHX - Angle_B1miB2/2.;
				// else cout <<  "beam1PvecBoosted.Pz() = 0?" << endl;
				if (Angle_B1ZHX > Angle_B2ZHX)
					delta = Angle_B2ZHX + Angle_B1miB2 / 2.;
				else if (Angle_B1ZHX < Angle_B2ZHX)
					delta = Angle_B1ZHX + Angle_B1miB2 / 2.;
				else
					cout << "beam1PvecBoosted.Pz() = 0?" << endl;

				// ******** Print out the angles ******** //
				// cout << endl;
				// cout << "angle between ZHX and b1: " << Angle_B1ZHX << endl;
				// cout << "angle between b1 and -b2: " << Angle_B1miB2 << " (half: " << (Angle_B1miB2)/2. << ")"<< endl;
				// cout << "angle between ZHX and ZCS: " << delta << " (" << delta*180./M_PI << "deg)" << endl;

				// ******** Rotate the coordinate along the y-axis by the angle between z_HX and z_CS ******** //
				TVector3 muplPvecBoostedCS(muplPvecBoosted.Px(), muplPvecBoosted.Py(), muplPvecBoosted.Pz());
				TVector3 mumiPvecBoostedCS(mumiPvecBoosted.Px(), mumiPvecBoosted.Py(), mumiPvecBoosted.Pz());

				ZHXunitVec.RotateY(delta);
				muplPvecBoostedCS.RotateY(delta);
				mumiPvecBoostedCS.RotateY(delta);

				// cout << endl;
				// cout << "Rotated unit Vec: (" << ZHXunitVec.Px() << ", " << ZHXunitVec.Py() << ", " << ZHXunitVec.Pz() << ")" << endl;
				// cout << "mupl CosTheta, mupl Phi: " << muplPvecBoostedCS.CosTheta() << ", " << muplPvecBoostedCS.Phi() << endl;

				// ******** Fill Ntuple with kinematics of upsilon and muons in the Lab, HX, and CS frames******** //
				float tuple[] = {
				  static_cast<float>(Reco_QQ_m),
				  static_cast<float>(Reco_QQ_y),
				  static_cast<float>(Reco_QQ_pt),
				  static_cast<float>(Reco_QQ_pz),
				  static_cast<float>(Reco_QQ_eta),
				  static_cast<float>(Reco_QQ_phi),
				  static_cast<float>(Reco_QQ_costheta),

				  static_cast<float>(Reco_mupl_pt),
				  static_cast<float>(Reco_mupl_pz),
				  static_cast<float>(Reco_mupl_eta),
				  static_cast<float>(Reco_mupl_phi),
				  static_cast<float>(Reco_mupl_costheta),
				  static_cast<float>(Reco_mumi_pt),
				  static_cast<float>(Reco_mumi_pz),
				  static_cast<float>(Reco_mumi_eta),
				  static_cast<float>(Reco_mumi_phi),

				  static_cast<float>(mupl4MomBoostedRot.Mag()),
				  static_cast<float>(muplPvecBoosted.CosTheta()),
				  static_cast<float>(muplPvecBoosted.Phi()),

				  static_cast<float>(muplPvecBoostedCS.CosTheta()),
				  static_cast<float>(muplPvecBoostedCS.Phi())};

				UpsMuNTuple->Fill(tuple);
			}
		}
	}

	// ******** Create a file and store the ntuples ******** //
	TFile* file = new TFile("Upsilon1S_Reference_Frames2_Reco.root", "RECREATE", "Upsilon 1S");

	UpsMuNTuple->Write();

	file->Close();

	// ******** End measuring time ******** //
	end = clock();
	cpu_time = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << cpu_time / 60. << "minutes" << endl;

	return;
}