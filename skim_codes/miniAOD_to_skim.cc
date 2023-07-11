// ##################################################### //
// The sizes of miniAOD samples are big so reduce file sizes
// Add inforamtion about transformation form Lab frame to HX and CS frame. 
// Input: miniAOD smaples of CMS run2 in 2018
// Output: tree
// ##################################################### //


// Add ROOT Headers for compile
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

// Add namespaces for compile
using namespace ROOT;
using std::cout; using std::endl;

// Refer to /Users/pjgwak/work/psi2S_Raa_Run2018_at2023/SkimTree_Event_Psi2S_Cent.C
void miniAOD_to_first_skim();
void infile_check(TFile* infile);


int main() {
    miniAOD_to_first_skim();
    return 0;
}


// Start main function
void miniAOD_to_first_skim() // include options
{
    auto input_path = "../skimmedFiles/OniaTree_miniAOD_HIDoubleMuonPD_addQVect_merged.root";
    auto infile = TFile::Open(input_path);
    auto intree = (TTree*)infile->Get("hionia/myTree");

    // Refer to "Eff_Acc/cutsAndBin.h:int kTrigJpsi = 12;"
    int kTrigSel = 12;

    // intree -> cpp 변수 준비
    int max_branch_size = 100; // Any enough big number

    UInt_t runNb;
    ULong64_t HLTriggers; // HLT triggers
    ULong64_t Reco_QQ_trig[max_branch_size]; // QQ trigger (?)
    Double_t px,py,pt;
	Double_t jpsi_mass;
    Short_t Reco_QQ_size; // 이벤트 별 엔트리 개수 (# of Jpsi ?)
    Short_t Reco_QQ_sign[max_branch_size];
    TLorentzVector* jpsi_reco = nullptr; // Will be connected with Reco_QQ_4mom
    TLorentzVector* mupl_reco = nullptr; // Will be connected with Reco_mu_4mom
    TLorentzVector* mumi_reco = nullptr; // Will be connected with Reco_mu_4mom
    TClonesArray *Reco_QQ_4mom = nullptr; // Initialization is mandatory for TClonesArray
    TClonesArray *Reco_mu_4mom = nullptr; 
    Float_t Reco_QQ_VtxProb[max_branch_size];
    Short_t Reco_QQ_mupl_idx[max_branch_size]; // positive muon id
    Short_t Reco_QQ_mumi_idx[max_branch_size]; // negative muon id
    Int_t Reco_mu_SelectionType[max_branch_size];
    Int_t Reco_mu_nTrkWMea[max_branch_size];
    Int_t Reco_mu_nPixWMea[max_branch_size];
    Float_t Reco_mu_dxy[max_branch_size];
    Float_t Reco_mu_dz[max_branch_size];

    // intree TBranch와 cpp 변수 연결
    intree->SetBranchAddress("runNb", &runNb);
    intree->SetBranchAddress("HLTriggers", &HLTriggers);
    intree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
    intree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
    intree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign);
    intree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
    intree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom);
    intree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb);
    intree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
    intree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
    intree->SetBranchAddress("Reco_mu_SelectionType", &Reco_mu_SelectionType);
    intree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea);
    intree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea);
    intree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy);
    intree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz);


    // outfile, outtree 생성 및 TBranch 연결
    auto outfile = new TFile("tree1.root", "recreate");
    outfile->cd();
    // int evt;
    // int runN;
    // int lumi;
    // int cBin;
    int nDimu;
    // float vz;
    float mass[max_branch_size];
    // float pt[max_branch_size];
    // float y[max_branch_size];
    // float phi[max_branch_size];
    // float eta[max_branch_size];
    // float eta1[max_branch_size];
    // float eta2[max_branch_size];
    // float phi1[max_branch_size];
    // float phi2[max_branch_size];
    // float pt1[max_branch_size];
    // float pt2[max_branch_size];
    // float weight0[max_branch_size];
    // float weight1[max_branch_size];
    int recoQQsign[max_branch_size];
    // float ctau3D[max_branch_size];
    // float ctau3DErr[max_branch_size];
    // float ctau3DRes[max_branch_size];
    // float ctau3D2S[max_branch_size];
    // float ctau3DErr2S[max_branch_size];
    // float ctau3DRes2S[max_branch_size];
    // double TnPweight[max_branch_size] = {1.};
    // double weight = 1;

    TTree *outtree = new TTree("outtree", "skim test");
    outtree->SetMaxTreeSize(100000000);
    // outtree->Branch("event", &evt, "event");
    // outtree->Branch("runN", &runN, "runN");
    // outtree->Branch("lumi", &lumi, "lumi");
    // outtree->Branch("cBin", &cBin, "cBin");
    // outtree->Branch("vz", &vz, "vz");
    outtree->Branch("nDimu", &nDimu, "nDimu/I");
    outtree->Branch("mass", mass, "mass[nDimu]/F");
    // outtree->Branch("y", y, "y[nDimu]");
    // outtree->Branch("pt", pt, "pt[nDimu]");
    // outtree->Branch("pt1", pt1, "pt1[nDimu]");
    // outtree->Branch("pt2", pt2, "pt2[nDimu]");
    // outtree->Branch("eta", eta, "eta[nDimu]");
    // outtree->Branch("eta1", eta1, "eta1[nDimu]");
    // outtree->Branch("eta2", eta2, "eta2[nDimu]");
    outtree->Branch("recoQQsign", recoQQsign, "recoQQsign[nDimu]/I");
    // outtree->Branch("ctau3D", ctau3D, "ctau3D[nDimu]");
    // outtree->Branch("ctau3DErr", ctau3DErr, "ctau3DErr[nDimu]");
    // outtree->Branch("ctau3DRes", ctau3DRes, "ctau3DRes[nDimu]");
    // outtree->Branch("ctau3D2S", ctau3D2S, "ctau3D2S[nDimu]");
    // outtree->Branch("ctau3DErr2S", ctau3DErr2S, "ctau3DErr2S[nDimu]");
    // outtree->Branch("ctau3DRes2S", ctau3DRes2S, "ctau3DRes2S[nDimu]");
    // outtree->Branch("weight", &weight, "weight");
    // outtree->Branch("TnPweight", TnPweight, "TnPweight[nDimu]");

    // 트리 안의 이벤트 읽기
    long n_evt = 1400;
    n_evt = intree->GetEntries();
    long total_dimuon = 0;
    for (long idx_evt = 0; idx_evt < n_evt; idx_evt++)
    {
        if (idx_evt % 100000 == 0)
            cout << ">>>>> EVENT " << idx_evt << " / " << intree->GetEntries() << " (" << (int)(100. * idx_evt / intree->GetEntries()) << "%)" << endl;

        intree->GetEntry(idx_evt);

        if (runNb >= 327123)
        {
            // Ignore not interesting CMS runs
            continue;
        }

        if (!(
                (HLTriggers & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel))))
        {
            // Exclude events which not matched with Jpsi HLT trigger
            // 주의! && (논리) 연산자가 아니라 & (비트) 연산자 사용하는 게 맞다. HLTriggers에 저장된 값을 읽는 거라서.
            continue;
        }

        // 한 이벤트 안에 있는 Jpsi 수 확인하고 작업하기-> Reco_QQ_size 이용
        // 이 코드에서 가장 중요한 부분이다.
        nDimu = 0;
        for (int idx_qq = 0; idx_qq < Reco_QQ_size; idx_qq++)
        {
            jpsi_reco = (TLorentzVector *)Reco_QQ_4mom->At(idx_qq); // Notice it's Reco_"QQ"_4mom
            mupl_reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mupl_idx[idx_qq]);
            mumi_reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mumi_idx[idx_qq]);
            // cout << "Jpsi(dimuon) invariant mass (evt: " << idx_evt << ", entry: " << idx_qq << "): " << jpsi_mass << endl;

            if ( !(
                    (Reco_QQ_trig[idx_qq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel))))
            {
                // Exclude events which not matched with trigger
                continue;
            }

            // Check muons' signs
            bool passMuonTypePl = true;
            passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[idx_qq]] & ((int)pow(2, 1)));
            passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[idx_qq]] & ((int)pow(2, 3)));

            bool passMuonTypeMi = true;
            passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[idx_qq]] & ((int)pow(2, 1)));
            passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[idx_qq]] & ((int)pow(2, 3)));

            // Check soft muons
            bool muplSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[idx_qq]]==true) &&
                (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[idx_qq]] > 5) &&
                (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[idx_qq]] > 0) &&
                (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[idx_qq]]) < 0.3) &&
                (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[idx_qq]]) < 20.) &&
                passMuonTypePl //                       &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[idx_qq]]==true)
            );

            bool mumiSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[idx_qq]]==true) &&
                (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[idx_qq]] > 5) &&
                (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[idx_qq]] > 0) &&
                (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[idx_qq]]) < 0.3) &&
                (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[idx_qq]]) < 20.) &&
                passMuonTypeMi //                        &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[idx_qq]]==true)
            );

            if (!(muplSoft && mumiSoft))
            {
                // Exclude if there is non-soft muon
                continue;
            }

            if ((Reco_QQ_VtxProb[idx_qq]) < 0.01)
            {
                // Exclude low probability for vertex fit
                continue;
            }
            
            if ( !(jpsi_reco->M() > 2.6 && jpsi_reco->M() < 3.4)) {
                // Allow Jpsi. Only for test.
                continue;
            }

            // 트리에 저장할 준비
            recoQQsign[idx_qq] = Reco_QQ_sign[idx_qq];
            mass[idx_qq] = jpsi_reco->M();
            //cout << mass[idx_qq] << endl;

            // phi[idx_qq] = JP_Reco->Phi();
            // phi1[idx_qq] = mupl_Reco->Phi();
            // phi2[idx_qq] = mumi_Reco->Phi();
            // eta[idx_qq] = JP_Reco->Eta();
            // y[idx_qq] = JP_Reco->Rapidity();
            // pt[idx_qq] = JP_Reco->Pt();
            // pt1[idx_qq] = mupl_Reco->Pt();
            // pt2[idx_qq] = mumi_Reco->Pt();
            // eta1[idx_qq] = mupl_Reco->Eta();
            // eta2[idx_qq] = mumi_Reco->Eta();
            // ctau3D[idx_qq] = Reco_QQ_ctau3D[idx_qq];
            // ctau3DErr[idx_qq] = Reco_QQ_ctauErr3D[idx_qq];
            // ctau3DRes[idx_qq] = (Reco_QQ_ctau3D[idx_qq]) / (Reco_QQ_ctauErr3D[idx_qq]);
            // ctau3D2S[idx_qq] = ctau3D[idx_qq] * (pdgMass.Psi2S / pdgMass.JPsi);
            // ctau3DErr2S[idx_qq] = ctau3DErr[idx_qq] * (pdgMass.Psi2S / pdgMass.JPsi);
            // ctau3DRes2S[idx_qq] = ctau3DRes[idx_qq] * (pdgMass.Psi2S / pdgMass.JPsi);
            nDimu++;
        } 
        // End of one event
        // Fill the tree at this point
        if (nDimu > 0) {
            // Fill the tree only when there is at least one muon
            outtree->Fill();
            total_dimuon += nDimu;
        }
        if (total_dimuon > 10000) {
            // Only for test
            break;
        }
    }
    // Save output
    outtree->Write();
    outfile->Close();
}
    
void infile_check(TFile* infile)
{
    if (infile->IsOpen()) {
        cout << "== File is opend ==" << endl;
    }
    else {
        cout << "== File isn't opend ==" << endl;
        exit(1);
    }
};