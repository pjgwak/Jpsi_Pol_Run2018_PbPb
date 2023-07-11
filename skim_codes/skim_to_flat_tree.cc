// ##################################################### //
// Extract indispensable variables for study
// Apply weight if needed
// Input: a tree you got by running "miniAOD_to_first_skim.cc"
// Output: tree for polarization study
// ##################################################### //

// Add ROOT Headers for compile


// Add namespaces for compile


// Refering the /Users/pjgwak/work/psi2S_Raa_Run2018_at2023/makeRooDataSet_psi2S_Data.C

// Start main function
void skim_to_flat_tree() // include options
{

    auto infile = TFile::Open("data_skim_tuple.root");
    auto intree = (TTree*) infile->Get("outtree");


    // ========== Prepare branches for input tree ========== //
    // int evt;
    // int runN;
    // int lumi;
    int cBin;
    int nDimu;
    // float vz;
    long const max_branch_size = 100;
    float mass[max_branch_size];
    float y[max_branch_size];
    float pt[max_branch_size];
    float pt1[max_branch_size];
    float pt2[max_branch_size];
    float eta[max_branch_size];
    float eta1[max_branch_size];
    float eta2[max_branch_size];
    // float phi[max_branch_size];
    // float phi1[max_branch_size];
    // float phi2[max_branch_size];
    // float weight0[max_branch_size];
    // float weight1[max_branch_size];
    int recoQQsign[max_branch_size];
    float ctau3D[max_branch_size];
    // float ctau3DErr[max_branch_size];
    float ctau3DRes[max_branch_size];
    // float ctau3D2S[max_branch_size];
    // float ctau3DErr2S[max_branch_size];
    // float ctau3DRes2S[max_branch_size];
    // double TnPweight[max_branch_size] = {1.};
    // double weight = 1;
    double cosTheta_cs[max_branch_size];

    TBranch *b_cBin;
    TBranch *b_nDimu;
    TBranch *b_mass;
    TBranch *b_y;
    TBranch *b_pt;
    TBranch *b_pt1;
    TBranch *b_pt2;
    TBranch *b_eta;
    TBranch *b_eta1;
    TBranch *b_eta2;
    TBranch *b_recoQQsign;
    TBranch *b_ctau3D;
    TBranch *b_ctau3DRes;
    TBranch *b_cosTheta_cs;

    // intree->SetBranchAddress("event", &evt, &event);
    // intree->SetBranchAddress("runN", &runN, &runN);
    // intree->SetBranchAddress("lumi", &lumi, &lumi);
    intree->SetBranchAddress("cBin", &cBin, &b_cBin);
    // intree->SetBranchAddress("vz", &vz, &vz);
    intree->SetBranchAddress("nDimu", &nDimu, &b_nDimu);
    intree->SetBranchAddress("mass", mass, &b_mass);
    intree->SetBranchAddress("y", y, &b_y);
    intree->SetBranchAddress("pt", pt, &b_pt);
    intree->SetBranchAddress("pt1", pt1, &b_pt1);
    intree->SetBranchAddress("pt2", pt2, &b_pt2);
    intree->SetBranchAddress("eta", eta, &b_eta);
    intree->SetBranchAddress("eta1", eta1, &b_eta1);
    intree->SetBranchAddress("eta2", eta2, &b_eta2);
    intree->SetBranchAddress("recoQQsign", recoQQsign, &b_recoQQsign);
    intree->SetBranchAddress("ctau3D", ctau3D, &b_ctau3D);
    // intree->SetBranchAddress("ctau3DErr", ctau3DErr, &ctau3DErr);
    intree->SetBranchAddress("ctau3DRes", ctau3DRes, &b_ctau3DRes);
    // intree->SetBranchAddress("ctau3D2S", ctau3D2S, &ctau3D2S);
    // intree->SetBranchAddress("ctau3DErr2S", ctau3DErr2S, &ctau3DErr2S);
    // intree->SetBranchAddress("ctau3DRes2S", ctau3DRes2S, &ctau3DRes2S);
    // intree->SetBranchAddress("weight", &weight, &weight);
    // intree->SetBranchAddress("TnPweight", TnPweight, &TnPweight);
    intree->SetBranchAddress("cosTheta_cs", cosTheta_cs, &b_cosTheta_cs);



    // ========== Prepare branches for output tree ========== //
    auto outfile = new TFile("skim_flat.root", "recreate");
    outfile->cd();

    int out_cBin;
    int out_nDimu;
    // float vz;
    float out_mass;
    float out_y;
    float out_pt;
    float out_pt1;
    float out_pt2;
    float out_eta;
    float out_eta1;
    float out_eta2;
    // float phi;
    // float phi1;
    // float phi2;
    // float weight0;
    // float weight1;
    int out_recoQQsign;
    float out_ctau3D;
    // float ctau3DErr;
    float out_ctau3DRes;
    // float ctau3D2S;
    // float ctau3DErr2S;
    // float ctau3DRes2S;
    // double TnPweight = {1.};
    // double weight = 1;
    double out_cosTheta_cs;

    TTree *outtree = new TTree("outtree", "skim test");
    outtree->SetMaxTreeSize(100000000);
    // outtree->Branch("event", &evt, "event");
    // outtree->Branch("runN", &runN, "runN");
    // outtree->Branch("lumi", &lumi, "lumi");
    outtree->Branch("cBin", &out_cBin, "cBin/I");
    // outtree->Branch("vz", &vz, "vz");
    outtree->Branch("nDimu", &out_nDimu, "nDimu/I");
    outtree->Branch("mass", &out_mass, "mass");
    outtree->Branch("y", &out_y, "y");
    outtree->Branch("pt", &out_pt, "pt");
    outtree->Branch("pt1", &out_pt1, "pt1");
    outtree->Branch("pt2", &out_pt2, "pt2");
    outtree->Branch("eta", &out_eta, "eta");
    outtree->Branch("eta1", &out_eta1, "eta1");
    outtree->Branch("eta2", &out_eta2, "eta2");
    outtree->Branch("recoQQsign", &out_recoQQsign, "recoQQsign/I");
    outtree->Branch("ctau3D", &out_ctau3D, "ctau3D");
    // outtree->Branch("ctau3DErr", ctau3DErr, "ctau3DErr");
    outtree->Branch("ctau3DRes", &out_ctau3DRes, "ctau3DRes");
    // outtree->Branch("ctau3D2S", ctau3D2S, "ctau3D2S");
    // outtree->Branch("ctau3DErr2S", ctau3DErr2S, "ctau3DErr2S");
    // outtree->Branch("ctau3DRes2S", ctau3DRes2S, "ctau3DRes2S");
    // outtree->Branch("weight", &weight, "weight");
    // outtree->Branch("TnPweight", TnPweight, "TnPweight");
    outtree->Branch("cosTheta_cs", &out_cosTheta_cs, "cosTheta_cs/D");

    long const n_evt = intree->GetEntries();
    cout << "n_evt : " << n_evt << endl;

    int nDimu_all=0;
    int nDimuPass=0;
    int nDimu_one=0;
    int nDimu_more=0;
    //Begin Loop
    for(long evt_idx = 0; evt_idx < n_evt; evt_idx++) {
        intree->GetEntry(evt_idx);

        // Fill Dimuon Loop
        for(int dimu_idx = 0; dimu_idx < nDimu; dimu_idx++) {
            out_cBin = cBin;
            out_nDimu = nDimu;
            out_mass = mass[dimu_idx];
            out_y = y[dimu_idx];
            out_pt = pt[dimu_idx];
            out_pt1 = pt1[dimu_idx];
            out_pt2 = pt2[dimu_idx];
            out_eta = eta[dimu_idx]; 
            out_eta1 = eta1[dimu_idx]; 
            out_eta2 = eta2[dimu_idx];
            out_recoQQsign = recoQQsign[dimu_idx];
            out_ctau3D = ctau3D[dimu_idx];
            out_ctau3DRes = ctau3DRes[dimu_idx];
            out_cosTheta_cs = cosTheta_cs[dimu_idx];
            outtree->Fill();
        }
    }
    outfile->Write();
    outfile->Close();
}
