void read_tree()
{
    auto infile = TFile::Open("tree1.root");
    auto intree = (TTree*) infile->Get("outtree");
    float mass[1000];
    int nDimu;
    intree->SetBranchAddress("mass", mass);
    intree->SetBranchAddress("nDimu", &nDimu);
    
    for (long idx_evt = 0; idx_evt < intree->GetEntries(); idx_evt++) {
        intree->GetEntry(idx_evt);
        for (long idx_dimuon = 0; idx_dimuon < nDimu; idx_dimuon++)
        {
            cout << mass[idx_dimuon] << endl;
        }
    }
}