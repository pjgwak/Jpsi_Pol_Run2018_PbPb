#include <ctime>
#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBin.h"
#include "tnp_weight_lowptPbPb.h"

static const long MAXTREESIZE = 1000000000000;
double getAccWeight(TH1D *h = 0, double pt = 0);
double getEffWeight(TH1D *h = 0, double pt = 0);

void SkimTree_Event_Psi2S_Cent(int nevt = -1, bool isMC = false, int kTrigSel = kTrigJpsi, int hiHFBinEdge = 0, int PDtype = 1)
{
    for (int iev = 0; iev < nevt; ++iev)
    {
        // cout << "Reco_QQ_size : " << Reco_QQ_size << endl;
        for (Int_t irqq = 0; irqq < Reco_QQ_size; ++irqq)
        {
            if (isMC)
            {
                tnp_weight = 1;
                tnp_trig_weight_mupl = -1;
                tnp_trig_weight_mumi = -1;
                tnp_weight = tnp_weight * tnp_weight_muid_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0) * tnp_weight_muid_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0); // mu id
                tnp_weight = tnp_weight * tnp_weight_trk_pbpb(mupl_Reco->Eta(), 0) * tnp_weight_trk_pbpb(mumi_Reco->Eta(), 0);                                     // inner tracker

                // Trigger part
                if (!((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))))
                {
                    //         cout << "irqq : " << irqq << " - iev : " << iev << endl;
                    //         cout << "TnP ERROR !!!! ::: No matched L2 filter1 " << endl;
                    continue;
                }
                bool mupl_L2Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))) ? true : false;
                bool mupl_L3Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter))) ? true : false;
                bool mumi_L2Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))) ? true : false;
                bool mumi_L3Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter))) ? true : false;
                if (mupl_L2Filter == false || mumi_L2Filter == false)
                {
                    cout << "TnP ERROR !!!! ::: No matched L2 filter2 " << endl;
                    cout << endl;
                }

                bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter) ? true : false;
                bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter) ? true : false;
                bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter) ? true : false;
                bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter) ? true : false;
                bool SelDone = false;

                if (mupl_isL2 && mumi_isL3)
                {
                    tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0);
                    tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0);
                    SelDone = true;
                }
                else if (mupl_isL3 && mumi_isL2)
                {
                    tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0);
                    tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0);
                    SelDone = true;
                }
                else if (mupl_isL3 && mumi_isL3)
                {
                    int t[2] = {-1, 1}; // mupl, mumi
                    int l = rand() % (2);
                    // pick up what will be L2
                    if (t[l] == -1)
                    {
                        tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0);
                        tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0);
                    }
                    else if (t[l] == 1)
                    {
                        tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0);
                        tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0);
                    }
                    else
                    {
                        cout << "ERROR :: No random selection done !!!!" << endl;
                        continue;
                    }
                    SelDone = true;
                }
                if (SelDone == false || (tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1))
                {
                    cout << "ERROR :: No muon filter combination selected !!!!" << endl;
                    continue;
                }
                tnp_weight = tnp_weight * tnp_trig_weight_mupl * tnp_trig_weight_mumi;
                counttnp++;
            }
        } // end of dimuon loop
    } // end of event loop
}
