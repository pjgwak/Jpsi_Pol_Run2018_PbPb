#include <iostream>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include <TMath.h>
#include "commonUtility.h"
#include "cutsAndBin.h"
#include "HiEvtPlaneList.h"
#include "Style.h"
#include "tdrstyle.C"
#include "CMS_lumi_v2mass.C"
#include "rootFitHeaders.h"

using namespace std;
using namespace RooFit;
using namespace hi;

static const long MAXTREESIZE = 1000000000000;

double getAccWeight(TH1D* h = 0, double pt = 0);
double getEffWeight(TH1D* h = 0, double pt = 0);
void GetHistSqrt(TH1D* h1 =0, TH1D* h2=0);
void Get2DHistSqrt(TH2D* h1 =0, TH2D* h2=0);
void GetHistBkg(TH1D* h1 =0, TH1D* h2=0);

void makeRooDataSet_psi2S_Data(
    int cLow = 0, int cHigh = 200,
    float massLow = 3.3, float massHigh =4.1, 
    bool dimusign=true, bool isMC = false, 
    bool fAccW = false, bool fEffW = false,
    bool isTnP = false, bool isPtW = false,
    int hiHFBinEdge = 0,
    int dtype = 1, 
    int weight_PR = 0
    ) 
{
  //Basic Setting
  gStyle->SetOptStat(0);

  TString DATE;
  //if(ptLow==6.5&&ptHigh==50) DATE=Form("%i_%i",0,180);
  DATE=Form("%i_%i",cLow/2,cHigh/2);
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir(Form("figs/nominal_Psi2S/v2mass_hist/%s",DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/nominal_Psi2S/q_vector/%s",DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/nominal_Psi2S/mass_dist/%s",DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/nominal_Psi2S/decayL/%s",DATE.Data()), kTRUE);
  gSystem->mkdir(Form("roots/nominal_Psi2S/v2mass_hist/%s",DATE.Data()), kTRUE);
  gStyle->SetOptStat(000000000);
  gROOT->ForceStyle();
  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  TH1::SetDefaultSumw2();
  //TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh, 1) ;
  //TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0, cLow, cHigh) ;
  TString dimusignString;
  if(dimusign) dimusignString = "OS";
  else if(!dimusign) dimusignString = "SS";

  TString sample;
  TString wName;
  TString cutName;

  if(weight_PR==0) wName="prompt";
  else if(weight_PR==1) wName="nprompt";

  TChain *tree = new TChain("mmepevt");
  if(!isMC){
    TString f1 = "./skimmedFiles/OniaFlowSkim_JpsiTrig_DBAllPD_miniADO_isMC0_HFNom_230513.root";
    tree->Add(f1.Data());
    sample="RD";}
  else if(dtype==1){
    TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/";
    tree->Add(f1.Data());
    sample="MC_PR";}
  else if(dtype==2){
    TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/";
    tree->Add(f1.Data());
    sample="MC_NP";}

  TFile *fEff;
  //if(cLow==0&&cHigh==20) fEff= new TFile(Form("./Eff_Acc/mc_eff_vs_pt_cent_0_to_20_rap_prompt_pbpb_psi2s_PtW%d_tnp%d_20220125.root",isPtW,isTnP),"read");
  //else if(cLow==20&&cHigh==120) fEff= new TFile(Form("./Eff_Acc/mc_eff_vs_pt_cent_20_to_120_rap_prompt_pbpb_psi2s_PtW%d_tnp%d_20220125.root",isPtW,isTnP),"read");
  //else if(cLow==0&&cHigh==200) fEff= new TFile(Form("./Eff_Acc/mc_eff_vs_pt_cent_0_to_200_rap_prompt_pbpb_psi2s_PtW%d_tnp%d_20221116.root",isPtW,isTnP),"read");
  fEff= new TFile(Form("./Eff_Acc/mc_eff_vs_pt_cent_0_to_200_rap_prompt_pbpb_psi2s_PtW%d_tnp%d_20221116.root",isPtW,isTnP),"read");

  TH1D* hEffPt[4];
  /*
  if(cLow==0&&cHigh==20){
    hEffPt[0] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy0_1p2",isTnP,isPtW));
    hEffPt[1] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy1p2_1p6",isTnP,isPtW));
    hEffPt[2] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy1p6_1p8",isTnP,isPtW));
    hEffPt[3] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_20_absy1p8_2p4",isTnP,isPtW));}
  else if(cLow==20&&cHigh==120){
    hEffPt[0] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy0_1p2",isTnP,isPtW));
    hEffPt[1] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p2_1p6",isTnP,isPtW));
    hEffPt[2] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p6_1p8",isTnP,isPtW));
    hEffPt[3] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p8_2p4",isTnP,isPtW));}
  else if(cLow==0&&cHigh==200){
    hEffPt[0] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_200_absy0_1p2",isTnP,isPtW));
    hEffPt[1] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_200_absy1p2_1p6",isTnP,isPtW));
    hEffPt[2] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_200_absy1p6_1p8",isTnP,isPtW));
    hEffPt[3] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_200_absy1p8_2p4",isTnP,isPtW));}
  */
    hEffPt[0] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_200_absy0_1p2",isTnP,isPtW));
    hEffPt[1] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_200_absy1p2_1p6",isTnP,isPtW));
    hEffPt[2] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_200_absy1p6_1p8",isTnP,isPtW));
    hEffPt[3] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_200_absy1p8_2p4",isTnP,isPtW));

  TFile *fAcc = new TFile(Form("./Eff_Acc/acceptance_Prompt_Psi2s_GenOnly_wgt%d_20220125.root",isPtW),"read");
  TH1D* hAccPt[3];
  hAccPt[0] = (TH1D*) fAcc -> Get("hAccPt_2021_ally");
  hAccPt[1] = (TH1D*) fAcc -> Get("hAccPt_2021_midy");
  hAccPt[2] = (TH1D*) fAcc -> Get("hAccPt_2021_Fory");

  //SetBranchAddress
  Int_t event; 
  Int_t cBin;
  Int_t nDimu; 
  const int nMaxDimu = 1000;
  int recoQQsign[nMaxDimu];
  float vz;
  float mass[nMaxDimu];
  float pt[nMaxDimu];
  float y[nMaxDimu]; 
  float pt1[nMaxDimu];
  float pt2[nMaxDimu];
  float eta[nMaxDimu];
  float eta1[nMaxDimu];
  float eta2[nMaxDimu];
  float ctau3D[nMaxDimu];
  float ctau3DErr[nMaxDimu];
  float ctau3D2S[nMaxDimu];
  float ctau3DErr2S[nMaxDimu];
  float ctau3DRes2S[nMaxDimu];
  double weight;

  TBranch *b_event;
  TBranch *b_cBin;
  TBranch *b_nDimu;
  TBranch *b_recoQQsign;
  TBranch *b_vz;
  TBranch *b_mass;
  TBranch *b_pt;
  TBranch *b_y;
  TBranch *b_eta;
  TBranch *b_eta1;
  TBranch *b_eta2;
  TBranch *b_ctau3D;
  TBranch *b_ctau3DErr;
  TBranch *b_ctau3D2S;
  TBranch *b_ctau3DErr2S;
  TBranch *b_ctau3DRes2S;
  TBranch *b_pt1;
  TBranch *b_pt2;
  TBranch *b_weight;

  tree->SetBranchAddress("event", &event, &b_event);
  tree->SetBranchAddress("cBin", &cBin, &b_cBin);
  tree->SetBranchAddress("nDimu", &nDimu, &b_nDimu);
  tree->SetBranchAddress("recoQQsign", recoQQsign, &b_recoQQsign);
  tree->SetBranchAddress("vz", &vz, &b_vz);
  tree->SetBranchAddress("mass", mass, &b_mass);
  tree->SetBranchAddress("y", y, &b_y);
  tree->SetBranchAddress("pt", pt, &b_pt);
  tree->SetBranchAddress("pt1", pt1, &b_pt1);
  tree->SetBranchAddress("pt2", pt2, &b_pt2);
  tree->SetBranchAddress("eta", eta, &b_eta);
  tree->SetBranchAddress("eta1", eta1, &b_eta1);
  tree->SetBranchAddress("eta2", eta2, &b_eta2);
  tree->SetBranchAddress("ctau3D", ctau3D, &b_ctau3D);
  tree->SetBranchAddress("ctau3DErr", ctau3DErr, &b_ctau3DErr);
  tree->SetBranchAddress("ctau3D2S", ctau3D2S, &b_ctau3D2S);
  tree->SetBranchAddress("ctau3DErr2S", ctau3DErr2S, &b_ctau3DErr2S);
  tree->SetBranchAddress("ctau3DRes2S", ctau3DRes2S, &b_ctau3DRes2S);
  tree->SetBranchAddress("weight", &weight, &b_weight);

  const int nV2Bin = 6;

  float massBinDiff[nV2Bin+1]={-3,-2,-1,0,1,2,3};
  float massBin_[nV2Bin+1];

  //kineLabel = kineLabel + Form("_m%.1f-%.1f",massLow,massHigh) + "_" + dimusignString;
  //cout<<kineLabel<<endl;



  //Define drawing histogram
  TH1D* h_ljpsi;

  double nMassBin = 45;
  double nQBin = 100;
  double Q_avg_low = -10000;
  double Q_avg_high = 40000;
  double Q_avg_low_dimu = -200;
  double Q_avg_high_dimu = 200;
  h_ljpsi = new TH1D("h_ljBkg",";l_{J/#psi(mm)};Counts/(0.03 mm)",1050,-1.5,2);

  const static int countMax = 1000000;


  ////////////////////////////////////////////////////////////////////////
  //////////////////  RooDataSet
  ////////////////////////////////////////////////////////////////////////
  RooRealVar* massVar  = new RooRealVar("mass","mass variable",1.0, 6.0,"GeV/c^{2}");
  RooRealVar* ptVar    = new RooRealVar("pt","pt variable", 0,100,"GeV/c");
  RooRealVar* yVar     = new RooRealVar("y","rapidity of the dimuon pair", -5,5,"");
  RooRealVar* pt1Var   = new RooRealVar("pt1","pt of muon+", 0,500,"GeV/c");
  RooRealVar* eta1Var  = new RooRealVar("eta1","eta of muon+", -4,4,"");
  RooRealVar* pt2Var   = (RooRealVar*)pt1Var->Clone("pt2");
  RooRealVar* eta2Var  = (RooRealVar*)eta1Var->Clone("eta2");
  RooRealVar* cBinVar   = new RooRealVar("cBin","Centrality bin", -100,500,"");
  RooRealVar* ep2Var   = new RooRealVar("ep2","2nd order event plane", -100,100,"");
  RooRealVar* evtWeight = new RooRealVar("weight","corr weight", 0, 10000,"");
  RooRealVar* recoQQ = new RooRealVar("recoQQsign","qq sign",-1,3,"");
  RooRealVar* ctau3DVar    = new RooRealVar("ctau3D","c_{#tau}", -100000.0, 100000.0, "mm");
  RooRealVar* ctau3DErrVar = new RooRealVar("ctau3DErr","#sigma_{c#tau}", -100000.0, 100000.0, "mm");
  RooRealVar* ctau3DResVar = new RooRealVar("ctau3DRes","c_{#tau}", -100000.0, 100000.0, "");
  RooRealVar* ctau3D2SVar    = new RooRealVar("ctau3D2S","c_{#tau}", -100000.0, 100000.0, "mm");
  RooRealVar* ctau3DErr2SVar = new RooRealVar("ctau3DErr2S","#sigma_{c#tau}", -100000.0, 100000.0, "mm");
  RooRealVar* ctau3DRes2SVar = new RooRealVar("ctau3DRes2S","c_{#tau}", -100000.0, 100000.0, "");
  RooRealVar* NumDimu = new RooRealVar("NumDimu","number of dimuon",0,100,"");
  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *eta1Var, *eta2Var,*evtWeight);
  argSet->add(*cBinVar); argSet->add(*recoQQ); argSet->add(*NumDimu); argSet->add(*ctau3DVar); argSet->add(*ctau3DErrVar); argSet->add(*ctau3DResVar);
  argSet->add(*ctau3D2SVar); argSet->add(*ctau3DErr2SVar); argSet->add(*ctau3DRes2SVar);
  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);

  ////////////////////////////////////////////////////////////////////////
  //////////////////  TreeDataSet
  ////////////////////////////////////////////////////////////////////////
  TString fCentSelHF = "HFNom";
  if(hiHFBinEdge==1) fCentSelHF = "HFUp";
  else if(hiHFBinEdge==-1) fCentSelHF = "HFDown";

  TFile* newfile;
  newfile = new TFile(Form("./skimmedFiles/OniaFlowSkim_isMC%d_%s_230119.root", isMC,fCentSelHF.Data()),"recreate");

  const static int nMaxDimu_ = 1000;
  int nDimu_;
  float mass_[nMaxDimu_];
  float ctau3D_[nMaxDimu_];
  float ctau3DErr_[nMaxDimu_];

  TTree* mmevttree = new TTree("mmepevt","dimuonAndEventPlanes in event based");
  mmevttree->SetMaxTreeSize(MAXTREESIZE);
  mmevttree->Branch("nDimu_",&nDimu_,"nDimu_/I");
  mmevttree->Branch("mass_",mass_,"mass_[nDimu_]/F");
  mmevttree->Branch("ctau3D_",ctau3D_,"ctau3D_[nDimu_]/F");
  mmevttree->Branch("ctau3DErr_",ctau3DErr_,"ctau3DErr_[nDimu_]/F");

  double weight_acc = 1;
  double weight_eff = 1;

  Int_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;

  int nDimu_all=0;
  int nDimuPass=0;
  int nDimu_one=0;
  int nDimu_more=0;
  //Begin Loop
    for(int i=0; i<nEvt; i++){
      tree->GetEntry(i);
      if(fabs(vz) > 15) continue;
      nDimuPass=0;
      //if(i==1000000) break;
      if(cBin>=cLow&&cBin<cHigh){ 
        for(int j=0; j<nDimu; j++){
          //cout<<"Evt: "<<i<<", Cent: "<<cBin<<", mass: "<<mass[j]<<", pt: "<<pt[j]<<", pt1:"<<pt1[j]<<", pt2: "<<pt2[j]<<", eta1: "<<eta1[j]<<", eta2: "<<eta2[j]<<", vz: "<<vz<<endl;
          if(! ((double)pt[j]<50 && recoQQsign[j]==0 && abs((double)y[j])<2.4
                && IsAcceptanceQQ(pt1[j],eta1[j]) && IsAcceptanceQQ(pt2[j],eta2[j])) ) continue;
          nDimuPass++;
        }
        nDimu_all++;
        if(nDimuPass>1) {nDimu_more++; continue;}
        if(nDimuPass==1) nDimu_one++;
        // Fill Dimuon Loop
        nDimu_=0;
        for(int j=0; j<nDimu; j++){
          if((double)pt[j]<50&&recoQQsign[j]==0&&mass[j]>massLow&&mass[j]<massHigh&&abs(y[j])<2.4
              && IsAcceptanceQQ(pt1[j],eta1[j]) && IsAcceptanceQQ(pt2[j],eta2[j])){
            weight_acc=1;
            weight_eff=1;
		//	cout << "pt : " << pt[j] << " y : " << y[j] << endl;
            if(fAccW){
              if ( abs((double)y[j])<1.6 ) {weight_acc = getAccWeight(hAccPt[1], pt[j]);}
              else if ( abs((double)y[j])>1.6 && abs((double)y[j])<2.4 ) { weight_acc = getAccWeight(hAccPt[2], pt[j]);}
            }
            if(fEffW){
              if ( abs((double)y[j])>=0.0 && abs((double)y[j])<1.2 ) { weight_eff = getEffWeight(hEffPt[0], pt[j]);  }
              else if ( abs((double)y[j])>=1.2 && abs((double)y[j])<1.6 ) { weight_eff = getEffWeight(hEffPt[1], pt[j]); }
              else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<1.8 ) { weight_eff = getEffWeight(hEffPt[2], pt[j]); }
              else if ( abs((double)y[j])>=1.8 && abs((double)y[j])<2.4 ) { weight_eff = getEffWeight(hEffPt[3], pt[j]); }
            }

            double weight_ = weight * weight_eff * weight_acc;

            recoQQ->setVal((int)recoQQsign[j]);
            massVar->setVal( (double)mass[j] );
            ptVar->setVal(   (double)pt[j]   );
            yVar->setVal(    (double)y[j]    );
            pt1Var->setVal(  (double)pt1[j]  );
            eta1Var->setVal( (double)eta1[j] );
            pt2Var->setVal(  (double)pt2[j]  );
            eta2Var->setVal( (double)eta2[j] );
            cBinVar->setVal( (double)cBin ) ;
            ctau3DVar->setVal( (double)ctau3D[j] ) ;
            ctau3DErrVar->setVal( (double)ctau3DErr[j] ) ;
            ctau3DResVar->setVal( (double)ctau3D[j]/ctau3DErr[j] ) ;
            ctau3D2SVar->setVal( (double)ctau3D[j] ) ;
            ctau3DErr2SVar->setVal( (double)ctau3DErr[j] ) ;
            ctau3DRes2SVar->setVal( (double)ctau3D2S[j] ) ;
            evtWeight->setVal( (double)weight_ ) ;
            NumDimu->setVal((int)nDimu);
            //cout<<"Evt: "<<j<<", Cent: "<<cBin<<", mass: "<<mass[j]<<", pt: "<<pt[j]<<", pt1:"<<pt1[j]<<", pt2: "<<pt2[j]<<", eta1: "<<eta1[j]<<", eta2: "<<eta2[j]<<", vz: "<<vz<<endl;
            dataSet->add( *argSet);

            mass_[j] = mass[j];
            ctau3D_[j] = ctau3D[j];
            ctau3DErr_[j] = ctau3DErr[j];
            nDimu_++;
          }
          if(nDimu_>0) mmevttree->Fill();
        }
      }
    }

    newfile->cd();
    mmevttree->Write();
    newfile->Close();

    TFile *wf = new TFile(Form("skimmedFiles/OniaRooDataSet_miniAOD_isMC%d_Psi2S_cent%i_%i_Effw%d_Accw%d_PtW%d_TnP%d_230513.root", 
          isMC,cLow,cHigh,fEffW,fAccW,isPtW,isTnP),"recreate");  wf->cd();
    dataSet->Write();

    cout<<"How many Jpsi??: "<<nDimu_one<<endl;


  }

  void GetHistSqrt(TH1D* h1, TH1D* h2){
    if(h1->GetNbinsX() != h2->GetNbinsX()){ cout << "Inconsistent # of bins b/w histograms !! " << endl;}
    double content;
    double err;
    for(int i=1; i<=h1->GetNbinsX(); i++){
      content=0;err=0;
      content = h1->GetBinContent(i);
      err = h1->GetBinError(i);
      err = 0.5*err*TMath::Power(content,-0.5);
      h2->SetBinContent(i,TMath::Sqrt(content));
      h2->SetBinError(i,err);
    }
  } 

  void Get2DHistSqrt(TH2D* h1, TH2D* h2){
    if(h1->GetNbinsX() != h2->GetNbinsX()){ cout << "Inconsistent # of bins b/w histograms !! " << endl;}
    double content;
    double err;
    for(int i=1; i<=h1->GetNbinsX(); i++){
      content=0;err=0;
      content = h1->GetBinContent(i);
      err = h1->GetBinError(i);
      err = 0.5*err*TMath::Power(content,-0.5);
      h2->SetBinContent(i,TMath::Sqrt(content));
      h2->SetBinError(i,err);
    }
  } 

  void GetHistBkg(TH1D* h1, TH1D* h2){
    double Nh1;
    double Nh2;
    double Nh12;
    Nh1 = h1->GetEntries();
    Nh2 = h2->GetEntries();
    Nh12 = Nh1+Nh2;
    return Nh1/(Nh1+Nh2);

  }

  double getAccWeight(TH1D* h, double pt){
    //cout<<"pt: "<<pt<<endl;
    double binN = h->FindBin(pt);
    double weight_ = 1./(h->GetBinContent(binN));
    return weight_;
  }

  double getEffWeight(TH1D *h, double pt){
    double binN = h->FindBin(pt);
    double weight_ = 1./(h->GetBinContent(binN));
    return weight_;
  }
