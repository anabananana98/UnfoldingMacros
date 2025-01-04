#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <iterator>
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include "TPaveLabel.h"
#include "TColor.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include <TRandom3.h>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "TRandom.h"

#endif

Double_t delR(double phi1, double phi2, double eta1, double eta2)
{
    double dphi = abs(phi1-phi2);
    if(dphi>TMath::Pi()){dphi = (2.*TMath::Pi() - dphi);}
    
    double deta = std::fabs(eta1 - eta2);
    Double_t delR = std::sqrt(dphi*dphi + deta*deta);
    
    return delR;
}

// Template function to check for NaNs or Infs in 1D, 2D, and 3D histograms
template <typename TH>
bool hasNaNsOrInfs(TH* hist) {
    if (!hist) {
        std::cerr << "Histogram is null!" << std::endl;
        return true;
    }
    
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for (int j = 1; j <= hist->GetNbinsY(); ++j) {
            for (int k = 1; k <= hist->GetNbinsZ(); ++k) {
                if (std::isnan(hist->GetBinContent(i, j, k)) || std::isinf(hist->GetBinContent(i, j, k))) {
                    return true;
                }
            }
        }
    }
    return false;
}

// Specialization for 2D histograms
template <>
bool hasNaNsOrInfs<TH2>(TH2* hist) {
    if (!hist) {
        std::cerr << "Histogram is null!" << std::endl;
        return true;
    }

    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for (int j = 1; j <= hist->GetNbinsY(); ++j) {
            if (std::isnan(hist->GetBinContent(i, j)) || std::isinf(hist->GetBinContent(i, j))) {
                return true;
            }
        }
    }
    return false;
}

// Specialization for 1D histograms
template <>
bool hasNaNsOrInfs<TH1>(TH1* hist) {
    if (!hist) {
        std::cerr << "Histogram is null!" << std::endl;
        return true;
    }

    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        if (std::isnan(hist->GetBinContent(i)) || std::isinf(hist->GetBinContent(i))) {
            return true;
        }
    }
    return false;
}


//For fastsim
double paireffwt(TH2D* h2, TH2D* h2_tru, double chpt, double dRval) {
    // Find the bin corresponding to chpt on the Y-axis
    int ybin = h2->GetYaxis()->FindBin(chpt);
    
    // Create projections of the 2D histograms on the X-axis
    TH1D* h1_det = h2->ProjectionX("h1_det", ybin, ybin);
    TH1D* h1_tru = h2_tru->ProjectionX("h1_tru", ybin, ybin);

    // Normalize the histograms
    h1_det->Scale(1.0 / h1_det->Integral(), "width");
    h1_tru->Scale(1.0 / h1_tru->Integral(), "width");
    
    // Check for empty histograms (avoid division by zero)
    if (h1_tru->Integral() == 0 || h1_det->Integral() == 0) {
        std::cerr << "Warning: One of the histograms is empty!" << std::endl;
        return 0.0;
    }

    // Divide the histograms (det / tru)
    h1_det->Divide(h1_tru);
    
    // Find the bin corresponding to dRval in the X-axis of h1_det
    double paireffwt = h1_det->GetBinContent(h1_det->GetXaxis()->FindBin(dRval));
    
    // Return the weight
    return paireffwt;
}


void Unfolding3dAnalysisOptNew(std::string tag = "", int pthardbin = -1, std::string fileName ="") {
    
    // tag = "rw";

    //***************************************************
    TH1D* h1_true; TH1D* h1_trueSplit;
    TH1D* h1_reco;
    TH1D* h1_unfolded; TH1D* h1_unfolded_rw;  TH1D* h1_unfoldedSplit;
    TH1D* h1_raw; TH1D* h1_rawSplit;
    TH1D* h1_smeared;
    TH1D* h1_fulleff;
    TH1D* h1_eff_match;
    TH1D* h1_fullreco;
    TH1D* h1_rw;
    TH1* h1_fold;
    TH1* h1_triv;
    
    TH3D* h3_true_eec; TH3D* h3_trueSplit;
    TH3D* h3_reco_eec;
    TH3D* h3_true_eecTriv;
    TH3D* h3_reco_eecTriv;
    TH3D* h3_unfolded; TH3D* h3_unfolded_rw;  TH3D* h3_unfoldedSplit;
    TH3D* h3_raw; TH3D* h3_rawSplit;
    TH3D* h3_smeared;
    TH3D* h3_fulleff;
    TH3D* h3_eff_match;
    TH3D* h3_fullreco;
    TH3D* h3_rw;
    TH3* h3_fold;
    TH3* h3_triv;
    
    TH2D* eec_pt_hist;
    TH2D* eec_pt_hist_det;
    TH2D* eec_pt_histTriv;
    TH2D* eec_pt_hist_detTriv;
    TH2D* eec_pt_histSplit;
    TH2D* eec_pt_hist_detSplit;
    
    TFile *fmc;

    TRandom3* rand = new TRandom3(0);
    
    RooUnfoldResponse response1D;
    RooUnfoldResponse response1DSplit;
    
    RooUnfoldResponse response3D;
    RooUnfoldResponse response3DTriv;
    RooUnfoldResponse response3DSplit;
   
    float pThard_val[21] = {5., 7., 9., 12., 16., 21., 28., 36., 45., 57., 70., 85., 99., 115., 132., 150., 169., 190., 212., 235., 1000.};
    
   
    //   //reco jet pT bins for 1D unfolding 
    // Double_t from_const_reco_jet = 20;
    // Double_t to_const_reco_jet = 80;
    // Int_t bins_const_reco_jet = 3;
    // Double_t width_const_reco_jet = (to_const_reco_jet-from_const_reco_jet)/bins_const_reco_jet;
    // Double_t xbins_jet[5] = {};
    // for (Int_t i = 0; i <= bins_const_reco_jet; i++)
    // {
    //     xbins_jet[i] = (from_const_reco_jet + i * width_const_reco_jet);
    // }
    // xbins_jet[4] = 120;
    
    
         //reco jet pT bins for 1D unfolding -- NOV 20 TRIAL
    Double_t from_const_reco_jet = 0;
    Double_t to_const_reco_jet = 80;
    Int_t bins_const_reco_jet = 4;
    Double_t width_const_reco_jet = (to_const_reco_jet-from_const_reco_jet)/bins_const_reco_jet;
    Double_t xbins_jet[6] = {};
    for (Int_t i = 0; i <= bins_const_reco_jet; i++)
    {
        xbins_jet[i] = (from_const_reco_jet + i * width_const_reco_jet);
    }
    xbins_jet[5] = 120;

    
    
      //true jet pT bins for 1D unfolding 
    Double_t from_const_jet = 0;
    Double_t to_const_jet= 80;
    Int_t bins_const_jet = 4;
    Double_t width_const_jet = (to_const_jet-from_const_jet)/bins_const_jet;
    Double_t tbins_jet[8] = {};
    for (Int_t i = 0; i <= bins_const_jet; i++)
    {
        tbins_jet[i] = (from_const_jet + i * width_const_jet);
    }
    tbins_jet[5]=120;
    tbins_jet[6]=160;
    tbins_jet[7]=200;
    
    
 
    // //true jet pT bins for 3D unfolding 
    // Double_t from_const = 0;
    // Double_t to_const = 200;
    // Int_t bins_const = 20;
    // //  Double_t to_const = 350;
    // // Int_t bins_const = 35;
    // Double_t width_const = (to_const-from_const)/bins_const;
    // Double_t tbins[21] = {};
    // for (Int_t i = 0; i <= bins_const; i++)
    // {
    //     tbins[i] = (from_const + i * width_const);
    // }
    
    //  //jet pT bins for 3D unfolding 
    // Double_t from_const_reco = 10;
    // Double_t to_const_reco = 120;
    // Int_t bins_const_reco = 11;
    // Double_t width_const_reco = (to_const_reco-from_const_reco)/bins_const_reco;
    // Double_t xbins[12] = {};
    // for (Int_t i = 0; i <= bins_const_reco; i++)
    // {
    //     xbins[i] = (from_const_reco + i * width_const_reco);
    // }
    
    
    Double_t from_const = 0;
    Double_t to_const= 80;
    Int_t bins_const = 4;
    Double_t width_const = (to_const-from_const)/bins_const;
    Double_t tbins[8] = {};
    for (Int_t i = 0; i <= bins_const; i++)
    {
        tbins[i] = (from_const + i * width_const_jet);
    }
    tbins[5]=120;
    tbins[6]=160;
    tbins[7]=200;
    
    ////Main one 
    // Double_t from_const_reco = 20;
    // Double_t to_const_reco = 80;
    // Int_t bins_const_reco = 3;
    // Double_t width_const_reco = (to_const_reco-from_const_reco)/bins_const_reco;
    // Double_t xbins[5] = {};
    // for (Int_t i = 0; i <= bins_const_reco; i++)
    // {
    //     xbins[i] = (from_const_reco + i * width_const_reco);
    // }
    // xbins[4] = 120;
    
      //jet pT bins for 3D unfolding -- NOV 20 TRIAL
    Double_t from_const_reco = 0;
    Double_t to_const_reco = 80;
    Int_t bins_const_reco = 4;
    Double_t width_const_reco = (to_const_reco-from_const_reco)/bins_const_reco;
    Double_t xbins[6] = {};
    for (Int_t i = 0; i <= bins_const_reco; i++)
    {
        xbins[i] = (from_const_reco + i * width_const_reco);
    }
    xbins[5] = 120;
    
    
    
    
    
    
    // reco dR bins --same for eec and e3c
    // double dRbins[] = {
    //     0.0001,
    //     0.01,
    //     0.0158489, 0.0190546, 0.0229087,
    //     0.0275423, 0.0331131,
    //     0.0398107, 0.0436516, 0.047863, 0.0524807, 0.057544,
    //     0.0630957, 0.0691831, 0.0758578, 0.0831764, 0.0912011,
    //     0.1, 0.109648, 0.120226, 0.131826, 0.144544,
    //     0.158489, 0.17378, 0.190546, 0.20893, 0.229087,
    //     0.251189, 0.275423, 0.301995, 0.331131, 0.363078,
    //     0.398107, 0.524807,
    //     1.0
    // };
    
//       //Nov3 trial
   double dRbins[] = {0.0001,0.01,0.0144613,0.0173904,0.0209128,0.0251487,0.0302425,
 0.0363681,0.0437345,0.0525929,0.0632456,0.0760559,0.091461,0.109986,0.132264,0.159054,0.191271,0.230012,0.276601,0.332627,0.4,0.524807,1};
    
    
    //  double dRbins[] = {
    //     0.0001,
    //     0.01,
    //     0.0158489, 0.0190546, 0.0229087,
    //     0.0275423, 0.0331131,
    //     0.0398107, 0.0436516, 0.047863, 0.0524807, 0.057544,
    //     0.0630957, 0.0691831, 0.0758578, 0.0831764, 0.0912011,
    //     0.1, 0.109648, 0.120226, 0.131826, 0.144544,
    //     0.158489, 0.17378, 0.190546, 0.20893, 0.229087,
    //     0.251189, 0.275423, 0.301995, 0.331131, 0.363078, 0.398107, 0.524807,
    //     1.0
    // };
    
    // //true dR bins --same for eec and e3c
    // double dRbinst[] = {
    //     0.0001,
    //     0.01,
    //     0.0158489, 0.0190546, 0.0229087,
    //     0.0275423, 0.0331131,
    //     0.0398107, 0.0436516, 0.047863, 0.0524807, 0.057544,
    //     0.0630957, 0.0691831, 0.0758578, 0.0831764, 0.0912011,
    //     0.1, 0.109648, 0.120226, 0.131826, 0.144544,
    //     0.158489, 0.17378, 0.190546, 0.20893, 0.229087,
    //     0.251189, 0.275423, 0.301995, 0.331131, 0.363078,
    //     0.398107, 0.524807,
    //     1.0
    // };

     double dRbinst[] = {0.0001,0.01,0.0144613,0.0173904,0.0209128,0.0251487,0.0302425,
 0.0363681,0.0437345,0.0525929,0.0632456,0.0760559,0.091461,0.109986,0.132264,0.159054,0.191271,0.230012,0.276601,0.332627,0.4,0.524807,1};
    
    
    // //reco wt bins eec
    // double wtbins[] = {
    //     0.001,0.00346737,0.00912011, 0.0138038,0.0165959,0.0190546,
    //     0.0288403,
    //     0.0346737,
    //     0.0457088, 0.060256,
    //     0.0794328,
    //     0.104713, 0.138038,
    //     0.18197, 0.239883,
    //     0.316228, 0.363078,
    //     0.416869, 0.5};
    
    // //true wt bins eec
    // double wtbinst[] = {
    //     0.001,0.00346737,0.00912011, 0.0138038,0.0165959,0.0190546,
    //     0.0288403,
    //     0.0346737,
    //     0.0457088, 0.060256,
    //     0.0794328,
    //     0.104713, 0.138038,
    //     0.18197, 0.239883,
    //     0.316228, 0.363078,
    //     0.416869, 0.5};
        
        //reco wt bins eec
    double wtbins[] = {
        0.0005,0.001,0.00346737,0.00912011, 0.0138038,0.0165959,0.0190546,
        0.0288403,
        0.0346737,
        0.0457088, 0.060256,
        0.0794328,
        0.104713,0.138038, 0.18197,0.26};
    
    //true wt bins eec
    double wtbinst[] = {
        0.0001,0.0005,0.001,0.00346737,0.00912011, 0.0138038,0.0165959,0.0190546,
        0.0288403,
        0.0346737,
        0.0457088, 0.060256,
        0.0794328,
        0.104713,0.138038, 0.18197,0.26
        };
        
   
    
    
    // //reco wt bins eec
    // double wtbins[] = {
    //     0.0001, 0.00301995, 0.00501187, 0.01, 0.0158489, 0.020893, 0.0301995, 0.0501187, 0.060256, 0.1,
    //     0.151356, 0.758578,1.0};
    
    // //true wt bins eec
    // double wtbinst[] = {
    //     0.0001, 0.00301995, 0.00501187, 0.01, 0.0158489, 0.020893, 0.0301995, 0.0501187, 0.060256, 0.1,
    //     0.151356, 0.758578,1.0};
    
    // Number of bins is one less than the number of edges
    int nWtbins = sizeof(wtbins) / sizeof(wtbins[0]) - 1;
    int ndRbins = sizeof(dRbins) / sizeof(dRbins[0]) - 1;
    int nJetPtbins = sizeof(xbins) / sizeof(xbins[0]) - 1;
     
    int nWtbinst = sizeof(wtbinst) / sizeof(wtbinst[0]) - 1;
    int ndRbinst = sizeof(dRbinst) / sizeof(dRbinst[0]) - 1;
    int nJetPtbinst = sizeof(tbins) / sizeof(tbins[0]) - 1;
    
    int nJetPtbins1Dreco = sizeof(xbins_jet) / sizeof(xbins_jet[0]) - 1;
    int nJetPtbins1Dtrue = sizeof(tbins_jet) / sizeof(tbins_jet[0]) - 1;
    
    
    eec_pt_hist = new TH2D("eec_pt_hist", "EEC and jet_pt 2D",  ndRbinst, dRbinst, nJetPtbinst, tbins);
    eec_pt_hist_det = new TH2D("eec_pt_hist_det", "EEC and jet_pt 2D det", ndRbins, dRbins, nJetPtbins,xbins);
    
    eec_pt_histTriv = new TH2D("eec_pt_histTriv", "EEC and jet_pt 2D Triv",  ndRbinst, dRbinst, nJetPtbinst, tbins);
    eec_pt_hist_detTriv = new TH2D("eec_pt_hist_detTriv", "EEC and jet_pt 2D det Triv", ndRbins, dRbins, nJetPtbins,xbins);
    
    eec_pt_histSplit = new TH2D("eec_pt_histSplit", "EEC and jet_pt 2D Split",  ndRbinst, dRbinst, nJetPtbinst, tbins);
    eec_pt_hist_detSplit = new TH2D("eec_pt_hist_detSplit", "EEC and jet_pt 2D det Split", ndRbins, dRbins, nJetPtbins,xbins);
    
    
    h1_true = new TH1D("h1_true", "h1_true", nJetPtbins1Dtrue, tbins_jet);
    h1_fulleff = new TH1D("h1_fulleff", "h1_fulleff", nJetPtbins1Dtrue, tbins_jet);
    h1_reco = new TH1D("h1_reco", "h1_reco", nJetPtbins1Dreco, xbins_jet);

    
   cout<<"nWtbinst "<<nWtbinst<<" and nWtbins "<<nWtbins<<endl;
   cout<<"ndRbinst "<<ndRbinst<<" and ndRbins "<< ndRbins <<endl;
   cout<<"nJetPtbinst "<<nJetPtbinst<<" and nJetPtbins "<<nJetPtbins<<endl;    
   
   cout<<"nJetPtbins1Dtrue "<<nJetPtbins1Dtrue<<" and nJetPtbins1Dreco "<<nJetPtbins1Dreco<<endl;
   
  
    
    h3_true_eec = new TH3D("h3_true_eec", "h3_true_eec", ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_reco_eec = new TH3D("h3_reco_eec", "h3_reco_eec",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
    h3_smeared = new TH3D("h3_smeared", "h3_smeared",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
    
    h3_true_eecTriv = new TH3D("h3_true_eecTriv", "h3_true_eec_triv", ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_reco_eecTriv = new TH3D("h3_reco_eecTriv", "h3_reco_eec_triv",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
    
    h3_trueSplit = new TH3D("h3_trueSplit", "h3_true_split", ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_rawSplit = new TH3D("h3_rawSplit", "h3_reco_split",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
   
    h3_fulleff = new TH3D("h3_fulleff","h3_fulleff",ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_eff_match = new TH3D("h3_eff_match","h3_eff_match",ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_fullreco = new TH3D("h3_fullreco","h3_fullreco",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
    
    
    h3_true_eec->Sumw2();
    h3_reco_eec->Sumw2();
    
    h3_true_eecTriv->Sumw2();
    h3_reco_eecTriv->Sumw2();
    
    h3_rawSplit->Sumw2();
    h3_trueSplit->Sumw2();

    cout<<h3_rawSplit->GetNbinsX()<<endl;
  
    h3_eff_match->Sumw2();
    h3_fulleff->Sumw2();
    h3_fullreco->Sumw2();
  
   cout<<"HERE"<<endl;
  
    cout<<"setting up response declarations"<<endl;
    response3D.Setup(h3_reco_eec,h3_true_eec);
    response3DTriv.Setup(h3_reco_eecTriv,h3_true_eecTriv);
    response3DSplit.Setup(h3_rawSplit,h3_trueSplit);
    response1D.Setup(h1_reco,h1_true);

    cout<<"ending declarations"<<endl;
    

      //Get the mc files
    fmc = TFile::Open(Form("%s",fileName.c_str()));
    TString file_name = Form("%s",fileName.c_str());

    cout<<"file name is: "<<file_name.Data()<<endl;
    
    TH1D *fHistEventsAfterSel = (TH1D*)fmc->Get(TString::Format("fHistEventsAfterSel"));
    Int_t nTrials = fHistEventsAfterSel->GetEntries();
    TProfile *fHistXsection = (TProfile*)fmc->Get(TString::Format("fHistXsection"));
    Double_t xSection = fHistXsection->GetBinContent(pthardbin+1); // we need i+1 to get the correct bin
    double mc = xSection/nTrials;
    //
    std::cout << "# of Trials: " << nTrials << "\n";
    std::cout << "X section: " << xSection << "\n";
    std::cout << "weight: " << mc << "\n\n\n"<<endl;
    
    //get the histogram for reweighting the prior
  TFile* f_rw = new TFile("/home/ar2545/Trial/RwHistEEC.root");
  h3_rw = (TH3D*)f_rw->Get("hResult");
     
 

    TTree *tree=(TTree*)fmc->Get("MatchTracksTree");
     //In this tree I only have matched jets, so tree entries where either jet pT is 0 is only for missed(true level) or fake tracks(det level)
     Long64_t nEntries = tree->GetEntries();
    cout<<nEntries<<endl;
    
    if (nEntries != 0) {
    tree->Show(0);
  
    Int_t nEv=tree->GetEntries();
   
    Double_t fJet_pt_det, fJet_pt_tru;
    Double_t fTrack_pt_det, fTrack_pt_tru, fTrack_pt_miss, fTrack_pt_fake;
    Double_t fTrack_eta_det, fTrack_eta_tru, fTrack_eta_miss, fTrack_eta_fake;
    Double_t fTrack_phi_det, fTrack_phi_tru, fTrack_phi_miss, fTrack_phi_fake;
    Double_t mc_weight;
    Double_t fTrack_choverpt_det;
    Double_t fTrack_choverpt_tru;
    Double_t fTrack_choverpt_miss;
    Double_t fTrack_choverpt_fake;
    
    tree->SetBranchAddress("fJet_pt_det", &fJet_pt_det);
    tree->SetBranchAddress("fJet_pt_tru", &fJet_pt_tru);
    tree->SetBranchAddress("fTrack_pt_det", &fTrack_pt_det);
    tree->SetBranchAddress("fTrack_pt_tru", &fTrack_pt_tru);
    tree->SetBranchAddress("fTrack_pt_miss", &fTrack_pt_miss);
    tree->SetBranchAddress("fTrack_pt_fake", &fTrack_pt_fake);
    tree->SetBranchAddress("fTrack_eta_det", &fTrack_eta_det);
    tree->SetBranchAddress("fTrack_eta_tru", &fTrack_eta_tru);
    tree->SetBranchAddress("fTrack_eta_miss", &fTrack_eta_miss);
    tree->SetBranchAddress("fTrack_eta_fake", &fTrack_eta_fake);
    tree->SetBranchAddress("fTrack_phi_det", &fTrack_phi_det);
    tree->SetBranchAddress("fTrack_phi_tru", &fTrack_phi_tru);
    tree->SetBranchAddress("fTrack_phi_miss", &fTrack_phi_miss);
    tree->SetBranchAddress("fTrack_phi_fake", &fTrack_phi_fake);

    tree->Branch("fTrack_choverpt_det", &fTrack_choverpt_det,"fTrack_choverpt_det/D");
    tree->Branch("fTrack_choverpt_tru", &fTrack_choverpt_tru,"fTrack_choverpt_tru/D");
    tree->Branch("fTrack_choverpt_miss", &fTrack_choverpt_miss,"fTrack_choverpt_miss/D");
    tree->Branch("fTrack_choverpt_fake", &fTrack_choverpt_fake,"fTrack_choverpt_fake/D");
    // tree->GetEntry(1);
    // double mc = mc_weight; //same for one run in 1 pTHat bin
    // double weightresp = mc;
    
    // double ptJetmin_det = 10.; //min det level pt -- NOV 20 TRIAL
     double ptJetmin_det = 15.; //min det level pt -- NOV 22 TRIAL
    // double ptJetmin_det = 20.; //min det level pt
    double ptJetmax_det = 120.; //max det level pt
    double ptJetmin_tru3D = 0.; //min tru level pt for 3D unfolding
    double ptJetmax_tru3D = 200; //max tru level pt for 3D unfolding
    
    
    //for 1D unfolding
    double ptJetmin_tru = 0.; //min tru level pt
    double ptJetmax_tru = 200.; //max tru level pt
    
    
    double jetpT_tru = -1;
    double jetpT_tru_fulleff = -1;
    double jetpT_det = -1;
    
   

    double pt_tru;
    double pt_det;
    double pt_tru2;
    double pt_det2;
    double pt_miss;
    double pt_fake;

    double eta_tru;
    double eta_det;
    double eta_tru2;
    double eta_det2;
    double eta_miss;
    double eta_fake;
    

    double phi_tru;
    double phi_det;
    double phi_tru2;
    double phi_det2;
    double phi_miss;
    double phi_fake;

    bool mi = false; bool mk = false; bool fi = false; bool fk = false;
    
    bool ifData = true;
    bool ifPairCut = true;
    bool ifTrivial = false;
    bool ifSplit = false;
    bool ifFastSim = false;
    

if(ifData && !ifPairCut){
  for (Int_t i = 0; i < nEntries; ++i) {
    
        tree->GetEntry(i);
        
    // Skip if either pT is zero (fake or missed track)
    if (fJet_pt_det == 0 || fJet_pt_tru == 0) continue;
    
    
        // Fill full efficiency histogram if the true jet pT has changed
    if (jetpT_tru_fulleff != fJet_pt_tru) {
        h1_fulleff->Fill(fJet_pt_tru, mc);
    }
    jetpT_tru_fulleff = fJet_pt_tru;
    
   
    if (fJet_pt_tru > ptJetmax_tru || fJet_pt_tru < 0) continue;
    // Fill histograms for unique jets
    if (jetpT_det != fJet_pt_det && jetpT_tru != fJet_pt_tru) {
        if (fJet_pt_det > ptJetmin_det && fJet_pt_det < ptJetmax_det) 
        {
            h1_true->Fill(fJet_pt_tru, mc);
            h1_reco->Fill(fJet_pt_det, mc);
            response1D.Fill(fJet_pt_det, fJet_pt_tru, mc);
            
        }
    }
    
     jetpT_tru = fJet_pt_tru;
     jetpT_det = fJet_pt_det;
  }

    jetpT_tru = -1;
    jetpT_det = -1;
    

  for (Int_t i = 0; i < nEntries; ++i) {
    
        tree->GetEntry(i);
        
        jetpT_tru = fJet_pt_tru;
        jetpT_det = fJet_pt_det;
        
    // Apply really wide cuts on true jet pT///THIS LINE IS NEW -- CHECK THIS WITH LAURA
       if (fJet_pt_tru > (ptJetmax_tru) || fJet_pt_tru < 0) continue;
     
        if(fTrack_pt_tru!=0){pt_tru = fTrack_pt_tru; mi = false;}
        else{pt_tru = fTrack_pt_miss; mi = true;}
        
        if(fTrack_pt_det!=0){pt_det = fTrack_pt_det; fi = false;}
        else{pt_det = fTrack_pt_fake; fi = true;}
        
        if(fTrack_eta_tru!=0){eta_tru = fTrack_eta_tru;}
        else{eta_tru = fTrack_eta_miss;}
        
        if(fTrack_eta_det!=0){eta_det = fTrack_eta_det;}
        else{eta_det = fTrack_eta_fake;}
        
        if(fTrack_phi_tru!=0){phi_tru = fTrack_phi_tru;}
        else{phi_tru = fTrack_phi_miss;}
        
        if(fTrack_phi_det!=0){phi_det = fTrack_phi_det;}
        else{phi_det = fTrack_phi_fake;}
       
    //   cout<<i<<" det: "<<jetpT_det<<" and true "<<jetpT_tru << endl;
        //     //2-pt corr
            for (Int_t k = i+1; k < nEntries; ++k) {
                tree->GetEntry(k);
                
                if(fTrack_pt_tru!=0){pt_tru2 = fTrack_pt_tru; mk = false;}
                else{pt_tru2 = fTrack_pt_miss; mk = true;}
        
                if(fTrack_pt_det!=0){pt_det2 = fTrack_pt_det; fk = false;}
                else{pt_det2 = fTrack_pt_fake; fk = true;}
        
                if(fTrack_eta_tru!=0){eta_tru2 = fTrack_eta_tru;}
                else{eta_tru2 = fTrack_eta_miss;}
        
                if(fTrack_eta_det!=0){eta_det2 = fTrack_eta_det;}
                else{eta_det2 = fTrack_eta_fake;}
         
                if(fTrack_phi_tru!=0){phi_tru2 = fTrack_phi_tru;}
                else{phi_tru2 = fTrack_phi_miss;}
        
                if(fTrack_phi_det!=0){phi_det2 = fTrack_phi_det;}
                else{phi_det2 = fTrack_phi_fake;}
               
            //   cout<<i<<" "<<k <<" jetpT_det "<< jetpT_det<< " and fJet_pt_det "<<fJet_pt_det<<endl;
            //   cout<<i<<" "<<k<< " jetpT_tru "<<jetpT_tru << " and fJet_pt_tru "<<fJet_pt_tru<<endl;
            
                //if both true and detector jet pT change then we have moved on to a new jet
                //because if one of them changes then that just means that its a missed or a fake track in the same jet
                if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                break;} 
                
                
                double dR_tru = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                double dR_det = delR(phi_det,phi_det2,eta_det,eta_det2);
                
                //fake correlations//////////////////////////////////
                if(fi || fk){ 
                //     cout<<"fake: "<<endl;
                //  tree->Show(k);
                 double wt_det =(2*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                 double wt_det_3D = (pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                 if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                 if(dR_det < 0.01 || dR_det > 0.4) continue;
                 if(wt_det_3D < 0.0005 || wt_det_3D > 0.26 ) continue;
                 h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                 h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                 
                 eec_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det*mc);
                }
                //missed correlations//////////////////////////////////
                else if(mi || mk){ 
                //  cout<<"miss: "<<endl;
                // tree->Show(k);
                double wt_tru= (2*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru);
                double wt_tru_3D = (pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru); //for filling 3D histogram
                // if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue; //only wide cuts needed for full eff
                if(dR_tru < 0.005 || dR_tru >0.4) continue;
                if(wt_tru_3D < 0.0001 || wt_tru_3D > 0.26) continue;
                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                
                eec_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru*mc);
                }
                //matched correlations//////////////////////////////////
                else if(!mi && !fi && !mk && !fk){
                    // cout<<"match: "<<endl; 
                    // tree->Show(k); 
                double wt_det =(2*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                double wt_det_3D = (pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                
                double wt_tru= (2*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru);
                double wt_tru_3D = (pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru); //for filling 3D histogram
                
                //3D KE histogram TRUTH level -- applied before any truth or reco level cuts //only bin edge cuts and super wide cuts for jet pT
                if(dR_tru > 0.01 && dR_tru < 0.4 && wt_tru_3D > 0.0001 && wt_tru_3D < 0.26){
                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                
                eec_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru*mc); //because this is how i would fill it otherwise.
                } 
                
                
                //3D KE histogram DET level -- applied before any truth level cuts 
                 if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                 if(dR_det < 0.01 || dR_det > 0.4) continue;
                 if(wt_det_3D < 0.0005 || wt_det_3D > 0.26 ) continue;
                 
                 h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                 h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
             
                // eec_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru*mc);
                eec_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det*mc);
                if (hasNaNsOrInfs(eec_pt_hist) || hasNaNsOrInfs(eec_pt_hist_det)) {
                    std::cerr << "Histogram eec contains NaNs or Infs!" << std::endl;
                    return;
                    }
               
                 h3_reco_eec->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                 h3_reco_eec->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
             
                     if (hasNaNsOrInfs(h3_reco_eec) || hasNaNsOrInfs(h3_true_eec)) {
                         std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                         return;
                        }
                    
                
                //Apply truth level cuts before filling matched truth histogram -- this mainly avoids underflow and overflow
                if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue;
                if(dR_tru < 0.005 || dR_tru >0.4) continue;
                if(wt_tru_3D < 0.0001 || wt_tru_3D > 0.26) continue;
                
                h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D, mc);
                h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D, mc);
                
                h3_true_eec->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                h3_true_eec->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                   
                  
                response3D.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
                response3D.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
                
                
                }
            else{cerr<<"something is wrong"<<endl;}
                
             
            }
            
        } 
        
}

if(ifData && ifPairCut){
  for (Int_t i = 0; i < nEntries; ++i) {
    
        tree->GetEntry(i);
        
       // Fill full efficiency histogram if the true jet pT has changed
    if (jetpT_tru_fulleff != fJet_pt_tru) {
        h1_fulleff->Fill(fJet_pt_tru, mc);
    }
    jetpT_tru_fulleff = fJet_pt_tru;
    
    // Skip if either pT is zero (fake or missed track)
    if (fJet_pt_det == 0 || fJet_pt_tru == 0) continue;
    
      // Apply cuts on true jet pT
    if (fJet_pt_tru > ptJetmax_tru || fJet_pt_tru < 0) continue;

    // Fill histograms for unique jets
    if (jetpT_det != fJet_pt_det && jetpT_tru != fJet_pt_tru) {
        if (fJet_pt_det > ptJetmin_det && fJet_pt_det < ptJetmax_det) {
            h1_true->Fill(fJet_pt_tru, mc); 
            h1_reco->Fill(fJet_pt_det, mc);
            response1D.Fill(fJet_pt_det, fJet_pt_tru, mc);
            
        }
    }
    
     jetpT_tru = fJet_pt_tru;
     jetpT_det = fJet_pt_det;
  }

    jetpT_tru = -1;
    jetpT_det = -1;
    

  for (Int_t i = 0; i < nEntries; ++i) {
    
        tree->GetEntry(i);
        
        jetpT_tru = fJet_pt_tru;
        jetpT_det = fJet_pt_det;
     
        if(fTrack_pt_tru!=0){pt_tru = fTrack_pt_tru; mi = false;}
        else{pt_tru = fTrack_pt_miss; mi = true;}
        
        if(fTrack_pt_det!=0){pt_det = fTrack_pt_det; fi = false;}
        else{pt_det = fTrack_pt_fake; fi = true;}
        
        if(fTrack_eta_tru!=0){eta_tru = fTrack_eta_tru;}
        else{eta_tru = fTrack_eta_miss;}
        
        if(fTrack_eta_det!=0){eta_det = fTrack_eta_det;}
        else{eta_det = fTrack_eta_fake;}
        
        if(fTrack_phi_tru!=0){phi_tru = fTrack_phi_tru;}
        else{phi_tru = fTrack_phi_miss;}
        
        if(fTrack_phi_det!=0){phi_det = fTrack_phi_det;}
        else{phi_det = fTrack_phi_fake;}
       
    //   cout<<i<<" det: "<<jetpT_det<<" and true "<<jetpT_tru << endl;
        //     //2-pt corr
            for (Int_t k = i+1; k < nEntries; ++k) {
                tree->GetEntry(k);
                
                if(fTrack_pt_tru!=0){pt_tru2 = fTrack_pt_tru; mk = false;}
                else{pt_tru2 = fTrack_pt_miss; mk = true;}
        
                if(fTrack_pt_det!=0){pt_det2 = fTrack_pt_det; fk = false;}
                else{pt_det2 = fTrack_pt_fake; fk = true;}
        
                if(fTrack_eta_tru!=0){eta_tru2 = fTrack_eta_tru;}
                else{eta_tru2 = fTrack_eta_miss;}
        
                if(fTrack_eta_det!=0){eta_det2 = fTrack_eta_det;}
                else{eta_det2 = fTrack_eta_fake;}
         
                if(fTrack_phi_tru!=0){phi_tru2 = fTrack_phi_tru;}
                else{phi_tru2 = fTrack_phi_miss;}
        
                if(fTrack_phi_det!=0){phi_det2 = fTrack_phi_det;}
                else{phi_det2 = fTrack_phi_fake;}
               
            //   cout<<i<<" "<<k <<" jetpT_det "<< jetpT_det<< " and fJet_pt_det "<<fJet_pt_det<<endl;
            //   cout<<i<<" "<<k<< " jetpT_tru "<<jetpT_tru << " and fJet_pt_tru "<<fJet_pt_tru<<endl;
            
                //if both true and detector jet pT change then we have moved on to a new jet
                //because if one of them changes then that just means that its a missed or a fake track in the same jet
                if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                break;} 
                
                
                double dR_tru = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                double dR_det = delR(phi_det,phi_det2,eta_det,eta_det2);
                
                //fake correlations//////////////////////////////////
                if(fi || fk){ 
                //     cout<<"fake: "<<endl;
                //  tree->Show(k);
                 double wt_det =(2*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                 double wt_det_3D = (pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                 if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                 if(dR_det < 0.01 || dR_det > 0.4) continue;
                 if(wt_det_3D < 0.0005 || wt_det_3D > 0.26 ) continue;
                 if(abs(eta_det - eta_det2)<0.008) continue;
                 h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                 h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                 
                 eec_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det*mc);
                }
                //missed correlations//////////////////////////////////
                else if(mi || mk){ 
                //  cout<<"miss: "<<endl;
                // tree->Show(k);
                double wt_tru= (2*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru);
                double wt_tru_3D = (pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru); //for filling 3D histogram
                if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue;
                if(dR_tru < 0.005 || dR_tru >0.4) continue;
                if(wt_tru_3D < 0.0001 || wt_tru_3D > 0.26) continue;
                if(abs(eta_tru - eta_tru2)<0.008) continue;
                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                
                eec_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru*mc);
                }
                //matched correlations//////////////////////////////////
                else if(!mi && !fi && !mk && !fk){
                    // cout<<"match: "<<endl; 
                    // tree->Show(k); 
                double wt_det =(2*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                double wt_det_3D = (pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                
                double wt_tru= (2*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru);
                double wt_tru_3D = (pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru); //for filling 3D histogram
                
                //3D KE histogram TRUTH level -- applied before any truth or reco level cuts 
                if(dR_tru > 0.005 && dR_tru < 0.4 && wt_tru_3D > 0.0001 && wt_tru_3D < 0.26 && 
                (abs(eta_tru - eta_tru2)>0.008)){
                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru,mc);
                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru,mc);
                } 
                
                
                //3D KE histogram DET level -- applied before any truth level cuts 
                 if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                 if(dR_det < 0.01 || dR_det > 0.4) continue;
                 if(wt_det_3D < 0.0005 || wt_det_3D > 0.26 ) continue;
                 if(abs(eta_det - eta_det2)<0.008) continue;
                 
                 h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                 h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
             
                eec_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru*mc);
                eec_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det*mc);
                if (hasNaNsOrInfs(eec_pt_hist) || hasNaNsOrInfs(eec_pt_hist_det)) {
                    std::cerr << "Histogram eec contains NaNs or Infs!" << std::endl;
                    return;
                    }
               
                 h3_reco_eec->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                 h3_reco_eec->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                
             
                     if (hasNaNsOrInfs(h3_reco_eec) || hasNaNsOrInfs(h3_true_eec)) {
                         std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                         return;
                        }
                    
                
                //Apply truth level cuts before filling matched truth histogram -- this mainly avoids underflow and overflow
                if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue;
                if(dR_tru < 0.005 || dR_tru >0.4) continue;
                if(wt_tru_3D < 0.0001 || wt_tru_3D > 0.26) continue;
                if(abs(eta_tru - eta_tru2)<0.008) continue;
                
                h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D, mc);
                h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D, mc);
                
                h3_true_eec->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                h3_true_eec->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                   
                  
                response3D.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
                response3D.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
                
                
                }
            else{cerr<<"something is wrong"<<endl;}
                
             
            }
            
        } 
        
}
        
if (tag == "rw") {cout<<"reweight"<<endl;}
        
        if(ifTrivial){
             for (Int_t i = 0; i < nEntries; ++i) {
    
        tree->GetEntry(i);
        
        jetpT_tru = fJet_pt_tru;
        jetpT_det = fJet_pt_det;
     
         if(fTrack_pt_tru!=0){pt_tru = fTrack_pt_tru; mi = false;}
        else{pt_tru = fTrack_pt_miss; mi = true;}
        
        if(fTrack_pt_det!=0){pt_det = fTrack_pt_det; fi = false;}
        else{pt_det = fTrack_pt_fake; fi = true;}
        
        if(fTrack_eta_tru!=0){eta_tru = fTrack_eta_tru;}
        else{eta_tru = fTrack_eta_miss;}
        
        if(fTrack_eta_det!=0){eta_det = fTrack_eta_det;}
        else{eta_det = fTrack_eta_fake;}
        
        if(fTrack_phi_tru!=0){phi_tru = fTrack_phi_tru;}
        else{phi_tru = fTrack_phi_miss;}
        
        if(fTrack_phi_det!=0){phi_det = fTrack_phi_det;}
        else{phi_det = fTrack_phi_fake;}
    //   cout<<i<<" det: "<<jetpT_det<<" and true "<<jetpT_tru << endl;
        //     //2-pt corr
            for (Int_t k = i+1; k < nEntries; ++k) {
                tree->GetEntry(k);
                
                if(fTrack_pt_tru!=0){pt_tru2 = fTrack_pt_tru; mk = false;}
                else{pt_tru2 = fTrack_pt_miss; mk = true;}
        
                if(fTrack_pt_det!=0){pt_det2 = fTrack_pt_det; fk = false;}
                else{pt_det2 = fTrack_pt_fake; fk = true;}
        
                if(fTrack_eta_tru!=0){eta_tru2 = fTrack_eta_tru;}
                else{eta_tru2 = fTrack_eta_miss;}
        
                if(fTrack_eta_det!=0){eta_det2 = fTrack_eta_det;}
                else{eta_det2 = fTrack_eta_fake;}
         
                if(fTrack_phi_tru!=0){phi_tru2 = fTrack_phi_tru;}
                else{phi_tru2 = fTrack_phi_miss;}
        
                if(fTrack_phi_det!=0){phi_det2 = fTrack_phi_det;}
                else{phi_det2 = fTrack_phi_fake;}
          
                //if both true and detector jet pT change then we have moved on to a new jet
                //because if one of them changes then that just means that its a missed or a fake track in the same jet
                if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                break;} 
                
                double dR_tru = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                double dR_det = delR(phi_det,phi_det2,eta_det,eta_det2);
                
                //matched correlations//////////////////////////////////
                if (!mi && !fi && !mk && !fk){
                    // cout<<"match: "<<endl; 
                    // tree->Show(k); 
                double wt_det =(2*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                double wt_det_3D = (pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                
                double wt_tru= (2*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru);
                double wt_tru_3D = (pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru); //for filling 3D histogram
              
                
                // //3D KE histogram DET level -- applied before any truth level cuts 
                //  if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                //  if(dR_det < 0.01 || dR_det > 0.4) continue;
                //  if(wt_det_3D < 0.0005 || wt_det_3D > 0.26 ) continue;
                 
              
                //  h3_reco_eecTriv->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                //  h3_reco_eecTriv->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
             
                //      if (hasNaNsOrInfs(h3_reco_eec) || hasNaNsOrInfs(h3_true_eec)) {
                //          std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                //          return;
                //         }
                    
                // if(fJet_pt_tru < ptJetmin_det || fJet_pt_tru >ptJetmax_det) continue;
                // if(dR_tru < 0.01 || dR_tru >0.4) continue;
                // if(wt_tru_3D < 0.0005 || wt_tru_3D > 0.26) continue;
               
                // h3_true_eecTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                // h3_true_eecTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                
                //reweighting the response to the prior
                int binx = -1;
                if (dR_tru < 0.01) binx = 1;
                else if (dR_tru >0.4) binx = h3_rw->GetNbinsX();
                else binx = h3_rw->GetXaxis()->FindBin(dR_tru);
                int biny = -1;
                if (fJet_pt_tru < 0) biny = 1;
                else if ( fJet_pt_tru > ptJetmax_det) biny = h3_rw->GetNbinsY();
                else biny = h3_rw->GetYaxis()->FindBin(fJet_pt_tru);
                int binz = -1;
                if (wt_tru_3D < 0.0005) binz = 1;
                else if (wt_tru_3D > 0.26) binz = h3_rw->GetNbinsZ();
                else binz = h3_rw->GetZaxis()->FindBin(wt_tru_3D);
                double rw_fac = h3_rw->GetBinContent(binx, biny, binz);
             
                // If rw_fac is zero, search for a nearby non-zero bin
            if (rw_fac == 0) {
                const int search_radius = 5;  // Radius for neighboring bins
                bool found_nonzero = false;

                for (int dx = -search_radius; dx <= search_radius && !found_nonzero; ++dx) {
                    for (int dy = -search_radius; dy <= search_radius && !found_nonzero; ++dy) {
                        for (int dz = -search_radius; dz <= search_radius && !found_nonzero; ++dz) {
                            int nx = binx + dx;
                            int ny = biny + dy;
                            int nz = binz + dz;

                            // Ensure indices are within valid ranges
                            if (nx >= 1 && nx <= h3_rw->GetNbinsX() &&
                                ny >= 1 && ny <= h3_rw->GetNbinsY() &&
                                nz >= 1 && nz <= h3_rw->GetNbinsZ()) {

                                rw_fac = h3_rw->GetBinContent(nx, ny, nz);
                                if (rw_fac > 0) {
                                    found_nonzero = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
                
                                
                
                
                 //3D KE histogram DET level -- applied before any truth level cuts 
                 if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det || fJet_pt_tru < ptJetmin_det || fJet_pt_tru >ptJetmax_det) continue;
                 if(dR_det < 0.01 || dR_det > 0.4 || dR_tru < 0.01 || dR_tru >0.4) continue;
                 if(wt_det_3D < 0.0005 || wt_det_3D > 0.26 || wt_tru_3D < 0.0005 || wt_tru_3D > 0.26 ) continue;
                 
                 
                // Apply reweighting if `rw` tag is set
                if (tag == "rw") mc *= rw_fac;
                // cout<<mc<<endl;
                // Check if mc is zero after reweighting and output a message if it is
                if (mc == 0) {
                    std::cout << "Warning: mc value is 0 after reweighting." << std::endl;
                    cout<<binx<<" "<<biny<<" "<<binz<<endl;
                }
              
                h3_reco_eecTriv->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                h3_reco_eecTriv->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
             
                h3_true_eecTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                h3_true_eecTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                
                eec_pt_histTriv->Fill(dR_tru,fJet_pt_tru,wt_tru*mc);
                eec_pt_hist_detTriv->Fill(dR_det,fJet_pt_det,wt_det*mc);
                  
                response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
                response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
                
                mc = xSection/nTrials;
                }

                
                
            }
          }    
        }



        if(ifSplit){
             for (Int_t i = 0; i < nEntries; ++i) {
    
        tree->GetEntry(i);
        
        jetpT_tru = fJet_pt_tru;
        jetpT_det = fJet_pt_det;
     
        if(fTrack_pt_tru!=0){pt_tru = fTrack_pt_tru; mi = false;}
        else{pt_tru = fTrack_pt_miss; mi = true;}
        
        if(fTrack_pt_det!=0){pt_det = fTrack_pt_det; fi = false;}
        else{pt_det = fTrack_pt_fake; fi = true;}
        
        if(fTrack_eta_tru!=0){eta_tru = fTrack_eta_tru;}
        else{eta_tru = fTrack_eta_miss;}
        
        if(fTrack_eta_det!=0){eta_det = fTrack_eta_det;}
        else{eta_det = fTrack_eta_fake;}
        
        if(fTrack_phi_tru!=0){phi_tru = fTrack_phi_tru;}
        else{phi_tru = fTrack_phi_miss;}
        
        if(fTrack_phi_det!=0){phi_det = fTrack_phi_det;}
        else{phi_det = fTrack_phi_fake;}
       
    //   cout<<i<<" det: "<<jetpT_det<<" and true "<<jetpT_tru << endl;
        //     //2-pt corr
            for (Int_t k = i+1; k < nEntries; ++k) {
                tree->GetEntry(k);
                
                if(fTrack_pt_tru!=0){pt_tru2 = fTrack_pt_tru; mk = false;}
                else{pt_tru2 = fTrack_pt_miss; mk = true;}
        
                if(fTrack_pt_det!=0){pt_det2 = fTrack_pt_det; fk = false;}
                else{pt_det2 = fTrack_pt_fake; fk = true;}
        
                if(fTrack_eta_tru!=0){eta_tru2 = fTrack_eta_tru;}
                else{eta_tru2 = fTrack_eta_miss;}
        
                if(fTrack_eta_det!=0){eta_det2 = fTrack_eta_det;}
                else{eta_det2 = fTrack_eta_fake;}
         
                if(fTrack_phi_tru!=0){phi_tru2 = fTrack_phi_tru;}
                else{phi_tru2 = fTrack_phi_miss;}
        
                if(fTrack_phi_det!=0){phi_det2 = fTrack_phi_det;}
                else{phi_det2 = fTrack_phi_fake;}
          
                //if both true and detector jet pT change then we have moved on to a new jet
                //because if one of them changes then that just means that its a missed or a fake track in the same jet
                if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                break;} 
                
                double dR_tru = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                double dR_det = delR(phi_det,phi_det2,eta_det,eta_det2);
                
                //matched correlations//////////////////////////////////
                if (!mi && !fi && !mk && !fk){
                    // cout<<"match: "<<endl; 
                    // tree->Show(k); 
                double wt_det =(2*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                double wt_det_3D = (pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                
                double wt_tru= (2*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru);
                double wt_tru_3D = (pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru); //for filling 3D histogram
               
                
                
                //3D KE histogram DET level -- applied before any truth level cuts 
                 if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det || fJet_pt_tru < ptJetmin_det || fJet_pt_tru >ptJetmax_det) continue;
                 if(dR_det < 0.01 || dR_det > 0.4 || dR_tru < 0.01 || dR_tru >0.4) continue;
                 if(wt_det_3D < 0.0005 || wt_det_3D > 0.26 || wt_tru_3D < 0.0005 || wt_tru_3D > 0.26 ) continue;
                 
                
                
                double split = rand->Rndm();
	                if (split < 0.5)
	                {
	          
	                response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
	                response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
	            
	                }
	             else
	                 {
	                    
	                h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
	                h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                    h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                    h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                    
                    eec_pt_histSplit->Fill(dR_tru,fJet_pt_tru,wt_tru*mc);
                    eec_pt_hist_detSplit->Fill(dR_det,fJet_pt_det,wt_det*mc);
	                }
                
                
                }
          

            }
          }    
        }

 
 if(ifFastSim){ 
    
   TFile* fpaireff = new TFile("/home/ar2545/ALICEmc/AnalysisResultsMC_May15_JESplot.root","OPEN");
   TH2D*    _h_n2_qpt_det = (TH2D*)fpaireff->Get("qpt_det");
   TH2D*    _h_n2_qpt_tru = (TH2D*)fpaireff->Get("qpt_tru");
    
      
        tree->Branch("fTrack_choverpt_det", &fTrack_choverpt_det,"fTrack_choverpt_det/D");
        tree->Branch("fTrack_choverpt_tru", &fTrack_choverpt_tru,"fTrack_choverpt_tru/D");
        tree->Branch("fTrack_choverpt_miss", &fTrack_choverpt_miss,"fTrack_choverpt_miss/D");
        tree->Branch("fTrack_choverpt_fake", &fTrack_choverpt_fake,"fTrack_choverpt_fake/D");
        
    double chpt_tru = 0;
    double chpt_det = 0;
    double chpt_tru2 = 0;
    double chpt_det2 = 0;
  
  for (Int_t i = 0; i < nEntries; ++i) {
    
        tree->GetEntry(i);
        
    // Skip if either pT is zero (fake or missed track)
    if (fJet_pt_det == 0 || fJet_pt_tru == 0) continue;
    
    
        // Fill full efficiency histogram if the true jet pT has changed
    if (jetpT_tru_fulleff != fJet_pt_tru) {
        h1_fulleff->Fill(fJet_pt_tru, mc);
    }
    jetpT_tru_fulleff = fJet_pt_tru;
    


    if (fJet_pt_tru > ptJetmax_tru || fJet_pt_tru < 0) continue;
    // Fill histograms for unique jets
    if (jetpT_det != fJet_pt_det && jetpT_tru != fJet_pt_tru) {
        if (fJet_pt_det > ptJetmin_det && fJet_pt_det < ptJetmax_det) 
        {
            h1_true->Fill(fJet_pt_tru, mc);
            h1_reco->Fill(fJet_pt_det, mc);
            response1D.Fill(fJet_pt_det, fJet_pt_tru, mc);
            
        }
    }
    
     jetpT_tru = fJet_pt_tru;
     jetpT_det = fJet_pt_det;
  }

    jetpT_tru = -1;
    jetpT_det = -1;
    

  for (Int_t i = 0; i < nEntries; ++i) {
    
        tree->GetEntry(i);
        
        jetpT_tru = fJet_pt_tru;
        jetpT_det = fJet_pt_det;
      
        
    // Apply really wide cuts on true jet pT///THIS LINE IS NEW -- CHECK THIS WITH LAURA
       if (fJet_pt_tru > (ptJetmax_tru) || fJet_pt_tru < 0) continue;
     
        if(fTrack_pt_tru!=0){pt_tru = fTrack_pt_tru; chpt_tru = fTrack_choverpt_tru; mi = false; }
        else{pt_tru = fTrack_pt_miss; chpt_tru = fTrack_choverpt_miss; mi = true;}
        
        if(fTrack_pt_det!=0){pt_det = fTrack_pt_det; chpt_det = fTrack_choverpt_det; fi = false;}
        else{pt_det = fTrack_pt_fake; chpt_det = fTrack_choverpt_fake; fi = true;}
        
        if(fTrack_eta_tru!=0){eta_tru = fTrack_eta_tru;}
        else{eta_tru = fTrack_eta_miss;}
        
        if(fTrack_eta_det!=0){eta_det = fTrack_eta_det;}
        else{eta_det = fTrack_eta_fake;}
        
        if(fTrack_phi_tru!=0){phi_tru = fTrack_phi_tru;}
        else{phi_tru = fTrack_phi_miss;}
        
        if(fTrack_phi_det!=0){phi_det = fTrack_phi_det;}
        else{phi_det = fTrack_phi_fake;}
       
    //   cout<<i<<" det: "<<jetpT_det<<" and true "<<jetpT_tru << endl;
        //     //2-pt corr
            for (Int_t k = i+1; k < nEntries; ++k) {
                tree->GetEntry(k);
                
                if(fTrack_pt_tru!=0){pt_tru2 = fTrack_pt_tru; chpt_tru2 = fTrack_choverpt_tru; mk = false;}
                else{pt_tru2 = fTrack_pt_miss; chpt_tru2 = fTrack_choverpt_miss; mk = true;}
        
                if(fTrack_pt_det!=0){pt_det2 = fTrack_pt_det; chpt_det2 = fTrack_choverpt_det; fk = false;}
                else{pt_det2 = fTrack_pt_fake; chpt_det2 = fTrack_choverpt_fake; fk = true;}
        
                if(fTrack_eta_tru!=0){eta_tru2 = fTrack_eta_tru;}
                else{eta_tru2 = fTrack_eta_miss;}
        
                if(fTrack_eta_det!=0){eta_det2 = fTrack_eta_det;}
                else{eta_det2 = fTrack_eta_fake;}
         
                if(fTrack_phi_tru!=0){phi_tru2 = fTrack_phi_tru;}
                else{phi_tru2 = fTrack_phi_miss;}
        
                if(fTrack_phi_det!=0){phi_det2 = fTrack_phi_det;}
                else{phi_det2 = fTrack_phi_fake;}
               
            //   cout<<i<<" "<<k <<" jetpT_det "<< jetpT_det<< " and fJet_pt_det "<<fJet_pt_det<<endl;
            //   cout<<i<<" "<<k<< " jetpT_tru "<<jetpT_tru << " and fJet_pt_tru "<<fJet_pt_tru<<endl;
            
                //if both true and detector jet pT change then we have moved on to a new jet
                //because if one of them changes then that just means that its a missed or a fake track in the same jet
                if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                break;} 
                
                
                double dR_tru = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                double dR_det = delR(phi_det,phi_det2,eta_det,eta_det2);
                
                double chptdiff = abs(chpt_det - chpt_det2);
                double pairwt = paireffwt(_h_n2_qpt_det, _h_n2_qpt_tru, chptdiff, dR_det);
                
                //Multiply the pair efficiency with the MC weight
                mc*=pairwt;
                //fake correlations//////////////////////////////////
                if(fi || fk){ 
                //     cout<<"fake: "<<endl;
                //  tree->Show(k);
                 double wt_det =(2*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                 double wt_det_3D = (pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                 if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                 if(dR_det < 0.01 || dR_det > 0.4) continue;
                 if(wt_det_3D < 0.0005 || wt_det_3D > 0.26 ) continue;
                 h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                 h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                 
                  eec_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det*mc);
                }
                //missed correlations//////////////////////////////////
                else if(mi || mk){ 
                //  cout<<"miss: "<<endl;
                // tree->Show(k);
                double wt_tru= (2*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru);
                double wt_tru_3D = (pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru); //for filling 3D histogram
                // if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue; //only wide cuts needed for full eff
                if(dR_tru < 0.005 || dR_tru >0.4) continue;
                if(wt_tru_3D < 0.0001 || wt_tru_3D > 0.26) continue;
                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                
                eec_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru*mc);
                }
                //matched correlations//////////////////////////////////
                else if(!mi && !fi && !mk && !fk){
                    // cout<<"match: "<<endl; 
                    // tree->Show(k); 
                double wt_det =(2*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                double wt_det_3D = (pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det);
                
                double wt_tru= (2*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru);
                double wt_tru_3D = (pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru); //for filling 3D histogram
                
                //3D KE histogram TRUTH level -- applied before any truth or reco level cuts //only bin edge cuts and super wide cuts for jet pT
                if(dR_tru > 0.01 && dR_tru < 0.4 && wt_tru_3D > 0.0001 && wt_tru_3D < 0.26){
                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                
                eec_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru*mc);
                } 
                
                 
                
                //3D KE histogram DET level -- applied before any truth level cuts 
                 if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                 if(dR_det < 0.01 || dR_det > 0.4) continue;
                 if(wt_det_3D < 0.0005 || wt_det_3D > 0.26 ) continue;
                 
                 h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                 h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
             
               
                eec_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det*mc);
                if (hasNaNsOrInfs(eec_pt_hist) || hasNaNsOrInfs(eec_pt_hist_det)) {
                    std::cerr << "Histogram eec contains NaNs or Infs!" << std::endl;
                    return;
                    }
               
                 h3_reco_eec->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                 h3_reco_eec->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
             
                     if (hasNaNsOrInfs(h3_reco_eec) || hasNaNsOrInfs(h3_true_eec)) {
                         std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                         return;
                        }
                    
                
                //Apply truth level cuts before filling matched truth histogram -- this mainly avoids underflow and overflow
                if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue;
                if(dR_tru < 0.005 || dR_tru >0.4) continue;
                if(wt_tru_3D < 0.0001 || wt_tru_3D > 0.26) continue;
                
                h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D, mc);
                h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D, mc);
                
                h3_true_eec->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                h3_true_eec->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                   
                  
                response3D.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
                response3D.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
                
                
                }
            else{cerr<<"something is wrong"<<endl;}
                
             
            }
            
        } 
        
    }
 
    // Find the last '/' in the string
    size_t pos = fileName.find_last_of("/\\");
    
    // Extract the substring after the last '/' (this is the filename)
    std::string outfileName = fileName.substr(pos + 1);
    
  
    TFile *fout=new TFile (Form("/home/ar2545/Trial/Unfold3DData_%i_%s",pthardbin,outfileName.c_str()),"RECREATE");
    
    TString foutname = Form("/home/ar2545/Trial/Unfold3DData_%i_%s.root", pthardbin,outfileName.c_str());
    cout<<"writing to output file named "<< foutname<<endl;
    
    
    cout<<"write histograms"<<endl;
    
    eec_pt_hist->Write();
    eec_pt_hist_det->Write();
    eec_pt_histTriv->Write();
    eec_pt_hist_detTriv->Write();
    eec_pt_histSplit->Write();
    eec_pt_hist_detSplit->Write();

    h1_true->Write();
    h1_reco->Write();
    h1_fulleff->Write();
    h3_true_eec->Write();
    h3_reco_eec->Write();
    h3_true_eecTriv->Write();
    h3_reco_eecTriv->Write();
    h3_rawSplit->Write();
    h3_trueSplit->Write();
    h3_fulleff->Write();
    h3_fullreco->Write();
    h3_eff_match->Write();
    h3_smeared->Write();
   

    cout<<"Writing response matrices"<<endl;
    response3D.Write();
    response3DSplit.Write();
    response3DTriv.Write();
    response1D.Write();
    
    fout->Close();
     cout<<"Response matrices written. Success!"<<endl;
  }  
  else{
      cout<<"0 ENTRIES!"<<endl;
    
  }  
    fmc->Close();
     
 

}
#ifndef __CINT__
int main () { Unfolding3dAnalysisOptNew(); return 0; }  // Main program when run stand-alone
#endif
