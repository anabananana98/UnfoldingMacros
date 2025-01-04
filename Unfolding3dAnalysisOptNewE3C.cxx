//Create response matrix for E3C 
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

void Unfolding3dAnalysisOptNewE3C(std::string tag = "", int pthardbin = -1, std::string fileName ="") {

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
    
    TH3D* h3_true_e3c; TH3D* h3_trueSplit;
    TH3D* h3_reco_e3c;
    TH3D* h3_true_e3cTriv;
    TH3D* h3_reco_e3cTriv;
    TH3D* h3_unfolded; TH3D* h3_unfolded_rw;  TH3D* h3_unfoldedSplit;
    TH3D* h3_raw; TH3D* h3_rawSplit;
    TH3D* h3_smeared;
    TH3D* h3_fulleff;
    TH3D* h3_eff_match;
    TH3D* h3_fullreco;
    TH3D* h3_rw;
    TH3* h3_fold;
    TH3* h3_triv;
    
    TH2D* e3c_pt_hist;
    TH2D* e3c_pt_hist_det;
    TH2D* e3c_pt_histTriv;
    TH2D* e3c_pt_hist_detTriv;
    TH2D* e3c_pt_histSplit;
    TH2D* e3c_pt_hist_detSplit;
    
    TFile *fmc;

    TRandom3* rand = new TRandom3(0);
    
    RooUnfoldResponse response1D;
    RooUnfoldResponse response1DSplit;
    
    RooUnfoldResponse response3D;
    RooUnfoldResponse response3DTriv;
    RooUnfoldResponse response3DSplit;
   
    float pThard_val[21] = {5., 7., 9., 12., 16., 21., 28., 36., 45., 57., 70., 85., 99., 115., 132., 150., 169., 190., 212., 235., 1000.};
    
    
    //  //reco jet pT bins for 1D unfolding 
    // Double_t from_const_reco_jet = 10;
    // Double_t to_const_reco_jet = 80;
    // Int_t bins_const_reco_jet = 14;
    // Double_t width_const_reco_jet = (to_const_reco_jet-from_const_reco_jet)/bins_const_reco_jet;
    // Double_t xbins_jet[16] = {};
    // for (Int_t i = 0; i <= bins_const_reco_jet; i++)
    // {
    //     xbins_jet[i] = (from_const_reco_jet + i * width_const_reco_jet);
    // }
    // xbins_jet[15] = 120;
    
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
    
    //  //true jet pT bins for 1D unfolding 
    // Double_t from_const_jet = 0;
    // Double_t to_const_jet= 200;
    // Int_t bins_const_jet = 40;
    // Double_t width_const_jet = (to_const_jet-from_const_jet)/bins_const_jet;
    // Double_t tbins_jet[41] = {};
    // for (Int_t i = 0; i <= bins_const_jet; i++)
    // {
    //     tbins_jet[i] = (from_const_jet + i * width_const_jet);
    // }
    
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
    
    
    
      //Nov3 trial
   double dRbins[] = {0.0001,0.01,0.0144613,0.0173904,0.0209128,0.0251487,0.0302425,
 0.0363681,0.0437345,0.0525929,0.0632456,0.0760559,0.091461,0.109986,0.132264,0.159054,0.191271,0.230012,0.276601,0.332627,0.4,0.524807,1};
 
 double dRbinst[] = {0.0001,0.01,0.0144613,0.0173904,0.0209128,0.0251487,0.0302425,
 0.0363681,0.0437345,0.0525929,0.0632456,0.0760559,0.091461,0.109986,0.132264,0.159054,0.191271,0.230012,0.276601,0.332627,0.4,0.524807,1};
    
    // //reco dR bins --same for eec and e3c
    // double dRbins[] = {
    //     0.000102329,
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
    
    //commented out on dec 4 
    // //true dR bins --same for eec and e3c
    // double dRbinst[] = {
    //     0.000102329,
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
    
        
    //     //reco wt bins e3c
    // double wtbins[] = {
    //   0.00010232902329,
    //     0.00151356, 0.00301995, 0.00506991, 0.00716143, 0.00912011, 0.0153109,
    //     0.0231739, 0.0305492,
    //     0.0462381,
    //     0.0922571,
    //     0.160325,
    //     1.0};
    
    // //true wt bins e3c
    // double wtbinst[] = {
    //   0.00010232902329,
    //     0.00151356, 0.00301995, 0.00506991, 0.00716143, 0.00912011, 0.0153109,
    //     0.0231739, 0.0305492,
    //     0.0462381,
    //     0.0922571,
    //     0.160325,
    //     1.0
    //     };
   
    // double wtbins[] = {     
    //     0.000102329,0.000121619,0.000251189, 0.000501187,
    //     0.00151356, 0.00301995, 0.00506991, 0.00716143, 0.00912011, 0.0101158, 0.0104713, 0.0108393, 0.0128825, 0.0133352, 0.0138038, 0.0142889, 0.0153109,
    //     0.0164059, 0.018197, 0.0201837, 0.025704, 0.0305492,
    //     0.0363078, 0.0402717, 0.0462381, 0.0512861, 0.0568853,0.0653131, 0.0724436,
    //     0.0803526, 0.0988553,
    //     0.15
    // };
    //  double wtbinst[] = {
    //   0.00001,0.000102329,0.000121619,0.000251189, 0.000501187,
    //     0.00151356, 0.00301995, 0.00506991, 0.00716143, 0.00912011, 0.0101158, 0.0104713, 0.0108393, 0.0128825, 0.0133352, 0.0138038, 0.0142889, 0.0153109,
    //     0.0164059, 0.018197, 0.0201837, 0.025704, 0.0305492,
    //     0.0363078, 0.0402717, 0.0462381, 0.0512861, 0.0568853,0.0653131, 0.0724436,
    //     0.0803526, 0.0988553,
    //     0.15
    //     };
    
    
    //Dec 3 trial
    double wtbins[] = {
        0.000102329,0.000121619,0.000251189, 0.000501187,
        0.00151356, 0.00301995, 0.00506991, 0.00716143, 0.00912011, 0.0101158,  0.0108393, 0.0128825,  0.0153109,
        0.0175792, 0.0201837,
        0.0363078, 0.0462381, 0.0512861,0.0724436,
        0.0988553,
        0.15
    };
    
      double wtbinst[] = {
      0.00001,0.000102329,0.000121619,0.000251189, 0.000501187,
        0.00151356, 0.00301995, 0.00506991, 0.00716143, 0.00912011, 0.0101158,  0.0108393, 0.0128825,  0.0153109,
        0.0175792, 0.0201837,
        0.0363078, 0.0462381, 0.0512861,0.0724436,
        0.0988553,
        0.15
        };
    
    
    // //reco wt bins eec
    // double wtbins[] = {
    //     0.000102329, 0.00301995, 0.00501187, 0.01, 0.0158489, 0.020893, 0.0301995, 0.0501187, 0.060256, 0.1,
    //     0.151356, 0.758578,1.0};
    
    // //true wt bins eec
    // double wtbinst[] = {
    //     0.000102329, 0.00301995, 0.00501187, 0.01, 0.0158489, 0.020893, 0.0301995, 0.0501187, 0.060256, 0.1,
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
    
    
    e3c_pt_hist = new TH2D("e3c_pt_hist", "E3C and jet_pt 2D",  ndRbinst, dRbinst, nJetPtbinst, tbins);
    e3c_pt_hist_det = new TH2D("e3c_pt_hist_det", "E3C and jet_pt 2D det", ndRbins, dRbins, nJetPtbins,xbins);
    
    e3c_pt_histTriv = new TH2D("e3c_pt_histTriv", "E3C and jet_pt 2D Triv",  ndRbinst, dRbinst, nJetPtbinst, tbins);
    e3c_pt_hist_detTriv = new TH2D("e3c_pt_hist_detTriv", "E3C and jet_pt 2D det Triv", ndRbins, dRbins, nJetPtbins,xbins);
    
    e3c_pt_histSplit = new TH2D("e3c_pt_histSplit", "E3C and jet_pt 2D Split",  ndRbinst, dRbinst, nJetPtbinst, tbins);
    e3c_pt_hist_detSplit = new TH2D("e3c_pt_hist_detSplit", "E3C and jet_pt 2D det Split", ndRbins, dRbins, nJetPtbins,xbins);
    
    
    h1_true = new TH1D("h1_true", "h1_true", nJetPtbins1Dtrue, tbins_jet);
    h1_fulleff = new TH1D("h1_fulleff", "h1_fulleff", nJetPtbins1Dtrue, tbins_jet);
    h1_reco = new TH1D("h1_reco", "h1_reco", nJetPtbins1Dreco, xbins_jet);

    
   cout<<"nWtbinst "<<nWtbinst<<" and nWtbins "<<nWtbins<<endl;
   cout<<"ndRbinst "<<ndRbinst<<" and ndRbins "<< ndRbins <<endl;
   cout<<"nJetPtbinst "<<nJetPtbinst<<" and nJetPtbins "<<nJetPtbins<<endl;    
   
   cout<<"nJetPtbins1Dtrue "<<nJetPtbins1Dtrue<<" and nJetPtbins1Dreco "<<nJetPtbins1Dreco<<endl;
   
  
    
    h3_true_e3c = new TH3D("h3_true_e3c", "h3_true_e3c", ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_reco_e3c = new TH3D("h3_reco_e3c", "h3_reco_e3c",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
    h3_smeared = new TH3D("h3_smeared", "h3_smeared",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
    
    h3_true_e3cTriv = new TH3D("h3_true_e3cTriv", "h3_true_e3c_triv", ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_reco_e3cTriv = new TH3D("h3_reco_e3cTriv", "h3_reco_e3c_triv",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
    
    h3_trueSplit = new TH3D("h3_trueSplit", "h3_true_split", ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_rawSplit = new TH3D("h3_rawSplit", "h3_reco_split",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
   
    h3_fulleff = new TH3D("h3_fulleff","h3_fulleff",ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_eff_match = new TH3D("h3_eff_match","h3_eff_match",ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_fullreco = new TH3D("h3_fullreco","h3_fullreco",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
    
    
    h3_true_e3c->Sumw2();
    h3_reco_e3c->Sumw2();
    
    h3_reco_e3cTriv->Sumw2();
    h3_reco_e3cTriv->Sumw2();
    
    h3_rawSplit->Sumw2();
    h3_rawSplit->Sumw2();
  
    h3_eff_match->Sumw2();
    h3_fulleff->Sumw2();
    h3_fullreco->Sumw2();
  
   cout<<"HERE"<<endl;
  
    cout<<"setting up response declarations"<<endl;
    response3D.Setup(h3_reco_e3c,h3_true_e3c);
    response3DTriv.Setup(h3_reco_e3cTriv,h3_true_e3cTriv);
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
        
        
        // tree->GetEntry(1);
        // double mc = mc_weight; //same for one run in 1 pTHat bin
        // double weightresp = mc;
        
        // double ptJetmin_det = 10.; //min det level pt
        double ptJetmin_det = 15.; //min det level pt
        double ptJetmax_det = 120.; //max det level pt
        double ptJetmin_tru3D = 0.; //min tru level pt for 3D unfolding
        double ptJetmax_tru3D = 200; //max tru level pt for 3D unfolding
        
        
        //for 1D unfolding
        double ptJetmin_tru = 0.; //min tru level pt
        double ptJetmax_tru = 200.; //max tru level pt
        
        
        double jetpT_tru = -1;
        double jetpT_tru_fulleff = -1;
        double jetpT_det = -1;
        
        
        
        double pt_tru=-1;double pt_det=-1;
        double pt_tru2=-1;double pt_det2=-1;
        double pt_miss=-1;double pt_fake=-1;
        
        double eta_tru=-1;double eta_det=-1;
        double eta_tru2=-1;double eta_det2=-1;
        double eta_miss=-1;double eta_fake=-1;
        
        
        double phi_tru=-1;double phi_det=-1;
        double phi_tru2=-1;double phi_det2=-1;
        double phi_miss=-1;double phi_fake=-1;
        
        double pt_tru3=-1;double pt_det3=-1; double eta_tru3=-1;double eta_det3=-1;double phi_tru3=-1;double phi_det3=-1;
        
        bool fi = false; bool fj = false; bool fk = false;
        bool mi = false; bool mj = false; bool mk =false;
        
        jetpT_tru = -1;
        jetpT_det = -1;
        
        bool ifData = true;
        bool ifTrivial = false;
        bool ifSplit = false;
        bool ifPaircut = false;
        bool ifFastSim = false;
        
        if(ifData && !ifPaircut){
            cout<<"data loop "<<endl;
            
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
                
                 // Apply really wide cuts on true jet pT///THIS LINE IS NEW -- CHECK THIS WITH LAURA
                 if (fJet_pt_tru > (ptJetmax_tru) || fJet_pt_tru < 0) continue;
                
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
                
                
                for (Int_t k = i+1; k < nEntries; ++k) {
                    tree->GetEntry(k);
                    
                    if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                        // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                        break;}
                        
                    if(fTrack_pt_tru!=0){pt_tru2 = fTrack_pt_tru;mk = false;}
                    else{pt_tru2 = fTrack_pt_miss; mk = true;}
                    
                    if(fTrack_pt_det!=0){pt_det2 = fTrack_pt_det;fk = false;}
                    else{pt_det2 = fTrack_pt_fake; fk = true;}
                    
                    if(fTrack_eta_tru!=0){eta_tru2 = fTrack_eta_tru;}
                    else{eta_tru2 = fTrack_eta_miss;}
                    
                    if(fTrack_eta_det!=0){eta_det2 = fTrack_eta_det;}
                    else{eta_det2 = fTrack_eta_fake;}
                    
                    if(fTrack_phi_tru!=0){phi_tru2 = fTrack_phi_tru;}
                    else{phi_tru2 = fTrack_phi_miss;}
                    
                    if(fTrack_phi_det!=0){phi_det2 = fTrack_phi_det;}
                    else{phi_det2 = fTrack_phi_fake;}
                    
                    
                    double dR_tru = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                    double dR_det = delR(phi_det,phi_det2,eta_det,eta_det2);
                    
                    //fake correlations//////////////////////////////////
                    if(fi || fk){
                        //     cout<<"fake: "<<endl;
                        //  tree->Show(k);
                        double wt_det_iik =(3*pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_ikk =(3*pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        double wt_det_3D_iik = (pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_3D_ikk = (pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                        if(dR_det > 0.01 && dR_det < 0.4){
                            if(wt_det_3D_iik > 0.000102329 && wt_det_3D_iik < 0.15 ){
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                
                                e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_iik*mc);
                            }
                            if(wt_det_3D_ikk > 0.000102329 && wt_det_3D_ikk < 0.15 ){
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                
                                e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ikk*mc);
                            }
                        }
                    }
                    //missed correlations//////////////////////////////////
                    else if(mi || mk){
                        //  cout<<"miss: "<<endl;
                        // tree->Show(k);
                        double wt_tru_iik= (3*pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_ikk= (3*pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        double wt_tru_3D_iik = (pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_3D_ikk = (pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue;
                        //was dR_det --fixed on Dec 18
                        if(dR_tru > 0.01 && dR_tru < 0.4){
                            if(wt_tru_3D_iik > 0.00001 && wt_tru_3D_iik < 0.15 ){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                
                                e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_iik*mc);
                            }
                            if(wt_tru_3D_ikk > 0.00001 && wt_tru_3D_ikk < 0.15 ){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                
                                e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ikk*mc);
                            }
                        }
                    }
                    //matched correlations//////////////////////////////////
                    else if(!mi && !mk && !fk && !fi){
                        
                        double wt_tru_iik= (3*pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_ikk= (3*pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        double wt_tru_3D_iik = (pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_3D_ikk = (pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        double wt_det_iik =(3*pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_ikk =(3*pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        double wt_det_3D_iik = (pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_3D_ikk = (pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        
                        //3D KE histogram TRUTH level -- applied before any truth or reco level cuts
                        if(dR_tru > 0.005 && dR_tru < 0.4){
                            if(wt_tru_3D_iik > 0.00001 && wt_tru_3D_iik < 0.15){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                
                                e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_iik*mc);
                            }
                            
                            if(wt_tru_3D_ikk > 0.00001 && wt_tru_3D_ikk < 0.15){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                
                                e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ikk*mc);
                            }
                        }
                        
                        //3D KE histogram DET level -- applied before any truth level cuts
                        if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                        
                        //because even if this condition is not met for ik , it might be met for ijk
                        if(dR_det > 0.01 && dR_det < 0.4){
                            if(wt_det_3D_iik > 0.000102329 && wt_det_3D_iik < 0.15 ){
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                
                                
                                e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_iik*mc);
                                if (hasNaNsOrInfs(e3c_pt_hist) || hasNaNsOrInfs(e3c_pt_hist_det)) {
                                    std::cerr << "Histogram e3c contains NaNs or Infs!" << std::endl;
                                    return;
                                }
                                
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                
                                if (hasNaNsOrInfs(h3_reco_e3c) || hasNaNsOrInfs(h3_true_e3c)) {
                                    std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                                    return;
                                }
                                
                               
                                //Apply truth level cuts before filling matched truth histogram -- this mainly avoids underflow and overflow
                                if(fJet_pt_tru > ptJetmin_tru && fJet_pt_tru < ptJetmax_tru && dR_tru > 0.005 && dR_tru < 0.4){
                                    if(wt_tru_3D_iik > 0.00001 && wt_tru_3D_iik < 0.15) {
                                        
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_iik, mc);
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_iik, mc);
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_iik, mc);
                                        
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                        
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                    }
                                }
                            }
                            if(wt_det_3D_ikk > 0.000102329 && wt_det_3D_ikk < 0.15 ){
                                
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                
                                
                                e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ikk*mc);
                                if (hasNaNsOrInfs(e3c_pt_hist) || hasNaNsOrInfs(e3c_pt_hist_det)) {
                                    std::cerr << "Histogram eec contains NaNs or Infs!" << std::endl;
                                    return;
                                }
                                
                                
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                
                                if (hasNaNsOrInfs(h3_reco_e3c) || hasNaNsOrInfs(h3_true_e3c)) {
                                    std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                                    return;
                                }
                                
                                if(fJet_pt_tru > ptJetmin_tru && fJet_pt_tru < ptJetmax_tru && dR_tru > 0.005 && dR_tru < 0.4){
                                    if(wt_tru_3D_ikk > 0.00001 && wt_tru_3D_ikk < 0.15){
                                        
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, mc);
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, mc);
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, mc);
                                        
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                        
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                        
                                    }
                                }
                            }
                        }
                    }
                    else{cerr<<"something is wrong"<<endl;
                    }
                    
                    for (Int_t j = k+1; j < nEntries; ++j) {
                        tree->GetEntry(j);
                        
                        if(fTrack_pt_tru!=0){pt_tru3 = fTrack_pt_tru;mj = false;}
                        else{pt_tru3 = fTrack_pt_miss; mj = true;}
                        
                        if(fTrack_pt_det!=0){pt_det3 = fTrack_pt_det;fj = false;}
                        else{pt_det3 = fTrack_pt_fake; fj = true;}
                        
                        if(fTrack_eta_tru!=0){eta_tru3 = fTrack_eta_tru;}
                        else{eta_tru3 = fTrack_eta_miss;}
                        
                        if(fTrack_eta_det!=0){eta_det3 = fTrack_eta_det;}
                        else{eta_det3 = fTrack_eta_fake;}
                        
                        if(fTrack_phi_tru!=0){phi_tru3 = fTrack_phi_tru;}
                        else{phi_tru3 = fTrack_phi_miss;}
                        
                        if(fTrack_phi_det!=0){phi_det3 = fTrack_phi_det;}
                        else{phi_det3 = fTrack_phi_fake;}
                        
                        
                        if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                            // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                            break;}
                        
                        double dR_tru = -1;
                        double dR_det = -1;
                        
                        double dR_tru_ik = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                        double dR_tru_jk = delR(phi_tru2,phi_tru3,eta_tru2,eta_tru3);
                        double dR_tru_ij = delR(phi_tru,phi_tru3,eta_tru,eta_tru3);
                        
                        double dR_det_ik = delR(phi_det,phi_det2,eta_det,eta_det2);
                        double dR_det_jk = delR(phi_det2,phi_det3,eta_det2,eta_det3);
                        double dR_det_ij = delR(phi_det,phi_det3,eta_det,eta_det3);
                        
                        if(dR_tru_ik > dR_tru_jk && dR_tru_ik > dR_tru_ij){dR_tru = dR_tru_ik;}
                        else if (dR_tru_jk > dR_tru_ik && dR_tru_jk > dR_tru_ij){dR_tru = dR_tru_jk;}
                        else{dR_tru = dR_tru_ij;}
                        
                        if(dR_det_ik > dR_det_jk && dR_det_ik > dR_det_ij){dR_det = dR_det_ik;}
                        else if (dR_det_jk > dR_det_ik && dR_det_jk > dR_det_ij){dR_det = dR_det_jk;}
                        else{dR_det = dR_det_ij;}
                        
                        if(fi || fk || fj){
                            double wt_det_ijk =(6*pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            double wt_det_3D_ijk = (pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            
                            if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                            if(dR_det < 0.01 || dR_det > 0.4) continue;
                            if(wt_det_3D_ijk < 0.000102329 || wt_det_3D_ijk > 0.15 ) continue;
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                             e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ijk*mc);
                        }
                        //missed correlations//////////////////////////////////
                        else if(mi || mk || mj){
                            //  cout<<"miss: "<<endl;
                            // tree->Show(k);
                            double wt_tru_ijk= (6*pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            double wt_tru_3D_ijk = (pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            
                            if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue;
                            if(dR_tru < 0.005 || dR_tru > 0.4) continue;
                            if(wt_tru_3D_ijk < 0.00001 || wt_tru_3D_ijk > 0.15) continue;
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            
                            e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ijk*mc);
                        }
                        //matched correlations//////////////////////////////////
                        else if(!mi && !mk && !mj && !fi && !fk && !fj){
                            
                            double wt_tru_ijk= (6*pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            double wt_tru_3D_ijk = (pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            
                            double wt_det_ijk =(6*pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            double wt_det_3D_ijk = (pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            
                            
                            //3D KE histogram TRUTH level -- applied before any truth or reco level cuts
                            if(dR_tru > 0.005 && dR_tru < 0.4 && wt_tru_3D_ijk > 0.00001 && wt_tru_3D_ijk < 0.15){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                
                                 e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ijk*mc);
                            }
                            
                            
                            //3D KE histogram DET level -- applied before any truth level cuts
                            if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                            if(dR_det < 0.01 || dR_det > 0.4) continue;
                            if(wt_det_3D_ijk < 0.000102329 || wt_det_3D_ijk > 0.15 ) continue;
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                           
                            e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ijk*mc);
                            if (hasNaNsOrInfs(e3c_pt_hist) || hasNaNsOrInfs(e3c_pt_hist_det)) {
                                std::cerr << "Histogram eec contains NaNs or Infs!" << std::endl;
                                return;
                            }
                            
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                            if (hasNaNsOrInfs(h3_reco_e3c) || hasNaNsOrInfs(h3_true_e3c)) {
                                std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                                return;
                            }
                            
                            
                            //Apply truth level cuts before filling matched truth histogram -- this mainly avoids underflow and overflow
                            if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue;
                            if(dR_tru < 0.005 || dR_tru >0.4) continue;
                            if(wt_tru_3D_ijk < 0.00001 || wt_tru_3D_ijk > 0.15) continue;
                            
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            
                            
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            
                        }
                        else{cerr<<"something is wrong"<<endl;}
                        
                    }
                }
                
            }
        }
        
        
 if(ifData && ifPaircut){
            cout<<"data loop w paircut "<<endl;
            
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
                
                 // Apply really wide cuts on true jet pT///THIS LINE IS NEW -- CHECK THIS WITH LAURA
                 if (fJet_pt_tru > (ptJetmax_tru) || fJet_pt_tru < 0) continue;
                
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
                
                
                for (Int_t k = i+1; k < nEntries; ++k) {
                    tree->GetEntry(k);
                    
                    if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                        // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                        break;}
                        
                    if(fTrack_pt_tru!=0){pt_tru2 = fTrack_pt_tru;mk = false;}
                    else{pt_tru2 = fTrack_pt_miss; mk = true;}
                    
                    if(fTrack_pt_det!=0){pt_det2 = fTrack_pt_det;fk = false;}
                    else{pt_det2 = fTrack_pt_fake; fk = true;}
                    
                    if(fTrack_eta_tru!=0){eta_tru2 = fTrack_eta_tru;}
                    else{eta_tru2 = fTrack_eta_miss;}
                    
                    if(fTrack_eta_det!=0){eta_det2 = fTrack_eta_det;}
                    else{eta_det2 = fTrack_eta_fake;}
                    
                    if(fTrack_phi_tru!=0){phi_tru2 = fTrack_phi_tru;}
                    else{phi_tru2 = fTrack_phi_miss;}
                    
                    if(fTrack_phi_det!=0){phi_det2 = fTrack_phi_det;}
                    else{phi_det2 = fTrack_phi_fake;}
                    
                    
                    double dR_tru = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                    double dR_det = delR(phi_det,phi_det2,eta_det,eta_det2);
                    
                    //fake correlations//////////////////////////////////
                    if(fi || fk){
                        //     cout<<"fake: "<<endl;
                        //  tree->Show(k);
                        double wt_det_iik =(3*pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_ikk =(3*pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        double wt_det_3D_iik = (pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_3D_ikk = (pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                        if(abs(eta_det - eta_det2)<0.008) continue;
                        if(dR_det > 0.01 && dR_det < 0.4){
                            if(wt_det_3D_iik > 0.000102329 && wt_det_3D_iik < 0.15 ){
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                
                                e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_iik*mc);
                            }
                            if(wt_det_3D_ikk > 0.000102329 && wt_det_3D_ikk < 0.15 ){
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                
                                e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ikk*mc);
                            }
                        }
                    }
                    //missed correlations//////////////////////////////////
                    else if(mi || mk){
                        //  cout<<"miss: "<<endl;
                        // tree->Show(k);
                        double wt_tru_iik= (3*pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_ikk= (3*pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        double wt_tru_3D_iik = (pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_3D_ikk = (pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue;
                        if(abs(eta_tru - eta_tru2)<0.008) continue;
                        if(dR_tru > 0.01 && dR_tru < 0.4){
                            if(wt_tru_3D_iik > 0.00001 && wt_tru_3D_iik < 0.15 ){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                
                                e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_iik*mc);
                            }
                            if(wt_tru_3D_ikk > 0.00001 && wt_tru_3D_ikk < 0.15 ){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                
                                e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ikk*mc);
                            }
                        }
                    }
                    //matched correlations//////////////////////////////////
                    else if(!mi && !mk && !fk && !fi){
                        
                        double wt_tru_iik= (3*pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_ikk= (3*pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        double wt_tru_3D_iik = (pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_3D_ikk = (pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        double wt_det_iik =(3*pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_ikk =(3*pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        double wt_det_3D_iik = (pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_3D_ikk = (pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        
                        //3D KE histogram TRUTH level -- applied before any truth or reco level cuts
                        if(dR_tru > 0.005 && dR_tru < 0.4){
                            if(abs(eta_tru - eta_tru2)<0.008) continue;
                            if(wt_tru_3D_iik > 0.00001 && wt_tru_3D_iik < 0.15){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                
                                e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_iik*mc);
                            }
                            
                            if(wt_tru_3D_ikk > 0.00001 && wt_tru_3D_ikk < 0.15){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                
                                e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ikk*mc);
                            }
                        }
                        
                        //3D KE histogram DET level -- applied before any truth level cuts
                        if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                        
                        //because even if this condition is not met for ik , it might be met for ijk
                        if(dR_det > 0.01 && dR_det < 0.4){
                            if(abs(eta_det - eta_det2)<0.008) continue;
                            if(wt_det_3D_iik > 0.000102329 && wt_det_3D_iik < 0.15 ){
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                
                                
                                e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_iik*mc);
                                if (hasNaNsOrInfs(e3c_pt_hist) || hasNaNsOrInfs(e3c_pt_hist_det)) {
                                    std::cerr << "Histogram e3c contains NaNs or Infs!" << std::endl;
                                    return;
                                }
                                
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                
                                if (hasNaNsOrInfs(h3_reco_e3c) || hasNaNsOrInfs(h3_true_e3c)) {
                                    std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                                    return;
                                }
                                
                               
                                //Apply truth level cuts before filling matched truth histogram -- this mainly avoids underflow and overflow
                                if(fJet_pt_tru > ptJetmin_tru && fJet_pt_tru < ptJetmax_tru && dR_tru > 0.005 && dR_tru < 0.4){
                                    if(wt_tru_3D_iik > 0.00001 && wt_tru_3D_iik < 0.15) {
                                        
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_iik, mc);
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_iik, mc);
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_iik, mc);
                                        
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                        
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                    }
                                }
                            }
                            if(wt_det_3D_ikk > 0.000102329 && wt_det_3D_ikk < 0.15 ){
                                
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                
                                
                                e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ikk*mc);
                                if (hasNaNsOrInfs(e3c_pt_hist) || hasNaNsOrInfs(e3c_pt_hist_det)) {
                                    std::cerr << "Histogram eec contains NaNs or Infs!" << std::endl;
                                    return;
                                }
                                
                                
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                
                                if (hasNaNsOrInfs(h3_reco_e3c) || hasNaNsOrInfs(h3_true_e3c)) {
                                    std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                                    return;
                                }
                                
                                if(fJet_pt_tru > ptJetmin_tru && fJet_pt_tru < ptJetmax_tru && dR_tru > 0.005 && dR_tru < 0.4){
                                    if(wt_tru_3D_ikk > 0.00001 && wt_tru_3D_ikk < 0.15){
                                        
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, mc);
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, mc);
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, mc);
                                        
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                        
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                        
                                    }
                                }
                            }
                        }
                    }
                    else{cerr<<"something is wrong"<<endl;
                    }
                    
                    for (Int_t j = k+1; j < nEntries; ++j) {
                        tree->GetEntry(j);
                        
                        if(fTrack_pt_tru!=0){pt_tru3 = fTrack_pt_tru;mj = false;}
                        else{pt_tru3 = fTrack_pt_miss; mj = true;}
                        
                        if(fTrack_pt_det!=0){pt_det3 = fTrack_pt_det;fj = false;}
                        else{pt_det3 = fTrack_pt_fake; fj = true;}
                        
                        if(fTrack_eta_tru!=0){eta_tru3 = fTrack_eta_tru;}
                        else{eta_tru3 = fTrack_eta_miss;}
                        
                        if(fTrack_eta_det!=0){eta_det3 = fTrack_eta_det;}
                        else{eta_det3 = fTrack_eta_fake;}
                        
                        if(fTrack_phi_tru!=0){phi_tru3 = fTrack_phi_tru;}
                        else{phi_tru3 = fTrack_phi_miss;}
                        
                        if(fTrack_phi_det!=0){phi_det3 = fTrack_phi_det;}
                        else{phi_det3 = fTrack_phi_fake;}
                        
                        
                        if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                            // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                            break;}
                        
                        double dR_tru = -1;
                        double dR_det = -1;
                        
                        double dR_tru_ik = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                        double dR_tru_jk = delR(phi_tru2,phi_tru3,eta_tru2,eta_tru3);
                        double dR_tru_ij = delR(phi_tru,phi_tru3,eta_tru,eta_tru3);
                        
                        double dR_det_ik = delR(phi_det,phi_det2,eta_det,eta_det2);
                        double dR_det_jk = delR(phi_det2,phi_det3,eta_det2,eta_det3);
                        double dR_det_ij = delR(phi_det,phi_det3,eta_det,eta_det3);
                        
                        if(dR_tru_ik > dR_tru_jk && dR_tru_ik > dR_tru_ij){dR_tru = dR_tru_ik;}
                        else if (dR_tru_jk > dR_tru_ik && dR_tru_jk > dR_tru_ij){dR_tru = dR_tru_jk;}
                        else{dR_tru = dR_tru_ij;}
                        
                        if(dR_det_ik > dR_det_jk && dR_det_ik > dR_det_ij){dR_det = dR_det_ik;}
                        else if (dR_det_jk > dR_det_ik && dR_det_jk > dR_det_ij){dR_det = dR_det_jk;}
                        else{dR_det = dR_det_ij;}
                        
                        if(fi || fk || fj){
                            double wt_det_ijk =(6*pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            double wt_det_3D_ijk = (pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            
                            if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                            if(dR_det < 0.01 || dR_det > 0.4) continue;
                            if(wt_det_3D_ijk < 0.000102329 || wt_det_3D_ijk > 0.15 ) continue;
                            if(abs(eta_det - eta_det2)<0.008) continue;
                            if(abs(eta_det2 - eta_det3)<0.008) continue;
                            if(abs(eta_det - eta_det3)<0.008) continue;
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                             e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ijk*mc);
                        }
                        //missed correlations//////////////////////////////////
                        else if(mi || mk || mj){
                            //  cout<<"miss: "<<endl;
                            // tree->Show(k);
                            double wt_tru_ijk= (6*pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            double wt_tru_3D_ijk = (pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            
                            if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue;
                            if(dR_tru < 0.005 || dR_tru > 0.4) continue;
                            if(wt_tru_3D_ijk < 0.00001 || wt_tru_3D_ijk > 0.15) continue;
                            if(abs(eta_tru - eta_tru2)<0.008) continue;
                            if(abs(eta_tru2 - eta_tru3)<0.008) continue;
                            if(abs(eta_tru - eta_tru3)<0.008) continue;
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            
                            e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ijk*mc);
                        }
                        //matched correlations//////////////////////////////////
                        else if(!mi && !mk && !mj && !fi && !fk && !fj){
                            
                            double wt_tru_ijk= (6*pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            double wt_tru_3D_ijk = (pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            
                            double wt_det_ijk =(6*pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            double wt_det_3D_ijk = (pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            
                            
                            //3D KE histogram TRUTH level -- applied before any truth or reco level cuts
                            if(dR_tru > 0.005 && dR_tru < 0.4 && wt_tru_3D_ijk > 0.00001 && wt_tru_3D_ijk < 0.15){
                                
                                if(abs(eta_tru - eta_tru2)<0.008) continue;
                                if(abs(eta_tru2 - eta_tru3)<0.008) continue;
                                if(abs(eta_tru - eta_tru3)<0.008) continue;
                                
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                
                                 e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ijk*mc);
                            }
                            
                            
                            //3D KE histogram DET level -- applied before any truth level cuts
                            if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                            if(dR_det < 0.01 || dR_det > 0.4) continue;
                            if(wt_det_3D_ijk < 0.000102329 || wt_det_3D_ijk > 0.15 ) continue;
                            if(abs(eta_det - eta_det2)<0.008) continue;
                            if(abs(eta_det2 - eta_det3)<0.008) continue;
                            if(abs(eta_det - eta_det3)<0.008) continue;
                            
                            
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                           
                            e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ijk*mc);
                            if (hasNaNsOrInfs(e3c_pt_hist) || hasNaNsOrInfs(e3c_pt_hist_det)) {
                                std::cerr << "Histogram eec contains NaNs or Infs!" << std::endl;
                                return;
                            }
                            
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                            if (hasNaNsOrInfs(h3_reco_e3c) || hasNaNsOrInfs(h3_true_e3c)) {
                                std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                                return;
                            }
                            
                            
                            //Apply truth level cuts before filling matched truth histogram -- this mainly avoids underflow and overflow
                            if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue;
                            if(dR_tru < 0.005 || dR_tru >0.4) continue;
                            if(wt_tru_3D_ijk < 0.00001 || wt_tru_3D_ijk > 0.15) continue;
                            
                            if(abs(eta_tru - eta_tru2)<0.008) continue;
                            if(abs(eta_tru2 - eta_tru3)<0.008) continue;
                            if(abs(eta_tru - eta_tru3)<0.008) continue;
                            
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            
                            
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            
                        }
                        else{cerr<<"something is wrong"<<endl;}
                        
                    }
                }
                
            }
        }
        

        
        if(ifTrivial){
            
            cout<<"trivial loop "<<endl;
            
            for (Int_t i = 0; i < nEntries; ++i) {
                
                tree->GetEntry(i);
                
                jetpT_tru = fJet_pt_tru;
                jetpT_det = fJet_pt_det;
                
                if(fTrack_pt_tru!=0){pt_tru = fTrack_pt_tru;  mi = false;}
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
                
                
                for (Int_t k = i+1; k < nEntries; ++k) {
                    tree->GetEntry(k);
                    
                    if(fTrack_pt_tru!=0){pt_tru2 = fTrack_pt_tru;mk = false;}
                    else{pt_tru2 = fTrack_pt_miss; mk = true;}
                    
                    if(fTrack_pt_det!=0){pt_det2 = fTrack_pt_det;fk = false;}
                    else{pt_det2 = fTrack_pt_fake; fk = true;}
                    
                    if(fTrack_eta_tru!=0){eta_tru2 = fTrack_eta_tru;}
                    else{eta_tru2 = fTrack_eta_miss;}
                    
                    if(fTrack_eta_det!=0){eta_det2 = fTrack_eta_det;}
                    else{eta_det2 = fTrack_eta_fake;}
                    
                    if(fTrack_phi_tru!=0){phi_tru2 = fTrack_phi_tru;}
                    else{phi_tru2 = fTrack_phi_miss;}
                    
                    if(fTrack_phi_det!=0){phi_det2 = fTrack_phi_det;}
                    else{phi_det2 = fTrack_phi_fake;}
                    
                    
                    if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                        // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                        break;}
                    
                    
                    //matched correlations//////////////////////////////////
                    if(!mi && !mk && !fk && !fi){
                        double dR_tru = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                        double dR_det = delR(phi_det,phi_det2,eta_det,eta_det2);
                        
                        double wt_tru_iik= (3*pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_ikk= (3*pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        double wt_tru_3D_iik = (pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_3D_ikk = (pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        double wt_det_iik =(3*pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_ikk =(3*pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        double wt_det_3D_iik = (pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_3D_ikk = (pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        
                        //                        //3D KE histogram TRUTH level -- applied before any truth or reco level cuts
                        //                        if(dR_tru > 0.01 && dR_tru < 0.4){
                        //                            if(wt_tru_3D_iik > 0.000102329 && wt_tru_3D_iik < 0.15){
                        //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                        //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                        //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                        //                            }
                        //
                        //                            if(wt_tru_3D_ikk > 0.000102329 && wt_tru_3D_ikk < 0.15){
                        //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                        //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                        //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                        //                            }
                        //                        }
                        //
                        //3D KE histogram DET level -- applied before any truth level cuts
                        if (fJet_pt_det > ptJetmin_det && fJet_pt_det < ptJetmax_det && fJet_pt_tru > ptJetmin_det && fJet_pt_tru < ptJetmax_det && dR_det > 0.01 && dR_det < 0.4 && dR_tru > 0.01 && dR_tru < 0.4){
                            if(wt_det_3D_iik > 0.000102329 && wt_det_3D_iik < 0.15 && wt_tru_3D_iik > 0.000102329 && wt_tru_3D_iik > 0.15)
                            {
                                h3_reco_e3cTriv->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_reco_e3cTriv->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_reco_e3cTriv->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                
                                h3_true_e3cTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_true_e3cTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_true_e3cTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                
                                e3c_pt_histTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_iik*mc);
                                e3c_pt_hist_detTriv->Fill(dR_det,fJet_pt_det,wt_det_iik*mc);
                                
                                response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                
                            }
                            if(wt_det_3D_ikk > 0.000102329 && wt_det_3D_ikk < 0.15 && wt_tru_3D_ikk > 0.000102329 && wt_tru_3D_ikk > 0.15){
                                h3_reco_e3cTriv->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_reco_e3cTriv->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_reco_e3cTriv->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                
                                h3_true_e3cTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_true_e3cTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_true_e3cTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                
                                e3c_pt_histTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_ikk*mc);
                                e3c_pt_hist_detTriv->Fill(dR_det,fJet_pt_det,wt_det_ikk*mc);
                                
                                response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                            }
                        }
                    }
                    
                    for (Int_t j = k+1; j < nEntries; ++j) {
                        tree->GetEntry(j);
                        
                        if(fTrack_pt_tru!=0){pt_tru3 = fTrack_pt_tru; mj = false;}
                        else{pt_tru3 = fTrack_pt_miss; mj = true;}
                        
                        if(fTrack_pt_det!=0){pt_det3 = fTrack_pt_det;fj = false;}
                        else{pt_det3 = fTrack_pt_fake; fj = true;}
                        
                        if(fTrack_eta_tru!=0){eta_tru3 = fTrack_eta_tru;}
                        else{eta_tru3 = fTrack_eta_miss;}
                        
                        if(fTrack_eta_det!=0){eta_det3 = fTrack_eta_det;}
                        else{eta_det3 = fTrack_eta_fake;}
                        
                        if(fTrack_phi_tru!=0){phi_tru3 = fTrack_phi_tru;}
                        else{phi_tru3 = fTrack_phi_miss;}
                        
                        if(fTrack_phi_det!=0){phi_det3 = fTrack_phi_det;}
                        else{phi_det3 = fTrack_phi_fake;}
                        
                        
                        if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                            // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                            break;}
                        double dR_tru = -1;
                        double dR_det = -1;
                        //matched correlations//////////////////////////////////
                        if(!mi && !mk && !mj && !fi && !fk && !fj){
                            
                            double dR_tru_ik = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                            double dR_tru_jk = delR(phi_tru2,phi_tru3,eta_tru2,eta_tru3);
                            double dR_tru_ij = delR(phi_tru,phi_tru3,eta_tru,eta_tru3);
                            
                            double dR_det_ik = delR(phi_det,phi_det2,eta_det,eta_det2);
                            double dR_det_jk = delR(phi_det2,phi_det3,eta_det2,eta_det3);
                            double dR_det_ij = delR(phi_det,phi_det3,eta_det,eta_det3);
                            
                            if(dR_tru_ik > dR_tru_jk && dR_tru_ik > dR_tru_ij){dR_tru = dR_tru_ik;}
                            else if (dR_tru_jk > dR_tru_ik && dR_tru_jk > dR_tru_ij){dR_tru = dR_tru_jk;}
                            else{dR_tru = dR_tru_ij;}
                            
                            if(dR_det_ik > dR_det_jk && dR_det_ik > dR_det_ij){dR_det = dR_det_ik;}
                            else if (dR_det_jk > dR_det_ik && dR_det_jk > dR_det_ij){dR_det = dR_det_jk;}
                            else{dR_det = dR_det_ij;}
                            
                            
                            double wt_tru_ijk= (6*pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            double wt_tru_3D_ijk = (pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            
                            double wt_det_ijk =(6*pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            double wt_det_3D_ijk = (pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            
                            
                            //                            //3D KE histogram TRUTH level -- applied before any truth or reco level cuts
                            //                            if(dR_tru > 0.01 && dR_tru < 0.4 && wt_tru_3D > 0.000102329 && wt_tru_3D < 0.15){
                            //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            //
                            //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            //                            }
                            
                            if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det || fJet_pt_tru < ptJetmin_det || fJet_pt_tru >ptJetmax_det) continue;
                            if(dR_det < 0.01 || dR_det > 0.4 || dR_tru < 0.01 || dR_tru >0.4) continue;
                            if(wt_det_3D_ijk < 0.000102329 || wt_det_3D_ijk > 0.15 || wt_tru_3D_ijk < 0.000102329 || wt_tru_3D_ijk > 0.15 ) continue;
                            
                            
                            e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ijk*mc);
                            e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ijk*mc);
                            if (hasNaNsOrInfs(e3c_pt_hist) || hasNaNsOrInfs(e3c_pt_hist_det)) {
                                std::cerr << "Histogram eec contains NaNs or Infs!" << std::endl;
                                return;
                            }
                            
                            h3_reco_e3cTriv->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3cTriv->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3cTriv->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                            h3_reco_e3cTriv->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3cTriv->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3cTriv->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                            
                            h3_true_e3cTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3cTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3cTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            
                            h3_true_e3cTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3cTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3cTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            
                            
                            response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            
                            response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            
                        }
                    }
                }
            }
        }
        
        
        
        if(ifSplit){
            for (Int_t i = 0; i < nEntries; ++i) {
                
                tree->GetEntry(i);
                
                jetpT_tru = fJet_pt_tru;
                jetpT_det = fJet_pt_det;
                
                if(fTrack_pt_tru!=0){pt_tru = fTrack_pt_tru;}
                else{pt_tru = fTrack_pt_miss; mi = true;}
                
                if(fTrack_pt_det!=0){pt_det = fTrack_pt_det;}
                else{pt_det = fTrack_pt_fake; fi = true;}
                
                if(fTrack_eta_tru!=0){eta_tru = fTrack_eta_tru;}
                else{eta_tru = fTrack_eta_miss;}
                
                if(fTrack_eta_det!=0){eta_det = fTrack_eta_det;}
                else{eta_det = fTrack_eta_fake;}
                
                if(fTrack_phi_tru!=0){phi_tru = fTrack_phi_tru;}
                else{phi_tru = fTrack_phi_miss;}
                
                if(fTrack_phi_det!=0){phi_det = fTrack_phi_det;}
                else{phi_det = fTrack_phi_fake;}
                
                
                for (Int_t k = i+1; k < nEntries; ++k) {
                    tree->GetEntry(k);
                    
                    if(fTrack_pt_tru!=0){pt_tru2 = fTrack_pt_tru;}
                    else{pt_tru2 = fTrack_pt_miss; mk = true;}
                    
                    if(fTrack_pt_det!=0){pt_det2 = fTrack_pt_det;}
                    else{pt_det2 = fTrack_pt_fake; fk = true;}
                    
                    if(fTrack_eta_tru!=0){eta_tru2 = fTrack_eta_tru;}
                    else{eta_tru2 = fTrack_eta_miss;}
                    
                    if(fTrack_eta_det!=0){eta_det2 = fTrack_eta_det;}
                    else{eta_det2 = fTrack_eta_fake;}
                    
                    if(fTrack_phi_tru!=0){phi_tru2 = fTrack_phi_tru;}
                    else{phi_tru2 = fTrack_phi_miss;}
                    
                    if(fTrack_phi_det!=0){phi_det2 = fTrack_phi_det;}
                    else{phi_det2 = fTrack_phi_fake;}
                    
                    
                    if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                        // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                        break;}
                    
                     double split = rand->Rndm();
                    //matched correlations//////////////////////////////////
                    if(!mi && !mk && !fk && !fi){
                        double dR_tru = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                        double dR_det = delR(phi_det,phi_det2,eta_det,eta_det2);
                        
                        double wt_tru_iik= (3*pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_ikk= (3*pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        double wt_tru_3D_iik = (pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_3D_ikk = (pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        double wt_det_iik =(3*pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_ikk =(3*pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        double wt_det_3D_iik = (pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_3D_ikk = (pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                       
                        
                        //                        //3D KE histogram TRUTH level -- applied before any truth or reco level cuts
                        //                        if(dR_tru > 0.01 && dR_tru < 0.4){
                        //                            if(wt_tru_3D_iik > 0.000102329 && wt_tru_3D_iik < 0.15){
                        //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                        //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                        //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                        //                            }
                        //
                        //                            if(wt_tru_3D_ikk > 0.000102329 && wt_tru_3D_ikk < 0.15){
                        //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                        //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                        //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                        //                            }
                        //                        }
                        //
                        //3D KE histogram DET level -- applied before any truth level cuts
                        if (fJet_pt_det > ptJetmin_det && fJet_pt_det < ptJetmax_det && fJet_pt_tru > ptJetmin_det && fJet_pt_tru < ptJetmax_det && dR_det > 0.01 && dR_det < 0.4 && dR_tru > 0.01 && dR_tru < 0.4){
                            if(wt_det_3D_iik > 0.000102329 && wt_det_3D_iik < 0.15 && wt_tru_3D_iik > 0.000102329 && wt_tru_3D_iik > 0.15)
                            {
                                
                                if (split < 0.5)
                                {
                                    h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                    h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                    h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                    
                                    h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                    h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                    h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                    
                                    e3c_pt_histTriv->Fill(dR_tru,fJet_pt_tru,wt_tru_iik*mc);
                                    e3c_pt_hist_detTriv->Fill(dR_det,fJet_pt_det,wt_det_iik*mc);
                                }
                                else{
                                    response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                    response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                    response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                }
                            }
                            if(wt_det_3D_ikk > 0.000102329 && wt_det_3D_ikk < 0.15 && wt_tru_3D_ikk > 0.000102329 && wt_tru_3D_ikk > 0.15){
                                
                                if (split < 0.5)
                                {
                                    h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                    h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                    h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                    
                                    h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                    h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                    h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                    
                                    e3c_pt_histSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_ikk*mc);
                                    e3c_pt_hist_detSplit->Fill(dR_det,fJet_pt_det,wt_det_ikk*mc);
                                }
                                else{
                                    response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                    response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                    response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                }
                            }
                        }
                    }
                    
                    for (Int_t j = k+1; j < nEntries; ++j) {
                        tree->GetEntry(j);
                        
                        if(fTrack_pt_tru!=0){pt_tru3 = fTrack_pt_tru; mj = false;}
                        else{pt_tru3 = fTrack_pt_miss; mj = true;}
                        
                        if(fTrack_pt_det!=0){pt_det3 = fTrack_pt_det; fj = true;}
                        else{pt_det3 = fTrack_pt_fake; fj = true;}
                        
                        if(fTrack_eta_tru!=0){eta_tru3 = fTrack_eta_tru;}
                        else{eta_tru3 = fTrack_eta_miss;}
                        
                        if(fTrack_eta_det!=0){eta_det3 = fTrack_eta_det;}
                        else{eta_det3 = fTrack_eta_fake;}
                        
                        if(fTrack_phi_tru!=0){phi_tru3 = fTrack_phi_tru;}
                        else{phi_tru3 = fTrack_phi_miss;}
                        
                        if(fTrack_phi_det!=0){phi_det3 = fTrack_phi_det;}
                        else{phi_det3 = fTrack_phi_fake;}
                        
                        
                        if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                            // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                            break;}
                        double dR_tru = -1;
                        double dR_det = -1;
                        //matched correlations//////////////////////////////////
                        if(!mi && !mk && !mj && !fi && !fk && !fj){
                            
                            double dR_tru_ik = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                            double dR_tru_jk = delR(phi_tru2,phi_tru3,eta_tru2,eta_tru3);
                            double dR_tru_ij = delR(phi_tru,phi_tru3,eta_tru,eta_tru3);
                            
                            double dR_det_ik = delR(phi_det,phi_det2,eta_det,eta_det2);
                            double dR_det_jk = delR(phi_det2,phi_det3,eta_det2,eta_det3);
                            double dR_det_ij = delR(phi_det,phi_det3,eta_det,eta_det3);
                            
                            if(dR_tru_ik > dR_tru_jk && dR_tru_ik > dR_tru_ij){dR_tru = dR_tru_ik;}
                            else if (dR_tru_jk > dR_tru_ik && dR_tru_jk > dR_tru_ij){dR_tru = dR_tru_jk;}
                            else{dR_tru = dR_tru_ij;}
                            
                            if(dR_det_ik > dR_det_jk && dR_det_ik > dR_det_ij){dR_det = dR_det_ik;}
                            else if (dR_det_jk > dR_det_ik && dR_det_jk > dR_det_ij){dR_det = dR_det_jk;}
                            else{dR_det = dR_det_ij;}
                            
                            
                            double wt_tru_ijk= (6*pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            double wt_tru_3D_ijk = (pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            
                            double wt_det_ijk =(6*pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            double wt_det_3D_ijk = (pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            
                            
                            //                            //3D KE histogram TRUTH level -- applied before any truth or reco level cuts
                            //                            if(dR_tru > 0.01 && dR_tru < 0.4 && wt_tru_3D > 0.000102329 && wt_tru_3D < 0.15){
                            //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            //
                            //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            //                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            //                            }
                            
                            if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det || fJet_pt_tru < ptJetmin_det || fJet_pt_tru >ptJetmax_det) continue;
                            if(dR_det < 0.01 || dR_det > 0.4 || dR_tru < 0.01 || dR_tru >0.4) continue;
                            if(wt_det_3D_ijk < 0.000102329 || wt_det_3D_ijk > 0.15 || wt_tru_3D_ijk < 0.000102329 || wt_tru_3D_ijk > 0.15 ) continue;
                            
                            
                            e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ijk*mc);
                            e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ijk*mc);
                            if (hasNaNsOrInfs(e3c_pt_hist) || hasNaNsOrInfs(e3c_pt_hist_det)) {
                                std::cerr << "Histogram eec contains NaNs or Infs!" << std::endl;
                                return;
                            }
                            
                            if(split<0.5){
                                h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                                h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                                h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                                
                                h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                                h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                                h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                                
                                
                                h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                
                                h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                
                                e3c_pt_histSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_ijk*mc);
                                e3c_pt_hist_detSplit->Fill(dR_det,fJet_pt_det,wt_det_ijk*mc);
                            }
                            else{
                                response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                                response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                                response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                                
                                response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                                response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                                response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            }
                        }
                        
                    }
                }
            }
        }
        
        
     if(ifFastSim){
            cout<<"fast sim loop"<<endl;
            
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
    double chpt_tru3 = 0;
    double chpt_det3 = 0;
            
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
                
                 // Apply really wide cuts on true jet pT///THIS LINE IS NEW -- CHECK THIS WITH LAURA
                 if (fJet_pt_tru > (ptJetmax_tru) || fJet_pt_tru < 0) continue;
                
                jetpT_tru = fJet_pt_tru;
                jetpT_det = fJet_pt_det;
                
                if(fTrack_pt_tru!=0){pt_tru = fTrack_pt_tru; chpt_tru = fTrack_choverpt_tru; mi = false;}
                else{pt_tru = fTrack_pt_miss; chpt_tru = fTrack_choverpt_miss; mi = true;}
                
                if(fTrack_pt_det!=0){pt_det = fTrack_pt_det; chpt_det = fTrack_choverpt_det; fi = false;}
                else{pt_det = fTrack_pt_fake;chpt_det = fTrack_choverpt_fake; fi = true;}
                
                if(fTrack_eta_tru!=0){eta_tru = fTrack_eta_tru;}
                else{eta_tru = fTrack_eta_miss;}
                
                if(fTrack_eta_det!=0){eta_det = fTrack_eta_det;}
                else{eta_det = fTrack_eta_fake;}
                
                if(fTrack_phi_tru!=0){phi_tru = fTrack_phi_tru;}
                else{phi_tru = fTrack_phi_miss;}
                
                if(fTrack_phi_det!=0){phi_det = fTrack_phi_det;}
                else{phi_det = fTrack_phi_fake;}
                
                
                for (Int_t k = i+1; k < nEntries; ++k) {
                    tree->GetEntry(k);
                    
                    if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                        // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                        break;}
                        
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
                    
                    
                    double dR_tru = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                    double dR_det = delR(phi_det,phi_det2,eta_det,eta_det2);
                    
                    double chptdiff_ik = abs(chpt_det - chpt_det2);
                    double pairwt_ik = paireffwt(_h_n2_qpt_det, _h_n2_qpt_tru, chptdiff_ik, dR_det);
                    
                    mc*=pairwt_ik;
                    //fake correlations//////////////////////////////////
                    if(fi || fk){
                        //     cout<<"fake: "<<endl;
                        //  tree->Show(k);
                        double wt_det_iik =(3*pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_ikk =(3*pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        double wt_det_3D_iik = (pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_3D_ikk = (pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                        if(dR_det > 0.01 && dR_det < 0.4){
                            if(wt_det_3D_iik > 0.000102329 && wt_det_3D_iik < 0.15 ){
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                
                                e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_iik*mc);
                            }
                            if(wt_det_3D_ikk > 0.000102329 && wt_det_3D_ikk < 0.15 ){
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                
                                e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ikk*mc);
                            }
                        }
                    }
                    //missed correlations//////////////////////////////////
                    else if(mi || mk){
                        //  cout<<"miss: "<<endl;
                        // tree->Show(k);
                        double wt_tru_iik= (3*pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_ikk= (3*pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        double wt_tru_3D_iik = (pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_3D_ikk = (pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue;
                        //was dR_det --fixed on Dec 18
                        if(dR_tru > 0.01 && dR_tru < 0.4){
                            if(wt_tru_3D_iik > 0.00001 && wt_tru_3D_iik < 0.15 ){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                
                                e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_iik*mc);
                            }
                            if(wt_tru_3D_ikk > 0.00001 && wt_tru_3D_ikk < 0.15 ){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                
                                e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ikk*mc);
                            }
                        }
                    }
                    //matched correlations//////////////////////////////////
                    else if(!mi && !mk && !fk && !fi){
                        
                        double wt_tru_iik= (3*pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_ikk= (3*pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        double wt_tru_3D_iik = (pt_tru*pt_tru*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        double wt_tru_3D_ikk = (pt_tru*pt_tru2*pt_tru2)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                        
                        double wt_det_iik =(3*pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_ikk =(3*pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        double wt_det_3D_iik = (pt_det*pt_det*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        double wt_det_3D_ikk = (pt_det*pt_det2*pt_det2)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                        
                        
                        //3D KE histogram TRUTH level -- applied before any truth or reco level cuts
                        if(dR_tru > 0.005 && dR_tru < 0.4){
                            if(wt_tru_3D_iik > 0.00001 && wt_tru_3D_iik < 0.15){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                
                                e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_iik*mc);
                            }
                            
                            if(wt_tru_3D_ikk > 0.00001 && wt_tru_3D_ikk < 0.15){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                
                                e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ikk*mc);
                            }
                        }
                        
                        //3D KE histogram DET level -- applied before any truth level cuts
                        if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                        
                        //because even if this condition is not met for ik , it might be met for ijk
                        if(dR_det > 0.01 && dR_det < 0.4){
                            if(wt_det_3D_iik > 0.000102329 && wt_det_3D_iik < 0.15 ){
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                
                                
                                e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_iik*mc);
                                if (hasNaNsOrInfs(e3c_pt_hist) || hasNaNsOrInfs(e3c_pt_hist_det)) {
                                    std::cerr << "Histogram e3c contains NaNs or Infs!" << std::endl;
                                    return;
                                }
                                
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_iik,mc);
                                
                                if (hasNaNsOrInfs(h3_reco_e3c) || hasNaNsOrInfs(h3_true_e3c)) {
                                    std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                                    return;
                                }
                                
                               
                                //Apply truth level cuts before filling matched truth histogram -- this mainly avoids underflow and overflow
                                if(fJet_pt_tru > ptJetmin_tru && fJet_pt_tru < ptJetmax_tru && dR_tru > 0.005 && dR_tru < 0.4){
                                    if(wt_tru_3D_iik > 0.00001 && wt_tru_3D_iik < 0.15) {
                                        
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_iik, mc);
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_iik, mc);
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_iik, mc);
                                        
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_iik,mc);
                                        
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_iik, dR_tru, fJet_pt_tru, wt_tru_3D_iik, mc);
                                    }
                                }
                            }
                            if(wt_det_3D_ikk > 0.000102329 && wt_det_3D_ikk < 0.15 ){
                                
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                
                                
                                e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ikk*mc);
                                if (hasNaNsOrInfs(e3c_pt_hist) || hasNaNsOrInfs(e3c_pt_hist_det)) {
                                    std::cerr << "Histogram eec contains NaNs or Infs!" << std::endl;
                                    return;
                                }
                                
                                
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ikk,mc);
                                
                                if (hasNaNsOrInfs(h3_reco_e3c) || hasNaNsOrInfs(h3_true_e3c)) {
                                    std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                                    return;
                                }
                                
                                if(fJet_pt_tru > ptJetmin_tru && fJet_pt_tru < ptJetmax_tru && dR_tru > 0.005 && dR_tru < 0.4){
                                    if(wt_tru_3D_ikk > 0.00001 && wt_tru_3D_ikk < 0.15){
                                        
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, mc);
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, mc);
                                        h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, mc);
                                        
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                        h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ikk,mc);
                                        
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                        response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ikk, dR_tru, fJet_pt_tru, wt_tru_3D_ikk, mc);
                                        
                                    }
                                }
                            }
                        }
                    }
                    else{cerr<<"something is wrong"<<endl;
                    }
                    
                    for (Int_t j = k+1; j < nEntries; ++j) {
                        tree->GetEntry(j);
                        
                        if(fTrack_pt_tru!=0){pt_tru3 = fTrack_pt_tru; chpt_tru3 = fTrack_choverpt_tru; mj = false;}
                        else{pt_tru3 = fTrack_pt_miss; chpt_tru3 = fTrack_choverpt_miss; mj = true;}
                        
                        if(fTrack_pt_det!=0){pt_det3 = fTrack_pt_det; chpt_det3 = fTrack_choverpt_det; fj = false;}
                        else{pt_det3 = fTrack_pt_fake; chpt_det3 = fTrack_choverpt_fake; fj = true;}
                        
                        if(fTrack_eta_tru!=0){eta_tru3 = fTrack_eta_tru;}
                        else{eta_tru3 = fTrack_eta_miss;}
                        
                        if(fTrack_eta_det!=0){eta_det3 = fTrack_eta_det;}
                        else{eta_det3 = fTrack_eta_fake;}
                        
                        if(fTrack_phi_tru!=0){phi_tru3 = fTrack_phi_tru;}
                        else{phi_tru3 = fTrack_phi_miss;}
                        
                        if(fTrack_phi_det!=0){phi_det3 = fTrack_phi_det;}
                        else{phi_det3 = fTrack_phi_fake;}
                        
                        
                        if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru){
                            // cout<<"DIFFERENT JET BEGINS, BREAK OUT OF LOOP"<<endl;
                            break;}
                        
                        double dR_tru = -1;
                        double dR_det = -1;
                        
                        double dR_tru_ik = delR(phi_tru,phi_tru2,eta_tru,eta_tru2);
                        double dR_tru_jk = delR(phi_tru2,phi_tru3,eta_tru2,eta_tru3);
                        double dR_tru_ij = delR(phi_tru,phi_tru3,eta_tru,eta_tru3);
                        
                        double dR_det_ik = delR(phi_det,phi_det2,eta_det,eta_det2);
                        double dR_det_jk = delR(phi_det2,phi_det3,eta_det2,eta_det3);
                        double dR_det_ij = delR(phi_det,phi_det3,eta_det,eta_det3);
                        
                        double chptdiff_ij = abs(chpt_det - chpt_det3);
                        double chptdiff_jk = abs(chpt_det2 - chpt_det3);
                        
                        double pairwt_ij = paireffwt(_h_n2_qpt_det, _h_n2_qpt_tru, chptdiff_ij, dR_det_ij);
                        double pairwt_jk = paireffwt(_h_n2_qpt_det, _h_n2_qpt_tru, chptdiff_jk, dR_det_jk);
                    
                         mc*=(pairwt_ij*pairwt_jk);
                        
                        if(dR_tru_ik > dR_tru_jk && dR_tru_ik > dR_tru_ij){dR_tru = dR_tru_ik;}
                        else if (dR_tru_jk > dR_tru_ik && dR_tru_jk > dR_tru_ij){dR_tru = dR_tru_jk;}
                        else{dR_tru = dR_tru_ij;}
                        
                        if(dR_det_ik > dR_det_jk && dR_det_ik > dR_det_ij){dR_det = dR_det_ik;}
                        else if (dR_det_jk > dR_det_ik && dR_det_jk > dR_det_ij){dR_det = dR_det_jk;}
                        else{dR_det = dR_det_ij;}
                        
                        if(fi || fk || fj){
                            double wt_det_ijk =(6*pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            double wt_det_3D_ijk = (pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            
                            if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                            if(dR_det < 0.01 || dR_det > 0.4) continue;
                            if(wt_det_3D_ijk < 0.000102329 || wt_det_3D_ijk > 0.15 ) continue;
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                             e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ijk*mc);
                        }
                        //missed correlations//////////////////////////////////
                        else if(mi || mk || mj){
                            //  cout<<"miss: "<<endl;
                            // tree->Show(k);
                            double wt_tru_ijk= (6*pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            double wt_tru_3D_ijk = (pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            
                            if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue;
                            if(dR_tru < 0.005 || dR_tru > 0.4) continue;
                            if(wt_tru_3D_ijk < 0.00001 || wt_tru_3D_ijk > 0.15) continue;
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            
                            e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ijk*mc);
                        }
                        //matched correlations//////////////////////////////////
                        else if(!mi && !mk && !mj && !fi && !fk && !fj){
                            
                            double wt_tru_ijk= (6*pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            double wt_tru_3D_ijk = (pt_tru*pt_tru2*pt_tru3)/(fJet_pt_tru*fJet_pt_tru*fJet_pt_tru);
                            
                            double wt_det_ijk =(6*pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            double wt_det_3D_ijk = (pt_det*pt_det2*pt_det3)/(fJet_pt_det*fJet_pt_det*fJet_pt_det);
                            
                            
                            //3D KE histogram TRUTH level -- applied before any truth or reco level cuts
                            if(dR_tru > 0.005 && dR_tru < 0.4 && wt_tru_3D_ijk > 0.00001 && wt_tru_3D_ijk < 0.15){
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                h3_fulleff->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                                
                                 e3c_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru_ijk*mc);
                            }
                            
                            
                            //3D KE histogram DET level -- applied before any truth level cuts
                            if (fJet_pt_det < ptJetmin_det || fJet_pt_det > ptJetmax_det) continue;
                            if(dR_det < 0.01 || dR_det > 0.4) continue;
                            if(wt_det_3D_ijk < 0.000102329 || wt_det_3D_ijk > 0.15 ) continue;
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                           
                            e3c_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det_ijk*mc);
                            if (hasNaNsOrInfs(e3c_pt_hist) || hasNaNsOrInfs(e3c_pt_hist_det)) {
                                std::cerr << "Histogram eec contains NaNs or Infs!" << std::endl;
                                return;
                            }
                            
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            h3_reco_e3c->Fill(dR_det,fJet_pt_det,wt_det_3D_ijk,mc);
                            
                            if (hasNaNsOrInfs(h3_reco_e3c) || hasNaNsOrInfs(h3_true_e3c)) {
                                std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                                return;
                            }
                            
                            
                            //Apply truth level cuts before filling matched truth histogram -- this mainly avoids underflow and overflow
                            if(fJet_pt_tru < ptJetmin_tru || fJet_pt_tru >ptJetmax_tru) continue;
                            if(dR_tru < 0.005 || dR_tru >0.4) continue;
                            if(wt_tru_3D_ijk < 0.00001 || wt_tru_3D_ijk > 0.15) continue;
                            
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            h3_smeared->Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, mc);
                            
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            h3_true_e3c->Fill(dR_tru,fJet_pt_tru,wt_tru_3D_ijk,mc);
                            
                            
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            response3D.Fill(dR_det, fJet_pt_det, wt_det_3D_ijk, dR_tru, fJet_pt_tru, wt_tru_3D_ijk, mc);
                            
                        }
                        else{cerr<<"something is wrong"<<endl;}
                        
                    }
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
    
    e3c_pt_hist->Write();
    e3c_pt_hist_det->Write();
    e3c_pt_histTriv->Write();
    e3c_pt_hist_detTriv->Write();
    e3c_pt_histSplit->Write();
    e3c_pt_hist_detSplit->Write();

    h1_true->Write();
    h1_reco->Write();
    h1_fulleff->Write();
    h3_true_e3c->Write();
    h3_reco_e3c->Write();
    h3_true_e3cTriv->Write();
    h3_reco_e3cTriv->Write();
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
    
    cout<<"Response matrices written. Success!"<<endl;
    fout->Close();
    }
   else{
      cout<<"0 ENTRIES!"<<endl;
    
  }  
    
     
   
 
    fmc->Close();
     
 

}
#ifndef __CINT__
int main () { Unfolding3dAnalysisOptNewE3C(); return 0; }  // Main program when run stand-alone
#endif
