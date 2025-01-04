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


void Unfolding3dAnalysisOpt(std::string tag = "", int pthardbin = -1, int startEntries = -1, int endEntries = -1) {

    //***************************************************
    TH1D* h1_true; TH1D* h1_trueSplit;
    TH1D* h1_reco;
    TH1D* h1_unfolded; TH1D* h1_unfolded_rw;  TH1D* h1_unfoldedSplit;
    TH1D* h1_raw; TH1D* h1_rawSplit;
    TH1D* h1_smeared;
    TH1D* h1_fulleff_match;
    TH1D* h1_eff_match;
    TH1D* h1_fullreco;
    TH1D* h1_rw;
    TH1* h1_fold;
    TH1* h1_triv;
    
    TH3D* h3_true_eec; TH3D* h3_trueSplit;
    TH3D* h3_reco_eec;
    TH3D* h3_true_eec_triv;
    TH3D* h3_reco_eec_triv;
    TH3D* h3_unfolded; TH3D* h3_unfolded_rw;  TH3D* h3_unfoldedSplit;
    TH3D* h3_raw; TH3D* h3_rawSplit;
    TH3D* h3_smeared;
    TH3D* h3_fulleff_match;
    TH3D* h3_eff_match;
    TH3D* h3_fullreco;
    TH3D* h3_rw;
    TH3* h3_fold;
    TH3* h3_triv;
    
    TH2D* eec_pt_hist;
    TH2D* eec_pt_hist_det;
    
    TFile *fmc;

    RooUnfoldResponse response1D;
    RooUnfoldResponse response1DSplit;
    
    RooUnfoldResponse response3D;
    RooUnfoldResponse response3DTriv;
    RooUnfoldResponse response3DSplit;
   
    float pThard_val[21] = {5., 7., 9., 12., 16., 21., 28., 36., 45., 57., 70., 85., 99., 115., 132., 150., 169., 190., 212., 235., 1000.};
    
    TRandom3* rand = new TRandom3(0);
    
     //reco jet pT bins for 1D unfolding 
    Double_t from_const_reco_jet = 10;
    Double_t to_const_reco_jet = 120;
    Int_t bins_const_reco_jet = 22;
    Double_t width_const_reco_jet = (to_const_reco_jet-from_const_reco_jet)/bins_const_reco_jet;
    Double_t xbins_jet[23] = {};
    for (Int_t i = 0; i <= bins_const_reco_jet; i++)
    {
        xbins_jet[i] = (from_const_reco_jet + i * width_const_reco_jet);
    }
    
     //true jet pT bins for 1D unfolding 
    Double_t from_const_jet = 0;
    Double_t to_const_jet= 200;
    Int_t bins_const_jet = 40;
    Double_t width_const_jet = (to_const_jet-from_const_jet)/bins_const_jet;
    Double_t tbins_jet[41] = {};
    for (Int_t i = 0; i <= bins_const_jet; i++)
    {
        tbins_jet[i] = (from_const_jet + i * width_const_jet);
    }
    
    //jet pT bins for 3D unfolding 
    Double_t from_const_reco = 10;
    Double_t to_const_reco = 120;
    Int_t bins_const_reco = 11;
    Double_t width_const_reco = (to_const_reco-from_const_reco)/bins_const_reco;
    Double_t xbins[12] = {};
    for (Int_t i = 0; i <= bins_const_reco; i++)
    {
        xbins[i] = (from_const_reco + i * width_const_reco);
    }
    
    
    //true jet pT bins for 3D unfolding 
    Double_t from_const = 0;
    Double_t to_const = 200;
    Int_t bins_const = 20;
    //  Double_t to_const = 350;
    // Int_t bins_const = 35;
    Double_t width_const = (to_const-from_const)/bins_const;
    Double_t tbins[21] = {};
    for (Int_t i = 0; i <= bins_const; i++)
    {
        tbins[i] = (from_const + i * width_const);
    }
    
    
    //reco dR bins --same for eec and e3c
    double dRbins[] = {
        0.0001,
        0.01,
        0.0158489, 0.0190546, 0.0229087,
        0.0275423, 0.0331131,
        0.0398107, 0.0436516, 0.047863, 0.0524807, 0.057544,
        0.0630957, 0.0691831, 0.0758578, 0.0831764, 0.0912011,
        0.1, 0.109648, 0.120226, 0.131826, 0.144544,
        0.158489, 0.17378, 0.190546, 0.20893, 0.229087,
        0.251189, 0.275423, 0.301995, 0.331131, 0.363078,
        0.398107, 0.524807,
        1.0
    };
    
    //true dR bins --same for eec and e3c
    double dRbinst[] = {
        0.0001,
        0.01,
        0.0158489, 0.0190546, 0.0229087,
        0.0275423, 0.0331131,
        0.0398107, 0.0436516, 0.047863, 0.0524807, 0.057544,
        0.0630957, 0.0691831, 0.0758578, 0.0831764, 0.0912011,
        0.1, 0.109648, 0.120226, 0.131826, 0.144544,
        0.158489, 0.17378, 0.190546, 0.20893, 0.229087,
        0.251189, 0.275423, 0.301995, 0.331131, 0.363078,
        0.398107, 0.524807,
        1.0
    };
    
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
    
    
    eec_pt_hist = new TH2D("eec_pt_hist", "EEC and jet_pt 2D",  33, dRbinst, nJetPtbinst, tbins);
    eec_pt_hist_det = new TH2D("eec_pt_hist_det", "EEC and jet_pt 2D det", 33, dRbins, nJetPtbins,xbins);
    
    
    h1_true = new TH1D("h1_true", "h1_true", bins_const_jet, tbins_jet);
    h1_reco = new TH1D("h1_reco", "h1_reco", bins_const_reco_jet, xbins_jet);

    
   cout<<"nWtbinst "<<nWtbinst<<" and nWtbins "<<nWtbins<<endl;
   cout<<"ndRbinst "<<ndRbinst<<" and ndRbins "<< ndRbins <<endl;
   cout<<"nJetPtbinst "<<nJetPtbinst<<" and nJetPtbins "<<nJetPtbins<<endl;    
   
   cout<<"nJetPtbins1Dtrue "<<nJetPtbins1Dtrue<<" and nJetPtbins1Dreco "<<nJetPtbins1Dreco<<endl;
   
  
    
    h3_true_eec = new TH3D("h3_true_eec", "h3_true_eec", ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_reco_eec = new TH3D("h3_reco_eec", "h3_reco_eec",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
    
    h3_true_eec_triv = new TH3D("h3_true_eecTriv", "h3_true_eec_triv", ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_reco_eec_triv = new TH3D("h3_reco_eecTriv", "h3_reco_eec_triv",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
    
    h3_trueSplit = new TH3D("h3_trueSplit", "h3_true_split", ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_rawSplit = new TH3D("h3_rawSplit", "h3_reco_split",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
   
    h3_fulleff_match = new TH3D("h3_fulleff_match","h3_fulleff_match",ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_eff_match = new TH3D("h3_eff_match","h3_eff_match",ndRbinst, dRbinst, nJetPtbinst, tbins, nWtbinst,wtbinst);
    h3_fullreco = new TH3D("h3_fullreco","h3_fullreco",ndRbins, dRbins, nJetPtbins, xbins, nWtbins,wtbins);
    
    
    h3_true_eec->Sumw2();
    h3_reco_eec->Sumw2();
    
    h3_true_eec_triv->Sumw2();
    h3_reco_eec_triv->Sumw2();
    
    h3_rawSplit->Sumw2();
    h3_rawSplit->Sumw2();
  
    h3_eff_match->Sumw2();
    h3_fulleff_match->Sumw2();
    h3_fullreco->Sumw2();
  
   cout<<"HERE"<<endl;
  
    cout<<"setting up response declarations"<<endl;
    response3D.Setup(h3_reco_eec,h3_true_eec);
    response3DTriv.Setup(h3_reco_eec_triv,h3_true_eec_triv);
    response3DSplit.Setup(h3_rawSplit,h3_trueSplit);
    response1D.Setup(h1_reco,h1_true);

    cout<<"ending declarations"<<endl;
    

    // pthard = i;
      //Get the mc files
     fmc = TFile::Open(Form("/home/ar2545/Trial/AnalysisResultsMC%i.root",pthardbin));
    
    TString file_name = Form("/home/ar2545/Trial/AnalysisResultsMC%i.root", pthardbin);

    cout<<"file name is: "<<file_name.Data()<<endl;

    TTree *tree=(TTree*)fmc->Get("MatchTracksTree");
    tree->Show(0);
  
    Int_t nEv=tree->GetEntries();
   
    Double_t fJet_pt_det, fJet_pt_tru;
    Double_t fTrack_pt_det, fTrack_pt_tru, fTrack_pt_miss, fTrack_pt_fake;
    Double_t fTrack_eta_det, fTrack_eta_tru, fTrack_eta_miss, fTrack_eta_fake;
    Double_t fTrack_phi_det, fTrack_phi_tru, fTrack_phi_miss, fTrack_phi_fake;
    Double_t mc_weight;
    
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
    tree->SetBranchAddress("mc_weight", &mc_weight);
    
    // tree->GetEntry(1);
    // double mc = mc_weight; //same for one run in 1 pTHat bin
    // double weightresp = mc;
    
    double ptJetmin_det = 10.; //min det level pt
    double ptJetmax_det = 90.; //max det level pt
    double ptJetmin_tru3D = 0.; //min tru level pt for 3D unfolding
    double ptJetmax_tru3D = to_const; //max tru level pt for 3D unfolding
    
    
    //for 1D unfolding
    double ptJetmin_tru = 10.; //min tru level pt
    double ptJetmax_tru = 120.; //max tru level pt
    
    
    double jetpT_tru = -1;
    double jetpT_det = -1;
    
   
    //In this tree I only have matched jets, so tree entries where either jet pT is 0 is only for missed(true level) or fake tracks(det level)
     Long64_t nEntries = tree->GetEntries();
    cout<<nEntries<<endl;


    double pt_tru;
    double pt_det;
    double pt_miss;
    double pt_fake;

    double eta_tru;
    double eta_det;
    double eta_miss;
    double eta_fake;

    double phi_tru;
    double phi_det;
    double phi_miss;
    double phi_fake;

    //hard-coding these from a crossx file I had written
    int hat = pthardbin-1;
    double mcarray[20] = {1.2035e-05,4.18184e-06,2.06124e-06,8.884e-07,3.0596e-07,1.22758e-07,3.98291e-08,1.59277e-08,6.27856e-09,2.41448e-09,9.90233e-10,3.64615e-10,2.03901e-10,1.0315e-10,5.21154e-11,2.95027e-11,1.60933e-11,9.73774e-12,5.62143e-12,9.658e-12};
    cout<<hat<<endl;
    cout<<"MC weight is: "<<mcarray[hat]<<setprecision(7)<< endl;
    // Loop over entries in the tree
    double mc = mcarray[hat];
    for (Int_t i = startEntries; i < endEntries; ++i) {
    
        tree->GetEntry(i);
        
        // matched correlations (if either is 0, then its a fake or missed track)
        //don't need to jet match here since this is the first loop
        if(fJet_pt_det==0 || fJet_pt_tru==0) continue;
        // cout<<"i "<<i<<endl;
        
        //want to fill the jet histogram only once for each jet
         if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru) {
         if (fJet_pt_det>ptJetmin_det && fJet_pt_det<ptJetmax_det && fJet_pt_tru<ptJetmax_tru && fJet_pt_tru>ptJetmin_tru){//putting jet pT cuts for unfolding data
                h1_reco->Fill(fJet_pt_det, mc);
                h1_true->Fill(fJet_pt_tru, mc); 
                response1D.Fill(fJet_pt_det, fJet_pt_tru, mc);
         }
         }
       

        // double mc = mc_weight;
        jetpT_tru = fJet_pt_tru;
        jetpT_det = fJet_pt_det;

        pt_tru = fTrack_pt_tru;
        pt_det = fTrack_pt_det;
        pt_miss = fTrack_pt_miss;
        pt_fake = fTrack_pt_fake;

        eta_tru = fTrack_eta_tru;
        eta_det = fTrack_eta_det;
        eta_miss = fTrack_eta_miss;
        eta_fake = fTrack_eta_fake;

        phi_tru = fTrack_phi_tru;
        phi_det = fTrack_phi_det;
        phi_miss = fTrack_phi_miss;
        phi_fake = fTrack_phi_fake;
        
      
        
        // cout<<"particle I "<<i<< " : true jet "<<fJet_pt_tru<<" and det jet "<<fJet_pt_det<<endl;
        
       
            //2-pt corr
            for (Int_t k = i+1; k < endEntries; ++k) {
                tree->GetEntry(k);
                // cout<<"particle k "<<k<< " : true jet "<<fJet_pt_tru<<" and det jet "<<fJet_pt_det<<endl;
               
                //if both true and detector jet pT change then we have moved on to a new jet
                //because if one of them changes then that just means that its a missed or a fake track in the same jet
                if(jetpT_det!=fJet_pt_det && jetpT_tru!=fJet_pt_tru) break; 
                
                //  cout<<"After IF statement on line 329-particle k: true jet "<<fJet_pt_tru<<" and det jet "<<fJet_pt_det<<endl;
                // cout<<" i k "<<i<< " "<<k<<endl;
                // cout<<" idet kdet "<<jetpT_det<< " "<<fJet_pt_det<<endl;
                
                 //this ensures matched correlations
                if(fTrack_pt_tru==0 || fTrack_pt_det==0) continue;
                
                    
                    double wt_tru = (2*pt_tru*fTrack_pt_tru)/(fJet_pt_tru*fJet_pt_tru);
                    double wt_det =(2*pt_det*fTrack_pt_det)/(fJet_pt_det*fJet_pt_det);
                    
                    double dR_tru = delR(phi_tru,fTrack_phi_tru,eta_tru,fTrack_eta_tru);
                    double dR_det = delR(phi_det,fTrack_phi_det,eta_det,fTrack_eta_det);
                    
                    double wt_tru_3D = (pt_tru*fTrack_pt_tru)/(fJet_pt_tru*fJet_pt_tru); //for filling 3D histogram
                    double wt_det_3D = (pt_det*fTrack_pt_det)/(fJet_pt_det*fJet_pt_det);
                    
                  if(dR_tru < 0.01 && dR_det<0.01) continue;
                  if(wt_det_3D<0.0005 && wt_tru_3D < 0.0001){ continue;} //if both are outside the range
                  if(fJet_pt_det>120 && fJet_pt_tru > 200){ continue;} //if both are outside the range
                  
                //   if(wt_det_3D<0.0005){ continue;} //if det is outside the range
                //   if(fJet_pt_det>120){ continue;} //if det is outside the range
                  
                  
                //   if(wt_tru_3D < 0.0001) continue; //if both are outside the range
                    
                    // if(wt_det_3D<0.0005){cout<<"less than 0.0005 filled"<<endl;}
                    
                    eec_pt_hist->Fill(dR_tru,fJet_pt_tru,wt_tru*mc);
                    eec_pt_hist_det->Fill(dR_det,fJet_pt_det,wt_det*mc);
                    if (hasNaNsOrInfs(eec_pt_hist) || hasNaNsOrInfs(eec_pt_hist_det)) {
                         std::cerr << "Histogram eec contains NaNs or Infs!" << std::endl;
                         return;
                        }
                        
                    // h3_reco_eec->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                    // h3_reco_eec->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                    // h3_true_eec->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                    // h3_true_eec->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                    
                     if (hasNaNsOrInfs(h3_reco_eec) || hasNaNsOrInfs(h3_true_eec)) {
                         std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                         return;
                        }
                    
                    // h3_reco_eec_triv->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                    // h3_reco_eec_triv->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                    // h3_true_eec_triv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                    // h3_true_eec_triv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                   
                    //       //For actual data
                    if (fJet_pt_det>=ptJetmin_det && fJet_pt_det<=ptJetmax_det && fJet_pt_tru<=ptJetmax_tru3D && fJet_pt_tru>ptJetmin_tru3D){ //putting jet pT cuts 
                    // if(dR_tru>0.4 || dR_det>0.4 || dR_tru<0.01 || dR_det<0.01) continue;
                    // if((wt_tru_3D)>0.25||(wt_det_3D)>0.25||(wt_tru_3D)<0.0001||(wt_det_3D)<0.0005) continue;
                    
                    if(dR_tru<=0.4 && dR_det<=0.4 && dR_tru>=0.01 && dR_det>=0.01 && wt_tru_3D<=0.25 && wt_det_3D<=0.25 && wt_tru_3D>=0.0001 && wt_det_3D>=0.0005 ){
                    // if(dR_tru<=0.4 && dR_det<=0.4 && dR_tru>=0.01 && dR_det>=0.01 && wt_tru_3D<=0.25 && wt_det_3D<=0.25 && wt_tru_3D>=0.0005 && wt_det_3D>=0.0005 ){
                 
                    h3_reco_eec->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                    h3_reco_eec->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                    h3_true_eec->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                    h3_true_eec->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                    //  if (hasNaNsOrInfs(h3_reco_eec) || hasNaNsOrInfs(h3_true_eec)) {
                    //      std::cerr << "Histogram h3 contains NaNs or Infs!" << std::endl;
                    //      return;
                    //     }
                        
                    h3_reco_eec_triv->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                    h3_reco_eec_triv->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                    h3_true_eec_triv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                    h3_true_eec_triv->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                    
                    
                    response3D.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
                    response3D.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
              
        //         //For Trivial Closure
                    response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
                    response3DTriv.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
                    
        //         //for MC split test
                    double split = rand->Rndm();
	                if (split < 0.5)
	                {
	          
	                response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
	                response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
	                response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
	                response3DSplit.Fill(dR_det, fJet_pt_det, wt_det_3D, dR_tru, fJet_pt_tru, wt_tru_3D, mc);
	            
	                }
	             else
	                 {
	                    
	                h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
	                h3_rawSplit->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                    h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                    h3_trueSplit->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
	                }
                    
                  }   
                }
                  
                  //Kinematic efficiency
                  //fill all true values for the given detector range
                   if (fJet_pt_det>=10 && fJet_pt_det<=90 && wt_det_3D>=0.0005 && wt_det_3D<=0.25 && dR_det<=0.4 && dR_det>=0.01) 
                    {
                        h3_eff_match->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                        h3_eff_match->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                    }
                    
                     if (fJet_pt_det>=10 && wt_det_3D>=0.0005 && dR_det>=0.01) 
                    {
                        h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                        h3_fullreco->Fill(dR_det,fJet_pt_det,wt_det_3D,mc);
                    }
                    
                     //fill all true values for the given detector range
                  if (fJet_pt_det>=10 && fJet_pt_det<=120 && wt_det_3D>=0.0005 && dR_det>=0.01)
                     {
                        
                        h3_fulleff_match->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                        h3_fulleff_match->Fill(dR_tru,fJet_pt_tru,wt_tru_3D,mc);
                     }
                
            }
            
        //   if (fJet_pt_det>=ptJetmin_det && fJet_pt_det<=ptJetmax_det && fJet_pt_tru<=to_const)//this is the correct thing to do
         
            // if (jetpT_det<ptJetmin_det || jetpT_det>ptJetmax_det || jetpT_tru<=ptJetmin_det ||  jetpT_tru>ptJetmax_det) continue;//putting jet pT cuts for Trivial closure
            //  if (jetpT_det<ptJetmin_det || jetpT_det>ptJetmax_det || jetpT_tru>ptJetmax_tru || jetpT_tru<ptJetmin_tru) continue;//putting jet pT cuts for unfolding data
            //     h1_reco->Fill(jetpT_det, mc);
            //     h1_true->Fill(jetpT_tru, mc); 
                
            //     response1D.Fill(jetpT_det, jetpT_tru, mc);

     
        } 
    

    
  
    TFile *fout=new TFile (Form("/home/ar2545/Trial/Unfold3DData_%i_%i.root",pthardbin,startEntries),"RECREATE");
    
    TString foutname = Form("/home/ar2545/Trial/Unfold3DData_%i_%i.root", pthardbin,startEntries);
    cout<<"writing to output file named "<< foutname<<endl;
    
    
    cout<<"write histograms"<<endl;

    h1_true->Write();
    h1_reco->Write();
    h3_true_eec->Write();
    h3_reco_eec->Write();
    h3_true_eec_triv->Write();
    h3_reco_eec_triv->Write();
    eec_pt_hist->Write();
    eec_pt_hist_det->Write();
    h3_rawSplit->Write();
    h3_trueSplit->Write();
    h3_fulleff_match->Write();
    h3_fullreco->Write();
    h3_eff_match->Write();

    cout<<"Writing response matrices"<<endl;
    response3D.Write();
    response3DSplit.Write();
    response3DTriv.Write();
    response1D.Write();
    
    
    fout->Close();
    fmc->Close();
     
  cout<<"Response matrices written. Success!"<<endl;
    

}
#ifndef __CINT__
int main () { Unfolding3dAnalysisOpt(); return 0; }  // Main program when run stand-alone
#endif
