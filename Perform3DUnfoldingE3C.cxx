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

#endif

// Function to check for NaNs or Infs in 1D, 2D, and 3D histograms
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

// Main function
void Perform3DUnfoldingE3C(const char* responseFile = "", const char* dataFile = "", int nIter = 4, std::string name = "") {
    
    TH1D* h1_triv;
    TH1D* h1_unfolded;
    
    TH3D* h3_unfolded; TH3D* h3_unfoldedTriv;  TH3D* h3_unfoldedSplit;
    
   // Load response matrix from file
    TFile* fResponse = TFile::Open(responseFile);
    if (!fResponse) {
        std::cerr << "Could not open response file!" << std::endl;
        return;
    }
    // Load measured data from file
    TFile* fData = TFile::Open(dataFile);
    if (!fData) {
        std::cerr << "Could not open data file!" << std::endl;
        return;
    }
    
    
    RooUnfoldResponse* response1D = (RooUnfoldResponse*)fResponse->Get("h1_reco_h1_true");
    if (!response1D) {
        std::cerr << "Could not retrieve 1Dresponse object!" << std::endl;
        return;
    }
    RooUnfoldResponse* response3D = (RooUnfoldResponse*)fResponse->Get("h3_reco_e3c_h3_true_e3c");
    if (!response3D) {
        std::cerr << "Could not retrieve nominal 3Dresponse object!" << std::endl;
        return;
    }
    
     RooUnfoldResponse* response3DTriv = (RooUnfoldResponse*)fResponse->Get("h3_reco_e3cTriv_h3_true_e3cTriv");
    if (!response3DTriv) {
        std::cerr << "Could not retrieve triv 3Dresponse object!" << std::endl;
        return;
    }
    
    RooUnfoldResponse* response3DSplit = (RooUnfoldResponse*)fResponse->Get("h3_rawSplit_h3_trueSplit");
    if (!response3DSplit) {
        std::cerr << "Could not retrieve split 3Dresponse object!" << std::endl;
        return;
    }
    
  
    TH1D* h1_reco = (TH1D*)fResponse->Get("h1_reco");
    TH1D* h1_true = (TH1D*)fResponse->Get("h1_true");
    TH1D* h1_fulleff = (TH1D*)fResponse->Get("h1_fulleff");
    
    TH3D* h3_reco = (TH3D*)fResponse->Get("h3_reco_e3c");
    TH3D* h3_true = (TH3D*)fResponse->Get("h3_true_e3c");
    
    TH3D* h3_fulleff= (TH3D*)fResponse->Get("h3_fulleff");
    TH3D* h3_eff_match = (TH3D*)fResponse->Get("h3_eff_match");
    TH3D* h3_fullreco = (TH3D*)fResponse->Get("h3_fullreco");
    TH3D* h3_smeared = (TH3D*)fResponse->Get("h3_smeared");
   
    TH3D* h3_recoTriv = (TH3D*)fResponse->Get("h3_reco_e3cTriv");
    TH3D* h3_trueTriv = (TH3D*)fResponse->Get("h3_true_e3cTriv");

    TH3D* h3_rawSplit = (TH3D*)fResponse->Get("h3_rawSplit");
    TH3D* h3_trueSplit = (TH3D*)fResponse->Get("h3_trueSplit");
    
    TH3D* h3_purity = (TH3D*)h3_reco->Clone("h3_purity");
    h3_purity->Divide(h3_fullreco);
  
     
    TH3D* h3_efficiency = (TH3D*)h3_true->Clone("h3_efficiency");
    h3_efficiency->Divide(h3_fulleff);
 
    

    cout<<"######Getting Data histograms#####"<<endl;
    TH1D* h1_raw = (TH1D*)fData->Get("jet_pt_hist");
    if (!h1_raw) {
        std::cerr << "Could not retrieve measured histogram!" << std::endl;
        return;
    }
    
    TFile* fUsefulHist = TFile::Open("/home/ar2545/Trial/DataFiles/AnalysisResults_DataCopied_Apr16.root");
    TH1F* h1_numEvents = (TH1F*)fUsefulHist->Get("fHistEventCount");
    // double scalefactor = 22.55e06;
    
    // h1_raw->Scale(1./scalefactor);
    
    TH3D* h3_raw = (TH3D*)fData->Get("Opt_Un_e3c_unf");
    h3_raw->SetName("h3_raw");
    if (!h3_raw) {
        std::cerr << "Could not retrieve measured histogram!" << std::endl;
        return;
    }

    // Check histograms for NaNs or Infs
    if (hasNaNsOrInfs(h3_raw)) {
        std::cerr << "Measured histogram contains NaNs or Infs!" << std::endl;
        return;
    }


    // 1D Unfolding using RooUnfold
    cout<<"######################################Performing 1D unfolding ###################################"<<endl;
//     // // Create an unfolding object using the Bayes method
//     RooUnfoldBayes unfold(response1D, h1_raw, nIter);
//     RooUnfoldBayes unfoldTriv(response1D, h1_reco, nIter);

//     // Perform unfolding
//     h1_unfolded = (TH1D*)unfold.Hreco();
//     h1_unfolded->SetName("h1_unfolded");
    
//     // h1_triv = (TH1D*)unfoldTriv.Hreco();//trivial unfolding 
//     // h1_triv->SetName("h1_triv");
   
//   //FOLD BACK
//     TH1* hfold1D = response1D->ApplyToTruth(h1_unfolded, "");
//     hfold1D->SetName("h1_fold");
//     // Check for NaNs or Infs in unfolded histogram
//     if (hasNaNsOrInfs(h1_unfolded) || hasNaNsOrInfs(h1_triv)) {
//         std::cerr << "Unfolded histogram contains NaNs or Infs!" << std::endl;
//         return;
//     }
 
 cout<<"######################################Performing 3D unfolding ###################################"<<endl;
//     // Create an unfolding object using the Bayes method
    
    //  RooUnfoldBayes unfoldSplit3D(response3DSplit, h3_raw, nIter);
   
//   cout<<"unfolding objects created without issue"<<endl;
//  //     // Perform unfolding
//     RooUnfoldBayes unfoldTriv3D(response3DTriv, h3_recoTriv, nIter);
//     h3_unfoldedTriv = (TH3D*)unfoldTriv3D.Hreco();//trivial unfolding 
//     h3_unfoldedTriv->SetName("h3_triv");
    
//     TH1* hfold3D = response3DTriv->ApplyToTruth(h3_unfoldedTriv, "");
//     hfold3D->SetName("h3_fold");
   
    
    //  // Check for NaNs or Infs in unfolded histogram
    // if (hasNaNsOrInfs(h3_unfoldedTriv)) {
    //     std::cerr << "Unfolded histogram contains NaNs or Infs!" << std::endl;
    //     return;
    // }
    
    cout<<"####Trivial unfolding finished####"<<endl;
    
    // RooUnfoldBayes unfoldSplit3D(response3DSplit, h3_rawSplit, nIter);
    // h3_unfoldedSplit = (TH3D*)unfoldSplit3D.Hreco();//mc split test unfolding 
    // h3_unfoldedSplit->SetName("h3_unfoldedSplit");
    
    //  // Check for NaNs or Infs in unfolded histogram
    // if (hasNaNsOrInfs(h3_unfoldedSplit)) {
    //     std::cerr << "Unfolded histogram contains NaNs or Infs!" << std::endl;
    //     return;
    // }
    
    cout<<"####Split unfolding finished####"<<endl;
      //correct the raw distribution for the purity
    TH3D* h3_raw_corr = (TH3D*)h3_raw->Clone("h3_raw_corr");
    h3_raw_corr->Multiply(h3_purity);
    
    RooUnfoldBayes unfold3D(response3D, h3_raw_corr, nIter);
    h3_unfolded = (TH3D*)unfold3D.Hreco();//data unfolding 
    h3_unfolded->SetName("h3_unfoldedData");
    
    TH1* hfold3D = response3D->ApplyToTruth(h3_unfolded, "");
    hfold3D->SetName("h3_fold");
    
    cout<<"####Nominal unfolding finished####"<<endl;
    
    // Save unfolded histogram to file
    TString foutname = Form("UnfoldingResultsE3CData_%i%s.root",nIter,name.c_str());
    TFile* fOutput = new TFile(foutname, "RECREATE");
    
    cout<<"1D HISTOS WRITING "<<endl;
    // h1_unfolded->Write();
    // h1_raw->Write();
    // h1_true->Write();
    // h1_reco->Write();
    // h1_fulleff->Write();
    // hfold1D->Write();
    
   
    //  h1_triv->Write();
    
    cout<<"3D HISTOS WRITING "<<endl;
    h3_raw->Write();
    h3_smeared->Write();
    h3_true->Write();
    h3_reco->Write();
    h3_fulleff->Write();
    h3_eff_match->Write();
    h3_fullreco->Write();
    
    h3_raw_corr->Write();
    h3_purity->Write();
    h3_efficiency->Write();
    
    h3_unfolded->Write();
    hfold3D->Write();
    
    
    cout<<"3D HISTOS TRIV WRITING "<<endl;
    // h3_raw->Write();
    // h3_trueTriv->Write();
    // h3_recoTriv->Write();
    // h3_unfoldedTriv->Write();
    // hfold3D->Write();
    
    cout<<"3D HISTOS SPLIT WRITING "<<endl;
    // h3_unfoldedSplit->Write();
    // h3_trueSplit->Write();
    // h3_rawSplit->Write();


    cout<<"WRITING 2D HISTOGRAMS"<<endl;
    
    TH2D* e3c_pt_hist = (TH2D*)fResponse->Get("e3c_pt_hist");
    TH2D* e3c_pt_hist_det = (TH2D*)fResponse->Get("e3c_pt_hist_det");
    TH2D* e3c_pt_histTriv = (TH2D*)fResponse->Get("e3c_pt_histTriv");
    TH2D* e3c_pt_hist_detTriv = (TH2D*)fResponse->Get("e3c_pt_hist_detTriv");
    TH2D* e3c_pt_histSplit = (TH2D*)fResponse->Get("e3c_pt_histSplit");
    TH2D* e3c_pt_hist_detSplit = (TH2D*)fResponse->Get("e3c_pt_hist_detSplit");
    
    e3c_pt_hist->Write();
    e3c_pt_hist_det->Write();
    e3c_pt_histTriv->Write();
    e3c_pt_hist_detTriv->Write();
    e3c_pt_histSplit->Write();
    e3c_pt_hist_detSplit->Write();
   
    
    fOutput->Close();

   // Clean up
    delete fResponse;
    delete fOutput;
    delete fData;

    std::cout << "Unfolding complete. Results saved to "<< foutname << std::endl;
}
