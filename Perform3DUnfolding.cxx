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
void Perform3DUnfolding(const char* responseFile = "", const char* dataFile = "", int nIter = 4) {
    
    TH1D* h1_triv;
    TH1D* h1_unfolded;
    
    
    TH3D* h3_true_eec; TH3D* h3_trueSplit;
    TH3D* h3_reco_eec;
    TH3D* h3_unfolded; TH3D* h3_unfolded_rw;  TH3D* h3_unfoldedSplit;
    TH3D* h3_raw; TH3D* h3_rawSplit;
    TH3D* h3_smeared;
    TH3D* h3_fulleff_match;
    TH3D* h3_eff_match;
    TH3D* h3_fullreco;
    TH3D* h3_rw;
    TH3D* h3_fold;
    TH3D* h3_triv;
    
    // RooUnfoldResponse* response1D;
    // RooUnfoldResponse* response1D_rw;
    // RooUnfoldResponse* response1DSplit;
    
    // RooUnfoldResponse* response3D;
    // RooUnfoldResponse* response3D_rw;
    // RooUnfoldResponse* response3DSplit;
  
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
    RooUnfoldResponse* response3D = (RooUnfoldResponse*)fResponse->Get("h3_reco_eec_h3_true_eec");
    if (!response1D) {
        std::cerr << "Could not retrieve 3Dresponse object!" << std::endl;
        return;
    }
    
    
    TH1D* h1_reco = (TH1D*)fResponse->Get("h1_reco");
    TH1D* h1_true = (TH1D*)fResponse->Get("h1_true");
    TH3D* h3_reco = (TH3D*)fResponse->Get("h3_reco_eec");
    TH3D* h3_true = (TH3D*)fResponse->Get("h3_true_eec");


    TH1* h1_raw = (TH1*) fData->Get("jet_pt_hist");
    if (!h1_raw) {
        std::cerr << "Could not retrieve measured histogram!" << std::endl;
        return;
    }

    // Check histograms for NaNs or Infs
    if (hasNaNsOrInfs(h1_raw)) {
        std::cerr << "Measured histogram contains NaNs or Infs!" << std::endl;
        return;
    }

    // 1D Unfolding using RooUnfold
        cout<<"Performing 1D unfolding"<<endl;
    // Create an unfolding object using the Bayes method
    RooUnfoldBayes unfold(response1D, h1_raw, nIter);
    RooUnfoldBayes unfoldTriv(response1D, h1_reco, nIter);

    // Perform unfolding
    h1_unfolded = (TH1D*)unfold.Hreco();
    h1_unfolded->SetName("h1_unfolded");
    
    h1_triv = (TH1D*)unfoldTriv.Hreco();//trivial unfolding 
    h1_triv->SetName("h1_triv");
   
   
    // Check for NaNs or Infs in unfolded histogram
    if (hasNaNsOrInfs(h1_unfolded) || hasNaNsOrInfs(h1_triv)) {
        std::cerr << "Unfolded histogram contains NaNs or Infs!" << std::endl;
        return;
    }
 
 cout<<"Performing 3D unfolding"<<endl;
// //     // Create an unfolding object using the Bayes method
    RooUnfoldBayes unfoldTriv3D(response3D, h3_reco, nIter);
   
   cout<<"this was okay"<<endl;
//  //     // Perform unfolding
    h3_triv = (TH3D*)unfoldTriv3D.Hreco();//trivial unfolding 
    h3_triv->SetName("h3_triv");
   
     // Check for NaNs or Infs in unfolded histogram
    if (hasNaNsOrInfs(h3_triv)) {
        std::cerr << "Unfolded histogram contains NaNs or Infs!" << std::endl;
        return;
    }
    
    
    // Save unfolded histogram to file
    TString foutname = Form("UnfoldingResults_%i.root",nIter);
    TFile* fOutput = new TFile(foutname, "RECREATE");
    
    h1_unfolded->Write();
    h1_raw->Write();
    h1_true->Write();
    h1_reco->Write();
    h1_triv->Write();
    
    h3_true->Write();
    h3_reco->Write();
    h3_triv->Write();
   

    
    fOutput->Close();

   // Clean up
    delete fResponse;
    delete fOutput;
    delete fData;

    std::cout << "Unfolding complete. Results saved to "<< foutname << std::endl;
}
