//Checking weight binning for EEC
//possibility to check trivial case or rebin a super fine histogram and then check
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TKey.h"
#include "TIterator.h"
#include <iostream>
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


void getUnfoldedEEC(int pt1, int pt2, int nIter) {
    
   const char* inputFile = Form("/Users/ar2545/Downloads/UnfoldingResultsDataEEC_%iEECnew6.root", nIter);
    std::string histName3D_unfolded = "h3_unfoldedData";
    std::string histName3D = "Opt_Un_eec";
    std::string histName2D = "eec_pt_hist_unf";
    std::string histName3D_true = "h3_true_eec";
    std::string histName3D_fulleff = "h3_fulleff";
    
 
    const char *str = "#it{p}_{T,jet}";
    const char *str1 = "#it{p}_{T,min}";
    const char *str2 = "anti-#it{k}_{T}";
    const char *str3 = "#sqrt{s} = 13 TeV";
    const char *str4 = "#it{R}_{L}";
    const char *str8 = "#it{#eta}_{jet}";
    const char *str11 = "#it{p}_{T}^{const}";
    
   // Open the input file
    TFile *f = TFile::Open(inputFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file: " << inputFile << std::endl;
        return;
    }

    TH3D* h3_unfolded = (TH3D*)f->Get(histName3D_unfolded.c_str());
    h3_unfolded->GetXaxis()->SetTitle(str4);
    h3_unfolded->GetYaxis()->SetTitle(str);
    h3_unfolded->GetZaxis()->SetTitle("wt");
    
    TH3D* h3_true = (TH3D*)f->Get(histName3D_true.c_str());
    TH3D* h3_fulleff = (TH3D*)f->Get(histName3D_fulleff.c_str());
    
    //Compute kinematic efficiency
    TH3D* h3_efficiency = (TH3D*)h3_true->Clone("h3_efficiency");
    h3_efficiency->Divide(h3_fulleff);
    
    
    TH1* h3_fulleff_x = (TH1*)h3_fulleff->Project3D("x");
    h3_fulleff_x->SetLineColor(kRed);
    h3_fulleff_x->Scale(1.,"width");
    TH1* h3_fulleff_y = (TH1*)h3_fulleff->Project3D("y");
    h3_fulleff_y->SetLineColor(kRed);
    h3_fulleff_y->Scale(1.,"width");
    TH1* h3_fulleff_z = (TH1*)h3_fulleff->Project3D("z");
    h3_fulleff_z->SetLineColor(kRed);
    h3_fulleff_z->Scale(1.,"width");
    
    TH1* h3_tru_x = (TH1*)h3_true->Project3D("x");
    h3_tru_x->Scale(1.,"width");
    TH1* h3_tru_y = (TH1*)h3_true->Project3D("y");
    h3_tru_y->Scale(1.,"width");
    TH1* h3_tru_z = (TH1*)h3_true->Project3D("z");
    h3_tru_z->Scale(1.,"width");
    
    //Correct unfolded histogram with kinematic efficiency
    h3_unfolded->Divide(h3_efficiency);
   
    if (!h3_unfolded) {
        std::cerr << "Error retrieving histogram: " << histName3D_unfolded << std::endl;
        f->Close();
        return;
    }
     h3_unfolded->GetYaxis()->SetRange(h3_unfolded->GetYaxis()->FindBin(pt1),h3_unfolded->GetYaxis()->FindBin(pt2));

//    // Get the 2D histogram
    TH2D *h2_unfolded = (TH2D*)h3_unfolded->Project3D("zx"); //Rl vs wt histogram
    TH1D *h1_unfolded = (TH1D*)h2_unfolded->ProjectionX("h1_unfolded");
    h1_unfolded->Reset();
    
        // Create the 1D histogram
    int nBinsX = h2_unfolded->GetNbinsX();
    double xMin = h2_unfolded->GetXaxis()->GetXmin();
    double xMax = h2_unfolded->GetXaxis()->GetXmax();
 

///############################SET UP WHICH KIND OF PROJECTION YOU WANT TO DO------------------------------------------------------------

    bool ifInterpolate = false;
    
 
    TCanvas* c = new TCanvas("c1Dhist","1dhist",800,650);
    //        h3->Draw();
    
    if(!ifInterpolate)
    {
        //Loop over the x bins to create a weighted average
       h1_unfolded->SetName(Form("EEC_unfolded_%i_%i",pt1,pt2));
        for (int i = 1; i <= h2_unfolded->GetNbinsX(); ++i) {
            int nBinsY = h2_unfolded->GetNbinsY();
            double weightedSum = 0;
            double sumOfWeights = 0;
            
            for (int j = 1; j <= nBinsY; ++j) {
                double binContent = h2_unfolded->GetBinContent(i, j);
                double binCenterY = h2_unfolded->GetYaxis()->GetBinCenter(j); // Bin center in Y direction
                
                // Compute weight as the bin center in Y direction
                double weight = binCenterY;
                
                // Accumulate weighted sum and sum of weights
                weightedSum += binContent * weight;
                sumOfWeights += weight;
            }
            // Compute weighted average, avoiding division by zero
            double weightedAverage = (sumOfWeights != 0) ? (weightedSum / sumOfWeights) : 0;
            double xCenter = h2_unfolded->GetXaxis()->GetBinCenter(i);
            h1_unfolded->SetBinContent(i, weightedSum);//this is closest to filling on the fly, not the avg
        }
    }
    else{
        cout<<"Interpolating"<<endl;
         h1_unfolded->SetName(Form("EEC_unfoldedInterpolated_%i_%i",pt1,pt2));
        for (int i = 1; i <= h2_unfolded->GetNbinsX(); ++i){
            int nBinsY = h2_unfolded->GetNbinsY();
            double weightedSum = 0;
            TH1D *h1_wt_unfolded = (TH1D*)h2_unfolded->ProjectionY(Form("h1_wt_unfolded_%i",i),i,i);
            for (int j = 1; j <= nBinsY; ++j) {
                double weight = h1_wt_unfolded->GetBinCenter(j);
                double bin_content = h1_wt_unfolded->Interpolate(weight);
                weightedSum += weight*bin_content;
               
            }
            h1_unfolded->SetBinContent(i, weightedSum);//this is closest to filling on the fly, not the avg
        }

    }
  
//    h1_unfolded->Scale(1.,"width");
    h1_unfolded->SetLineColor(kRed);
    
//    TH2D* h2_nom;
//    h2_nom = (TH2D*)f->Get(histName2D.c_str());
//    
//    TH1D* h1_nom = (TH1D*)h2_nom->ProjectionX(Form("h1nom"),h2_nom->GetYaxis()->FindBin(pt1),h2_nom->GetYaxis()->FindBin(pt2));
//    h1_nom->Scale(1.,"width");
//    //        h1->Divide(h1_nom);
//    h1_unfolded->Divide(h1_nom);
          
//    h1_unfolded->Divide(h1);
//     h1->Divide(h1_nom);

//    h1->Draw();
//    h1_nom->Draw("same");
    h1_unfolded->Draw("same");
   
    h3_tru_z->Draw();
    h3_fulleff_z->Draw("same");
    
    const char* fileOutName = Form("~/Desktop/Unfolding/ProjectedOutputEECNew6_%i.root",nIter);
    TFile* fileOut = TFile::Open(fileOutName,"UPDATE");
    h1_unfolded->Write();
    fileOut->Close();
    
  
}
