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


void gettingwt(int pt1, int pt2, const char* inputFile, const char* histName, const char* outputFile, std::string tag="") {
    
    int pt1jet = pt1;
    int pt2jet = pt2;
    
    const char *str = "#it{p}_{T,jet}";
    const char *str1 = "#it{p}_{T,min}";
    const char *str2 = "anti-#it{k}_{T}";
    const char *str3 = "#sqrt{s} = 13 TeV";
    const char *str4 = "#it{R}_{L}";
    const char *str8 = "#it{#eta}_{jet}";
    const char *str11 = "#it{p}_{T}^{const}";
    
    // Open the input file
    TFile *file = TFile::Open(inputFile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << inputFile << std::endl;
        return;
    }
    
    TH1D* hjet = (TH1D*)file->Get("jet_pt_hist");
    hjet->Scale(1.,"width");
    
     TCanvas* cjet;
     hjet->Draw();
     
    TH3D* h3 = (TH3D*)file->Get(histName);
    h3->GetXaxis()->SetTitle(str4);
    h3->GetYaxis()->SetTitle(str);
    h3->GetZaxis()->SetTitle("wt");
    

    if (!h3) {
        std::cerr << "Error retrieving histogram: " << histName << std::endl;
        file->Close();
        return;
    }
    h3->GetYaxis()->SetRange(pt1,pt2);
    // Get the 2D histogram
    TH2D *h2 = (TH2D*)h3->Project3D("zx"); //Rl vs wt histogram
    
    
    // Create the 1D histogram
    int nBinsX = h2->GetNbinsX();
    double xMin = h2->GetXaxis()->GetXmin();
    double xMax = h2->GetXaxis()->GetXmax();
    
    TH1D *h1 = (TH1D*)h2->ProjectionX("h1");
    h1->Reset();
   
    TCanvas* c;
    h3->Draw();
    // Loop over the x bins to create a simple average
//    for (int i = 1; i <= nBinsX; ++i) {
//        int nBinsY = h2->GetNbinsY();
//        double sum = 0;
//        for (int j = 1; j <= nBinsY; ++j) {
//            sum += h2->GetBinContent(i, j);
//        }
//        double average = sum / nBinsY;
//        double xCenter = h2->GetXaxis()->GetBinCenter(i);
//        h1->SetBinContent(i, average);
//    }
    
    
    //Loop over the x bins to create a weighted average
    for (int i = 1; i <= nBinsX; ++i) {
    int nBinsY = h2->GetNbinsY();
    double weightedSum = 0;
    double sumOfWeights = 0;

    for (int j = 1; j <= nBinsY; ++j) {
        double binContent = h2->GetBinContent(i, j);
        double binCenterY = h2->GetYaxis()->GetBinCenter(j); // Bin center in Y direction
        
        // Compute weight as the bin center in Y direction
        double weight = binCenterY;
        
        // Accumulate weighted sum and sum of weights
        weightedSum += binContent * weight;
        sumOfWeights += weight;
    }
    // Compute weighted average, avoiding division by zero
    double weightedAverage = (sumOfWeights != 0) ? (weightedSum / sumOfWeights) : 0;
    double xCenter = h2->GetXaxis()->GetBinCenter(i);
//    h1->SetBinContent(i, weightedAverage);
    h1->SetBinContent(i, weightedSum);//this is closest to filling on the fly, not the avg
}

//cout<<h1->Integral()<<endl;
//for very fine bins wt*n
// for (int i = 1; i <= nBinsX; ++i) {
//    int nBinsY = h2->GetNbinsY();
//    double weightedSum = 0;
//    double sumOfWeights = 0;
//    double weight =0;
//    for (int j = 1; j <= nBinsY; ++j) {
//        double binContent = h2->GetBinContent(i, j);
//        double binCenterY = h2->GetYaxis()->GetBinCenter(j); // Bin center in Y direction
//        
//        // Compute weight as the bin center in Y direction
////        double weight = binCenterY;
//        weight = binCenterY;
//        
//        // Accumulate weighted sum and sum of weights
////        weightedSum += binContent * weight;
//        weightedSum += weight;
//    }
//    double xCenter = h2->GetXaxis()->GetBinCenter(i);
//    h1->SetBinContent(i, weightedSum);
//}
  
    // Save the 1D histogram to a new file
    TFile *outFile = TFile::Open(outputFile, "UPDATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error opening output file: " << outputFile << std::endl;
        file->Close();
        return;
    }
    
    int low = h3->GetYaxis()->GetBinLowEdge(pt1);
    int high = h3->GetYaxis()->GetBinUpEdge(pt2);
    
   
    if(tag == "area"){
        h1->Scale(1./h1->Integral(),"width");
        h1->SetName(Form("h1_e3c_AreaNorm_%i_%i",low,high));
        h1->SetLineColor(kRed);
        h1->Write();
    }
    else if(tag == "njet"){
        h1->Scale(1./hjet->Integral(pt1,pt2),"width");
        h1->SetName(Form("h1_e3c_NjetNorm_%i_%i",low,high));
        h1->SetLineColor(kBlack);
        h1->Write();
    }
    else{
        h1->SetLineColor(kBlue);
        h1->Scale(1.,"width");
        h1->SetName(Form("h1_e3c_%i_%i",low,high));
        h1->Write();
    }
    
    
    //What you would normally get
    TH2D* h2_nom = (TH2D*)file->Get("eec_pt_hist");
    TH1D* h1_nom = (TH1D*)h2_nom->ProjectionX(Form("h1nom_eec_%i_%i",low,high),pt1,pt2);
    if(tag == "area"){
        h1_nom->Scale(1./h1_nom->Integral(),"width");
        h1_nom->SetName(Form("h1nom_eec_AreaNorm_%i_%i",low,high));
    }
    else if(tag == "njet")
    {
        
        h1_nom->Scale(1./hjet->Integral(pt1,pt2),"width");
        h1_nom->SetName(Form("h1nom_eec_NjetNorm_%i_%i",low,high));
        
    }
    else{
        h1_nom->Scale(1.,"width");
    }

    h1_nom->Write();
    
    TCanvas* c1= new TCanvas();
    h1->Draw();
    
//    TCanvas* c2 = new TCanvas();
    h1_nom->Draw("same");
//    outFile->Close();
    //    file->Close();
    
    std::cout << "1D histogram saved to " << outputFile << std::endl;
}


