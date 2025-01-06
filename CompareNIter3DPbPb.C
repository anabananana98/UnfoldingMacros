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



void FillWeightedSum(TH2* h2, TH1* h1) {
    if (!h2 || !h1) {
        std::cerr << "Error: Null histogram pointer provided!" << std::endl;
        return;
    }

    int nBinsX = h2->GetNbinsX();
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
        
        // Fill the 1D histogram bin
        h1->SetBinContent(i, weightedSum); // Using the weighted sum directly
    }
}


void CompareNIter3DPbPb(std::string tag, std::string projection) {
    bool ifFold = false;
    bool ifCheckProjected = false;
    bool ifbasic = false;
    bool ifeec = true;
    bool ifcheckActualENC = true;
    
    int num=4;
    int pt1 = 40;
    int pt2 = 100-1;
    double rlow = 0.01;
    double rhigh = 0.4;
    
    const char *str = "#it{p}_{T,jet}";
    const char *str1 = "#it{p}_{T,min}";
    const char *str2 = "anti-#it{k}_{T}";
    const char *str3 = "#sqrt{s} = 13 TeV";
    const char *str4 = "#it{R}_{L}";
    const char *str8 = "#it{#eta}_{jet}";
    const char *str11 = "#it{p}_{T}^{const}";
    
    std::string encstr = "EEC";
    if(!ifeec){
    encstr = "E3C";
    }
    
    TLatex latex;latex.SetNDC ();latex.SetTextSize(0.045);latex.SetTextFont(42);
    TLatex latex1;latex1.SetNDC ();latex1.SetTextSize(0.045);latex1.SetTextFont(42);
    TLatex latex2;latex2.SetNDC ();latex2.SetTextSize(0.045);latex2.SetTextFont(42);
    TLatex latex3;latex3.SetNDC ();latex3.SetTextSize(0.045);latex3.SetTextFont(42);
    TLatex latex4;latex4.SetNDC ();latex4.SetTextSize(0.045);latex4.SetTextFont(42);
    TLatex latex5;latex5.SetNDC ();latex5.SetTextSize(0.045);latex5.SetTextFont(42);
    
    // Create a legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.045);
    
        // Create a legend
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.88, 0.88);
    legend1->SetBorderSize(0);
    legend1->SetTextFont(42);
    legend1->SetTextSize(0.045);
    
     TLegend *legend2 = new TLegend(0.7, 0.7, 0.88, 0.88);
    legend2->SetBorderSize(0);
    legend2->SetTextFont(42);
    legend2->SetTextSize(0.045);

    
    const char* filenames[4] = {
    
    //PYTHIA TESTS FOR PBPB DATA SETUP
           "/Users/ar2545/Downloads/UnfoldingResultsPythiaEEC_1EECPythiaMin60.root",
                "/Users/ar2545/Downloads/UnfoldingResultsPythiaEEC_4EECPythiaMin60.root",
                "/Users/ar2545/Downloads/UnfoldingResultsPythiaEEC_7EECPythiaMin60.root",
                "/Users/ar2545/Downloads/UnfoldingResultsPythiaEEC_10EECPythiaMin60.root",
        
  
    };
    
    const char* base = filenames[1];
    
    
    // Array of unfolded histogram names
    const char* histNamesUnf3d[4] = {
        "iter13d",
        "iter43d",
        "iter73d",
        "iter103d",
    };
    
    const char* histNamesUnf[4] = {
        "iter1",
        "iter4",
        "iter7",
        "iter10",
    };
    
    //    // Array of true histogram names
    const char* histNamesTru3d[4] = {
        "true13d",
        "true23d",
        "true33d",
        "true43d",
    };
    
    const char* histNamesTru[4] = {
        "true1",
        "true2",
        "true3",
        "true4",
    };
    //    const char* histNamesTru = "h3_true";
    
    
    // Original 3D histogram name in the files
    const char* originalHistNameUnf;
    const char* originalHistNameTru;
    
    if(tag == "split")
    {
        originalHistNameUnf = "h3_unfoldedSplit";
        originalHistNameTru = "h3_trueSplit";
    }
    else if(tag == "triv")
    {
        originalHistNameUnf = "h3_triv";
        if(ifeec){
            originalHistNameTru = "h3_true_eecTriv";
        }
        else{
            originalHistNameTru = "h3_true_e3cTriv";
        }
    }
    else if(tag == "fold")
    {
        originalHistNameUnf = "h3_fold";
            originalHistNameTru = "h3_raw_corr";
        //     originalHistNameTru = "Opt_Un_eec_unf";
        //     originalHistNameTru = "h3_reco_eecTriv";
        //     originalHistNameTru = "h3_reco_eec";
    }
    else{
        originalHistNameUnf = "h3_unfoldedData";
        if(ifeec){
            originalHistNameTru = "h3_true_eec";
        }
        else{
            originalHistNameTru = "h3_true_e3c";
        }
    
    }
    

    // Array of colors for histograms
    Color_t colors[4] = {kRed, kBlue, kMagenta, kCyan};
    
    // Array of markers for histograms
    Color_t markers[4] = {20, 21, 29, 33};
    
    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Overlay Histograms", 800, 600);
//    c1->SetGrid();
     
    TFile *file1 = TFile::Open(base);
    TH3 *hist3DUnf_base = (TH3*)file1->Get(originalHistNameUnf);
    hist3DUnf_base->SetName("base");
    TH1 *hist1DUnf_base;
     if(projection=="x"){
            hist1DUnf_base = hist3DUnf_base->Project3D("x");
           
        }
        else if(projection=="y")
        {
            hist1DUnf_base = hist3DUnf_base->Project3D("y");

        }
        else{
            hist1DUnf_base = hist3DUnf_base->Project3D("z");
        }


    if(ifcheckActualENC){
//        TH3D* h3_fulleff_base = (TH3D*)file1->Get("h3_fulleff");
//        TH3 *hist3DTru_base = (TH3*)file1->Get("h3_true_eec");
//        
//        //Compute kinematic efficiency
//        TH3D* h3_efficiency_base = (TH3D*)hist3DTru_base->Clone("h3_efficiency_base");
//        h3_efficiency_base->Divide(h3_fulleff_base);
//        hist3DUnf_base->Divide(h3_efficiency_base);
//        
        hist3DUnf_base->GetXaxis()->SetRange(hist3DUnf_base->GetXaxis()->FindBin(pt1),hist3DUnf_base->GetXaxis()->FindBin(pt2));
        TH2D *h2D_base = (TH2D*)hist3DUnf_base->Project3D("zy"); //Rl vs wt histogra
        TH1D *h1_base = (TH1D*)h2D_base->ProjectionX("h1_base");
        h1_base->Reset();
        FillWeightedSum(h2D_base, h1_base);
        h1_base->Scale(1.,"width");
        h1_base->Draw();
        
        // Create a canvas
        TCanvas *c1 = new TCanvas("c12", " New Overlay Histograms", 600, 800);
        // Create the top and bottom pads
        TPad *topPad = new TPad("topPad", "Top Pad", 0, 0.3, 1, 1.0);
        TPad *bottomPad = new TPad("bottomPad", "Bottom Pad", 0, 0, 1, 0.3);
        // Adjust margins for the pads
        topPad->SetBottomMargin(0); // Remove the bottom margin for the top pad
        bottomPad->SetTopMargin(0); // Remove the top margin for the bottom pad
        bottomPad->SetBottomMargin(0.3); // Increase the bottom margin of the bottom pad for axis labels
        
        // Draw the pads on the canvas
        topPad->Draw();
        bottomPad->Draw();
        
        
        for (int i = 0; i < num; i++) {
            // Open the file
            TFile *file = TFile::Open(filenames[i]);
            if (!file || file->IsZombie()) {
                std::cerr << "Error opening file: " << filenames[i] << std::endl;
                continue;
            }
            
            topPad->cd();
            gPad->SetLogx();
            // Get the 3D histogram
            TH3 *hist3DUnf = (TH3*)file->Get(originalHistNameUnf);
            if (!hist3DUnf) {
                std::cerr << "Error retrieving histogram: " << originalHistNameUnf << " from file: " << filenames[i] << std::endl;
                file->Close();
                continue;
            }
            hist3DUnf->SetName(histNamesUnf3d[i]);
            hist3DUnf->GetXaxis()->SetTitle(str);
            hist3DUnf->GetYaxis()->SetTitle(str4);
            hist3DUnf->GetZaxis()->SetTitle("wt");
            hist3DUnf->SetMarkerStyle(markers[i]);
            hist3DUnf->SetMarkerColor(colors[i]);
            hist3DUnf->Sumw2();
            
            
            TH3 *hist3DTru = (TH3*)file->Get(originalHistNameTru);
            if (!hist3DTru) {
                std::cerr << "Error retrieving histogram: " << originalHistNameTru << " from file: " << filenames[i] << std::endl;
                file->Close();
                continue;
            }
            hist3DTru->SetName(histNamesTru3d[i]);
            
            
//            TH3D* h3_fulleff = (TH3D*)file->Get("h3_fulleff");
//            //Compute kinematic efficiency
//            TH3D* h3_efficiency = (TH3D*)hist3DTru->Clone(Form("h3_efficiency_%i",i));
//            h3_efficiency->Divide(h3_fulleff);
//            
//            hist3DUnf->Divide(h3_efficiency);
            
            cout<<"hher"<<endl;
            // Project the 3D histogram to a 1D histogram along the X-axis
            TH1 *hist1DUnf; //= hist3DUnf->ProjectionX(histNamesUnf[i],2,2,10,10);
            TH1 *hist1DTru; //= hist3DTru->ProjectionX(histNamesTru[i],2,2,10,10);
            
            
            hist3DTru->GetXaxis()->SetRange(hist3DTru->GetXaxis()->FindBin(pt1),hist3DTru->GetXaxis()->FindBin(pt2));
            TH2D *h2DTru = (TH2D*)hist3DTru->Project3D("zy"); //Rl vs wt histogra
            TH1D *h1Tru = (TH1D*)h2DTru->ProjectionX(Form("h1Tru_%i",i));
            h1Tru->Reset();
            FillWeightedSum(h2DTru, h1Tru);
            
            hist3DUnf->GetXaxis()->SetRange(hist3DUnf->GetXaxis()->FindBin(pt1),hist3DUnf->GetXaxis()->FindBin(pt2));
            TH2D *h2DUnf = (TH2D*)hist3DUnf->Project3D("zy"); //Rl vs wt histogra
            TH1D *h1Unf = (TH1D*)h2DUnf->ProjectionX(Form("h1Unf_%i",i));
            h1Unf->Reset();
            FillWeightedSum(h2DUnf, h1Unf);
            
            h1Unf->SetLineColor(colors[i]); // Set line color
            h1Unf->SetLineWidth(2);         // Set line thickness to 2
            h1Unf->SetMarkerColor(colors[i]);
            h1Unf->Scale(1.,"width");
            h1Unf->SetStats(0);
            h1Unf->GetXaxis()->SetRangeUser(0.01,0.4);
            h1Unf->SetTitle(Form("Unfolded %s iterations",encstr.c_str()));
            h1Unf->Draw("same");
            legend->AddEntry(h1Unf, histNamesUnf[i], "l");
            
            bottomPad->cd();
            gPad->SetLogx();
            
            TH1* h1Unf_clone = (TH1*)h1Unf->Clone(Form("h1UnfClone_%i",i));
            h1Unf_clone->Divide(h1_base);
            h1Unf_clone->SetStats(0);
            h1Unf_clone->SetTitle("");
            h1Unf_clone->GetYaxis()->SetLabelSize(0.075);
            h1Unf_clone->GetYaxis()->SetTitleSize(0.08);
            h1Unf_clone->GetYaxis()->SetTitleOffset(0.7);
            
            h1Unf_clone->GetXaxis()->SetRangeUser(0.01,0.4);
            h1Unf_clone->GetXaxis()->SetLabelSize(0.09);
            h1Unf_clone->GetXaxis()->SetTitleSize(0.09);
            h1Unf_clone->GetXaxis()->SetTitleOffset(1.07);
            h1Unf_clone->GetYaxis()->SetRangeUser(0.9,1.11);
            h1Unf_clone->GetYaxis()->SetTitle("Unfolded/iter4");
            h1Unf_clone->Draw("same");
            
        }
        topPad->cd();
        legend->Draw();
    }




    if(ifbasic){
        // Loop over the files, retrieve and rename histograms
        for (int i = 0; i < num; i++) {
            // Open the file
            TFile *file = TFile::Open(filenames[i]);
            if (!file || file->IsZombie()) {
                std::cerr << "Error opening file: " << filenames[i] << std::endl;
                continue;
            }
            
            // Get the 3D histogram
            TH3 *hist3DUnf = (TH3*)file->Get(originalHistNameUnf);
            if (!hist3DUnf) {
                std::cerr << "Error retrieving histogram: " << originalHistNameUnf << " from file: " << filenames[i] << std::endl;
                file->Close();
                continue;
            }
            hist3DUnf->SetName(histNamesUnf3d[i]);
            hist3DUnf->GetXaxis()->SetTitle(str4);
            hist3DUnf->GetYaxis()->SetTitle(str);
            hist3DUnf->GetZaxis()->SetTitle("wt");
            hist3DUnf->SetMarkerStyle(markers[i]);
            hist3DUnf->SetMarkerColor(colors[i]);
            hist3DUnf->Sumw2();
            
            
            TH3 *hist3DTru = (TH3*)file->Get(originalHistNameTru);
            if (!hist3DTru) {
                std::cerr << "Error retrieving histogram: " << originalHistNameTru << " from file: " << filenames[i] << std::endl;
                file->Close();
                continue;
            }
            hist3DTru->SetName(histNamesTru3d[i]);
            
            cout<<"hher"<<endl;
            // Project the 3D histogram to a 1D histogram along the X-axis
            TH1 *hist1DUnf; //= hist3DUnf->ProjectionX(histNamesUnf[i],2,2,10,10);
            TH1 *hist1DTru; //= hist3DTru->ProjectionX(histNamesTru[i],2,2,10,10);
            
            if(projection=="x"){
                //            hist1DUnf = hist3DUnf->ProjectionX(histNamesUnf[i],2,5,10,10);
                //            hist1DTru = hist3DTru->ProjectionX(histNamesTru[i],2,5,10,10);
                
                //            hist3DUnf->GetYaxis()->SetRangeUser(pt1,pt2);
                //            hist3DUnf->GetYaxis()->SetRangeUser(pt1,pt2);
                hist1DUnf = hist3DUnf->Project3D("x");
                hist1DTru = hist3DTru->Project3D("x");
            }
            else if(projection=="y")
            {
                //            hist1DUnf = hist3DUnf->ProjectionY(histNamesUnf[i],2,2,10,10);
                //            hist1DTru = hist3DTru->ProjectionY(histNamesTru[i],2,2,10,10);
                
                //            hist3DUnf->GetXaxis()->SetRangeUser(rlow,rhigh);
                //            hist3DUnf->GetXaxis()->SetRangeUser(rlow,rhigh);
                hist1DUnf = hist3DUnf->Project3D("y");
                hist1DTru = hist3DTru->Project3D("y");
            }
            else{
                //            hist1DUnf = hist3DUnf->ProjectionY(histNamesUnf[i],2,2,10,10);
                //            hist1DTru = hist3DTru->ProjectionY(histNamesTru[i],2,2,10,10);
                
                //            hist3DUnf->GetXaxis()->SetRangeUser(rlow,rhigh);
                //            hist3DUnf->GetXaxis()->SetRangeUser(rlow,rhigh);
                //            hist3DUnf->GetYaxis()->SetRangeUser(pt1,pt2);
                //            hist3DUnf->GetYaxis()->SetRangeUser(pt1,pt2);
                hist1DUnf = hist3DUnf->Project3D("z");
                hist1DTru = hist3DTru->Project3D("z");
            }
            
            hist1DUnf->SetName(histNamesUnf[i]);
            hist1DTru->SetName(histNamesTru[i]);
            //        hist1DUnf->GetYaxis()->SetRangeUser(0.5,1.5);
            //        hist1DTru->GetYaxis()->SetRangeUser(0.5,1.5);
            //        TH1 *hist1DUnf = hist3DUnf->Project3D(histNamesUnf[i],"x");
            //        TH1 *hist1DTru = hist3DTru->Project3D(histNamesTru[i],"x");
            
            // Set line color and thickness
            hist1DUnf->SetLineColor(colors[i]); // Set line color
            hist1DUnf->SetLineWidth(2);         // Set line thickness to 2
            hist1DUnf->SetMarkerColor(colors[i]);
            
            hist1DUnf->Divide(hist1DUnf_base);
            
            //        hist1DUnf->Divide(hist1DTru);
            // Draw the histogram
            
            hist1DUnf->GetYaxis()->SetTitle("Unfolded/iter4");
            hist1DUnf->SetStats(0);
            hist1DUnf->GetYaxis()->SetTitle("Unfolded/iter4");
            hist1DUnf->Draw("SAME");
            
            
            
            // Add histogram to legend
            legend->AddEntry(hist1DUnf, histNamesUnf[i], "l");
            
            // Close the file
            //        file->Close();
        }
        
        // Draw the legend
        legend->Draw();
        
    }
   
///##########################Studying effects of folding
 
    if(ifFold){
        cout<<"folded histogram "<<endl;
        // Loop over the files, retrieve and rename histograms
        for (int i = 0; i < num; i++) {
            // Open the file
            TFile *file = TFile::Open(filenames[i]);
            if (!file || file->IsZombie()) {
                std::cerr << "Error opening file: " << filenames[i] << std::endl;
                continue;
            }
            
            // Get the 3D histogram
            TH3 *hist3DUnf = (TH3*)file->Get(originalHistNameUnf);
            if (!hist3DUnf) {
                std::cerr << "Error retrieving histogram: " << originalHistNameUnf << " from file: " << filenames[i] << std::endl;
                file->Close();
                continue;
            }
            hist3DUnf->SetName(histNamesUnf3d[i]);
            hist3DUnf->GetXaxis()->SetTitle(str);
            hist3DUnf->GetYaxis()->SetTitle(str4);
            hist3DUnf->GetZaxis()->SetTitle("wt");
            hist3DUnf->SetMarkerStyle(markers[i]);
            hist3DUnf->SetMarkerColor(colors[i]);
            hist3DUnf->Sumw2();
            
            
            TH3 *hist3DTru = (TH3*)file->Get(originalHistNameTru);
            if (!hist3DTru) {
                std::cerr << "Error retrieving histogram: " << originalHistNameTru << " from file: " << filenames[i] << std::endl;
                file->Close();
                continue;
            }
            hist3DTru->SetName(histNamesTru3d[i]);
            
            cout<<"hher"<<endl;
            // Project the 3D histogram to a 1D histogram along the X-axis
            TH1 *hist1DUnf; //= hist3DUnf->ProjectionX(histNamesUnf[i],2,2,10,10);
            TH1 *hist1DTru; //= hist3DTru->ProjectionX(histNamesTru[i],2,2,10,10);
            
            if(projection=="x"){
                //            hist1DUnf = hist3DUnf->ProjectionX(histNamesUnf[i],2,5,10,10);
                //            hist1DTru = hist3DTru->ProjectionX(histNamesTru[i],2,5,10,10);
                
                //            hist3DUnf->GetYaxis()->SetRangeUser(pt1,pt2);
                //            hist3DUnf->GetYaxis()->SetRangeUser(pt1,pt2);
                hist1DUnf = hist3DUnf->Project3D("x");
                hist1DTru = hist3DTru->Project3D("x");
            }
            else if(projection=="y")
            {
                //            hist1DUnf = hist3DUnf->ProjectionY(histNamesUnf[i],2,2,10,10);
                //            hist1DTru = hist3DTru->ProjectionY(histNamesTru[i],2,2,10,10);
                
                //            hist3DUnf->GetXaxis()->SetRangeUser(rlow,rhigh);
                //            hist3DUnf->GetXaxis()->SetRangeUser(rlow,rhigh);
                hist1DUnf = hist3DUnf->Project3D("y");
                hist1DTru = hist3DTru->Project3D("y");
            }
            else{
                //            hist1DUnf = hist3DUnf->ProjectionY(histNamesUnf[i],2,2,10,10);
                //            hist1DTru = hist3DTru->ProjectionY(histNamesTru[i],2,2,10,10);
                
                //            hist3DUnf->GetXaxis()->SetRangeUser(rlow,rhigh);
                //            hist3DUnf->GetXaxis()->SetRangeUser(rlow,rhigh);
                //            hist3DUnf->GetYaxis()->SetRangeUser(pt1,pt2);
                //            hist3DUnf->GetYaxis()->SetRangeUser(pt1,pt2);
                hist1DUnf = hist3DUnf->Project3D("z");
                hist1DTru = hist3DTru->Project3D("z");
            }
            
            hist1DUnf->SetName(histNamesUnf[i]);
            hist1DTru->SetName(histNamesTru[i]);
            //        hist1DUnf->GetYaxis()->SetRangeUser(0.5,1.5);
            //        hist1DTru->GetYaxis()->SetRangeUser(0.5,1.5);
            //        TH1 *hist1DUnf = hist3DUnf->Project3D(histNamesUnf[i],"x");
            //        TH1 *hist1DTru = hist3DTru->Project3D(histNamesTru[i],"x");
            
            // Set line color and thickness
            hist1DUnf->SetLineColor(colors[i]); // Set line color
            hist1DUnf->SetLineWidth(2);         // Set line thickness to 2
            hist1DUnf->SetMarkerColor(colors[i]);
            
            hist1DUnf->Divide(hist1DTru);
            // Draw the histogram
            
            hist1DUnf->GetYaxis()->SetTitle("folded/raw");
            hist1DUnf->SetStats(0);
            hist1DUnf->GetYaxis()->SetTitle("folded/raw");
            hist1DUnf->Draw("SAME");
            
            // Add histogram to legend
            legend->AddEntry(hist1DUnf, histNamesUnf[i], "l");
            
            cout<<hist1DUnf->GetNbinsX()<<" and "<<hist1DTru->GetNbinsX()<<endl;
            // Close the file
            //        file->Close();
        }
        
        // Draw the legend
        legend->Draw();
        
    }
   
///##########################Studying effects of weights
  
    if(ifCheckProjected){
        
//        const char* inputFile = "/Users/ar2545/Downloads/E3C_dataOct24.root";
//        std::string histName3D_unfolded = "Opt_Un_e3c";
//        std::string histName3D = "Opt_Un_e3c";
//        std::string histName2D = "e3c_pt_hist_unf";
        
//        std::string histName3D_unfolded = "h3_triv";
//        std::string histName3D = "h3_true_eecTriv";
//        std::string histName2D = "eec_pt_histTriv";
        
        const char* inputFile = filenames[3];
//        std::string histName3D_unfolded = "h3_unfoldedSplit";
//        std::string histName3D = "h3_trueSplit";
//        std::string histName2D = "eec_pt_histSplit";
        
        std::string histName3D_unfolded = "h3_fold";
        std::string histName3D = "h3_rawSplit";
        std::string histName2D = "eec_pt_hist_detSplit";
        
 
        bool ifInterpolate = true;
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
        
//        TH1D* hjet = (TH1D*)f->Get("jet_pt_hist");
//        hjet->Scale(1.,"width");
//        
//        TCanvas* cjet;
//        hjet->Draw();
        
        TH3D* h3 = (TH3D*)f->Get(histName3D.c_str());
        h3->GetXaxis()->SetTitle(str4);
        h3->GetYaxis()->SetTitle(str);
        h3->GetZaxis()->SetTitle("wt");
        
        TH3D* h3_unfolded = (TH3D*)f->Get(histName3D_unfolded.c_str());
        h3_unfolded->GetXaxis()->SetTitle(str4);
        h3_unfolded->GetYaxis()->SetTitle(str);
        h3_unfolded->GetZaxis()->SetTitle("wt");
        
        if (!h3_unfolded) {
            std::cerr << "Error retrieving histogram: " << histName3D_unfolded << std::endl;
            f->Close();
            return;
        }
        
         
        int pt1jet = h3->GetYaxis()->FindBin(pt1);
        int pt2jet = h3->GetYaxis()->FindBin(pt2);
        
        cout<<"Bin edges are: "<<h3->GetYaxis()->GetBinLowEdge(pt1jet)<<" and "<<h3->GetYaxis()->GetBinUpEdge(pt2jet);
        
        h3->GetYaxis()->SetRange(pt1jet,pt2jet);
        h3_unfolded->GetYaxis()->SetRange(pt1jet,pt2jet);
        // Get the 2D histogram
        TH2D *h2 = (TH2D*)h3->Project3D("zx"); //Rl vs wt histogram
        TH2D *h2_unfolded = (TH2D*)h3_unfolded->Project3D("zx"); //Rl vs wt histogram
        
        // Create the 1D histogram
        int nBinsX = h3->GetNbinsX();
        double xMin = h3->GetXaxis()->GetXmin();
        double xMax = h3->GetXaxis()->GetXmax();
        
        TH1D *h1 = (TH1D*)h2->ProjectionX("h1");
        h1->Reset();
        
        TH1D *h1_unfolded = (TH1D*)h2_unfolded->ProjectionX("h1_unfolded");
        h1_unfolded->Reset();
        
        TCanvas* c;
//        h3->Draw();
        
       if(!ifInterpolate)
       {
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
           
           
           for (int i = 1; i <= nBinsX; ++i) {
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
         
           for (int i = 1; i <= nBinsX; ++i){
               int nBinsY = h2->GetNbinsY();
               double weightedSum = 0;
               TH1D *h1_wt = (TH1D*)h2->ProjectionY(Form("h1_wt_%i",i),i,i);
               //               h1_wt->Draw();
               
               for (int j = 1; j <= nBinsY; ++j) {
                   double weight = h1_wt->GetBinCenter(j);
                   double bin_content = h1_wt->Interpolate(weight);
                   weightedSum += weight*bin_content;
               }
               h1->SetBinContent(i, weightedSum);//this is closest to filling on the fly, not the avg
           }
           
           for (int i = 1; i <= nBinsX; ++i){
               int nBinsY = h2_unfolded->GetNbinsY();
               double weightedSum = 0;
               TH1D *h1_wt_unfolded = (TH1D*)h2_unfolded->ProjectionY(Form("h1_wt_unfolded_%i",i),i,i);
               for (int j = 1; j <= nBinsY; ++j) {
                   double weight = h1_wt_unfolded->GetBinCenter(j);
                   double bin_content = h1_wt_unfolded->Interpolate(weight);
                   weightedSum += weight*bin_content;
                   
                   if(bin_content!=0 && h1_wt_unfolded->GetBinContent(j)!=bin_content ){
                       cout<<"bin center: " <<setprecision(10)<<weight<<" and " <<h1_wt_unfolded->GetBinContent(j)<<" and interpolated "<<setprecision(10)<<bin_content<<endl;
                   }
               }
               h1_unfolded->SetBinContent(i, weightedSum);//this is closest to filling on the fly, not the avg
           }
       }
        
        h1->Scale(1.,"width");
        h1->SetLineColor(kGreen);
        h1_unfolded->Scale(1.,"width");
        h1_unfolded->SetLineColor(kRed);
        
        TH2D* h2_nom;
        h2_nom = (TH2D*)f->Get(histName2D.c_str());
    
        TH1D* h1_nom = (TH1D*)h2_nom->ProjectionX(Form("h1nom"),pt1jet,pt2jet);
        h1_nom->Scale(1.,"width");
         h1_nom->SetLineColor(kMagenta);
        cout<<h1_nom->GetNbinsX()<<" and "<<nBinsX<<endl;
//        h1->Divide(h1_nom);
       
        
//        TH1D* h1_clone_unfolded = (TH1D*)h1_unfolded->Clone("denominator");
//        h1_unfolded->Add(h1_nom,-1);
//        h1_unfolded->Divide(h1_clone_unfolded);
//
        h1_unfolded->GetYaxis()->SetTitle("percent diff");
         h1->GetYaxis()->SetTitle("3D projection/2D projection");
      
       h1->Divide(h1_nom);
       h1_unfolded->Divide(h1_nom);
        h1_nom->Draw("same,HIST");
        h1->Draw();
        h1_unfolded->Draw("same");

    legend2->AddEntry(h1, "true", "l");
    legend2->AddEntry(h1_unfolded, "unfolded", "l");
    
//    legend2->AddEntry(h1, "raw", "l");
//    legend2->AddEntry(h1_unfolded, "folded", "l");
//    legend2->AddEntry(h1_nom, "nom", "l");
//    legend2->Draw("same");

        
    }

    
    
}

