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


void gettingwtCheckEEC(int pt1, int pt2) {
    
        const char* inputFile = "/Users/ar2545/Downloads/EEC_dataNov4.root";
    std::string histName3D_unfolded = "Opt_Un_eec_unf";
    std::string histName3D = "Opt_Un_eec";
    std::string histName2D = "eec_pt_hist_unf";
    
    
//    std::string histName3D_unfolded = "h3_triv";
//    std::string histName3D = "h3_true_eecTriv";
//    std::string histName2D = "eec_pt_histTriv";
//    const char* inputFile = "/Users/ar2545/Downloads/UnfoldingResultsTrivEEC_4EEC.root";
    
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

    TH3D* h3 = (TH3D*)f->Get(histName3D.c_str());
    h3->GetXaxis()->SetTitle(str4);
    h3->GetYaxis()->SetTitle(str);
    h3->GetZaxis()->SetTitle("wt");
    
    int pt1jet = h3->GetYaxis()->FindBin(pt1);
    int pt2jet = h3->GetYaxis()->FindBin(pt2);
    
    h3->GetYaxis()->SetRange(pt1jet,pt2jet);
    
    TH2D *h2 = (TH2D*)h3->Project3D("zx"); //Rl vs wt histogram
    TH1D *h1 = (TH1D*)h2->ProjectionX("h1");
    h1->Reset();
 
    
       // Create the 1D histogram
    int nBinsX = h2->GetNbinsX();
    double xMin = h2->GetXaxis()->GetXmin();
    double xMax = h2->GetXaxis()->GetXmax();
 
    TH3D* h3_unfolded = (TH3D*)f->Get(histName3D_unfolded.c_str());
    h3_unfolded->GetXaxis()->SetTitle(str4);
    h3_unfolded->GetYaxis()->SetTitle(str);
    h3_unfolded->GetZaxis()->SetTitle("wt");
    
    if (!h3_unfolded) {
        std::cerr << "Error retrieving histogram: " << histName3D_unfolded << std::endl;
        f->Close();
        return;
    }
     h3_unfolded->GetYaxis()->SetRange(h3_unfolded->GetYaxis()->FindBin(pt1),h3_unfolded->GetYaxis()->FindBin(pt2));

//    // Get the 2D histogram
//    TH2D *h2_unfolded = (TH2D*)h3_unfolded->Project3D("zx"); //Rl vs wt histogram
//    TH1D *h1_unfolded = (TH1D*)h2_unfolded->ProjectionX("h1_unfolded");
//    h1_unfolded->Reset();
 
 
///############################SET UP WHICH KIND OF PROJECTION YOU WANT TO DO------------------------------------------------------------

    Double_t from = -2;
    Double_t to = -0.39794;
    Int_t bins = 20;
    Double_t width = (to-from)/bins;
    Double_t new_bins[21] = {};
    for (Int_t i = 0; i <= bins; i++)
    {
        new_bins[i] = TMath::Power(10.0, from + i * width);
        cout<<new_bins[i]<<endl;
    }
    
// double RvaluesReb[] = {0.0001,0.01,0.0120255,0.0144613,0.0173904,0.0209128,0.0251487,0.0302425,
// 0.0363681,0.0437345,0.0525929,0.0632456,0.0760559,0.091461,0.109986,0.132264,0.159054,0.191271,0.230012,0.276601,0.332627,0.4,0.524807,1};
// 
 double RvaluesReb[] = {0.0001,0.01,0.0144613,0.0173904,0.0209128,0.0251487,0.0302425,
 0.0363681,0.0437345,0.0525929,0.0632456,0.0760559,0.091461,0.109986,0.132264,0.159054,0.191271,0.230012,0.276601,0.332627,0.4,0.524807,1};


//    constexpr double Rebvalues_e3c[] = {
//        0.0001, 0.000524807, 0.001, 0.00151356,0.0020893, 0.00251189, 0.00301995, 0.00346737,
//        0.00416869,0.00457088, 0.00501187, 0.00549541,
//        0.0060256, 0.00660693, 0.00758578, 0.00831764,
//        0.00912011, 0.0104713, 0.0138038, 0.018197, 0.0263027
//    };
    
     constexpr double Rebvalues_e3c[] = {
        0.0005,0.001,0.00346737,0.00912011, 0.0138038,0.0165959,0.0190546,
        0.0288403,
        0.0346737,
        0.0457088, 0.060256,
        0.0794328,
        0.104713, 0.138038, 0.18197,
        0.26};
    
//constexpr double Rebvalues_e3c[] = {
//        0.0001, 0.000104713, 0.000109648, 0.000114815, 0.000120226, 0.000125893, 0.000131826,
//        0.000138038, 0.000144544, 0.000151356, 0.000158489, 0.000165959, 0.00017378, 0.00018197,
//        0.000190546, 0.000199526, 0.00020893, 0.000218776, 0.000229087, 0.000239883, 0.000251189,
//        0.000263027, 0.000275423, 0.000288403, 0.000301995, 0.000316228, 0.000331131, 0.000346737,
//        0.000363078, 0.000380189, 0.000398107, 0.000416869, 0.000436516, 0.000457088, 0.00047863,
//        0.000501187, 0.000524807, 0.000549541, 0.00057544, 0.00060256, 0.000630957, 0.000660693,
//        0.000691831, 0.000724436, 0.000758578, 0.000794328, 0.000831764, 0.000870964, 0.000912011,
//        0.000954993, 0.001, 0.00104713, 0.00109648, 0.00114815, 0.00120226, 0.00125893, 0.00131826,
//        0.00138038, 0.00144544, 0.00151356, 0.00158489, 0.00165959, 0.0017378, 0.0018197, 0.00190546,
//        0.00199526, 0.0020893, 0.00218776, 0.00229087, 0.00239883, 0.00251189, 0.00263027, 0.00275423,
//        0.00288403, 0.00301995, 0.00316228, 0.00331131, 0.00346737, 0.00363078, 0.00380189, 0.00398107,
//        0.00416869, 0.00436516, 0.00457088, 0.0047863, 0.00501187, 0.00524807, 0.00549541, 0.0057544,
//        0.0060256, 0.00630957, 0.00660693, 0.00691831, 0.00724436, 0.00758578, 0.00794328, 0.00831764,
//        0.00870964, 0.00912011, 0.00954993, 0.01, 0.0104713, 0.0109648, 0.0114815, 0.0120226, 0.0125893,
//        0.0131826, 0.0138038, 0.0144544, 0.0151356, 0.0158489, 0.0165959, 0.017378, 0.018197, 0.0190546,
//        0.0199526, 0.020893, 0.0218776, 0.0229087, 0.0239883, 0.0251189, 0.0263027, 0.0275423, 0.0288403,
//        0.0301995, 0.0316228, 0.0331131, 0.0346737, 0.0363078, 0.0380189, 0.0398107, 0.0416869, 0.0436516,
//        0.0457088, 0.047863, 0.0501187, 0.0524807, 0.0549541, 0.057544, 0.060256, 0.0630957, 0.0660693,
//        0.0691831, 0.0724436, 0.0758578, 0.0794328, 0.0831764, 0.0870964, 0.0912011, 0.0954993, 0.1,
//        0.104713, 0.109648, 0.114815, 0.120226, 0.125893, 0.131826, 0.138038, 0.144544, 0.151356,
//        0.158489, 0.165959, 0.17378, 0.18197, 0.190546, 0.199526, 0.20893, 0.218776, 0.229087, 0.239883,
//        0.251189, 0.263027, 0.275423, 0.288403, 0.301995, 0.316228, 0.331131, 0.346737, 0.363078,
//        0.380189, 0.398107, 0.416869, 0.436516, 0.457088, 0.47863, 0.501187, 0.524807, 0.549541,
//        0.57544, 0.60256, 0.630957, 0.660693, 0.691831, 0.724436, 0.758578, 0.794328, 0.831764,
//        0.870964, 0.912011, 0.954993, 1.0
//    };

//
//    constexpr double Rebvalues_e3c[] = {
//        0.0001, 0.000524807, 0.001, 0.00151356,0.0020893, 0.00251189, 0.00301995, 0.00346737,
//        0.00416869,0.00457088, 0.00501187, 0.00549541,
//        0.0060256, 0.00660693, 0.00724436, 0.00758578, 0.00831764,
//        0.00912011, 0.00954993, 0.01, 0.0104713, 0.0138038, 0.018197, 0.0263027
//    };
//  double Rebvalues_e3c[] = {
//        0.0005,0.001,0.00346737,0.00912011, 0.0138038,0.0165959,0.0190546,
//        0.0288403,
//        0.0346737,
//        0.0457088, 0.060256,
//        0.0794328,
//        0.104713, 0.138038, 0.18197,
//        0.26};


//

//    constexpr double Rebvalues_e3c[] = {
//        0.000102329,0.000121619,0.000251189, 0.000501187,
//        0.00151356, 0.00301995, 0.00506991, 0.00716143, 0.00912011, 0.0101158, 0.0104713, 0.0108393, 0.0128825, 0.0133352, 0.0138038, 0.0142889, 0.0147911, 0.0153109, 0.0158489,
//        0.0164059, 0.0169824, 0.0175792, 0.018197, 0.0188365, 0.0194984, 0.0201837, 0.025704, 0.0305492,
//        0.0363078, 0.0402717, 0.0462381, 0.0512861, 0.0568853,0.0653131, 0.0724436, 0.0749894,
//        0.0803526, 0.0860994, 0.0954993, 0.0988553,
//        0.154882
//    };
 

    double numbins_x_e3c = sizeof(RvaluesReb)/sizeof(RvaluesReb[0]) - 1;
    double numbins_y_e3c = sizeof(Rebvalues_e3c)/sizeof(Rebvalues_e3c[0]) - 1;
    
    TH2D* h2_unfolded = new TH2D("h2_reb", "Rebinned 2D Histogram (Both Axes)", numbins_x_e3c, RvaluesReb, numbins_y_e3c, Rebvalues_e3c);
    
    for (int i = 1; i <= h2->GetNbinsX(); ++i) {
        for (int j = 1; j <= h2->GetNbinsY(); ++j) {
            double content = h2->GetBinContent(i, j);
            double error = h2->GetBinError(i, j);
            double x = h2->GetXaxis()->GetBinCenter(i);
            double y = h2->GetYaxis()->GetBinCenter(j);
            h2_unfolded->Fill(x, y, content);
        }
    }
    TH1D *h1_unfolded = (TH1D*)h2_unfolded->ProjectionX("h1_unfolded");
    h1_unfolded->Reset();
    
    TCanvas *canvas1 = new TCanvas();
    canvas1->SetLogx();
    canvas1->SetLogy();
    h2_unfolded->SetStats(0);
    h2_unfolded->GetYaxis()->SetTitle("wt");
    h2_unfolded->GetXaxis()->SetTitle(str4);
    h2_unfolded->Draw("text");
    h2_unfolded->Draw("colz");
    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.04);
    latex1.SetTextFont(42);
    latex1.DrawLatex(0.15,0.80 ,Form("%.1d < %s < %.1d GeV/c, EEC", pt1, str, pt2));

///############################SET UP WHICH KIND OF PROJECTION YOU WANT TO DO------------------------------------------------------------

    bool ifInterpolate = false;
    
 
    TCanvas* c = new TCanvas("c1Dhist","1dhist",800,650);
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
    
    h1->Scale(1.,"width");
    h1->SetLineColor(kGreen);
    
    h1_unfolded->Scale(1.,"width");
    h1_unfolded->SetLineColor(kRed);
    
    TH2D* h2_nom;
    h2_nom = (TH2D*)f->Get(histName2D.c_str());
    
    TH1D* h1_nom = (TH1D*)h2_nom->ProjectionX(Form("h1nom"),h2_nom->GetYaxis()->FindBin(pt1),h2_nom->GetYaxis()->FindBin(pt2));
    h1_nom->Scale(1.,"width");
    //        h1->Divide(h1_nom);
    h1_unfolded->Divide(h1_nom);
          
//    h1_unfolded->Divide(h1);
//     h1->Divide(h1_nom);

//    h1->Draw();
//    h1_nom->Draw("same");
    h1_unfolded->Draw("same");
    
    TFile* fileOut = TFile::Open("ProjectedOutputEEC.root")
  
}
