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


void CompareResultsNIter(int pt1, int pt2) {
    std::string unfResults[4] = {
        "/Users/ar2545/Downloads/UnfoldingResultsDataEEC_1EECnew6.root",
        "/Users/ar2545/Downloads/UnfoldingResultsDataEEC_4EECnew6.root",
        "/Users/ar2545/Downloads/UnfoldingResultsDataEEC_7EECnew6.root",
        "/Users/ar2545/Downloads/UnfoldingResultsDataEEC_10EECnew6.root"
        
//        "/Users/ar2545/Downloads/UnfoldingResultsDataEEC_1EECnew2.root",
//        "/Users/ar2545/Downloads/UnfoldingResultsDataEEC_4EECnew2.root",
//        "/Users/ar2545/Downloads/UnfoldingResultsDataEEC_7EECnew2.root",
//        "/Users/ar2545/Downloads/UnfoldingResultsDataEEC_10EECnew2.root"
    };
    
    std::string rawDat = "/Users/ar2545/Downloads/EEC_dataNov22.root";
    
    std::string unfDat[4] = {
        "/Users/ar2545/Desktop/Unfolding/ProjectedOutputEECNew6_1.root",
        "/Users/ar2545/Desktop/Unfolding/ProjectedOutputEECNew6_4.root",
        "/Users/ar2545/Desktop/Unfolding/ProjectedOutputEECNew6_7.root",
        "/Users/ar2545/Desktop/Unfolding/ProjectedOutputEECNew6_10.root"
        
//         "/Users/ar2545/Desktop/Unfolding/ProjectedOutputEECNew2_1.root",
//        "/Users/ar2545/Desktop/Unfolding/ProjectedOutputEECNew2_4.root",
//        "/Users/ar2545/Desktop/Unfolding/ProjectedOutputEECNew2_7.root",
//        "/Users/ar2545/Desktop/Unfolding/ProjectedOutputEECNew2_10.root"
    };
    
    bool ifarea = false;
    int num = 4;
    double rlow = 0.01;
    double rhigh = 0.4;
    const char *str = "#it{p}_{T,jet}^{ch}";
    const char *str1 = "#it{p}_{T,jet-min}^{ch}";
    const char *str2 = "anti-#it{k}_{T}";
    const char *str3 = "#sqrt{s} = 13 TeV";
    const char *str4 = "#it{R}_{L}";
    const char *str8 = "#it{#eta}_{jet}";
    const char *str11 = "#it{p}_{T}^{const}";
    
    TLatex latex;latex.SetNDC ();latex.SetTextSize(0.045);latex.SetTextFont(42);
    TLatex latex1;latex1.SetNDC ();latex1.SetTextSize(0.045);latex1.SetTextFont(42);
    TLatex latex2;latex2.SetNDC ();latex2.SetTextSize(0.045);latex2.SetTextFont(42);
    TLatex latex3;latex3.SetNDC ();latex3.SetTextSize(0.045);latex3.SetTextFont(42);
    TLatex latex4;latex4.SetNDC ();latex4.SetTextSize(0.045);latex4.SetTextFont(42);
    TLatex latex5;latex5.SetNDC ();latex5.SetTextSize(0.045);latex5.SetTextFont(42);

    Color_t colors[4] = {kRed, kBlue, kMagenta, kGreen+2};
    Color_t markers[4] = {20, 25, 29, 33};
    
      // Create a legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.045);
    
    
//###############Get jet histograms
    TFile *fjetmc = TFile::Open(Form("%s",unfResults[1].c_str()));
    TH1 *jet_den = (TH1*)fjetmc->Get("h1_true");
    jet_den->SetName("h1_den");
    TH1 *jet_num = (TH1*)fjetmc->Get("h1_reco");
    TH1 *jet_dat = (TH1*)fjetmc->Get("jet_pt_hist");
    
//###############Bin-by-bin processing (similar to original code)
    TFile *fdat = TFile::Open(Form("%s",rawDat.c_str()));
    TH2 *h2_dat = (TH2*)fdat->Get("eec_pt_hist_unf");
    h2_dat->SetName("eec_pt_hist_dat");
    
    TFile *fmc = TFile::Open(Form("%s",unfResults[1].c_str()));
    TH2 *h2_den = (TH2*)fmc->Get("eec_pt_hist");
    TH2 *h2_num = (TH2*)fmc->Get("eec_pt_hist_det");
    
    int pt1jet = h2_den->GetYaxis()->FindBin(pt1);
    int pt2jet = h2_den->GetYaxis()->FindBin(pt2);
    
    int pt1jet_num = h2_num->GetYaxis()->FindBin(pt1);
    int pt2jet_num = h2_num->GetYaxis()->FindBin(pt2);
    
    int pt1jet_dat = h2_dat->GetYaxis()->FindBin(pt1);
    int pt2jet_dat = h2_dat->GetYaxis()->FindBin(pt2);
    
    TH1* h1_dat = (TH1*)h2_dat->ProjectionX("h1_dat", pt1jet_dat, pt2jet_dat);
    TH1* h1_den = (TH1*)h2_den->ProjectionX("h1_den",pt1jet,pt2jet);
    TH1* h1_num = (TH1*)h2_num->ProjectionX("h1_num",pt1jet_num,pt2jet_num);

    if(ifarea){
        h1_den->Scale(1./h1_den->Integral(),"width");
        h1_num->Scale(1./h1_num->Integral(),"width");
        h1_dat->Scale(1./h1_dat->Integral(),"width");
    }
    else{
        h1_den->Scale(1./jet_den->Integral(jet_den->GetXaxis()->FindBin(pt1),jet_den->GetXaxis()->FindBin(pt2)),"width");
        h1_num->Scale(1./jet_num->Integral(jet_num->GetXaxis()->FindBin(pt1),jet_num->GetXaxis()->FindBin(pt2)),"width");
        h1_dat->Scale(1./jet_dat->Integral(jet_dat->GetXaxis()->FindBin(pt1),jet_dat->GetXaxis()->FindBin(pt2)),"width");
    }
    
    h1_num->Divide(h1_den);
    h1_dat->Divide(h1_num);
//################################################

     // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Overlay Histograms", 600, 800);
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
    topPad->cd();
    gPad->SetLogx();
    h1_dat->SetStats(0);
    h1_dat->SetTitle(Form("Excluding jets with %s < 15 GeV/#it{c}",str));
    h1_dat->GetYaxis()->SetRangeUser(0.5,6);
    h1_dat->GetYaxis()->SetTitle("Njet Normalized EEC");
    h1_dat->Draw();
    h1_dat->GetXaxis()->SetRangeUser(0.01,0.4);
    h1_dat->GetXaxis()->SetLabelSize(0.09);
    h1_dat->GetXaxis()->SetTitleSize(0.09);
    h1_dat->GetXaxis()->SetTitleOffset(1.07);
    legend->AddEntry(h1_dat, "bin by bin");
    
    
    for (int i = 0; i < num; i++) {
    
        int iterNum = -1;
        if(i==0) {iterNum = 1;}
        else if(i==1){iterNum = 4; }
        else if(i==2){ iterNum = 7;}
        else { iterNum = 10;}
        
        
        TFile *fjetmc = TFile::Open(unfResults[i].c_str());
        if (!fjetmc || fjetmc->IsZombie()) {
            std::cerr << "Error opening file: " << unfResults[i] << std::endl;
            continue;
        }

        TH1 *jet_unf = (TH1*)fjetmc->Get("h1_unfolded");
        if (!jet_den || !jet_num || !jet_dat || !jet_unf) {
            std::cerr << "Error retrieving histograms from file: " << unfResults[i] << std::endl;
            fjetmc->Close();
            continue;
        }

        // Apply kinematic efficiency
        TH1* h1_true = (TH1*)fjetmc->Get("h1_true");
        TH1* h1_fulleff = (TH1*)fjetmc->Get("h1_fulleff");
        if (!h1_true || !h1_fulleff) {
            std::cerr << "Error retrieving efficiency histograms from file: " << unfResults[i] << std::endl;
            fjetmc->Close();
            continue;
        }
        TH1* h1_eff = (TH1*)h1_true->Clone("h1_eff");
        h1_eff->Divide(h1_fulleff);
        jet_unf->Divide(h1_eff);

        // Open corresponding unfDat file
        TFile *fUnfolded = TFile::Open(unfDat[i].c_str());
        if (!fUnfolded || fUnfolded->IsZombie()) {
            std::cerr << "Error opening file: " << unfDat[i] << std::endl;
            fjetmc->Close();
            continue;
        }

        TH1 *h1_unf = (TH1*)fUnfolded->Get(Form("EEC_unfolded_%i_%i", pt1, pt2));
        if (!h1_unf) {
            std::cerr << "Error retrieving unfolded histogram from file: " << unfDat[i] << std::endl;
            fjetmc->Close();
            fUnfolded->Close();
            continue;
        }

     
        if (ifarea) {
            
            h1_unf->Scale(1./h1_unf->Integral(),"width");
        } else {
            
            h1_unf->Scale(1. / jet_unf->Integral(jet_unf->GetXaxis()->FindBin(pt1), jet_unf->GetXaxis()->FindBin(pt2)), "width");
        }

        


        topPad->cd();
        h1_unf->SetLineColor(colors[i]);
        h1_unf->SetMarkerColor(colors[i]);
        h1_unf->SetMarkerStyle(markers[i]);
//        h1_unf->SetLineWidth(2);
        h1_unf->Draw("same");
        
        legend->AddEntry(h1_unf, Form("iter = %i",iterNum),"l");

        bottomPad->cd();
        gPad->SetLogx();
        TH1* h1_ratio = (TH1*)h1_unf->Clone(Form("h1_ratio_%d", i));
        h1_ratio->Divide(h1_dat);
        h1_ratio->Draw("same");
        h1_ratio->SetStats(0);
        h1_ratio->SetTitle("");
        h1_ratio->SetLineColor(colors[i]);
        h1_ratio->GetYaxis()->SetRangeUser(0.8,1.22);
        h1_ratio->GetYaxis()->SetLabelSize(0.06);
        h1_ratio->GetYaxis()->SetTitleSize(0.08);
        h1_ratio->GetYaxis()->SetTitleOffset(0.53);
        
        h1_ratio->GetXaxis()->SetRangeUser(0.01,0.4);
        h1_ratio->GetXaxis()->SetLabelSize(0.09);
        h1_ratio->GetXaxis()->SetTitleSize(0.09);
        h1_ratio->GetXaxis()->SetTitleOffset(1.07);
        
        h1_ratio->GetYaxis()->SetTitle("Unf/BbB");
       
       
        // Cleanup
//        fjetmc->Close();
//        fUnfolded->Close();
    }
    topPad->cd();
    legend->Draw();
    latex1.DrawLatex(0.15,0.80 ,Form("%.1d < %s < %.1d GeV/c", pt1, str, pt2));
}
