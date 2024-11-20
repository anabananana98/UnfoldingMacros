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


void CompareResults(int pt1, int pt2) {
    
    std::string unfResults = "/Users/ar2545/Downloads/UnfoldingResultsDataEEC_4EECnew2.root";
    std::string rawDat = "/Users/ar2545/Downloads/EEC_dataNov18.root";
    std::string UnfDat = "/Users/ar2545/Desktop/Unfolding/ProjectedOutputEECNew2.root";
//    std::string jetFileMC = "/Users/ar2545/Downloads/Unfolding/Unfold3DData_Nov18.root";
    
    bool ifarea = false;
    int num=4;
    double rlow = 0.01;
    double rhigh = 0.4;
    const char *str = "#it{p}_{T,jet} GeV/#it{c}";
    const char *str1 = "#it{p}_{T,min}";
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

    Color_t colors[4] = {kRed, kBlue, kMagenta, kCyan};
    
    // Array of markers for histograms
    Color_t markers[4] = {20, 25, 29, 33};
    
   
   
    TFile *fjetmc = TFile::Open(Form("%s",unfResults.c_str()));
    TH1 *jet_den = (TH1*)fjetmc->Get("h1_true");
    jet_den->SetName("h1_den");
    TH1 *jet_num = (TH1*)fjetmc->Get("h1_reco");
    TH1 *jet_dat = (TH1*)fjetmc->Get("jet_pt_hist");
    TH1 *jet_unf = (TH1*)fjetmc->Get("h1_unfolded");
    
    //Apply kinematic efficiency to unfolded jet pT hist
    TH1* h1_true = (TH1*)fjetmc->Get("h1_true");
    TH1* h1_fulleff = (TH1*)fjetmc->Get("h1_fulleff");
    TH1* h1_eff = (TH1*)h1_true->Clone("h1_eff");
    h1_eff->Divide(h1_fulleff);
    jet_unf->Divide(h1_eff);
    
    TFile *fmc = TFile::Open(Form("%s",unfResults.c_str()));
    TH2 *h2_den = (TH2*)fmc->Get("eec_pt_hist");
    TH2 *h2_num = (TH2*)fmc->Get("eec_pt_hist_det");
    
    TFile *fdat = TFile::Open(Form("%s",rawDat.c_str()));
    TH2 *h2_dat = (TH2*)fdat->Get("eec_pt_hist_unf");
    h2_dat->SetName("eec_pt_hist_dat");
    
    TFile *fUnfolded = TFile::Open(Form("%s",UnfDat.c_str()));
    TH1 *h1_unf = (TH1*)fUnfolded->Get(Form("EEC_unfolded_%i_%i",pt1,pt2));
//    TH1 *h1_unf = (TH1*)fUnfolded->Get(Form("EEC_unfoldedInterpolated_%i_%i",pt1,pt2));
//    h2_unf->SetName("eec_pt_hist_dat");
    
    int pt1jet = h2_den->GetYaxis()->FindBin(pt1);
    int pt2jet = h2_den->GetYaxis()->FindBin(pt2);
    
    int pt1jet_num = h2_num->GetYaxis()->FindBin(pt1);
    int pt2jet_num = h2_num->GetYaxis()->FindBin(pt2);
    
    int pt1jet_dat = h2_dat->GetYaxis()->FindBin(pt1);
    int pt2jet_dat = h2_dat->GetYaxis()->FindBin(pt2);
    
    cout<<pt1jet<<" and "<<pt2jet<<endl;
    cout<<pt1jet_num<<" and "<<pt2jet_num<<endl;
    cout<<pt1jet_dat<<" and "<<pt2jet_dat<<endl;
    
    TH1* h1_den = (TH1*)h2_den->ProjectionX("h1_den",pt1jet,pt2jet);
    TH1* h1_num = (TH1*)h2_num->ProjectionX("h1_num",pt1jet_num,pt2jet_num);
    TH1* h1_dat = (TH1*)h2_dat->ProjectionX("h1_dat",pt1jet_dat,pt2jet_dat);
    
    if(ifarea){
        h1_den->Scale(1./h1_den->Integral(),"width");
        h1_num->Scale(1./h1_num->Integral(),"width");
        h1_dat->Scale(1./h1_dat->Integral(),"width");
        h1_unf->Scale(1./h1_unf->Integral(),"width");
    }
    else{
        h1_den->Scale(1./jet_den->Integral(jet_den->GetXaxis()->FindBin(pt1),jet_den->GetXaxis()->FindBin(pt2)),"width");
        h1_num->Scale(1./jet_num->Integral(jet_num->GetXaxis()->FindBin(pt1),jet_num->GetXaxis()->FindBin(pt2)),"width");
        h1_dat->Scale(1./jet_dat->Integral(jet_dat->GetXaxis()->FindBin(pt1),jet_dat->GetXaxis()->FindBin(pt2)),"width");
        h1_unf->Scale(1./jet_unf->Integral(jet_unf->GetXaxis()->FindBin(pt1),jet_unf->GetXaxis()->FindBin(pt2)),"width");
    }
    
    h1_num->Divide(h1_den);
    h1_dat->Divide(h1_num);


    TCanvas *c2 = new TCanvas("c2", "Overlay Histograms", 600, 800);
    jet_dat->Draw();
    jet_unf->Draw("same");
    h1_eff->Draw();


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
    h1_dat->SetTitle("");
    h1_dat->GetYaxis()->SetRangeUser(0.5,6);
    h1_dat->GetYaxis()->SetTitle("Normalized EEC");
    h1_dat->Draw();
    h1_unf->Draw("same");
    h1_dat->GetXaxis()->SetRangeUser(0.01,0.4);
    
    bottomPad->cd();
    gPad->SetLogx();
    TH1* h1_datClone = (TH1*)h1_dat->Clone("h1_datClone");
    TH1* h1_unfClone = (TH1*)h1_unf->Clone("h1_unfClone");
    
    h1_unfClone->SetStats(0);
    h1_unfClone->SetTitle("");
    
    h1_unfClone->GetYaxis()->SetRangeUser(0.8,1.22);
    h1_unfClone->GetYaxis()->SetLabelSize(0.06);
    h1_unfClone->GetYaxis()->SetTitleSize(0.08);
    h1_unfClone->GetYaxis()->SetTitleOffset(0.53);
    
    h1_unfClone->GetXaxis()->SetRangeUser(0.01,0.4);
    h1_unfClone->GetXaxis()->SetLabelSize(0.09);
    h1_unfClone->GetXaxis()->SetTitleSize(0.09);
    h1_unfClone->GetXaxis()->SetTitleOffset(1.07);
    
    h1_unfClone->GetYaxis()->SetTitle("Unf/BbB");
    
    h1_unfClone->Divide(h1_datClone);
    h1_unfClone->Draw();
    

        
    }

