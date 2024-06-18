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

void CompareUnfolding(std::string inputFileName, std::string inputFileName2) {
    
    TFile* f = TFile::Open(inputFileName.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening input file " << inputFileName << std::endl;
        return;
    }
    
    
    TFile* f2 = TFile::Open(inputFileName2.c_str());
    if (!f2 || f2->IsZombie()) {
        std::cerr << "Error opening input file " << inputFileName2 << std::endl;
        return;
    }
    
    TH1D* h1_raw = (TH1D*)f->Get("h1_raw");
    
    TH1D* h1_fold = (TH1D*)f->Get("h1_fold");
    TH1D* h1_fold2 = (TH1D*)f2->Get("h1_fold");
    h1_fold2->SetName("h1_fold2");
    
    TH1D* h1_unfold = (TH1D*)f->Get("h1_unfolded");
    TH1D* h1_unfold2 = (TH1D*)f2->Get("h1_unfolded");
    h1_fold2->SetName("h1_unfolded2");
    
    h1_raw->SetLineColor(kRed);
    h1_fold->SetLineColor(kBlue);
    h1_fold2->SetLineColor(kGreen+4);
    h1_unfold->SetLineColor(kMagenta);
    h1_unfold2->SetLineColor(kCyan);
    
    
    TCanvas* c = new TCanvas();
    c->cd();
    h1_fold->SetStats(0);
    h1_fold->GetXaxis()->SetTitle("jet_pt");
    h1_fold->GetYaxis()->SetTitle("fold/raw");
    h1_fold->Divide(h1_raw);
    h1_fold2->Divide(h1_raw);
    h1_fold->Draw();
    h1_fold2->Draw("same");
    
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    //    legend->AddEntry(h1_raw, "Raw", "l");
    legend->AddEntry(h1_fold, "iter = 4", "l");
    legend->AddEntry(h1_fold2, "iter = 10", "l");
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);
    legend->Draw();
    
    c->Update();
    
    //
//    TCanvas* c1 = new TCanvas();
//    h1_unfold2->SetStats(0);
//    h1_unfold2->Divide(h1_unfold);
//    h1_unfold2->GetXaxis()->SetTitle("jet_pt");
//    h1_unfold2->GetYaxis()->SetTitle("iter=10/iter=4");
    
    //    TLegend* legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    //    legend1->AddEntry(h1_unfold, "iter = 4", "l");
    //    legend1->AddEntry(h1_unfold2, "iter = 10", "l");
    //    legend1->Draw();
    
    
    
    
    //    // Close the files
    //    f->Close();
    //    f2->Close();
    
}

