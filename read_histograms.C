#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"

void read_histograms() {
    // Open the ROOT file
    TFile *file = TFile::Open("Unfold3Dtest_4.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file." << std::endl;
        return;
    }
//histogram names: h1_true, h1_reco, eec_pt_hist, eec_pt_hist_det


    // Read the 1D histogram
    TH1D *hist1D = nullptr;
    file->GetObject("h1_true", hist1D);
    if (!hist1D) {
        std::cerr << "Error reading 1D histogram." << std::endl;
        file->Close();
        return;
    }

    // Read the 2D histogram
    TH2D *hist2D = nullptr;
    file->GetObject("eec_pt_hist", hist2D);
    if (!hist2D) {
        std::cerr << "Error reading 2D histogram." << std::endl;
        file->Close();
        return;
    }
    
    int pt1 = 10;
    int pt2 = 13;
    
    TH1D* h1 = hist2D->ProjectionX("EEC",pt1,pt2,"e");
    h1->Scale(1./hist1D->Integral(pt1,pt2),"width");
    
    
    // Create a canvas to draw the histograms
    TCanvas *c1 = new TCanvas("c1", "Histograms", 800, 600);
//    c1->Divide(3, 1);

    // Draw the 1D histogram
//    c1->cd(1);
//    gStyle->SetOptStat(0);
//    hist1D->Draw();
//
//    // Draw the 2D histogram
//    c1->cd(2);
//    gStyle->SetOptStat(0);
//    hist2D->Draw("COLZ");
    
    c1->cd(3);
    gStyle->SetOptStat(0);
    gPad->SetLogx();
    h1->Draw();


    // Save the canvas as an image
    c1->SaveAs("histograms.pdf");

    // Clean up
//    file->Close();
//    delete c1;
}

