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


void CompareNIter3D(std::string tag, std::string projection) {
    
    int num=4;
    int pt1 =20;
    int pt2 = 40;
    double rlow = 0.01;
    double rhigh = 0.4;
    const char *str = "#it{p}_{T,jet}";
    const char *str1 = "#it{p}_{T,min}";
    const char *str2 = "anti-#it{k}_{T}";
    const char *str3 = "#sqrt{s} = 13 TeV";
    const char *str4 = "#it{R}_{L}";
    const char *str8 = "#it{#eta}_{jet}";
    const char *str11 = "#it{p}_{T}^{const}";
    
    // Create a legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.88, 0.88);
    //    legend->SetHeader("Histograms", "C");
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.045);
    // Array of file names
    const char* filenames[4] = {
        "UnfoldingResultsSplitTest_1.root",
        "UnfoldingResultsSplitTest_4.root",
        "UnfoldingResultsSplitTest_7.root",
        "UnfoldingResultsSplitTest_10.root",
    };
    
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
        originalHistNameTru = "h3_true_eec";
    }
    else{
        originalHistNameUnf = "h3_unfolded";
        originalHistNameTru = "h3_true_eec";
    
    }
    
    // Array of colors for histograms
    Color_t colors[4] = {kRed, kBlue, kMagenta, kCyan};
    
    // Array of markers for histograms
    Color_t markers[4] = {20, 21, 29, 33};
    
    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Overlay Histograms", 800, 600);
//    c1->SetGrid();
    
    
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
        
        hist3DUnf->GetXaxis()->SetRangeUser(0.005,0.5);

        
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
        hist1DUnf->GetYaxis()->SetRangeUser(0.5,1.5);
        hist1DTru->GetYaxis()->SetRangeUser(0.5,1.5);
//        TH1 *hist1DUnf = hist3DUnf->Project3D(histNamesUnf[i],"x");
//        TH1 *hist1DTru = hist3DTru->Project3D(histNamesTru[i],"x");
       
        // Set line color and thickness
        hist1DUnf->SetLineColor(colors[i]); // Set line color
        hist1DUnf->SetLineWidth(2);         // Set line thickness to 2
        
        hist1DUnf->Divide(hist1DTru);
        // Draw the histogram
        if (i == 1) {
            hist1DUnf->GetYaxis()->SetTitle("Unfolded/True");
            hist1DUnf->Draw();
            hist1DUnf->SetStats(0);
        } else {
            hist1DUnf->GetYaxis()->SetTitle("Unfolded/True");
            hist1DUnf->Draw("SAME");
        }
        
        // Add histogram to legend
        legend->AddEntry(hist1DUnf, histNamesUnf[i], "l");
        
        // Close the file
        //        file->Close();
    }
    
    // Draw the legend
    legend->Draw();
    
    // Update the canvas
//    c1->Update();
}
