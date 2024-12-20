//This script compares error methods from RooUnfold
//Various comparisons between kCov and kCovToy, kErrors and kNoErrors can be made 

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraphErrors.h>

void plot_histogram_errors() {
    // File paths for two sets of four files
//    const char* set1_files[4] = {"UnfoldingResultsDataEEC_1ErrorTestJetPt.root", "UnfoldingResultsDataEEC_4ErrorTestJetPt.root", "UnfoldingResultsDataEEC_7ErrorTestJetPt.root", "UnfoldingResultsDataEEC_10ErrorTestJetPt.root"};
//    const char* set2_files[4] = {"UnfoldingResultsDataEEC_1ErrorTestMCToyJetPt.root", "UnfoldingResultsDataEEC_4ErrorTestMCToyJetPt.root", "UnfoldingResultsDataEEC_7ErrorTestMCToyJetPt.root", "UnfoldingResultsDataEEC_10ErrorTesMCToytJetPt.root"};
    
     const char* set1_files[2] = {"UnfoldingResultsDataEEC_1kNoError.root", "UnfoldingResultsDataEEC_10kCov.root"};
    const char* set2_files[2] = {"UnfoldingResultsDataEEC_1kErrors.root", "UnfoldingResultsDataEEC_10kErrors.root"};

    // Histogram name to be retrieved
    const char* hist_name = "h1_unfolded";

    // Colors for histograms
    int colors[4] = {kRed, kBlue, kGreen+2, kMagenta};

    // Arrays to store histograms and graphs
    TH1D* hist_set1[4];
    TH1D* hist_set2[4];
    TGraphErrors* errors_set1[4];
    TGraphErrors* errors_set2[4];
    TGraphErrors* diff[4];

    // Canvas and legend
    TCanvas* c1 = new TCanvas("c1", "Histogram Errors", 800, 600);
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    gStyle->SetOptStat(0); // Disable stats box
    
    // Canvas and legend
    TCanvas* c2 = new TCanvas("c2", "Histogram Errors Diff", 800, 600);
    TLegend* legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    gStyle->SetOptStat(0); // Disable stats box

    // Retrieve histograms and plot errors
    for (int i = 0; i < 1; i++) {
        // Open files
        TFile* file1 = TFile::Open(set1_files[i]);
        TFile* file2 = TFile::Open(set2_files[i]);

        if (!file1 || !file2) {
            printf("Error: Could not open file %d\n", i);
            return;
        }

        // Retrieve histograms
        hist_set1[i] = (TH1D*)file1->Get(hist_name);
        hist_set1[i]->SetName(Form("set1_%i",i));
        hist_set1[i]->Scale(1.,"width");
        hist_set2[i] = (TH1D*)file2->Get(hist_name);
        hist_set2[i]->SetName(Form("set2_%i",i));
        hist_set2[i]->Scale(1.,"width");

        if (!hist_set1[i] || !hist_set2[i]) {
            printf("Error: Histogram %s not found in file %d\n", hist_name, i);
            return;
        }

        // Create TGraphErrors for Set 1 (errors only with markers)
        int nBins = hist_set1[i]->GetNbinsX();
        errors_set1[i] = new TGraphErrors(nBins);
        errors_set1[i]->SetMarkerStyle(24); // Open circle markers
        errors_set1[i]->SetLineColor(colors[i]);
        errors_set1[i]->SetMarkerColor(colors[i]);
        errors_set1[i]->SetLineWidth(2);
        
       
        // Fill the graph with bin errors
        for (int bin = 1; bin <= nBins; bin++) {
            double x = hist_set1[i]->GetBinCenter(bin);
            double y = 0; // Content is zero for errors-only plot
            double ex = 0;
//            double ey = hist_set1[i]->GetBinError(bin);
            double ey = hist_set1[i]->GetStdDev(bin);
            cout<<ey<<endl;
//            cout<<hist_set1[i]->GetBinError(bin)<<endl;
//            double ey = hist_set1[i]->GetStdErrorY(bin);
            errors_set1[i]->SetPoint(bin - 1, x,ey);
            errors_set1[i]->SetPointError(bin - 1, ex, 0);


        }

        // Create TGraphErrors for Set 2 (errors as lines)
        errors_set2[i] = new TGraphErrors(nBins);
        errors_set2[i]->SetLineColor(colors[i]);
        errors_set2[i]->SetMarkerStyle(21);
        errors_set2[i]->SetMarkerColor(colors[i]);
        errors_set2[i]->SetLineWidth(2);

        for (int bin = 1; bin <= nBins; bin++) {
            double x = hist_set2[i]->GetBinCenter(bin);
            double y = 0; // Content is zero for errors-only plot
            double ex = 0;
            double ey = hist_set2[i]->GetStdDev(bin);
//            double ey = hist_set2[i]->GetBinError(bin);
            double ey_up = hist_set2[i]->GetBinErrorUp(bin);
            double ey_low = hist_set2[i]->GetBinErrorLow(bin);
             cout<<"toy "<<ey<<endl;
//             cout<<"toy up "<<ey_up<<endl;
//             cout<<"toy low "<<ey_low<<endl;
//            double ey = hist_set2[i]->GetErrorY(bin);
            errors_set2[i]->SetPoint(bin - 1, x, ey);
            errors_set2[i]->SetPointError(bin - 1, ex, 0);
        }
         
         c1->cd();
        // Draw the errors
        if (i == 0) { // Draw the first histogram
            errors_set1[i]->Draw("AP"); // Set 1: Points with errors
            errors_set2[i]->Draw("L SAME"); // Set 2: Lines with errors
        } else {
            errors_set1[i]->Draw("P SAME");
            errors_set2[i]->Draw("L SAME");
        }
        
        int iter=-1;
        if(i+1==1){
        iter=1;
        }
        else if(i+1==2){iter=10;}
        else if(i+1==3){iter=7;}
        else if(i+1==4){iter=10;}

        // Add to legend
//        legend->AddEntry(errors_set1[i], Form("kCov iter%d", iter), "P");
        legend->AddEntry(errors_set1[i], Form("kNoError iter%d", iter), "P");
//        legend->AddEntry(errors_set2[i], Form("kCovToy iter%d", iter), "P");
        legend->AddEntry(errors_set2[i], Form("kErrors iter%d", iter), "L");
        
        
        c2->cd();
        diff[i] = new TGraphErrors(nBins);
        diff[i]->SetLineColor(colors[i]);
        diff[i]->SetMarkerStyle(24);
        diff[i]->SetMarkerColor(colors[i]);
        for (int bin = 1; bin <= nBins; bin++) {
            double x = hist_set2[i]->GetBinCenter(bin);
            double y = 0; // Content is zero for errors-only plot
            double ex = 0;
            double ey = hist_set2[i]->GetBinError(bin);
            double ey1 = hist_set1[i]->GetBinError(bin);
            double diffy = (ey - ey1)/ey;
            diff[i]->SetPoint(bin - 1, x, diffy);
            diff[i]->SetPointError(bin - 1, ex, 0);
        }
       if (i == 0) { // Draw the first histogram
            diff[i]->Draw("AP"); // Set 1: Points with errors
          
        } else {
           diff[i]->Draw("P SAME");
        }
         // Add to legend
        legend2->AddEntry(diff[i], Form("(kCovToy - kCov)/kCovToy iter%d", iter), "P");
        
    }

    // Finalize the canvas
    c1->cd();
    legend->Draw();
    c1->Update();
    c2->cd();
    legend2->Draw();
    c2->Update();
    c1->SaveAs("histogram_errors.png"); // Save the output
}

