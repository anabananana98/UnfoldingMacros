//This script correctly propagates errors from 3D unfolded histogram to
//1D projection
//kCovToy is used to construct the covariance matrix from unfolding
//Bins are mapped and the covariance errors are summed for projections


#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMatrixD.h>
#include <iostream>

void PropagateErrorsTo1D() {
    // Open the ROOT file containing the 3D histogram and covariance matrix
    TFile* file = TFile::Open("UnfoldingResultsDataEEC_TestCovMatrix3DToy.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return;
    }
    
    // Load the 3D histogram
    TH3D* h3D = (TH3D*)file->Get("h3_unfoldedData");
    if (!h3D) {
        std::cerr << "Error: 3D histogram not found!" << std::endl;
        return;
    }
    h3D->Sumw2();
    
    // Load the covariance matrix
    TMatrixD* covarianceMatrix = (TMatrixD*)file->Get("covMatrix3D");
    if (!covarianceMatrix) {
        std::cerr << "Error: Covariance matrix not found!" << std::endl;
        return;
    }
    
    // Define the range in Y
    int yMin = 2; // Start bin in Y (inclusive)
    int yMax = 2; // End bin in Y (inclusive)
    
    // Get bin information
    int nBinsX = h3D->GetNbinsX();
    int nBinsY = h3D->GetNbinsY();
    int nBinsZ = h3D->GetNbinsZ();
    
    // Create the 2D histogram for the Z-X projection
    TH2D* h2D = new TH2D("h2D", "Projected Histogram;X;Z", nBinsX, h3D->GetXaxis()->GetXbins()->GetArray(),
                         nBinsZ, h3D->GetZaxis()->GetXbins()->GetArray());
    
    // Loop over the bins in Z and X
    for (int x = 1; x <= nBinsX; ++x) {
        for (int z = 1; z <= nBinsZ; ++z) {
            double binContent = 0.0;
            double errorSquared = 0.0;
            
            // Sum the content and propagate the error across the Y range
            for (int y1 = yMin; y1 <= yMax; ++y1) {
                for (int y2 = yMin; y2 <= yMax; ++y2) {
                    int bin1 = (x - 1) + (y1 - 1) * nBinsX + (z - 1) * nBinsX * nBinsY;
                    int bin2 = (x - 1) + (y2 - 1) * nBinsX + (z - 1) * nBinsX * nBinsY;
                    errorSquared += (*covarianceMatrix)(bin1, bin2);
                }
                int bin3D = h3D->GetBin(x, y1, z);
                binContent += h3D->GetBinContent(bin3D);
            }
            
            // Set the content and error for the 2D bin
            h2D->SetBinContent(x, z, binContent);
            h2D->SetBinError(x, z, std::sqrt(errorSquared));
        }
    }
    TCanvas* c1 = new TCanvas("c1", "Normalized Covariance Matrix", 800, 600);
    h2D->Draw("COLZ");
    
    TH2D* h2 = (TH2D*)h2D->Clone("h2");
    TH1D *h1 = (TH1D*)h2->ProjectionX("h1");
    h1->Reset();
    
    for (int i = 1; i <= nBinsX; ++i) {
        int nBinsY2D = h2->GetNbinsY();
        double weightedSum = 0;
        double sumOfWeights = 0;
        double errorSumSquared = 0;
        
        for (int j = 1; j <= nBinsY2D; ++j) {
            double binContent = h2->GetBinContent(i, j);
            double binError = h2->GetBinError(i, j); // Error for (i, j)
            double binCenterY = h2->GetYaxis()->GetBinCenter(j); // Bin center in Y direction
            cout<<"here"<<endl;
            // Weight based on Y-bin center
            double weight = binCenterY;
            
            // Accumulate weighted sum and sum of weights
            weightedSum += binContent * weight;
            
            // Propagate error: Add (weight^2 * error^2) to the sum
            errorSumSquared += std::pow(weight * binError, 2);
        }
        
      
        // Fill the 1D histogram
        h1->SetBinContent(i, weightedSum); // Weighted sum
       
        h1->SetBinError(i, std::sqrt(errorSumSquared)); // Propagated error
        
    }
    
    TCanvas* c2 = new TCanvas("c2", "Normalized Covariance Matrix", 800, 600);
    h1->Scale(1.,"width");
    h1->Draw();
    
    
    //######### Now compare this with how you do this generally #######
    TH3D* h3_old = (TH3D*)h3D->Clone("h3_old");
    h3_old->Sumw2();
    h3_old->GetYaxis()->SetRange(yMin,yMax);
    TH2D *h2_old = (TH2D*)h3_old->Project3D("zx"); //Rl vs wt histogram
    TH1D *h1_old = (TH1D*)h2_old->ProjectionX("h1_old");
    h1_old->Reset();
    
    for (int i = 1; i <= nBinsX; ++i) {
        int nBinsY = h2_old->GetNbinsY();
        double weightedSum = 0;
        double sumOfWeights = 0;
        
        for (int j = 1; j <= nBinsY; ++j) {
            double binContent = h2_old->GetBinContent(i, j);
            double binCenterY = h2_old->GetYaxis()->GetBinCenter(j); // Bin center in Y direction
            
            // Compute weight as the bin center in Y direction
            double weight = binCenterY;
            
            // Accumulate weighted sum and sum of weights
            weightedSum += binContent * weight;
            
        }
        
        h1_old->SetBinContent(i, weightedSum);//this is closest to filling on the fly,not the avg
    }
    h1_old->Scale(1.,"width");
    h1_old->SetLineColor(kRed);
    h1_old->Draw("same");
    
    
    //Errors are not propagated in this case
    for (int i = 1; i <= nBinsX; ++i) {
    cout<<h1_old->GetBinError(i)<<setprecision(10)<<" and "<<h1->GetBinError(i)<<endl;
}


    
//    int nBins = h1->GetNbinsX();
//    TGraphErrors* graph = new TGraphErrors(nBins);
//    TGraphErrors* graph1 = new TGraphErrors(nBins);
//    
//    for (int i = 1; i <= nBins; ++i) {
//        double x = h1->GetBinCenter(i);  // X value (bin center)
//        double y = h1->GetBinContent(i); // Y value (bin content)
//        double ex = h1->GetBinWidth(i) / 2.0;  // X error (half of bin width)
//        double ey = h1->GetBinError(i);        // Y error (bin error)
//        graph->SetPoint(i-1, x, 0);          // Set X, Y values
//        graph->SetPointError(i-1, ex, ey);   // Set error values
//
//        
//        double x1 = h1_old->GetBinCenter(i);  // X value (bin center)
//        double y1 = h1_old->GetBinContent(i); // Y value (bin content)
//        double ex1 = h1_old->GetBinWidth(i) / 2.0;  // X error (half of bin width)
//        double ey1 = h1_old->GetBinError(i);        // Y error (bin error)
//        graph1->SetPoint(i-1, x1, 0);          // Set X, Y values
//        graph1->SetPointError(i-1, ex1, ey1);   // Set error values
//    }
//    
//    graph->SetTitle("1D Histogram with Errors");
//    graph->SetLineColor(kRed);
//     graph->SetMarkerStyle(21);
//    graph->Draw("AP");
//
//    graph1->Draw("SAME");
//    
    
    
    
    
   
}

