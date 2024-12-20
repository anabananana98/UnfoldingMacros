//This script plots covariance matrix
//Normalizes the matrix by bin width
//Then computes the correlation matrix
//ALL FOR 1D inputs 

#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMatrixD.h>
#include <TCanvas.h>
#include <iostream>

void PlotCovarianceMatrix() {
    // Open the ROOT file containing the unfolded histogram and covariance matrix
    TFile* file = TFile::Open("UnfoldingResultsDataEEC_TestCovMatrixToy.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return;
    }

    // Load the unfolded histogram
    TH1D* hUnfolded = (TH1D*)file->Get("h1_unfolded");
    if (!hUnfolded) {
        std::cerr << "Error: Unfolded histogram not found!" << std::endl;
        return;
    }

    // Load the covariance matrix
    TMatrixD* covarianceMatrixPtr = (TMatrixD*)file->Get("covMatrix1D");
    if (!covarianceMatrixPtr) {
        std::cerr << "Error: Covariance matrix not found!" << std::endl;
        return;
    }
    TMatrixD covarianceMatrix = *covarianceMatrixPtr;

    cout<<std::pow(hUnfolded->GetBinError(1), 2)<<endl;
    hUnfolded->Scale(1.,"width");
    cout<<std::pow(hUnfolded->GetBinError(1), 2)<<endl;
    
    
    //scale covariance matrix by bin width
    int nBins = covarianceMatrix.GetNrows();
    for (int i = 0; i < nBins; ++i) {
        double binWidthI = hUnfolded->GetBinWidth(i + 1); // +1 because ROOT bins start at 1
        for (int j = 0; j < nBins; ++j) {
            double binWidthJ = hUnfolded->GetBinWidth(j + 1);
            covarianceMatrix(i, j) /= (binWidthI * binWidthJ);
        }
    }
   
    cout<<covarianceMatrix(0, 0)<<endl;
    
//    // Step 2: Normalize the covariance matrix

 TMatrixD corrMatrix(nBins, nBins);
    for (int i = 0; i < nBins; ++i) {
        for (int j = 0; j < nBins; ++j) {
            double variance_i = covarianceMatrix(i, i);  // Variance of bin i (diagonal element)
            double variance_j = covarianceMatrix(j, j);  // Variance of bin j (diagonal element)
            if (variance_i > 0 && variance_j > 0) {
                corrMatrix(i, j) = covarianceMatrix(i, j) / (sqrt(variance_i) * sqrt(variance_j));
            } else {
                corrMatrix(i, j) = 0;  // Set correlations involving zero variance to 0
            }
        }
    }

    // Step 3: Verify the normalized covariance matrix
    std::cout << "Normalized Covariance Matrix:" << std::endl;
    covarianceMatrix.Print();

    // Verify diagonal elements match squared errors
    for (int i = 0; i < nBins; ++i) {
        double expectedVariance = std::pow(hUnfolded->GetBinError(i + 1), 2);
        double actualVariance = covarianceMatrix(i, i);
        std::cout << "Bin " << i + 1 << ": Expected Variance = " << expectedVariance
                  << ", Actual Variance = " << actualVariance << std::endl;
    }

    // Step 4: Plot the normalized covariance matrix
    TH2D* hCovNorm = new TH2D("hCovNorm", "Normalized Covariance Matrix;Bin i;Bin j",
                               nBins, 0, nBins , nBins, 0, nBins );

    for (int i = 0; i < nBins; ++i) {
        for (int j = 0; j < nBins; ++j) {
            hCovNorm->SetBinContent(i + 1, j + 1, corrMatrix(i, j));
        }
    }

    TCanvas* c1 = new TCanvas("c1", "Normalized Covariance Matrix", 800, 600);
    gPad->SetLogz();
    hCovNorm->SetStats(0); // Disable statistics box
    hCovNorm->Draw("COLZ");
    c1->SaveAs("NormalizedCovarianceMatrix.png");

    // Step 5: Save the normalized histogram and covariance matrix
    TFile* outputFile = TFile::Open("normalized_results.root", "RECREATE");
    hUnfolded->Write("normalizedHistogram");
    hCovNorm->Write("normalizedCovarianceMatrix");
    outputFile->Close();

    std::cout << "Normalized histogram and covariance matrix saved to normalized_results.root." << std::endl;

    // Clean up
//    delete c1;
//    delete hCovNorm;
//    file->Close();
//    delete file;
}

