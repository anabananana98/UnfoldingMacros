//This script accesses the covariance matrix from RooUnfold
//The covariance matrix is for a 3D hitsogram
//Hence bins are flattened and a mapping is required
//The option to create a correlation matrix is also written out
//Should not scale by widths until after projections are made 

#include <TROOT.h>
#include <TFile.h>
#include <TH3D.h>
#include <TMatrixD.h>
#include <TCanvas.h>
#include <iostream>

void PlotCovarianceMatrix3D() {
    // Open the ROOT file containing the unfolded histogram and covariance matrix
    TFile* file = TFile::Open("UnfoldingResultsDataEEC_TestCovMatrix3DToy.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return;
    }

    // Load the unfolded 3D histogram
    TH3D* hUnfolded = (TH3D*)file->Get("h3_unfoldedData");
    if (!hUnfolded) {
        std::cerr << "Error: Unfolded histogram not found!" << std::endl;
        return;
    }

    // Load the covariance matrix
    TMatrixD* covarianceMatrixPtr = (TMatrixD*)file->Get("covMatrix3D");
    if (!covarianceMatrixPtr) {
        std::cerr << "Error: Covariance matrix not found!" << std::endl;
        return;
    }
    TMatrixD covarianceMatrix = *covarianceMatrixPtr;
    cout<<covarianceMatrix.GetNrows()<<endl;
    // Scaling the covariance matrix by bin width
    int nBinsX = hUnfolded->GetNbinsX();
    int nBinsY = hUnfolded->GetNbinsY();
    int nBinsZ = hUnfolded->GetNbinsZ();
    int nBins = nBinsX * nBinsY * nBinsZ;  // Total number of bins in 3D

cout<<nBinsX<<" "<<nBinsY<<" "<<nBinsZ<<endl;
    
    cout<<hUnfolded->GetBinContent(10,5,3)<<endl;
    cout<<std::pow(hUnfolded->GetBinError(10,5,3), 2)<<endl;
    cout<<std::pow(covarianceMatrix(405,405), 1)<<endl;
//    cout<<std::pow(covarianceMatrix(102,102), 2)<<endl;
//    hUnfolded->Scale(1.,"width");

//    // Rescale covariance matrix based on bin widths -- you dont scale until after projecting
// //Map from 3D bins to flattened bins: (i+1,j+1,k+1) = (i + j * nBinsX + k * nBinsX * nBinsY)
////////////////////////////////////////////////////////////////////////////////////////////////
//    for (int i = 0; i < nBinsX; ++i) {
//        for (int j = 0; j < nBinsY; ++j) {
//            for (int k = 0; k < nBinsZ; ++k) {
//                double binWidthI = hUnfolded->GetXaxis()->GetBinWidth(i + 1);
//                double binWidthJ = hUnfolded->GetYaxis()->GetBinWidth(j + 1);
//                double binWidthK = hUnfolded->GetZaxis()->GetBinWidth(k + 1);
//                for (int l = 0; l < nBinsX; ++l) {
//                    for (int m = 0; m < nBinsY; ++m) {
//                        for (int n = 0; n < nBinsZ; ++n) {
//                            covarianceMatrix(i + j * nBinsX + k * nBinsX * nBinsY, 
//                                              l + m * nBinsX + n * nBinsX * nBinsY) /= 
//                                              (binWidthI * binWidthJ * binWidthK);
//                        }
//                    }
//                }
//            }
//        }
//    }
////////////////////////////////////////////////////////////////////////////////////
cout<<std::pow(covarianceMatrix(405,405), 1)<<endl;


//    // Normalize covariance matrix (if applicable)
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

    // Step 2: Plot the normalized covariance matrix
    TH2D* hCovNorm = new TH2D("hCovNorm", "Normalized Covariance Matrix;Bin i;Bin j",
                               nBins, 0, nBins , nBins, 0, nBins);

//    for (int i = 0; i < nBins; ++i) {
//        for (int j = 0; j < nBins; ++j) {
//            hCovNorm->SetBinContent(i + 1, j + 1, corrMatrix(i, j));
//        }
//    }
   
    TCanvas* c1 = new TCanvas("c1", "Normalized Covariance Matrix", 800, 600);
    gPad->SetLogz();
    hCovNorm->SetStats(0); // Disable statistics box
    hCovNorm->Draw("COLZ");
    c1->SaveAs("NormalizedCovarianceMatrix3D.png");

    // Step 3: Save the results
    TFile* outputFile = TFile::Open("normalized_results_3D.root", "RECREATE");
    hUnfolded->Write("normalizedHistogram");
    hCovNorm->Write("normalizedCovarianceMatrix");
    outputFile->Close();

    std::cout << "Normalized histogram and covariance matrix saved to normalized_results_3D.root." << std::endl;

    // Clean up
//    delete c1;
//    delete hCovNorm;
//    file->Close();
//    delete file;
}

