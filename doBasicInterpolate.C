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

//bins means is the plot infinite stats pythia or limited stats data bins

void doBasicInterpolate(int pt1, int pt2, const char* inputFile, const char* histName, const char* outputFile, std::string tag="", std::string type = "EEC",std::string bins = "data") {
    
    int pt1jet = pt1;
    int pt2jet = pt2;
    
    const char *str = "#it{p}_{T,jet}";
    const char *str1 = "#it{p}_{T,min}";
    const char *str2 = "anti-#it{k}_{T}";
    const char *str3 = "#sqrt{s} = 13 TeV";
    const char *str4 = "#it{R}_{L}";
    const char *str8 = "#it{#eta}_{jet}";
    const char *str11 = "#it{p}_{T}^{const}";
    
    // Open the input file
    TFile *file = TFile::Open(inputFile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << inputFile << std::endl;
        return;
    }
    
    TH1D* hjet = (TH1D*)file->Get("jet_pt_hist");
    hjet->Scale(1.,"width");
    
    TCanvas* cjet;
    hjet->Draw();
    
    TH3D* h3 = (TH3D*)file->Get(histName);
    h3->GetXaxis()->SetTitle(str4);
    h3->GetYaxis()->SetTitle(str);
    h3->GetZaxis()->SetTitle("wt");
    
    
    if (!h3) {
        std::cerr << "Error retrieving histogram: " << histName << std::endl;
        file->Close();
        return;
    }
    h3->GetYaxis()->SetRange(pt1,pt2);
    // Get the 2D histogram
    TH2D *h2 = (TH2D*)h3->Project3D("zx"); //Rl vs wt histogram
    
    
    // Create the 1D histogram
    int nBinsX = h2->GetNbinsX();
    int nBinsY = h2->GetNbinsY();
    double xMin = h2->GetXaxis()->GetXmin();
    double xMax = h2->GetXaxis()->GetXmax();
  
    
    TH1D *h1 = (TH1D*)h2->ProjectionX("h1");
    h1->Reset();
    
    TH1D *h1_old = (TH1D*)h2->ProjectionX("h1_old");
    h1->Reset();
    
    TH1D *h1_otherWay = (TH1D*)h2->ProjectionX("h1_otherWay");
    h1_otherWay->Reset();
    
//    int nBins_new = 30; // Finer number of bins
//    TH1D *h1_interpolated = new TH1D("h1_scaled", "Scaled Histogram",wt_bins,wt_new_bins);
//    TH1D *h1_try = (TH1D*)h2->ProjectionY(Form("h1_try_%i",33),1,2);
//    
//    
//    // Create a TGraph to store the original bin centers and contents
//    int nBins_original = h1_try->GetNbinsX();
//    double* x = new double[nBins_original];
//    double* y = new double[nBins_original];
//    
//    // Fill x and y arrays with the bin centers and contents from h1_try
//    for (int i = 1; i <= nBins_original; ++i) {
//        x[i-1] = h1_try->GetXaxis()->GetBinCenter(i);
//        y[i-1] = h1_try->GetBinContent(i);
//    }
    
//    // Create a TSpline3 object for cubic spline interpolation
//    TGraph* graph = new TGraph(nBins_original, x, y);
//    TSpline3* spline = new TSpline3("spline", graph);
//    spline->SetLineColor(kRed);
//    
//    
//    for (int i = 1; i <= wt_bins; ++i) {
//        double x = h1_interpolated->GetXaxis()->GetBinCenter(i);
//        //        double interpolatedValue = h1_try->Interpolate(x);
//        
//        double interpolatedValue = spline->Eval(x);
//        
//        cout<<x<<" and interpolated "<<interpolatedValue<<endl;
//        h1_interpolated->SetBinContent(i, interpolatedValue);
//    }
//    
    
    
    
    TCanvas* c;
    h3->Draw();
    h1->Draw();
//    h1_try->Draw();
//    h1_interpolated->SetLineColor(kRed);
//    h1_interpolated->Draw("same");
//    spline->Draw("same");
    
    
    TCanvas* ctry = new TCanvas;
    
    //Loop over and interpolate
//    cout<<nBinsX<<endl;
//    for (int i = 1; i <= nBinsX; ++i){
//        //    for (int i = 1; i <= 1; ++i){
//        int nBinsY = h2->GetNbinsY();
//        double weightedSum = 0;
//        TH1D *h1_wt = (TH1D*)h2->ProjectionY(Form("h1_wt_%i",i),i,i+1);
//        TH1D *h1_interpolatedn = new TH1D(Form("h1_interpolatedn_%i",i),"Interp",wt_bins,wt_new_bins);
//        h1_wt->Draw();
//        
//        
//        
//        // Create TGraph and TSpline3 for non-linear interpolation
//        int nBins_original = h1_wt->GetNbinsX();
//        double* x = new double[nBins_original];
//        double* y = new double[nBins_original];
//        
//        // Fill x and y arrays with bin centers and contents from h1_wt
//        for (int j = 1; j <= nBins_original; ++j) {
//            x[j - 1] = h1_wt->GetXaxis()->GetBinCenter(j);
//            y[j - 1] = h1_wt->GetBinContent(j);
//        }
//        
//        // Create TGraph and TSpline3 objects
//        TGraph* graph = new TGraph(nBins_original, x, y);
//        TSpline3* spline = new TSpline3("spline", graph);
//        
//        
        
//        for (int j = 1; j <= h1_interpolatedn->GetNbinsX(); ++j){
//                    for (int j = 1; j <= h1_wt->GetNbinsX(); ++j) {
            
//            double x = h1_interpolatedn->GetXaxis()->GetBinCenter(j);

//                        double interpolatedValue = h1_wt->Interpolate(x);
            
            
//            double interpolatedValue = spline->Eval(x);
            
            //The interpolated value is giving an estimate of the y counts based on bin centers
//                        cout<<x<<" and interpolated "<<interpolatedValue<<endl;
            
//            weightedSum = interpolatedValue*x;
//            h1_interpolatedn->SetBinContent(j, interpolatedValue);
//        }
        
//        int nBins = h1_interpolatedn->GetNbinsX();
//        for (int k = 1; k <= nBins; ++k) {
//            // Get the center of the current bin
//            double binCenter = h1_interpolatedn->GetXaxis()->GetBinCenter(k);
//            // Get the content of the current bin
//            double binContent = h1_interpolatedn->GetBinContent(k);
//            // Accumulate the product of binCenter and binContent
//            weightedSum += binCenter * binContent;
//        }
        
//        h1->SetBinContent(i, weightedSum);
//    }
    
//    for (int i = 1; i <= nBinsX; ++i){
//        //    for (int i = 1; i <= 1; ++i){
//        int nBinsY = h2->GetNbinsY();
//        double weightedSum = 0;
//        TH1D *h1_wt = (TH1D*)h2->ProjectionY(Form("h1_wt_%i",i),i,i+1);
//        TH1D *h1_interpolatedn = new TH1D(Form("h1_interpolatedn_%i",i),"Interp",wt_bins,wt_new_bins);
//        h1_wt->Draw();
//       
//        // Create TGraph and TSpline3 for non-linear interpolation
//        int nBins_original = h1_wt->GetNbinsX();
//        double* xOriginal = new double[nBins_original];
//        double* yOriginal = new double[nBins_original];
//        
//        // Fill xOriginal and yOriginal arrays with bin centers and contents from h1_wt
//        for (int j = 1; j <= nBins_original; ++j) {
//            xOriginal[j - 1] = h1_wt->GetXaxis()->GetBinCenter(j);
//            yOriginal[j - 1] = h1_wt->GetBinContent(j);
//        }
//        
//        TGraph* graph = new TGraph(nBins_original, xOriginal, yOriginal);
//        TSpline3* spline = new TSpline3(Form("spline_%i", i), graph);
//        
//        for (int j = 1; j <= h1_interpolatedn->GetNbinsX(); ++j) {
//            double xCenter = h1_interpolatedn->GetXaxis()->GetBinCenter(j);
//            double binMin = h1_interpolatedn->GetXaxis()->GetBinLowEdge(j);
//            double binMax = h1_interpolatedn->GetXaxis()->GetBinUpEdge(j);
//            int numIntegrationSteps = 100;
//            // Integrate the spline over the bin range
//              // Calculate integral using trapezoidal rule
//            double integral = 0.0;
//            double step = (binMax - binMin) / numIntegrationSteps;
//            for (int stepIndex = 0; stepIndex < numIntegrationSteps; ++stepIndex) {
//                double x1 = binMin + stepIndex * step;
//                double x2 = x1 + step;
//                integral += 0.5 * (spline->Eval(x1) + spline->Eval(x2)) * (x2 - x1);
//            }
//            
//            double originalContent = h1_wt->Integral(h1_wt->FindBin(binMin), h1_wt->FindBin(binMax), "width");
//
//            
//            // Scale the spline to match the original bin content
//            if (integral != 0) {
//                double scaleFactor = originalContent / integral;
//                for (int k = 0; k < spline->GetNp(); ++k) {
//                    double xSpline;
//                    double ySpline;
//                    spline->GetKnot(k, xSpline, ySpline);
//                    ySpline *= scaleFactor;
//                }
//            }
//            
//            double interpolatedValue = spline->Eval(xCenter);
//            h1_interpolatedn->SetBinContent(j, interpolatedValue);
//        }
//        
//        int nBins = h1_interpolatedn->GetNbinsX();
//        for (int k = 1; k <= nBins; ++k) {
//            double binCenter = h1_interpolatedn->GetXaxis()->GetBinCenter(k);
//            double binContent = h1_interpolatedn->GetBinContent(k);
//            weightedSum += binCenter * binContent;
//        }
//        
//        h1->SetBinContent(i, weightedSum);
//    }
    
//    //Just with bins
//    for (int i = 1; i <= nBinsX; ++i){
//        double weightedSum = 0;
//        for (int j = 1; j <= nBinsY; ++j) {
//            double binContent = h2->GetBinContent(i, j);
//            double binCenterY = h2->GetYaxis()->GetBinCenter(j); // Bin center in Y direction
//            
//            // Compute weight as the bin center in Y direction
//            double weight = binCenterY;
//            
//            // Accumulate weighted sum and sum of weights
//            weightedSum += binContent * weight;
//        }
//      
//        h1_old->SetBinContent(i, weightedSum);//this is closest to filling on the fly, not the avg
//    }
//    
    
    
    
    
    
//    // Save the 1D histogram to a new file
//    TFile *outFile = TFile::Open(outputFile, "UPDATE");
//    if (!outFile || outFile->IsZombie()) {
//        std::cerr << "Error opening output file: " << outputFile << std::endl;
//        file->Close();
//        return;
//    }
    
    int low = h3->GetYaxis()->GetBinLowEdge(pt1);
    int high = h3->GetYaxis()->GetBinUpEdge(pt2);
    cout<<"low "<<low<< " and high "<<high<<endl;
    
    if(tag == "area"){
        h1->Scale(1./h1->Integral(),"width");
        h1->SetName(Form("h1%s_%s_AreaNorm_%i_%i",bins.c_str(),type.c_str(),low,high));
        h1->SetLineColor(kRed);
        
    }
    else if(tag == "njet"){
        h1->Scale(1./hjet->Integral(pt1,pt2),"width");
        h1->SetName(Form("h1%s_%s_NjetNorm_%i_%i",bins.c_str(),type.c_str(),low,high));
        h1->SetLineColor(kBlack);
        
    }
    else{
        h1->SetLineColor(kBlue);
        h1->Scale(1.,"width");
        h1->SetName(Form("h1%s_%s_%i_%i",bins.c_str(),type.c_str(),low,high));
        
    }
    
    
    
    if(tag == "area"){
        h1_old->Scale(1./h1->Integral(),"width");
        h1_old->SetName(Form("h1%s_%s_AreaNorm_%i_%i",bins.c_str(),type.c_str(),low,high));
        h1_old->SetLineColor(kRed);
        
    }
    else if(tag == "njet"){
        h1_old->Scale(1./hjet->Integral(pt1,pt2),"width");
        h1_old->SetName(Form("h1%s_%s_NjetNorm_%i_%i",bins.c_str(),type.c_str(),low,high));
        h1_old->SetLineColor(kBlack);
        
    }
    else{
        h1_old->SetLineColor(kBlack);
        h1_old->Scale(1.,"width");
        h1_old->SetName(Form("h1%s_%s_%i_%i",bins.c_str(),type.c_str(),low,high));
        
    }
    
    
    if(tag == "area"){
        h1_otherWay->Scale(1./h1->Integral(),"width");
        h1_otherWay->SetName(Form("h1%s_%s_AreaNorm_%i_%i",bins.c_str(),type.c_str(),low,high));
        h1_otherWay->SetLineColor(kRed);
        
    }
    else if(tag == "njet"){
        h1_otherWay->Scale(1./hjet->Integral(pt1,pt2),"width");
        h1_otherWay->SetName(Form("h1%s_%s_NjetNorm_%i_%i",bins.c_str(),type.c_str(),low,high));
        h1_otherWay->SetLineColor(kBlack);
        
    }
    else{
        h1_otherWay->SetLineColor(kGreen);
        h1_otherWay->Scale(1.,"width");
        h1_otherWay->SetName(Form("h1%s_%s_%i_%i",bins.c_str(),type.c_str(),low,high));
        
    }
    
    
    
    
    //What you would normally get
    TH2D* h2_nom;
    if(type =="EEC")
    {
        h2_nom = (TH2D*)file->Get("eec_pt_hist_unf");
    }
    else{
        h2_nom = (TH2D*)file->Get("e3c_pt_hist_unf");
    }
    TH1D* h1_nom = (TH1D*)h2_nom->ProjectionX(Form("h1nom_%s_%i_%i",type.c_str(),low,high),pt1,pt2);
    if(tag == "area"){
        h1_nom->Scale(1./h1_nom->Integral(),"width");
        h1_nom->SetName(Form("h1nom_%s_AreaNorm_%i_%i",type.c_str(),low,high));
    }
    else if(tag == "njet")
    {
        
        h1_nom->Scale(1./hjet->Integral(pt1,pt2),"width");
        h1_nom->SetName(Form("h1nom_%s_NjetNorm_%i_%i",type.c_str(),low,high));
        
    }
    else{
        h1_nom->Scale(1.,"width");
        h1_nom->SetLineColor(kRed);
    }
    
    //    h1_nom->Write();
    //    h1_wt->Write();
    //
    TCanvas* c1= new TCanvas();
//    h1->Draw();
//    h1->SetStats(0);
    h1_nom->Draw("same");
//    h1_old->Draw("same");
    h1_otherWay->Draw("Same");
    
    TLegend *legend = new TLegend(0.12, 0.25, 0.30, 0.50);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);
    legend->AddEntry(h1_nom,"Nominal");
//    legend->AddEntry(h1,"Spline and summed new 1D curve");
    legend->AddEntry(h1_otherWay,"Interpolate");
//    legend->AddEntry(h1_old,"Basic bin weights summation");
    legend->Draw();
    
    
    //    outFile->Close();
    //    file->Close();
    
    std::cout << "1D histogram saved to " << outputFile << std::endl;
}


