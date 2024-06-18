#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;


#include "TFile.h"
#include "TH1D.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#endif



void UnfoldingAnalysis(std::string tag = "") {

    //***************************************************
    TH1D* h1_true;
    TH1D* h1_reco;
    TH1D* h1_unfolded;
    TH1D* h1_raw;
    TH1D* h1_smeared;
    TH1D* h1_fulleff_match;
    TH1D* h1_fullreco;
    TH1D* h1_rw;
    TH1* h1_fold;
    
    TFile *fmc;
    TFile *file;
    TFile* f_rw;
    RooUnfoldResponse response1D;
    
    float pThard_val[21] = {5., 7., 9., 12., 16., 21., 28., 36., 45., 57., 70., 85., 99., 115., 132., 150., 169., 190., 212., 235., 1000.};
    
    //true & reco jet pT bins
    Double_t from_const = 15;
    Double_t to_const = 120;
    Int_t bins_const = 21;
    Double_t width_const = (to_const-from_const)/bins_const;
    Double_t tbins[22] = {};
    Double_t xbins[22] = {};
    for (Int_t i = 0; i <= bins_const; i++)
    {
        tbins[i] = (from_const + i * width_const);
        xbins[i] = (from_const + i * width_const);
    }
    
    
    //reco dR bins --same for eec and e3c
    double ybins[] = {
        0.0001,
        0.01,
        0.0158489, 0.0190546, 0.0229087,
        0.0275423, 0.0331131,
        0.0398107, 0.0436516, 0.047863, 0.0524807, 0.057544,
        0.0630957, 0.0691831, 0.0758578, 0.0831764, 0.0912011,
        0.1, 0.109648, 0.120226, 0.131826, 0.144544,
        0.158489, 0.17378, 0.190546, 0.20893, 0.229087,
        0.251189, 0.275423, 0.301995, 0.331131, 0.363078,
        0.398107, 0.524807,
        1.0
    };
    
    //true dR bins --same for eec and e3c
    double ybinst[] = {
        0.0001,
        0.01,
        0.0158489, 0.0190546, 0.0229087,
        0.0275423, 0.0331131,
        0.0398107, 0.0436516, 0.047863, 0.0524807, 0.057544,
        0.0630957, 0.0691831, 0.0758578, 0.0831764, 0.0912011,
        0.1, 0.109648, 0.120226, 0.131826, 0.144544,
        0.158489, 0.17378, 0.190546, 0.20893, 0.229087,
        0.251189, 0.275423, 0.301995, 0.331131, 0.363078,
        0.398107, 0.524807,
        1.0
    };
    
    
    //reco wt bins eec
    double wtbins[] = {
        0.0001, 0.00301995, 0.00501187, 0.01, 0.0158489, 0.020893, 0.0301995, 0.0501187, 0.060256, 0.1,
        0.151356, 0.758578,1.0};
    
    //true wt bins eec
    double wtbinst[] = {
        0.0001, 0.00301995, 0.00501187, 0.01, 0.0158489, 0.020893, 0.0301995, 0.0501187, 0.060256, 0.1,
        0.151356, 0.758578,1.0};
    
    //reco wt bins e3c
    double wtbins_e3c[] = {
        0.0001, 0.00301995, 0.00501187, 0.01, 0.0158489, 0.020893, 0.0301995, 0.0501187, 0.060256, 0.1,
        0.151356, 0.758578,1.0};
    
    //true wt bins e3c
    double wtbinst_e3c[] = {
        0.0001, 0.00301995, 0.00501187, 0.01, 0.0158489, 0.020893, 0.0301995, 0.0501187, 0.060256, 0.1,
        0.151356, 0.758578,1.0};
    
    // Create a ROOT file to read histograms
    file = TFile::Open("/home/ar2545/Trial/Data.root");
    // f_rw = new TFile("../../Lund/Reweight_Nov23.root"); //this file to reweight prior comes from herwig & herwig fast sim
    
    h1_raw = (TH1D*)file->Get("jet_pt_hist");
    h1_raw->SetName("h1_raw");
    // h1_rw = (TH1D*)f_rw->Get("h1_smooth");
   
    int pthard = -1;
    
    
      // TH1D* h1_true(0);
    h1_true = new TH1D("h1_true", "h1_true", 21, tbins);
    h1_fulleff_match = new TH1D("h1_fulleff_match", "h1_fulleff_match", 21, tbins);
     
    h1_reco = new TH1D("h1_reco", "h1_reco", 21, xbins);
    h1_fullreco = new TH1D("h1_fullreco", "h1_fullreco", 21, xbins);
    
    h1_smeared = new TH1D("h1_smeared", "h1_smeared", 21, xbins);
   
    h1_raw->Sumw2();
    h1_smeared->Sumw2();
    h1_true->Sumw2();
    h1_fullreco->Sumw2();
    
    response1D.Setup(h1_reco,h1_true);
    
    for(int i = 1; i<21; i++)
    {
    pthard = i;
      //Get the mc files
     fmc = TFile::Open(Form("/home/ar2545/Trial/AnalysisResultsMC%i.root",pthard));
    

    TTree *tree=(TTree*)fmc->Get("MatchTracksTree");
    tree->Show(0);
  
    Int_t nEv=tree->GetEntries();
   
    Double_t fJet_pt_det, fJet_pt_tru;
    Double_t fTrack_pt_det, fTrack_pt_tru, fTrack_pt_miss, fTrack_pt_fake;
    Double_t fTrack_eta_det, fTrack_eta_tru, fTrack_eta_miss, fTrack_eta_fake;
    Double_t fTrack_phi_det, fTrack_phi_tru, fTrack_phi_miss, fTrack_phi_fake;
    Double_t mc_weight;
    
    tree->SetBranchAddress("fJet_pt_det", &fJet_pt_det);
    tree->SetBranchAddress("fJet_pt_tru", &fJet_pt_tru);
    tree->SetBranchAddress("fTrack_pt_det", &fTrack_pt_det);
    tree->SetBranchAddress("fTrack_pt_tru", &fTrack_pt_tru);
    tree->SetBranchAddress("fTrack_pt_miss", &fTrack_pt_miss);
    tree->SetBranchAddress("fTrack_pt_fake", &fTrack_pt_fake);
    tree->SetBranchAddress("fTrack_eta_det", &fTrack_eta_det);
    tree->SetBranchAddress("fTrack_eta_tru", &fTrack_eta_tru);
    tree->SetBranchAddress("fTrack_eta_miss", &fTrack_eta_miss);
    tree->SetBranchAddress("fTrack_eta_fake", &fTrack_eta_fake);
    tree->SetBranchAddress("fTrack_phi_det", &fTrack_phi_det);
    tree->SetBranchAddress("fTrack_phi_tru", &fTrack_phi_tru);
    tree->SetBranchAddress("fTrack_phi_miss", &fTrack_phi_miss);
    tree->SetBranchAddress("fTrack_phi_fake", &fTrack_phi_fake);
    tree->SetBranchAddress("mc_weight", &mc_weight);
    
    tree->GetEntry(1);
    double mc = mc_weight; //same for one run in 1 pTHat bin
    
    double ptJetconst = 0.;
    double jetpT_tru = 0;
    double jetpT_det = 0;
    
     
    
    for(int iEntry=0; iEntry< nEv; iEntry++){
        tree->GetEntry(iEntry);
        
        if(jetpT_tru == fJet_pt_tru || jetpT_det == fJet_pt_det) continue;
        //fill the 1D response for unfolding the jet pT spectra for the normalization
        float dpt = fJet_pt_det - fJet_pt_tru;
        dpt = std::sqrt(dpt*dpt);
        if ((dpt > 0.) && (jetpT_det<=120) && (jetpT_det>=15))
        {
            // ptJetconst = ptJet;
            h1_reco->Fill(fJet_pt_det, mc);
            h1_true->Fill(fJet_pt_tru, mc);
            h1_smeared->Fill(fJet_pt_det, mc);
            response1D.Fill(fJet_pt_det, fJet_pt_tru, mc);
        }
        
        jetpT_tru = fJet_pt_tru;
        jetpT_det = fJet_pt_det;
        
    }
        
    //     //reweighting the response to the prior
    //  double weightresp = mc;
    //  int binx = -1;
    //   if (lnr_part < 0) binx = 1;
    //   else if (lnr_part > 1.4) binx = h3_rw->GetNbinsX();
    //   else binx = h3_rw->GetXaxis()->FindBin(lnr_part);
    //   int biny = -1;
    //   if (lnz_part < -1.0) biny = 1;
    //   else if (lnz_part > 1.5) biny = h3_rw->GetNbinsY();
    //   else biny = h3_rw->GetYaxis()->FindBin(lnz_part);
    //   int binz = -1;
    //   if (ptJetMatch < pTmin) binz = 1;
    //   else if (ptJetMatch > 120) binz = h3_rw->GetNbinsZ();
    //   else binz = h3_rw->GetZaxis()->FindBin(ptJetMatch);
    //   double rw_fac = h3_rw->GetBinContent(binx, biny, binz);
    //   if (tag == "rw") weightresp*=rw_fac;

    }
    
    int nIter = 10; // Number of iterations
    TFile *fout=new TFile (Form("/home/ar2545/Trial/Unfold1DJetPt_%i.root", nIter),"RECREATE");
    fout->cd();
    
    h1_raw->Write();
    cout<<"write"<<endl;
    h1_true->Write();
    cout<<"write"<<endl;
    h1_reco->Write();
    cout<<"write"<<endl;
    response1D.Write();
    
    
   
    cout<<"Performing unfolding"<<endl;
    // Create an unfolding object using the Bayes method
   
    RooUnfoldBayes unfold(&response1D, h1_raw, nIter);

    // Perform unfolding
    h1_unfolded = (TH1D*)unfold.Hreco();
    h1_unfolded->SetName("h1_unfolded");
    
   
    //FOLD BACK
    h1_fold = response1D.ApplyToTruth(h1_unfolded, "h1_fold");
    
    // Draw the unfolded histogram
    h1_unfolded->Write();
    h1_fold->Write();

    fout->Close();
    // Close the file
    file->Close();
    fmc->Close();
    

}
#ifndef __CINT__
int main () { UnfoldingAnalysis(); return 0; }  // Main program when run stand-alone
#endif

