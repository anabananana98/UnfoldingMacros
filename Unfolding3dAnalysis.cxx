#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;
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
#include <TRandom3.h>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#endif

Double_t delR(double phi1, double phi2, double eta1, double eta2)
{
    double dphi = abs(phi1-phi2);
    if(dphi>TMath::Pi()){dphi = (2.*TMath::Pi() - dphi);}
    
    double deta = std::fabs(eta1 - eta2);
    Double_t delR = std::sqrt(dphi*dphi + deta*deta);
    
    return delR;
}

void Unfolding3dAnalysis(std::string tag = "", int pthardbin = -1) {

    //***************************************************
    TH1D* h1_true; TH1D* h1_trueSplit;
    TH1D* h1_reco;
    TH1D* h1_unfolded; TH1D* h1_unfolded_rw;  TH1D* h1_unfoldedSplit;
    TH1D* h1_raw; TH1D* h1_rawSplit;
    TH1D* h1_smeared;
    TH1D* h1_fulleff_match;
    TH1D* h1_eff_match;
    TH1D* h1_fullreco;
    TH1D* h1_rw;
    TH1* h1_fold;
    TH1* h1_triv;
    
    TH3D* h3_true_eec; TH3D* h3_trueSplit;
    TH3D* h3_reco_eec;
    TH3D* h3_unfolded; TH3D* h3_unfolded_rw;  TH3D* h3_unfoldedSplit;
    TH3D* h3_raw; TH3D* h3_rawSplit;
    TH3D* h3_smeared;
    TH3D* h3_fulleff_match;
    TH3D* h3_eff_match;
    TH3D* h3_fullreco;
    TH3D* h3_rw;
    TH3* h3_fold;
    TH3* h3_triv;
    
    TFile *fmc;
    TFile *file;
    TFile* f_rw;
    RooUnfoldResponse response1D;
    RooUnfoldResponse response1D_rw;
    RooUnfoldResponse response1DSplit;
    
    RooUnfoldResponse response3D;
    RooUnfoldResponse response3D_rw;
    RooUnfoldResponse response3DSplit;
   
    float pThard_val[21] = {5., 7., 9., 12., 16., 21., 28., 36., 45., 57., 70., 85., 99., 115., 132., 150., 169., 190., 212., 235., 1000.};
    
    //reco jet pT bins
    Double_t from_const_reco = 15;
    Double_t to_const_reco = 120;
    Int_t bins_const_reco = 21;
    Double_t width_const_reco = (to_const_reco-from_const_reco)/bins_const_reco;
    Double_t xbins[22] = {};
    for (Int_t i = 0; i <= bins_const_reco; i++)
    {
        xbins[i] = (from_const_reco + i * width_const_reco);
    }
    
    //true jet pT bins
    Double_t from_const = 0;
    Double_t to_const = 400;
    // Double_t to_const = 800;
    // Int_t bins_const = 37;
    // Int_t bins_const = 160;
    Int_t bins_const = 10;
    // Int_t bins_const = 22;
    Double_t width_const = (to_const-from_const)/bins_const;
    Double_t tbins[11] = {};
    for (Int_t i = 0; i <= bins_const; i++)
    {
        tbins[i] = (from_const + i * width_const);
    }
    
    
    //reco dR bins --same for eec and e3c
    double dRbins[] = {
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
    double dRbinst[] = {
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
    
    // h3_raw = (TH3D*)file->Get("jet_pt_hist");
    // h3_raw->SetName("h3_raw");
    // h1_rw = (TH1D*)f_rw->Get("h1_smooth");
   
    int pthard = -1;
    
    
      // TH1D* h1_true(0);
    h1_true = new TH1D("h1_true", "h1_true", bins_const, tbins);
    h1_trueSplit = new TH1D("h1_trueSplit", "h1_trueSplit", bins_const, tbins);
    h1_eff_match = new TH1D("h1_eff_match", "h1_eff_match", bins_const, tbins);
    h1_fulleff_match = new TH1D("h1_fulleff_match", "h1_fulleff_match", bins_const, tbins);
     
    h1_reco = new TH1D("h1_reco", "h1_reco", 21, xbins);
    h1_fullreco = new TH1D("h1_fullreco", "h1_fullreco", 21, xbins);
    
    h1_rawSplit = new TH1D("h1_rawSplit", "h1_rawSplit", 21, xbins);
    
    h1_smeared = new TH1D("h1_smeared", "h1_smeared", 21, xbins);
    h1_rw = (TH1D*)h1_raw->Clone("h1_rw");  // Clone h1_raw to preserve it
    
    
    h1_raw->Sumw2();
    h1_smeared->Sumw2();
    h1_true->Sumw2();
    h1_fullreco->Sumw2();
    h1_rw->Sumw2();
    
    response1D.Setup(h1_reco,h1_true);
    response1D_rw.Setup(h1_reco,h1_true);
    response1DSplit.Setup(h1_reco,h1_true);
    
    
    h3_true_eec = new TH3D("h3_true_eec", "h3_true_eec", 33, dRbinst, bins_const, tbins, 12,wtbinst);
    h3_true_eec->RebinY(4);
    h3_reco_eec = new TH3D("h3_reco_eec", "h3_reco",33, dRbins, bins_const_reco, xbins, 12,wtbins);
    h3_reco_eec->RebinY(2);

    h3_true_eec->Sumw2();
    h3_reco_eec->Sumw2();
  
    cout<<"setting up response declarations"<<endl;
    response3D.Setup(h3_reco_eec,h3_true_eec);

    cout<<"ending declarations"<<endl;
    
    // for(int i = 1; i<21; i++)
    // for(int i = 2; i<5; i++)
    // {
    // pthard = i;
      //Get the mc files
     fmc = TFile::Open(Form("/home/ar2545/Trial/AnalysisResultsMC%i.root",pthardbin));
    
    TString file_name = Form("/home/ar2545/Trial/AnalysisResultsMC%i.root", pthardbin);

    cout<<"file name is: "<<file_name.Data()<<endl;

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
    
    // tree->GetEntry(1);
    // double mc = mc_weight; //same for one run in 1 pTHat bin
    // double weightresp = mc;
    
    double ptJetmin_det = 20.; //min det level pt
    double ptJetmax_det = 80.; //max det level pt
    double jetpT_tru = -1;
    double jetpT_det = -1;
    
   
    //In this tree I only have matched jets, so tree entries where either jet pT is 0 is only for missed(true level) or fake tracks(det level)
     Long64_t nEntries = tree->GetEntries();
    cout<<nEntries<<endl;


    double pt_tru;
    double pt_det;
    double pt_miss;
    double pt_fake;

    double eta_tru;
    double eta_det;
    double eta_miss;
    double eta_fake;

    double phi_tru;
    double phi_det;
    double phi_miss;
    double phi_fake;

    // cout<<"MC weight is: "<< mc_weight<< endl;
    // Loop over entries in the tree
    for (Int_t i = 0; i < nEntries; ++i) {
        if(i==nEntries-1) break;
        tree->GetEntry(i);
        // cout<<"i "<<i<<endl;
        // chain.GetEntry(i);
        double mc = mc_weight;
        jetpT_tru = fJet_pt_tru;
        jetpT_det = fJet_pt_det;

        pt_tru = fTrack_pt_tru;
        pt_det = fTrack_pt_det;
        pt_miss = fTrack_pt_miss;
        pt_fake = fTrack_pt_fake;

        eta_tru = fTrack_eta_tru;
        eta_det = fTrack_eta_det;
        eta_miss = fTrack_eta_miss;
        eta_fake = fTrack_eta_fake;

        phi_tru = fTrack_phi_tru;
        phi_det = fTrack_phi_det;
        phi_miss = fTrack_phi_miss;
        phi_fake = fTrack_phi_fake;

        
        
        
        
       // matched correlations
        if(jetpT_det==0 || jetpT_tru==0) continue;
       
        h1_reco->Fill(jetpT_det, mc);
        h1_true->Fill(jetpT_tru, mc); 
            //2-pt corr
            for (Int_t k = i+1; k < nEntries; ++k) {
                tree->GetEntry(k);
                // chain.GetEntry(k);
                
                if(jetpT_det!=fJet_pt_det || jetpT_tru!=fJet_pt_tru) break; //make sure they are the same jet
                
                // cout<<" i k "<<i<< " "<<k<<endl;
                // cout<<" idet kdet "<<jetpT_det<< " "<<fJet_pt_det<<endl;
                
                if(fTrack_pt_tru!=0 && fTrack_pt_det!=0)
                {
                    double wt_tru = (2*pt_tru*fTrack_pt_tru)/(jetpT_tru*jetpT_tru);
                    double wt_det =(pt_det*fTrack_pt_det)/(jetpT_det*jetpT_det);
                    
                    double dR_tru = delR(phi_tru,fTrack_phi_tru,eta_tru,fTrack_eta_tru);
                    double dR_det = delR(phi_det,fTrack_phi_det,eta_det,fTrack_eta_det);
                    
                    // eec_pt_hist->Fill(dR_tru,jetpT_tru,wt_tru*mc);
                    // eec_pt_hist_det->Fill(dR_det,jetpT_det,wt_det*mc);
                   
                    h3_reco_eec->Fill(dR_det,jetpT_det,wt_det*mc);
                    h3_true_eec->Fill(dR_tru,jetpT_tru,wt_tru*mc);
                    
                    response3D.Fill(dR_det, dR_tru, jetpT_det, jetpT_tru, wt_det, wt_tru, mc);
                }
            }
            
           if (fJet_pt_det>=ptJetmin_det && fJet_pt_det<=ptJetmax_det && fJet_pt_tru<=to_const)//this is the correct thing to do
            { 
                h1_smeared->Fill(fJet_pt_det, mc);
                response1D.Fill(fJet_pt_det, fJet_pt_tru, mc);
            
           }
        
        jetpT_tru = fJet_pt_tru;
        jetpT_det = fJet_pt_det;
    
        } 
    
    // }
    
    
    
    int nIter = 4; // Number of iterations
    TFile *fout=new TFile (Form("/home/ar2545/Trial/Unfold3DTrivClosure_%i_%i.root", nIter,pthardbin),"RECREATE");
    // TFile *fout=new TFile (Form("/home/ar2545/Trial/testUnf.root", nIter),"RECREATE");
    fout->cd();
    
    cout<<"write histograms"<<endl;
    h1_raw->Write();
    h1_true->Write();
    h1_reco->Write();
    h1_fullreco->Write();
    response1D.Write();
  
     
    h1_eff_match->Write();
    h1_fulleff_match->Write();
    
   
    cout<<"Performing 1D unfolding"<<endl;
    // Create an unfolding object using the Bayes method
    RooUnfoldBayes unfold(&response1D, h1_raw, nIter);
    RooUnfoldBayes unfoldTriv(&response1D, h1_reco, nIter);
    

    // Perform unfolding
    h1_unfolded = (TH1D*)unfold.Hreco();
    h1_unfolded->SetName("h1_unfolded");
    
    h1_triv = (TH1D*)unfoldTriv.Hreco();//trivial unfolding 
    h1_triv->SetName("h1_triv");
    h1_triv->Write();
    
    
    // Draw the unfolded histogram
    cout<<"write unfolded & folded histograms"<<endl;
  //For 3D unfolding 
    h1_unfolded->Write();
    // h1_fold->Write();
   
///////////////////////////#####################################//////////////////////////////
    cout<<"Performing 3D unfolding"<<endl;
    response3D.Write();

    // Create an unfolding object using the Bayes method
    // RooUnfoldBayes unfold3D(&response3D, h3_raw, nIter);
    RooUnfoldBayes unfoldTriv3D(&response3D, h3_reco_eec, nIter);
   
    // Perform unfolding
    // h3_unfolded = (TH3D*)unfold3D.Hreco();
    // h3_unfolded->SetName("h3_unfolded");
    
    h3_triv = (TH3D*)unfoldTriv3D.Hreco();//trivial unfolding 
    h3_triv->SetName("h3_triv");
    h3_triv->Write();
  
    //FOLD BACK
    // h3_fold = response3D.ApplyToTruth(h3_unfolded, "h3_fold");
    
    // Draw the unfolded histogram
    cout<<"write unfolded & folded histograms for 3D"<<endl;
    
    h3_true_eec->Write();
    h3_reco_eec->Write();
    h3_triv->Write();


 
    // Close the files
    fout->Close();
    file->Close();
    fmc->Close();
   
 
    

}
#ifndef __CINT__
int main () { Unfolding3dAnalysis(); return 0; }  // Main program when run stand-alone
#endif
