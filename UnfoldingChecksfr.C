#include <iostream>
#include <fstream>
//#include "PlottingHeaderUC.h"

///Checks for unfolding
void UnfoldingChecksfr(int pt1, int pt2, std::string hist_name, bool ifres, bool ifeec, bool ifR, bool ifwt,bool ifwtnojet, bool ifsave)
{
    
    //void UnfoldingChecksfr(int pt1, int pt2)
    //{
    
    //    std::string fname = "/Users/ar2545/AnalysisDownloads/AnalysisResultsMC_Apr19_CorrectBinsUC.root";
            std::string fname = "/Users/ar2545/AnalysisDownloads/AnalysisResultsMC_May15_JESplot.root";
    
//    std::string fname = "/Users/ar2545/AnalysisDownloads/AnalysisResults_Data_Apr16.root";
//    std::string fname = "/Users/ar2545/Downloads/TestE3COutput.root";
//     std::string fname = "/Users/ar2545/Downloads/TestOutputMCE3C_new.root";
    
    TFile *f = new TFile(fname.c_str(), "READ");
    
    //-------------------------------------------------------------------
    //    //-------------Defining strings------------------------------
    //-------------------------------------------------------------------
    const char *str = "p^{ch,det}_{T,jet}";
    const char *str3 = "#sqrt{s} = 13 TeV";
    const char *str2 = "anti-#it{k}_{T}";
    const char *str1 = "#it{p}^{ch}_{T,min}";
    const char *str4 = "#it{R}_{L}";
    const char *str5 = "#it{p}^{ch,true}_{T,jet}";
    const char *str6 = "#it{p}^{ch}_{T,jet}";
    const char *str13 = "#it{R}_{L}^{det}";
    const char *str14 = "#it{R}_{L}^{tru}";
    const char *str15 = "#it{wt}^{det}";
    const char *str16 = "#it{wt}^{tru}";
    const char *str17 = "#it{R}_{L}";
    const char *str18 = "#it{z}^{true}_{T}";
    const char *str19 = "#it{z}^{det}_{T}";
    const char *str20 = "#it{p}^{true}_{T}";
    const char *str21 = "#it{p}^{det}_{T}";
    const char *str22 = "#it{wt}^{det_nj}";
    const char *str23 = "#it{wt}^{tru_nj}";
    const char *str24 = "#it{wt}^{nj}";
    //-------------------------------------------------------------------
    Double_t textSizeRel = 0.04;
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(textSizeRel);
    latex.SetTextFont(42);
    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(textSizeRel);
    latex1.SetTextFont(42);
    TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextSize(textSizeRel);
    latex2.SetTextFont(42);
    
    latex.SetTextColor(kBlack);
    latex2.SetTextColor(kBlack);
    //-------------------------------------------------------------------
    
    string jetpt = str6;
    string plot_title ="";
    
   
    //    TH1D *constituentId;
    //    TH3D *R_res_eec_tru_debug;
    //    TH2D *track_pt_res_debug;
    //    TH3D *track_eta_debug;
    //    TH2D *track_phi_debug;
    //    TH3D *track_R_debug;
    //
    //    TH2D *track_pt_response_debug;
    //    TH2D *track_pt_wt_res_debug;
    //    TH2D *track_pt_wt_response_debug;
    //
    //    TH1D* fHistPtHard;
    //    TH1D* jet_pt_det;
    //    TH1D* jet_pt_tru;
    //
    //    TH2D* R_matrix_jetpt;
    ////    TH3D* R_match_eec_tru;
    //    TH3D* R_match_e3c_tru;
    //    TH3D* wtnojet_match_eec_tru;
    //    TH3D* wt_match_e3c_tru;
    //    TH3D* wt_match_eec_tru;
    //    TH3D* R_res_eec_tru;
    //
    //    TH2D* R_match_eec1;
    //    TH2D* R_match_e3c1;
    //    TH2D* R_match_debug;
    //
    //    TH3D* wt_res_eec;
    //    TH3D* wt_res_eec_tru;
    //    TH3D* wt_res_e3c_tru;
    //    TH3D* R_res_eec;
    //    TH3D* R_res_e3c;
    //
    //    TH3D* wtnojet_match_e3c_tru;
    //    TH3D* wtnojet_res_eec_tru;
    //    TH3D* wtnojet_res_e3c_tru;
    //
    //    TH2D* wt_res_eec_2d_tru;
    //    TH2D* wt_res_e3c_2d_tru;
    //    TH2D* wtnojet_res_eec_2d_tru;
    //    TH2D* wtnojet_res_eec_2d_tru;
    //    TH2D* wtnojet_res_e3c_2d_tru;
    //    TH2D* R_res_eec_2d;
    //    TH2D* R_res_e3c_2d;
    //    TH2D* wtnojet_match_eec_2d_tru;
    //    TH2D* wtnojet_match_e3c_2d_tru;
    
    //    TH2D* wtnojet_res_e3c_2d;
    //    TH2D* wt_match_e3c_2d_tru;
    //    TH2D* wt_match_eec_2d_tru;
    //
    //    TH1D* wt_res_eec_bin;
    //    TH1D* wt_res_e3c_bin;
    //    TH1D* R_res_eec_bin;
    //    TH1D* R_res_e3c_bin;
    //    TH1D* wtnojet_match_eec_bin;
    //    TH1D* wtnojet_match_e3c_bin;
    //    TH1D* wtnojet_res_eec_bin;
    //    TH1D* wtnojet_res_e3c_bin;
    
    //    jet_pt_det = (TH1D*)f->Get(TString::Format("jet_pt_hist"));
    //
    //    TCanvas *c1new = new TCanvas();
    //    jet_pt_det->Draw();
    
    
//        TH3D* jet_res_w_wt= (TH3D*)f->Get("jetpt_res_w_wt_e3c");
//    //
//        jet_res_w_wt->GetZaxis()->SetRange(pt1,pt2);
//        TH2D* jet_res_w_wt_yx= (TH2D*)jet_res_w_wt->Project3D("yx");
//    
//        TCanvas *canvas1 = new TCanvas();
//        canvas1->SetLogz();
//        jet_res_w_wt_yx->SetTitle("");
//        jet_res_w_wt_yx->GetXaxis()->SetRangeUser(1, 1600);
//        jet_res_w_wt_yx->GetYaxis()->SetRangeUser(-1, 1);
//        jet_res_w_wt_yx->Draw("colz");
//        jet_res_w_wt_yx->SetStats(0);
//        jet_res_w_wt_yx->GetXaxis()->SetTitle(str23);
//        jet_res_w_wt_yx->GetYaxis()->SetTitle(Form("%s Res",jetpt.c_str()));
//        float low_pt_bin_jet = jet_res_w_wt->GetZaxis()->GetBinLowEdge(pt1);
//        float high_pt_bin_jet = jet_res_w_wt->GetZaxis()->GetBinUpEdge(pt2);
//        TLatex latex_jet;latex_jet.SetNDC();latex_jet.SetTextSize(0.04);latex_jet.SetTextFont(42);
//        latex_jet.DrawLatex(0.60,0.70 ,Form("%.1f < %s < %.1f GeV/c", low_pt_bin_jet, jetpt.c_str(), high_pt_bin_jet));
//    
//        TH2D *jet_res_w_wt_yx_new = (TH2D*)jet_res_w_wt_yx->Clone(Form("newhisto_%s", jet_res_w_wt->GetName()));
//        jet_res_w_wt_yx_new->Reset();
//    
//        for (int binX = 1; binX <= jet_res_w_wt_yx->GetNbinsX(); ++binX) {
//            TH1D *yProjection = jet_res_w_wt_yx->ProjectionY("", binX, binX);
//    
//            for (int binY = 1; binY <= jet_res_w_wt_yx->GetNbinsY(); ++binY) {
//                if (yProjection->GetEntries() != 0) {
//                    double normalizedY = (jet_res_w_wt_yx->GetBinContent(binX, binY)) / yProjection->Integral();
//                    jet_res_w_wt_yx_new->SetBinContent(binX, binY, normalizedY);
//                }
//            }
//            yProjection->Reset();
//        }
//    
//        TCanvas *canvas2 = new TCanvas();
//        canvas2->SetLogz();
//        jet_res_w_wt_yx_new->SetTitle("");
//        jet_res_w_wt_yx_new->GetXaxis()->SetRangeUser(1, 1600);
//        // jet_res_w_wt_yx_new->GetYaxis()->SetRangeUser(0.09, 0.12);
//        jet_res_w_wt_yx_new->Draw("colz");
//        jet_res_w_wt_yx_new->SetStats(0);
//        jet_res_w_wt_yx_new->GetXaxis()->SetTitle(str23);
//        jet_res_w_wt_yx_new->GetYaxis()->SetTitle(Form("%s Res",jetpt.c_str()));
    //
    //    //Jet energy resolution with weight
    //    TCanvas *canssn1 = new TCanvas();
    //    canssn1->SetLogx();
    //    std::stringstream ssn1;
    //    std::stringstream ssn2;
    //    TH1D* h1 = (TH1D*)jet_res_w_wt_yx->ProjectionX("h1");
    //    TH1D* h2 = (TH1D*)jet_res_w_wt_yx->ProjectionX("h2");
    //    // h2->Reset();
    //    for (int i = 1; i <= jet_res_w_wt_yx->GetNbinsX(); i++)
    //    {
    //        ssn1 <<hist_name<<"_mean" << i;
    //        TH1D* hproj = (TH1D*)jet_res_w_wt_yx->ProjectionY(ssn1.str().c_str(), i, i);
    //        ssn1.str("");
    //        h1->SetBinContent(i, hproj->GetStdDev());
    //        h1->SetBinError(i, hproj->GetRMSError());
    //        h2->SetBinContent(i, hproj->GetMean());
    //    }
    //    h1->GetXaxis()->SetRangeUser(9,1600);
    //    h1->GetYaxis()->SetRangeUser(0.02,0.2);
    //    h1->SetStats(0);
    //    h1->SetTitle("");
    //    h1->GetYaxis()->SetTitle(Form("%s Res",jetpt.c_str()));
    //    h1->SetMarkerStyle(8);
    //    //        h1->GetYaxis()->SetTitle(label_yaxis.c_str());
    //    //        h1->GetXaxis()->SetTitle(label_xaxis.c_str());
    //    h1->Draw();
    //    latex_jet.DrawLatex(0.60,0.70 ,Form("%.1f < %s < %.1f GeV/c", low_pt_bin_jet, jetpt.c_str(), high_pt_bin_jet));
    //
    //    //Jet energy scale with weight
    //    TCanvas *canssn2 = new TCanvas();
    //    canssn2->SetLogx();
    //    h2->SetStats(0);
    //    h2->GetXaxis()->SetRangeUser(9,1600);
    //    h2->GetYaxis()->SetRangeUser(0.02,0.2);
    //    h2->SetTitle("");
    //    h2->GetYaxis()->SetTitle(Form("Mean %s",jetpt.c_str()));
    //    h2->SetMarkerStyle(8);
    //    h2->SetMarkerColor(kRed);
    //    h2->Draw();
    //    latex_jet.DrawLatex(0.60,0.70 ,Form("%.1f < %s < %.1f GeV/c", low_pt_bin_jet, jetpt.c_str(), high_pt_bin_jet));
    ////     std::stringstream ssjet;
    ////        ssjet<< "/Users/ar2545/Desktop/UnfoldingChecks/"<<"jetpt_res_w_wt"<<"_" <<  pt1  <<  "to"<< pt2<<"_"<<"Apr29.pdf";
    ////        canvas1->SaveAs(ssjet.str().c_str());
    //
    //  std::stringstream ssjet_width;
    //        ssjet_width<< "/Users/ar2545/Desktop/UnfoldingChecks/"<<"jetpt_res_w_wt"<<"_" <<  pt1  <<  "to"<< pt2<<"_"<<
    //        "_WIDTH_Apr29.pdf";
    //        canssn1->SaveAs(ssjet_width.str().c_str());
    //        std::stringstream ssjet_mean;
    //     ssjet_mean<< "/Users/ar2545/Desktop/UnfoldingChecks/"<<"jetpt_res_w_wt"<<"_" <<  pt1  <<  "to"<< pt2<<"_"<<"_MEAN_Apr29.pdf";
    //     canssn2->SaveAs(ssjet_mean.str().c_str());
    //
    //
    
    
    //     //--------------------------------------------------------------------------------
    //     //Draw response matrices
    //     //--------------------------------------------------------------------------------
    
    
         if(!ifres)
         {
         TCanvas *c1 = new TCanvas();
             string label_xaxis;
             string label_yaxis;
             TH3D *histo = (TH3D*)f->Get(hist_name.c_str());
    
             if (ifR)
             {
                 c1->SetLogx();
                 c1->SetLogy();
                 label_xaxis = str14;
                 label_yaxis = str13;
             }
             else if (ifwt)
             {
                 label_xaxis = str16;
                 label_yaxis = str15;
             }
             else
             {
     //        if(!ifeec)
     //        {
     //            histo->RebinX(10);
     //            histo->RebinY(10);
     //        }
                 label_xaxis = str23;
                 label_yaxis = str24;
             }
             c1->SetLogz();
    
    
             histo->GetZaxis()->SetRange(pt1, pt2);
             TH2D *h_yx = (TH2D*)histo->Project3D("yx");
    
             float low_pt_bin = histo->GetZaxis()->GetBinLowEdge(pt1);
             cout<<low_pt_bin<<endl;
             cout<<histo->GetName()<<endl;
             float high_pt_bin = histo->GetZaxis()->GetBinUpEdge(pt2);
    
             TH2D *h_yx_new = (TH2D*)h_yx->Clone(Form("newhisto_%s", histo->GetName()));
             h_yx_new->Reset();
    
             for (int binX = 1; binX <= h_yx->GetNbinsX(); ++binX) {
                 TH1D *yProjection = h_yx->ProjectionY("", binX, binX);
    
                 for (int binY = 1; binY <= h_yx->GetNbinsY(); ++binY) {
                     if (yProjection->GetEntries() != 0) {
                         double normalizedY = (h_yx->GetBinContent(binX, binY)) / yProjection->Integral();
                         h_yx_new->SetBinContent(binX, binY, normalizedY);
                     }
                 }
                 yProjection->Reset();
             }
    
             h_yx_new->SetTitle(plot_title.c_str());
             h_yx_new->GetYaxis()->SetTitle(label_yaxis.c_str());
             h_yx_new->GetXaxis()->SetTitle(label_xaxis.c_str());
             h_yx_new->SetStats(0);
             if (ifR)
             {
                 h_yx_new->GetXaxis()->SetRangeUser(0.005, 0.4);
                 h_yx_new->GetYaxis()->SetRangeUser(0.005, 0.4);
             }
             if(ifwt){
                 h_yx_new->GetXaxis()->SetRangeUser(0.0, 0.6);
                 h_yx_new->GetYaxis()->SetRangeUser(0.0, 0.6);
             }
             if(ifwtnojet){
                 if(ifeec){
                     if(pt1==2)
                     {
                         h_yx_new->GetXaxis()->SetRangeUser(0.0,1500);
                         h_yx_new->GetYaxis()->SetRangeUser(0.0,1500);
                     }
                     if(pt1==6)
                     {
                         h_yx_new->GetXaxis()->SetRangeUser(0.0,2000);
                         h_yx_new->GetYaxis()->SetRangeUser(0.0,2000);
                     }
                     if(pt1==10)
                     {
                         h_yx_new->GetXaxis()->SetRangeUser(0.0,3500);
                         h_yx_new->GetYaxis()->SetRangeUser(0.0,3500);
                     }
                 }
             }
             h_yx_new->Draw("colz");
    
              for (int binX = 1; binX <= h_yx_new->GetNbinsX(); ++binX) {
                 TH1D *yProjection = h_yx->ProjectionY("", binX, binX);
    
                 for (int binY = 1; binY <= h_yx_new->GetNbinsY(); ++binY) {
                     double mean = yProjection->GetMean();
                     double stddev = yProjection->GetStdDev();
    
                     }
                 }
    
    
    
    
    
    
             if(ifeec)
             {
                 if(ifR)
                 {
                     latex.DrawLatex(0.15, 0.75, Form("Response matrix for %s, EEC",str17));
                 }
                 else
                 {
                     latex.DrawLatex(0.15, 0.75, Form("Response matrix for weight, EEC"));
                 }
             }
             else
             {
                 if(ifR)
                 {
                     latex.DrawLatex(0.15, 0.75, Form("Response matrix for %s, E3C",str17));
                 }
                 else
                 {
                     latex.DrawLatex(0.15, 0.75, Form("Response matrix for weight, E3C"));
                 }
             }
             latex1.DrawLatex(0.15, 0.65, Form("%.1f < %s < %.1f GeV/c", low_pt_bin, jetpt.c_str(), high_pt_bin));
             latex2.DrawLatex(0.15, 0.70, Form("R = 0.4"));
    
    
    
             std::stringstream ss;
             if(ifeec)
             {
                 if(ifR)
                 {
                     ss<< "/Users/ar2545/Desktop/UnfoldingChecks/EEC_Rmatrix"<<"_" <<  pt1  <<  "to"<< pt2 << "Apr22.pdf";
                 }
                 else if(ifwt)
                 {
                     ss<< "/Users/ar2545/Desktop/UnfoldingChecks/Wt_EEC_Rmatrix"<<"_" <<  pt1  <<  "to"<< pt2 << "Apr22.pdf";
                 }
                 else {
                     ss<< "/Users/ar2545/Desktop/UnfoldingChecks/WtNJ_EEC_Rmatrix"<<"_" <<  pt1  <<  "to"<< pt2 << "Apr22.pdf";
                 }
             }
             else
             {
                 if(ifR)
                 {
                     ss<< "/Users/ar2545/Desktop/UnfoldingChecks/E3C_Rmatrix"<<"_" <<  pt1  <<  "to"<< pt2 << "Apr22.pdf";
                 }
                 else if(ifwt)
                 {
                     ss<< "/Users/ar2545/Desktop/UnfoldingChecks/Wt_E3C_Rmatrix"<<"_" <<  pt1  <<  "to"<< pt2 << "Apr22.pdf";
                 }
                 else {
                     ss<< "/Users/ar2545/Desktop/UnfoldingChecks/WtNJ_E3C_Rmatrix"<<"_" <<  pt1  <<  "to"<< pt2 << "Apr22.pdf";
                 }
             }
    
             if(ifsave)
             {
                 c1->SaveAs(ss.str().c_str());
             }
         }
         //--------------------------------------------------------------------------------
         //Draw resolution histograms
         //--------------------------------------------------------------------------------
        if(ifres)
        {
                TCanvas *c1 = new TCanvas();
            string label_xaxis;
                 string label_yaxis;
    //        string label_yaxis = "Res";
    //        string label_yaxis = R;
            TH3D *histo = (TH3D*)f->Get(hist_name.c_str());
            histo->GetZaxis()->SetRange(pt1, pt2);
            float low_pt_bin = histo->GetZaxis()->GetBinLowEdge(pt1);
            cout<<low_pt_bin<<endl;
            cout<<histo->GetName()<<endl;
            float high_pt_bin = histo->GetZaxis()->GetBinUpEdge(pt2);
    
            c1->cd();
            c1->SetLogz();
    
    
            TH2D *h_yx = (TH2D*)histo->Project3D("yx");
//            TCanvas *canssnabbd = new TCanvas();
//            h_yx->Draw();
    
    //        TH2D *h_yx = (TH2D*)f->Get(hist_name.c_str());
    //        float low_pt_bin = h_yx->GetZaxis()->GetBinLowEdge(pt1);
    //        cout<<low_pt_bin<<endl;
    //        float high_pt_bin = h_yx->GetZaxis()->GetBinUpEdge(pt2);
            h_yx->SetMarkerStyle(8);
    
            TH2D *h_yx_new = (TH2D*)h_yx->Clone(Form("newhisto_%s", histo->GetName()));
    //        TH2D *h_yx_new = (TH2D*)h_yx->Clone(Form("newhisto_%s", "JES"));
            h_yx_new->Reset();
            h_yx_new->Sumw2();
    
            for (int binX = 1; binX <= h_yx->GetNbinsX(); ++binX) {
                TH1D *yProjection = h_yx->ProjectionY("", binX, binX);
    
                for (int binY = 1; binY <= h_yx->GetNbinsY(); ++binY) {
                    if (yProjection->GetEntries() != 0) {
                        double normalizedY = (h_yx->GetBinContent(binX, binY)) / yProjection->Integral();
                        h_yx_new->SetBinContent(binX, binY, normalizedY);
                    }
                }
                yProjection->Reset();
            }
    
    
    
    
            if(ifR)
            {
                c1->SetLogx();
                label_xaxis = str14;
    
            }
            else if (ifwt)
            {
                label_xaxis = str16;
    
            }
            else
            {
    
                label_xaxis = str23;
                label_yaxis = str24;
            }
    
    cout<<"HERE "<<label_yaxis<<endl;
    
            if(ifR==true)
            {
                h_yx_new->GetXaxis()->SetRangeUser(0.01,0.4);
                h_yx_new->GetYaxis()->SetRangeUser(-0.1,0.1);
    
            }
            h_yx_new->SetTitle(plot_title.c_str());
            h_yx_new->GetYaxis()->SetTitle(label_yaxis.c_str());
            h_yx_new->GetXaxis()->SetTitle(label_xaxis.c_str());
            h_yx_new->SetStats(0);
            h_yx_new->Draw("colz");
    
            TLatex latex;latex.SetNDC();latex.SetTextSize(textSizeRel);latex.SetTextFont(42);
            TLatex latex1;latex1.SetNDC();latex1.SetTextSize(textSizeRel);latex1.SetTextFont(42);
    
    
            if(ifeec)
            {
                if(ifR)
                {
                    latex.DrawLatex(0.15, 0.85, Form("Resolution for %s, EEC",str17));
                }
                else
                {
                    latex.DrawLatex(0.15, 0.75, Form("Resolution for weight, EEC"));
                }
            }
            else
            {
                if(ifR)
                {
                    latex.DrawLatex(0.15, 0.75, Form("Resolution for %s, E3C",str17));
                }
                else
                {
                    latex.DrawLatex(0.15, 0.75, Form("Resolution for weight, E3C"));
                }
            }
            latex1.DrawLatex(0.15,0.80 ,Form("%.1f < %s < %.1f GeV/c", low_pt_bin, jetpt.c_str(), high_pt_bin));
    
            std::stringstream ss;
            ss<< "/Users/ar2545/Desktop/UnfoldingChecks/"<<hist_name<<"_" <<  pt1  <<  "to"<< pt2<<"_norm_"<<"Apr29.pdf";
    
            if(ifsave)
            {
                c1->SaveAs(ss.str().c_str());
            }
    
    
    //         //Plot mean values
            TCanvas *canssn = new TCanvas();
            canssn->SetLogx();
            std::stringstream ssn;
            TH1D* h1_JER = (TH1D*)h_yx->ProjectionX("h1_JER");
            TH1D* h1_JES = (TH1D*)h_yx->ProjectionX("h1_JES");
            for (int i = 1; i <= h_yx->GetNbinsX(); i++)
            {
                ssn <<hist_name<<"_mean" << i;
                TH1D* hproj = (TH1D*)h_yx->ProjectionY(ssn.str().c_str(), i, i);
                ssn.str("");
                h1_JER->SetBinContent(i, hproj->GetStdDev());
                h1_JER->SetBinError(i, hproj->GetRMSError());
                h1_JES->SetBinContent(i, hproj->GetMean());
            }
            h1_JER->GetXaxis()->SetRangeUser(0.005,0.5);
            h1_JER->GetYaxis()->SetRangeUser(0,0.25);
            h1_JER->SetStats(0);
            h1_JER->SetTitle("");
//            h1_JER->GetYaxis()->SetTitle(Form("%s Res", label_yaxis.c_str()));
            h1_JER->GetYaxis()->SetTitle("JER");
            //        h1_JER->GetXaxis()->SetTitle(label_xaxis.c_str());
            h1_JER->GetXaxis()->SetTitle(label_xaxis.c_str());
//            h1_JER->GetYaxis()->SetTitle(label_yaxis.c_str());
    //        h1_JER->GetXaxis()->SetTitle(str5);
    //         h1_JER->SetMaximum(1.8);
            // h1_JER->SetMinimum(0.0);
            // h1_JER->RebinX(2);
            h1_JER->SetMarkerStyle(8);
    //        h1_JER->SetMarkerSize(0.05);
            h1_JER->SetMarkerColor(kBlack);
            h1_JER->Draw();
            latex1.DrawLatex(0.15,0.80 ,Form("%.1f < %s < %.1f GeV/c, EEC", low_pt_bin, jetpt.c_str(), high_pt_bin));
    
            TCanvas *canssn1 = new TCanvas();
            canssn1->SetLogx();
            h1_JES->SetStats(0);
            h1_JES->SetTitle("");
//            h1_JES->GetXaxis()->SetRangeUser(-1,1);
            h1_JES->GetXaxis()->SetRangeUser(0.005,0.5);
            h1_JES->GetYaxis()->SetRangeUser(0.01,0.2);
//            h1_JES->GetYaxis()->SetTitle(Form("Mean %s",label_yaxis.c_str()));
            h1_JES->GetYaxis()->SetTitle("JES");
            //         h1_JES->GetYaxis()->SetTitle("JES");
            h1_JES->GetYaxis()->SetTitleOffset(1.2);
            h1_JES->GetXaxis()->SetTitle(label_xaxis.c_str());
//            h1_JES->GetYaxis()->SetTitle(label_yaxis.c_str());
    //        h1_JES->GetXaxis()->SetTitle(str5);
    //         h1_JES->GetXaxis()->SetTitle(str23);
             // h1_JES->SetMaximum(0.15);
             // h1_JES->SetMinimum(0.0);
             h1_JES->SetMarkerStyle(8);
    //         h1_JES->SetMarkerSize(0.02);
             h1_JES->SetMarkerColor(kRed);
             h1_JES->Draw();
    
//              if(ifeec)
//             {
//                 if(ifR)
//                 {
//                     latex.DrawLatex(0.15, 0.85, Form("Resolution for %s, EEC",str17));
//                 }
//                 else
//                 {
//                     latex.DrawLatex(0.15, 0.75, Form("Resolution for weight, EEC"));
//                 }
//             }
//             else
//             {
//                 if(ifR)
//                 {
//                     latex.DrawLatex(0.15, 0.75, Form("Resolution for %s, E3C",str17));
//                 }
//                 else
//                 {
//                     latex.DrawLatex(0.15, 0.75, Form("Resolution for weight, E3C"));
//                 }
//             }
    
    
            latex1.DrawLatex(0.15,0.80 ,Form("%.1f < %s < %.1f GeV/c, EEC", low_pt_bin, jetpt.c_str(), high_pt_bin));
             if(ifsave)
             {
                std::stringstream ssname;
                ssname<<"/Users/ar2545/Desktop/UnfoldingChecks/" <<hist_name<<"_WIDTH_"<<pt1<<"to"<<pt2<<"_"<<"May15.pdf";
                canssn->SaveAs(ssname.str().c_str());
    
                 std::stringstream ssname1;
                 ssname1<<"/Users/ar2545/Desktop/UnfoldingChecks/" <<hist_name<<"_MEAN_"<<pt1<<"to"<<pt2<<"_"<<"May15.pdf";
                 canssn1->SaveAs(ssname1.str().c_str());
             }
    
    
    
        }
    
}

