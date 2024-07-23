#include <iostream>
#include <fstream>
//#include "PlottingHeaderUC.h"

///Checks for unfolding
void OptBin(int pt1, int pt2, std::string hist_name, bool ifres, bool ifeec, bool ifR, bool ifwt,bool ifwtnojet, bool ifsave)
{

//Opt_Un_eec
//Opt_Un_e3c
    std::string fname = "/Users/ar2545/AnalysisDownloads/AnalysisResults_Data_Apr16.root";
    
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
    
    TList *hList = (TList*)f->Get(Form("JetsEEC_Jet_AKTChargedR040_tracks_pT0150_E_scheme_TCRaw_Data_NoSub_Incl;1"));
    
    TCanvas *ca;
    TH3D* h_un3d= (TH3D*)hList->FindObject("Opt_Un_eec");
    h_un3d->GetYaxis()->SetRange(pt1,pt2);
    TH2D* h_un2d= (TH2D*)h_un3d->Project3D("zx");
    h_un2d->Sumw2();
    cout<<h_un2d->GetNbinsX()<<endl;
    TCanvas *canvas = new TCanvas();
    h_un2d->Draw("text");
    
    TCanvas *canvas1 = new TCanvas();
    Double_t from = -4;
    Double_t to = 0;
    Int_t bins = 100;
    Double_t width = (to-from)/bins;
    Double_t new_bins[101] = {};
    for (Int_t i = 0; i <= bins; i++)
    {
        new_bins[i] = TMath::Power(10.0, from + i * width);
//        cout<<new_bins[i]<<endl;
    }

//    //RL bins
    Int_t mergedR = 31;
    Double_t fromR = -4;
//    Double_t toR = -3;
    Double_t toR = -2.80086;
    Int_t binsR = 1;
    Double_t widthR = (toR-fromR)/binsR;
    Int_t Rbins_reb_new = (bins-mergedR)+1;
    Double_t new_binsR[2] = {};
    Double_t knew_reb_binsR[70] = {};
    for (Int_t i = 0; i <= binsR; i++)
    {
        new_binsR[i] = TMath::Power(10.0, fromR + i * widthR);
        knew_reb_binsR[i] = new_binsR[i];
    }
    for (Int_t i = binsR+1; i <= Rbins_reb_new; i++)
    {
//    cout<<i<<" and "<<i + (mergedR-1)<<endl;
        knew_reb_binsR[i] = new_bins[i + (mergedR-1)];
//        cout<<"next bin : "<<knew_reb_binsR[i]<<endl;
    }
    
    
//    Int_t mergedR_again = 43;
//    Double_t fromR_again = -4;
//    Double_t toR_again = -2.319;
    Int_t mergedR_again = 12;
    Double_t toR_again = -2.279;
    Double_t fromR_again = -2.319;
    Int_t binsR_again = 1;
    Double_t widthR_again = (toR_again-fromR_again)/binsR_again;
    Int_t Rbins_reb_new_again = (bins-mergedR_again)+1;
    Double_t new_binsR_again[2] = {};
    Double_t knew_reb_binsR_again[59] = {};
    for (Int_t i = 0; i <= binsR; i++)
    {
        cout<<"first loop "<<i<<endl;
        new_binsR[i] = TMath::Power(10.0, fromR + i * widthR);
        knew_reb_binsR_again[i] = new_binsR[i];
        cout<<"next bin : "<<knew_reb_binsR_again[i]<<endl;
//        cout<<knew_reb_binsR_again[i]<<endl;
    }
//    for (Int_t i = binsR+1; i <= Rbins_reb_new; i++)
//    {
//    cout<<i<<" and "<<i + (mergedR-1)<<endl;
//        knew_reb_binsR_again[i] = new_bins[i + (mergedR-1)];
//        cout<<"next bin : "<<knew_reb_binsR_again[i]<<endl;
//    }
    for (Int_t i = 0 ; i <= binsR_again; i++)
    {
        cout<<"second loop "<<i<<endl;
        new_binsR_again[i] = TMath::Power(10.0, fromR_again + i * widthR_again);
//        cout<<new_binsR_again[i]<<endl;
        knew_reb_binsR_again[i+ binsR+1] = new_binsR_again[i];
        cout<<i+ binsR+1<<endl;
        cout<<"next bin : "<<knew_reb_binsR_again[i+binsR+1]<<endl;
    }
    for (Int_t i = binsR+3; i <= 59; i++)
    {
        cout<<"third loop "<<i<<endl;
        knew_reb_binsR_again[i] = new_bins[i + (mergedR + mergedR_again-2)];
        cout<<"next bin : "<<knew_reb_binsR_again[i]<<endl;
    }
    
    
   
    
    //eec weight bins
    Double_t wt_from = -4;
    Double_t wt_to = 0;
    Int_t wt_bins = 200;
    Double_t wt_width = (wt_to-wt_from)/wt_bins;
    Double_t wt_new_bins[201] = {};
    for (Int_t i = 0; i <= wt_bins; i++)
    {
        wt_new_bins[i] = TMath::Power(10.0, wt_from + i * wt_width);
    }
    

    // Define the array as a constexpr to ensure compile-time constant
    constexpr double values[] = {
        0.0001, 0.000104713, 0.000109648, 0.000114815, 0.000120226, 0.000125893, 0.000131826,
        0.000138038, 0.000144544, 0.000151356, 0.000158489, 0.000165959, 0.00017378, 0.00018197,
        0.000190546, 0.000199526, 0.00020893, 0.000218776, 0.000229087, 0.000239883, 0.000251189,
        0.000263027, 0.000275423, 0.000288403, 0.000301995, 0.000316228, 0.000331131, 0.000346737,
        0.000363078, 0.000380189, 0.000398107, 0.000416869, 0.000436516, 0.000457088, 0.00047863,
        0.000501187, 0.000524807, 0.000549541, 0.00057544, 0.00060256, 0.000630957, 0.000660693,
        0.000691831, 0.000724436, 0.000758578, 0.000794328, 0.000831764, 0.000870964, 0.000912011,
        0.000954993, 0.001, 0.00104713, 0.00109648, 0.00114815, 0.00120226, 0.00125893, 0.00131826,
        0.00138038, 0.00144544, 0.00151356, 0.00158489, 0.00165959, 0.0017378, 0.0018197, 0.00190546,
        0.00199526, 0.0020893, 0.00218776, 0.00229087, 0.00239883, 0.00251189, 0.00263027, 0.00275423,
        0.00288403, 0.00301995, 0.00316228, 0.00331131, 0.00346737, 0.00363078, 0.00380189, 0.00398107,
        0.00416869, 0.00436516, 0.00457088, 0.0047863, 0.00501187, 0.00524807, 0.00549541, 0.0057544,
        0.0060256, 0.00630957, 0.00660693, 0.00691831, 0.00724436, 0.00758578, 0.00794328, 0.00831764,
        0.00870964, 0.00912011, 0.00954993, 0.01, 0.0104713, 0.0109648, 0.0114815, 0.0120226, 0.0125893,
        0.0131826, 0.0138038, 0.0144544, 0.0151356, 0.0158489, 0.0165959, 0.017378, 0.018197, 0.0190546,
        0.0199526, 0.020893, 0.0218776, 0.0229087, 0.0239883, 0.0251189, 0.0263027, 0.0275423, 0.0288403,
        0.0301995, 0.0316228, 0.0331131, 0.0346737, 0.0363078, 0.0380189, 0.0398107, 0.0416869, 0.0436516,
        0.0457088, 0.047863, 0.0501187, 0.0524807, 0.0549541, 0.057544, 0.060256, 0.0630957, 0.0660693,
        0.0691831, 0.0724436, 0.0758578, 0.0794328, 0.0831764, 0.0870964, 0.0912011, 0.0954993, 0.1,
        0.104713, 0.109648, 0.114815, 0.120226, 0.125893, 0.131826, 0.138038, 0.144544, 0.151356,
        0.158489, 0.165959, 0.17378, 0.18197, 0.190546, 0.199526, 0.20893, 0.218776, 0.229087, 0.239883,
        0.251189, 0.263027, 0.275423, 0.288403, 0.301995, 0.316228, 0.331131, 0.346737, 0.363078,
        0.380189, 0.398107, 0.416869, 0.436516, 0.457088, 0.47863, 0.501187, 0.524807, 0.549541,
        0.57544, 0.60256, 0.630957, 0.660693, 0.691831, 0.724436, 0.758578, 0.794328, 0.831764,
        0.870964, 0.912011, 0.954993, 1.0
    };
    
    double Rebvalues[] = {
        0.0001, 0.000104713, 0.000109648, 0.000114815, 0.000120226, 0.000125893, 0.000131826,
        0.000138038, 0.000144544, 0.000151356, 0.000158489, 0.000165959, 0.00017378, 0.00018197,
        0.000190546, 0.000199526, 0.00020893, 0.000218776, 0.000229087, 0.000239883, 0.000251189,
        0.000263027, 0.000275423, 0.000288403, 0.000301995, 0.000316228, 0.000331131, 0.000346737,
        0.000363078, 0.000380189, 0.000398107, 0.000416869, 0.000436516, 0.000457088, 0.00047863,
        0.000501187, 0.000524807, 0.000549541, 0.00057544, 0.00060256, 0.000630957, 0.000660693,
        0.000691831, 0.000724436, 0.000758578, 0.000794328, 0.000831764, 0.000870964, 0.000912011,
        0.000954993, 0.001, 0.00104713, 0.00109648, 0.00114815, 0.00120226, 0.00125893, 0.00131826,
        0.00138038, 0.00144544, 0.00151356, 0.00158489, 0.00165959, 0.0017378, 0.0018197, 0.00190546,
        0.00199526, 0.0020893, 0.00218776, 0.00229087, 0.00239883, 0.00251189, 0.00263027, 0.00275423,
        0.00288403, 0.00301995, 0.00316228, 0.00331131, 0.00346737, 0.00363078, 0.00380189, 0.00398107,
        0.00416869, 0.00436516, 0.00457088, 0.0047863, 0.00501187, 0.00524807, 0.00549541, 0.0057544,
        0.0060256, 0.00630957, 0.00660693, 0.00691831, 0.00724436, 0.00758578, 0.00794328, 0.00831764,
        0.00870964, 0.00912011, 0.00954993, 0.01, 0.0104713, 0.0109648, 0.0114815, 0.0120226, 0.0125893,
        0.0131826, 0.0138038, 0.0144544, 0.0151356, 0.0158489, 0.0165959, 0.017378, 0.018197, 0.0190546,
        0.0199526, 0.020893, 0.0218776, 0.0229087, 0.0239883, 0.0251189, 0.0263027, 0.0275423, 0.0288403,
        0.0301995, 0.0316228, 0.0331131, 0.0346737, 0.0363078, 0.0380189, 0.0398107, 0.0416869, 0.0436516,
        0.0457088, 0.047863, 0.0501187, 0.0524807, 0.0549541, 0.057544, 0.060256, 0.0630957, 0.0660693,
        0.0691831, 0.0724436, 0.0758578, 0.0794328, 0.0831764, 0.0870964, 0.0912011, 0.0954993, 0.1,
        0.104713, 0.109648, 0.114815, 0.120226, 0.125893, 0.131826, 0.138038, 0.144544, 0.151356,
        0.158489, 0.165959, 0.17378, 0.18197, 0.190546, 0.199526, 0.20893, 0.218776, 0.229087, 0.239883,
        0.251189, 0.263027, 0.275423, 0.288403, 0.301995, 0.316228, 0.331131, 0.346737, 0.363078,
        0.380189, 0.398107, 0.416869, 0.436516, 0.457088, 0.47863, 0.501187, 0.524807, 0.549541,
        0.57544, 0.60256, 0.630957, 0.660693, 0.691831, 0.724436, 0.758578, 0.794328, 0.831764,
        0.870964, 0.912011, 0.954993, 1.0
    };

    //Rebinning low wt by a factor of 5 and keeping high wt the same
    //#no. of merged bins is 101
    Int_t merged = 151;
    Double_t from_reb = -4;
    Double_t to_reb = -2;
    Int_t bins_reb = 1;
    Int_t bins_reb_new = (wt_bins-merged)+1;
    Double_t width_reb = (to_reb-from_reb)/bins_reb;
    Double_t new_bins_reb[2] = {};
    Double_t kBinsNewReb[151] = {};
    for (Int_t i = 0; i <= bins_reb; i++)
    {
        new_bins_reb[i] = TMath::Power(10.0, from_reb + i * width_reb);
        kBinsNewReb[i] = new_bins_reb[i];
//        cout<<kBinsNewReb[i]<<endl;
    }
    
    for (Int_t i = bins_reb+1; i <= bins_reb_new; i++)
    {
        kBinsNewReb[i] = wt_new_bins[i + (merged-1)];
//        cout<<kBinsNewReb[i]<<endl;
    }
    

    bool isIncreasingR = true;
    for (int i = 0; i <= bins_reb_new; ++i) {
        if (kBinsNewReb[i] <= kBinsNewReb[i-1]) {
            isIncreasingR = false;
            cout<<i<< " and "<<i-1<<endl;
            cout<<kBinsNewReb[i]<<" and "<<kBinsNewReb[i-1]<<endl;
            break;
        }
    }
    
    if (isIncreasingR) {
        printf("The array is in increasing order.\n");
    } else {
        printf("The array is not in increasing order.\n");
    }

// TH2D* h_un2d_reb = (TH2D*)h_un2d->Clone("h_un2d_reb");
// h_un2d->Reset();
    TH2D* h_un2d_reb_X = new TH2D("h_un2d_reb_X", "Rebinned x (Both Axes)", 59, knew_reb_binsR_again, 200, wt_new_bins);
      TH2D* h_un2d_reb = new TH2D("h_un2d_reb", "Rebinned 2D Histogram (Both Axes)", 59, knew_reb_binsR_again, bins_reb_new, kBinsNewReb);
//      TH2D* texthist = new TH2D("texthist", "Rebinned 2D Histogram (Both Axes)", Rbins_reb_new_again, knew_reb_binsR_again, bins_reb_new, kBinsNewReb);
//    TH2D* h_un2d_reb_X = new TH2D("h_un2d_reb_X", "Rebinned x (Both Axes)", Rbins_reb_new, knew_reb_binsR, 200, wt_new_bins);
//    TH2D* h_un2d_reb = new TH2D("h_un2d_reb", "Rebinned 2D Histogram (Both Axes)", Rbins_reb_new, knew_reb_binsR, bins_reb_new, kBinsNewReb);

//
//    double summed_x = 0;
//    for (int ybin = 1; ybin <= h_un2d->GetNbinsY(); ++ybin) {
//        for (int xbin = 1; xbin <= mergedR; ++xbin) {
////        cout<<"x bin val "<<xbin<<endl;
//            summed_x += h_un2d->GetBinContent(xbin,ybin);
//        }
//////        cout<<"y bin val "<<ybin<<endl;
//       h_un2d_reb_X->SetBinContent(1,ybin,h_un2d_reb_X->GetBinContent(1,ybin) + summed_x);
//       summed_x = 0;
//    }
////}

//     for (int ybin = 1; ybin <= h_un2d_reb_X->GetNbinsY(); ++ybin) {
//        for (int xbin = 2; xbin <= h_un2d_reb_X->GetNbinsX(); ++xbin) {
//             h_un2d_reb_X->SetBinContent(xbin,ybin,h_un2d->GetBinContent(xbin,ybin + (mergedR-1)));
////             if(h_un2d->GetBinContent(xbin,ybin+101)!=0)
////             {
////            cout<<h_un2d->GetBinContent(xbin,ybin+101)<<endl;
////            }
//        }
//    }
////    
//    double summed = 0;
//    for (int xbin = 1; xbin <= h_un2d_reb_X->GetNbinsX(); ++xbin) {
//        for (int ybin = 1; ybin <= merged; ++ybin) {
//            summed += h_un2d_reb_X->GetBinContent(xbin,ybin);
//        }
//    
//       h_un2d_reb->SetBinContent(xbin,1,summed);
//       summed = 0;
//    }
//
for (int i = 1; i <= h_un2d->GetNbinsX(); ++i) {
        for (int j = 1; j <= h_un2d->GetNbinsY(); ++j) {
            double content = h_un2d->GetBinContent(i, j);
            double error = h_un2d->GetBinError(i, j);
            double x = h_un2d->GetXaxis()->GetBinCenter(i);
            double y = h_un2d->GetYaxis()->GetBinCenter(j);
            h_un2d_reb_X->Fill(x, y, content);
            }
        }
        
//        h_un2d_reb_X->Scale(1.,"width");
for (int i = 1; i <= h_un2d_reb_X->GetNbinsX(); ++i) {
        for (int j = 1; j <= h_un2d_reb_X->GetNbinsY(); ++j) {
            double content = h_un2d_reb_X->GetBinContent(i, j);
            double error = h_un2d_reb_X->GetBinError(i, j);
            double x = h_un2d_reb_X->GetXaxis()->GetBinCenter(i);
            double y = h_un2d_reb_X->GetYaxis()->GetBinCenter(j);
            h_un2d_reb->Fill(x, y, content);
//            if (content < 20) {
//                texthist->SetBinContent(i, j, content);
//            } else {
//                texthist->SetBinContent(i, j, 0); // Ensure other bins are empty
//            }
            }
        }
        
        
//    for (int xbin = 1; xbin <= h_un2d_reb->GetNbinsX(); ++xbin) {
//        for (int ybin = 2; ybin <= h_un2d_reb->GetNbinsY(); ++ybin) {
//             h_un2d_reb->SetBinContent(xbin,ybin,h_un2d_reb_X->GetBinContent(xbin,ybin + (merged-1)));
////             if(h_un2d->GetBinContent(xbin,ybin+101)!=0)
////             {
////            cout<<h_un2d->GetBinContent(xbin,ybin+101)<<endl;
////            }
//        }
//    }

//    for (int xbin = 1; xbin <= h_un2d_reb->GetNbinsX(); ++xbin) {
//        for (int ybin = 1; ybin <= h_un2d_reb->GetNbinsY(); ++ybin) {
//             h_un2d_reb->SetBinContent(xbin,ybin,ybin);
////             if(h_un2d->GetBinContent(xbin,ybin+101)!=0)
////             {
////            cout<<h_un2d->GetBinContent(xbin,ybin+101)<<endl;
////            }
//        }
//    }
//
////    
//    // Clean up intermediate histograms if needed
//    delete h_un2d_reb_x;
//h_un2d_reb->Scale(1.,"width");

h_un2d_reb->Draw("text");
//texthist->Draw("colz,text");
// Now, h_un2d_reb contains the rebinned 2D histogram with both axes rebinned

 
 
 
    ////    h_un->RebinX(2);
    //    h_un->RebinY(5);
//    h_un->SetTitle("");
//    
//    h_un->GetXaxis()->SetRangeUser(0.01, 0.4);
//    //    h_un->GetYaxis()->SetRangeUser(-1, 1);
//    h_un2d->Draw("text");
//    h_un2d_reb_X->Draw("text");
//    h_un->SetStats(0);
//    h_un->GetXaxis()->SetTitle(str17);
//    h_un->GetYaxis()->SetTitle("weight");

    
}

