//Checking weight binning for E3C
//possibility to check trivial case or rebin a super fine histogram and then check
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


void gettingwtCheck(int pt1, int pt2) {
    const char* inputFile = "/Users/ar2545/Downloads/E3C_dataOct24.root";
    std::string histName3D_unfolded = "Opt_Un_e3c_unf";
    std::string histName3D = "Opt_Un_e3c";
    std::string histName2D = "e3c_pt_hist_unf";
    
    
//    std::string histName3D_unfolded = "h3_triv";
//    std::string histName3D = "h3_true_eecTriv";
//    std::string histName2D = "eec_pt_histTriv";
//    const char* inputFile = "/Users/ar2545/Downloads/UnfoldingResultsTrivEEC_4EEC.root";
    
    const char *str = "#it{p}_{T,jet}";
    const char *str1 = "#it{p}_{T,min}";
    const char *str2 = "anti-#it{k}_{T}";
    const char *str3 = "#sqrt{s} = 13 TeV";
    const char *str4 = "#it{R}_{L}";
    const char *str8 = "#it{#eta}_{jet}";
    const char *str11 = "#it{p}_{T}^{const}";
    
   // Open the input file
    TFile *f = TFile::Open(inputFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file: " << inputFile << std::endl;
        return;
    }

    TH3D* h3 = (TH3D*)f->Get(histName3D.c_str());
    h3->GetXaxis()->SetTitle(str4);
    h3->GetYaxis()->SetTitle(str);
    h3->GetZaxis()->SetTitle("wt");
    
    int pt1jet = h3->GetYaxis()->FindBin(pt1);
    int pt2jet = h3->GetYaxis()->FindBin(pt2);
    
    h3->GetYaxis()->SetRange(pt1jet,pt2jet);
    
    TH2D *h2 = (TH2D*)h3->Project3D("zx"); //Rl vs wt histogram
    TH1D *h1 = (TH1D*)h2->ProjectionX("h1");
    h1->Reset();
 
    
       // Create the 1D histogram
    int nBinsX = h2->GetNbinsX();
    double xMin = h2->GetXaxis()->GetXmin();
    double xMax = h2->GetXaxis()->GetXmax();
 
    TH3D* h3_unfolded = (TH3D*)f->Get(histName3D_unfolded.c_str());
    h3_unfolded->GetXaxis()->SetTitle(str4);
    h3_unfolded->GetYaxis()->SetTitle(str);
    h3_unfolded->GetZaxis()->SetTitle("wt");
    
    if (!h3_unfolded) {
        std::cerr << "Error retrieving histogram: " << histName3D_unfolded << std::endl;
        f->Close();
        return;
    }
     h3_unfolded->GetYaxis()->SetRange(h3_unfolded->GetYaxis()->FindBin(pt1),h3_unfolded->GetYaxis()->FindBin(pt2));

//    // Get the 2D histogram
//    TH2D *h2_unfolded = (TH2D*)h3_unfolded->Project3D("zx"); //Rl vs wt histogram
//    TH1D *h1_unfolded = (TH1D*)h2_unfolded->ProjectionX("h1_unfolded");
//    h1_unfolded->Reset();
 
 
///############################SET UP WHICH KIND OF PROJECTION YOU WANT TO DO------------------------------------------------------------

   
     double RvaluesReb[] = {
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
//    
//    constexpr double values_e3c[] = {
//        1e-06, 1.03514e-06, 1.07152e-06, 1.10917e-06, 1.14815e-06, 1.1885e-06, 1.23027e-06, 1.2735e-06,
//        1.31826e-06, 1.36458e-06, 1.41254e-06, 1.46218e-06, 1.51356e-06, 1.56675e-06, 1.62181e-06, 1.6788e-06,
//        1.7378e-06, 1.79887e-06, 1.86209e-06, 1.92752e-06, 1.99526e-06, 2.06538e-06, 2.13796e-06, 2.21309e-06,
//        2.29087e-06, 2.37137e-06, 2.45471e-06, 2.54097e-06, 2.63027e-06, 2.7227e-06, 2.81838e-06, 2.91743e-06,
//        3.01995e-06, 3.12608e-06, 3.23594e-06, 3.34965e-06, 3.46737e-06, 3.58922e-06, 3.71535e-06, 3.84592e-06,
//        3.98107e-06, 4.12098e-06, 4.2658e-06, 4.4157e-06, 4.57088e-06, 4.73151e-06, 4.89779e-06, 5.06991e-06,
//        5.24807e-06, 5.4325e-06, 5.62341e-06, 5.82103e-06, 6.0256e-06, 6.23735e-06, 6.45654e-06, 6.68344e-06,
//        6.91831e-06, 7.16143e-06, 7.4131e-06, 7.67361e-06, 7.94328e-06, 8.22243e-06, 8.51138e-06, 8.81049e-06,
//        9.12011e-06, 9.44061e-06, 9.77237e-06, 1.01158e-05, 1.04713e-05, 1.08393e-05, 1.12202e-05, 1.16145e-05,
//        1.20226e-05, 1.24451e-05, 1.28825e-05, 1.33352e-05, 1.38038e-05, 1.42889e-05, 1.47911e-05, 1.53109e-05,
//        1.58489e-05, 1.64059e-05, 1.69824e-05, 1.75792e-05, 1.8197e-05, 1.88365e-05, 1.94984e-05, 2.01837e-05,
//        2.0893e-05, 2.16272e-05, 2.23872e-05, 2.31739e-05, 2.39883e-05, 2.48313e-05, 2.5704e-05, 2.66073e-05,
//        2.75423e-05, 2.85102e-05, 2.95121e-05, 3.05492e-05, 3.16228e-05, 3.27341e-05, 3.38844e-05, 3.50752e-05,
//        3.63078e-05, 3.75837e-05, 3.89045e-05, 4.02717e-05, 4.16869e-05, 4.31519e-05, 4.46684e-05, 4.62381e-05,
//        4.7863e-05, 4.9545e-05, 5.12861e-05, 5.30884e-05, 5.49541e-05, 5.68853e-05, 5.88844e-05, 6.09537e-05,
//        6.30957e-05, 6.53131e-05, 6.76083e-05, 6.99842e-05, 7.24436e-05, 7.49894e-05, 7.76247e-05, 8.03526e-05,
//        8.31764e-05, 8.60994e-05, 8.91251e-05, 9.22571e-05, 9.54993e-05, 9.88553e-05, 0.000102329, 0.000105925,
//        0.000109648, 0.000113501, 0.00011749, 0.000121619, 0.000125893, 0.000130317, 0.000134896, 0.000139637,
//        0.000144544, 0.000149624, 0.000154882, 0.000160325, 0.000165959, 0.000171791, 0.000177828, 0.000184077,
//        0.000190546, 0.000197242, 0.000204174, 0.000211349, 0.000218776, 0.000226464, 0.000234423, 0.000242661,
//        0.000251189, 0.000260016, 0.000269153, 0.000278612, 0.000288403, 0.000298538, 0.00030903, 0.00031989,
//        0.000331131, 0.000342768, 0.000354813, 0.000367282, 0.000380189, 0.00039355, 0.00040738, 0.000421697,
//        0.000436516, 0.000451856, 0.000467735, 0.000484172, 0.000501187, 0.0005188, 0.000537032, 0.000555904,
//        0.00057544, 0.000595662, 0.000616595, 0.000638263, 0.000660693, 0.000683912, 0.000707946, 0.000732825,
//        0.000758578, 0.000785236, 0.000812831, 0.000841395, 0.000870964, 0.000901571, 0.000933254, 0.000966051,
//        0.001, 0.00103514, 0.00107152, 0.00110917, 0.00114815, 0.0011885, 0.00123027, 0.0012735, 0.00131826,
//        0.00136458, 0.00141254, 0.00146218, 0.00151356, 0.00156675, 0.00162181, 0.0016788, 0.0017378, 0.00179887,
//        0.00186209, 0.00192752, 0.00199526, 0.00206538, 0.00213796, 0.00221309, 0.00229087, 0.00237137, 0.00245471,
//        0.00254097, 0.00263027, 0.0027227, 0.00281838, 0.00291743, 0.00301995, 0.00312608, 0.00323594, 0.00334965,
//        0.00346737, 0.00358922, 0.00371535, 0.00384592, 0.00398107, 0.00412098, 0.0042658, 0.0044157, 0.00457088,
//        0.00473151, 0.00489779, 0.00506991, 0.00524807, 0.0054325, 0.00562341, 0.00582103, 0.0060256, 0.00623735,
//        0.00645654, 0.00668344, 0.00691831, 0.00716143, 0.0074131, 0.00767361, 0.00794328, 0.00822243, 0.00851138,
//        0.00881049, 0.00912011, 0.00944061, 0.00977237, 0.0101158, 0.0104713, 0.0108393, 0.0112202, 0.0116145,
//        0.0120226, 0.0124451, 0.0128825, 0.0133352, 0.0138038, 0.0142889, 0.0147911, 0.0153109, 0.0158489,
//        0.0164059, 0.0169824, 0.0175792, 0.018197, 0.0188365, 0.0194984, 0.0201837, 0.020893, 0.0216272, 0.0223872,
//        0.0231739, 0.0239883, 0.0248313, 0.025704, 0.0266073, 0.0275423, 0.0285102, 0.0295121, 0.0305492,
//        0.0316228, 0.0327341, 0.0338844, 0.0350752, 0.0363078, 0.0375837, 0.0389045, 0.0402717, 0.0416869,
//        0.0431519, 0.0446684, 0.0462381, 0.047863, 0.049545, 0.0512861, 0.0530884, 0.0549541, 0.0568853,
//        0.0588844, 0.0609537, 0.0630957, 0.0653131, 0.0676083, 0.0699842, 0.0724436, 0.0749894, 0.0776247,
//        0.0803526, 0.0831764, 0.0860994, 0.0891251, 0.0922571, 0.0954993, 0.0988553, 0.102329, 0.105925,
//        0.109648, 0.113501, 0.11749, 0.121619, 0.125893, 0.130317, 0.134896, 0.139637, 0.144544, 0.149624,
//        0.154882, 0.160325, 0.165959, 0.171791, 0.177828, 0.184077, 0.190546, 0.197242, 0.204174, 0.211349,
//        0.218776, 0.226464, 0.234423, 0.242661, 0.251189, 0.260016, 0.269153, 0.278612, 0.288403, 0.298538,
//        0.30903, 0.31989, 0.331131, 0.342768, 0.354813, 0.367282, 0.380189, 0.39355, 0.40738, 0.421697,
//        0.436516, 0.451856, 0.467735, 0.484172, 0.501187, 0.5188, 0.537032, 0.555904, 0.57544, 0.595662,
//        0.616595, 0.638263, 0.660693, 0.683912, 0.707946, 0.732825, 0.758578, 0.785236, 0.812831, 0.841395,
//        0.870964, 0.901571, 0.933254, 0.966051, 1
//    };
//    
////    constexpr double Rebvalues_e3c[] = {
////        1e-06, 1.03514e-06, 1.07152e-06, 1.10917e-06, 1.14815e-06, 1.1885e-06, 1.23027e-06, 1.2735e-06,
////        1.31826e-06, 1.36458e-06, 1.41254e-06, 1.46218e-06, 1.51356e-06, 1.56675e-06, 1.62181e-06, 1.6788e-06,
////        1.7378e-06, 1.79887e-06, 1.86209e-06, 1.92752e-06, 1.99526e-06, 2.06538e-06, 2.13796e-06, 2.21309e-06,
////        2.29087e-06, 2.37137e-06, 2.45471e-06, 2.54097e-06, 2.63027e-06, 2.7227e-06, 2.81838e-06, 2.91743e-06,
////        3.01995e-06, 3.12608e-06, 3.23594e-06, 3.34965e-06, 3.46737e-06, 3.58922e-06, 3.71535e-06, 3.84592e-06,
////        3.98107e-06, 4.12098e-06, 4.2658e-06, 4.4157e-06, 4.57088e-06, 4.73151e-06, 4.89779e-06, 5.06991e-06,
////        5.24807e-06, 5.4325e-06, 5.62341e-06, 5.82103e-06, 6.0256e-06, 6.23735e-06, 6.45654e-06, 6.68344e-06,
////        6.91831e-06, 7.16143e-06, 7.4131e-06, 7.67361e-06, 7.94328e-06, 8.22243e-06, 8.51138e-06, 8.81049e-06,
////        9.12011e-06, 9.44061e-06, 9.77237e-06, 1.01158e-05, 1.04713e-05, 1.08393e-05, 1.12202e-05, 1.16145e-05,
////        1.20226e-05, 1.24451e-05, 1.28825e-05, 1.33352e-05, 1.38038e-05, 1.42889e-05, 1.47911e-05, 1.53109e-05,
////        1.58489e-05, 1.64059e-05, 1.69824e-05, 1.75792e-05, 1.8197e-05, 1.88365e-05, 1.94984e-05, 2.01837e-05,
////        2.0893e-05, 2.16272e-05, 2.23872e-05, 2.31739e-05, 2.39883e-05, 2.48313e-05, 2.5704e-05, 2.66073e-05,
////        2.75423e-05, 2.85102e-05, 2.95121e-05, 3.05492e-05, 3.16228e-05, 3.27341e-05, 3.38844e-05, 3.50752e-05,
////        3.63078e-05, 3.75837e-05, 3.89045e-05, 4.02717e-05, 4.16869e-05, 4.31519e-05, 4.46684e-05, 4.62381e-05,
////        4.7863e-05, 4.9545e-05, 5.12861e-05, 5.30884e-05, 5.49541e-05, 5.68853e-05, 5.88844e-05, 6.09537e-05,
////        6.30957e-05, 6.53131e-05, 6.76083e-05, 6.99842e-05, 7.24436e-05, 7.49894e-05, 7.76247e-05, 8.03526e-05,
////        8.31764e-05, 8.60994e-05, 8.91251e-05, 9.22571e-05, 9.54993e-05, 9.88553e-05, 0.000102329, 0.000105925,
////        0.000109648, 0.000113501, 0.00011749, 0.000121619, 0.000125893, 0.000130317, 0.000134896, 0.000139637,
////        0.000144544, 0.000149624, 0.000154882, 0.000160325, 0.000165959, 0.000171791, 0.000177828, 0.000184077,
////        0.000190546, 0.000197242, 0.000204174, 0.000211349, 0.000218776, 0.000226464, 0.000234423, 0.000242661,
////        0.000251189, 0.000260016, 0.000269153, 0.000278612, 0.000288403, 0.000298538, 0.00030903, 0.00031989,
////        0.000331131, 0.000342768, 0.000354813, 0.000367282, 0.000380189, 0.00039355, 0.00040738, 0.000421697,
////        0.000436516, 0.000451856, 0.000467735, 0.000484172, 0.000501187, 0.0005188, 0.000537032, 0.000555904,
////        0.00057544, 0.000595662, 0.000616595, 0.000638263, 0.000660693, 0.000683912, 0.000707946, 0.000732825,
////        0.000758578, 0.000785236, 0.000812831, 0.000841395, 0.000870964, 0.000901571, 0.000933254, 0.000966051,
////        0.001, 0.00103514, 0.00107152, 0.00110917, 0.00114815, 0.0011885, 0.00123027, 0.0012735, 0.00131826,
////        0.00136458, 0.00141254, 0.00146218, 0.00151356, 0.00156675, 0.00162181, 0.0016788, 0.0017378, 0.00179887,
////        0.00186209, 0.00192752, 0.00199526, 0.00206538, 0.00213796, 0.00221309, 0.00229087, 0.00237137, 0.00245471,
////        0.00254097, 0.00263027, 0.0027227, 0.00281838, 0.00291743, 0.00301995, 0.00312608, 0.00323594, 0.00334965,
////        0.00346737, 0.00358922, 0.00371535, 0.00384592, 0.00398107, 0.00412098, 0.0042658, 0.0044157, 0.00457088,
////        0.00473151, 0.00489779, 0.00506991, 0.00524807, 0.0054325, 0.00562341, 0.00582103, 0.0060256, 0.00623735,
////        0.00645654, 0.00668344, 0.00691831, 0.00716143, 0.0074131, 0.00767361, 0.00794328, 0.00822243, 0.00851138,
////        0.00881049, 0.00912011, 0.00944061, 0.00977237, 0.0101158, 0.0104713, 0.0108393, 0.0112202, 0.0116145,
////        0.0120226, 0.0124451, 0.0128825, 0.0133352, 0.0138038, 0.0142889, 0.0147911, 0.0153109, 0.0158489,
////        0.0164059, 0.0169824, 0.0175792, 0.018197, 0.0188365, 0.0194984, 0.0201837, 0.020893, 0.0216272, 0.0223872,
////        0.0231739, 0.0239883, 0.0248313, 0.025704, 0.0266073, 0.0275423, 0.0285102, 0.0295121, 0.0305492,
////        0.0316228, 0.0327341, 0.0338844, 0.0350752, 0.0363078, 0.0375837, 0.0389045, 0.0402717, 0.0416869,
////        0.0431519, 0.0446684, 0.0462381, 0.047863, 0.049545, 0.0512861, 0.0530884, 0.0549541, 0.0568853,
////        0.0588844, 0.0609537, 0.0630957, 0.0653131, 0.0676083, 0.0699842, 0.0724436, 0.0749894, 0.0776247,
////        0.0803526, 0.0831764, 0.0860994, 0.0891251, 0.0922571, 0.0954993, 0.0988553, 0.102329, 0.105925,
////        0.109648, 0.113501, 0.11749, 0.121619, 0.125893, 0.130317, 0.134896, 0.139637, 0.144544, 0.149624,
////        0.154882, 0.160325, 0.165959, 0.171791, 0.177828, 0.184077, 0.190546, 0.197242, 0.204174, 0.211349,
////        0.218776, 0.226464, 0.234423, 0.242661, 0.251189, 0.260016, 0.269153, 0.278612, 0.288403, 0.298538,
////        0.30903, 0.31989, 0.331131, 0.342768, 0.354813, 0.367282, 0.380189, 0.39355, 0.40738, 0.421697,
////        0.436516, 0.451856, 0.467735, 0.484172, 0.501187, 0.5188, 0.537032, 0.555904, 0.57544, 0.595662,
////        0.616595, 0.638263, 0.660693, 0.683912, 0.707946, 0.732825, 0.758578, 0.785236, 0.812831, 0.841395,
////        0.870964, 0.901571, 0.933254, 0.966051, 1
////    };
//    
////    constexpr double Rebvalues_e3c[] = {
////        0.000102329,0.000121619,0.000251189, 0.000501187,
////        0.00151356, 0.00301995, 0.00506991, 0.00716143, 0.00912011, 0.0101158, 0.0104713, 0.0108393, 0.0128825, 0.0133352, 0.0138038, 0.0142889, 0.0147911, 0.0153109, 0.0158489,
////        0.0164059, 0.0169824, 0.0175792, 0.018197, 0.0188365, 0.0194984, 0.0201837, 0.025704, 0.0305492,
////        0.0363078, 0.0402717, 0.0462381, 0.0512861, 0.0568853,0.0653131, 0.0724436, 0.0749894,
////        0.0803526, 0.0860994, 0.0954993, 0.0988553,
////        0.154882
////    };
////    
//    
//     constexpr double Rebvalues_e3c[] = {
//        0.000102329,0.000121619,0.000251189, 0.000501187,
//        0.00151356, 0.00301995, 0.00506991, 0.00716143, 0.00912011, 0.0101158, 0.0104713, 0.0108393, 0.0128825, 0.0133352, 0.0138038, 0.0142889, 0.0153109,
//        0.0164059, 0.0175792, 0.018197, 0.0201837, 0.025704, 0.0305492,
//        0.0363078, 0.0402717, 0.0462381, 0.0512861, 0.0568853,0.0653131, 0.0724436, 0.0749894,
//        0.0803526, 0.0860994, 0.0954993, 0.0988553,
//        0.154882
//    };
    
     constexpr double Rebvalues_e3c[] = {
        0.000102329,0.000121619,0.000251189, 0.000501187,
        0.00151356, 0.00301995, 0.00506991, 0.00716143, 0.00912011, 0.0101158, 0.0104713, 0.0108393, 0.0128825, 0.0133352, 0.0138038, 0.0142889, 0.0153109,
        0.0164059, 0.018197, 0.0201837, 0.025704, 0.0305492,
        0.0363078, 0.0402717, 0.0462381, 0.0512861, 0.0568853,0.0653131, 0.0724436,
        0.0803526, 0.0988553,
        0.154882
    };
//
////     constexpr double Rebvalues_e3c[] = {
////        0.000102329,0.000121619,0.000251189, 0.000501187,
////        0.00151356, 0.00301995, 0.00506991, 0.00716143, 0.00912011, 0.0101158, 0.0104713, 0.0108393, 0.0128825, 0.0133352, 0.0138038, 0.0142889, 0.0153109,
////        0.0164059, 0.0175792, 0.018197, 0.0201837, 0.025704, 0.0305492,
////        0.0363078, 0.0402717, 0.0462381, 0.0512861, 0.0568853,0.0653131, 0.0749894,
////       0.0860994, 0.0954993,
////        0.154882
////    };
////    
//      constexpr double Rebvalues_e3c[] = {
//        0.000102329,0.000121619,0.000251189, 0.000501187,
//        0.00151356, 0.00301995, 0.00506991, 0.00716143, 0.00912011, 0.0101158, 0.0104713, 0.0108393, 0.0128825, 0.0133352, 0.0138038, 0.0142889, 0.0153109,
//        0.0803526, 0.0860994, 0.0954993, 0.0988553,
//        0.154882
//    };
//    
//    
//    
    double numbins_x_e3c = sizeof(RvaluesReb)/sizeof(RvaluesReb[0]) - 1;
    double numbins_y_e3c = sizeof(Rebvalues_e3c)/sizeof(Rebvalues_e3c[0]) - 1;
    
    TH2D* h2_unfolded = new TH2D("h2_reb", "Rebinned 2D Histogram (Both Axes)", numbins_x_e3c, RvaluesReb, numbins_y_e3c, Rebvalues_e3c);
    
    for (int i = 1; i <= h2->GetNbinsX(); ++i) {
        for (int j = 1; j <= h2->GetNbinsY(); ++j) {
            double content = h2->GetBinContent(i, j);
            double error = h2->GetBinError(i, j);
            double x = h2->GetXaxis()->GetBinCenter(i);
            double y = h2->GetYaxis()->GetBinCenter(j);
            h2_unfolded->Fill(x, y, content);
        }
    }
    TH1D *h1_unfolded = (TH1D*)h2_unfolded->ProjectionX("h1_unfolded");
    h1_unfolded->Reset();
    
     TCanvas *canvas1 = new TCanvas();
    canvas1->SetLogx();
    canvas1->SetLogy();
    h2_unfolded->SetStats(0);
    h2_unfolded->GetYaxis()->SetTitle("wt");
    h2_unfolded->GetXaxis()->SetTitle(str4);
    h2_unfolded->Draw("text");
    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.04);
    latex1.SetTextFont(42);
    latex1.DrawLatex(0.15,0.80 ,Form("%.1d < %s < %.1d GeV/c, E3C", pt1, str, pt2));

///############################SET UP WHICH KIND OF PROJECTION YOU WANT TO DO------------------------------------------------------------

    bool ifInterpolate = true;
    
 
    TCanvas* c = new TCanvas("c1Dhist","1dhist",800,650);
    //        h3->Draw();
    
    if(!ifInterpolate)
    {
        //Loop over the x bins to create a weighted average
        for (int i = 1; i <= nBinsX; ++i) {
            int nBinsY = h2->GetNbinsY();
            double weightedSum = 0;
            double sumOfWeights = 0;
            
            for (int j = 1; j <= nBinsY; ++j) {
                double binContent = h2->GetBinContent(i, j);
                double binCenterY = h2->GetYaxis()->GetBinCenter(j); // Bin center in Y direction
                
                // Compute weight as the bin center in Y direction
                double weight = binCenterY;
                
                // Accumulate weighted sum and sum of weights
                weightedSum += binContent * weight;
                sumOfWeights += weight;
            }
            // Compute weighted average, avoiding division by zero
            double weightedAverage = (sumOfWeights != 0) ? (weightedSum / sumOfWeights) : 0;
            double xCenter = h2->GetXaxis()->GetBinCenter(i);
            //    h1->SetBinContent(i, weightedAverage);
            h1->SetBinContent(i, weightedSum);//this is closest to filling on the fly, not the avg
        }
        
        
        for (int i = 1; i <= h2_unfolded->GetNbinsX(); ++i) {
            int nBinsY = h2_unfolded->GetNbinsY();
            double weightedSum = 0;
            double sumOfWeights = 0;
            
            for (int j = 1; j <= nBinsY; ++j) {
                double binContent = h2_unfolded->GetBinContent(i, j);
                double binCenterY = h2_unfolded->GetYaxis()->GetBinCenter(j); // Bin center in Y direction
                
                // Compute weight as the bin center in Y direction
                double weight = binCenterY;
                
                // Accumulate weighted sum and sum of weights
                weightedSum += binContent * weight;
                sumOfWeights += weight;
            }
            // Compute weighted average, avoiding division by zero
            double weightedAverage = (sumOfWeights != 0) ? (weightedSum / sumOfWeights) : 0;
            double xCenter = h2_unfolded->GetXaxis()->GetBinCenter(i);
            h1_unfolded->SetBinContent(i, weightedSum);//this is closest to filling on the fly, not the avg
        }
    }
    else{
        cout<<"Interpolating"<<endl;
        for (int i = 1; i <= nBinsX; ++i){
            int nBinsY = h2->GetNbinsY();
            double weightedSum = 0;
            TH1D *h1_wt = (TH1D*)h2->ProjectionY(Form("h1_wt_%i",i),i,i);
            //               h1_wt->Draw();
            for (int j = 1; j <= nBinsY; ++j) {
                double weight = h1_wt->GetBinCenter(j);
                double bin_content = h1_wt->Interpolate(weight);
                weightedSum += weight*bin_content;
                //                if(bin_content!=0){
                //                    cout<<setprecision(7)<<h1_wt->GetBinContent(j)<<" and "<<setprecision(7)<<bin_content<<endl;
                //                }
                
            }
            h1->SetBinContent(i, weightedSum);//this is closest to filling on the fly, not the avg
        }
        
        for (int i = 1; i <= h2_unfolded->GetNbinsX(); ++i){
            int nBinsY = h2_unfolded->GetNbinsY();
            double weightedSum = 0;
            TH1D *h1_wt_unfolded = (TH1D*)h2_unfolded->ProjectionY(Form("h1_wt_unfolded_%i",i),i,i);
            for (int j = 1; j <= nBinsY; ++j) {
                double weight = h1_wt_unfolded->GetBinCenter(j);
                double bin_content = h1_wt_unfolded->Interpolate(weight);
                weightedSum += weight*bin_content;
                
                //                 if(bin_content!=h1_wt_unfolded->GetBinContent(j)){
                //                    cout<<setprecision(12)<<h1_wt_unfolded->GetBinContent(j)<<" and "<<setprecision(12)<<bin_content<<endl;
//                cout<<"Bin Content: j "<<j<<" "<<h1_wt_unfolded->GetBinContent(j-1)<<" and "<<h1_wt_unfolded->GetBinContent(j)<<" and "<<h1_wt_unfolded->GetBinContent(j+1)<<endl;
//                cout<<"Interpolated: j "<<j<<" "<<h1_wt_unfolded->Interpolate(h1_wt_unfolded->GetBinCenter(j-1))<<" and "<<h1_wt_unfolded->Interpolate(h1_wt_unfolded->GetBinCenter(j))<<" and "<<h1_wt_unfolded->Interpolate(h1_wt_unfolded->GetBinCenter(j+1))<<endl;
//                //                }
                
            }
            h1_unfolded->SetBinContent(i, weightedSum);//this is closest to filling on the fly, not the avg
        }
        
        //        for (int i = 1; i <= h2_unfolded->GetNbinsX(); ++i){
        //            int nBinsY = h2_unfolded->GetNbinsY();
        //            double weightedSum = 0;
        //            double weightedSum_unf = 0;
        //
        //            TH1D *h1_wt_unfolded = (TH1D*)h2_unfolded->ProjectionY(Form("h1_wt_unfolded_%i",i),i,i);
        //            TH1D *h1_wt = (TH1D*)h2->ProjectionY(Form("h1_wt_%i",i),i,i);
        //
        //            cout<<"h1_unf "<<h1_wt_unfolded->GetXaxis()->GetXmin()<<endl;
        //            cout<<"h2_unf "<<h2_unfolded->GetYaxis()->GetXmin()<<endl;
        //
        //            cout<<h2->GetYaxis()->GetXmin()<<endl;
        //            cout<<h1_wt->GetXaxis()->GetXmin()<<endl;
        //
        //            cout<<"h1 "<<h1_wt->GetXaxis()->GetXmax()<<endl;
        //            cout<<"h2 "<<h2->GetYaxis()->GetXmax()<<endl;
        //
        //            cout<<"h1_unf "<<h1_wt_unfolded->GetXaxis()->GetXmax()<<endl;
        //            cout<<"h2_unf "<<h2_unfolded->GetYaxis()->GetXmax()<<endl;
        //
        //
        ////            cout<<h2_unfolded->GetNbinsY()<<endl;
        ////            cout<<h2->GetNbinsY()<<endl;
        ////
        //            cout<<h1_wt_unfolded->GetNbinsX()<<endl;
        //            cout<<h1_wt->GetNbinsX()<<endl;
        //
        //
        //
        //            for (int j = 1; j <= nBinsY; ++j) {
        //
        ////                cout<<h2_unfolded->GetBinContent(i,j)<<" and "<<h2->GetBinContent(i,j)<<endl;
        ////                cout<<h1_wt_unfolded->GetBinContent(j)<<" and "<<h1_wt->GetBinContent(j)<<endl;
        //
        //                double weight = h1_wt->GetBinCenter(j);
        //                double bin_content = h1_wt->Interpolate(weight);
        //
        //                double weight_unf = h1_wt_unfolded->GetBinCenter(j);
        //                double bin_content_unf = h1_wt_unfolded->Interpolate(weight_unf);
        //
        //                weightedSum += weight*bin_content;
        //                weightedSum_unf += weight_unf*bin_content_unf;
        //
        ////                cout<<"weight "<<weight<<" and weight_unf "<<weight_unf<<endl;
        ////                cout<<"bin+content "<<bin_content<<" and bin+content_unf "<<bin_content_unf<<endl;
        ////
        ////
        //
        //            }
        //            h1_unfolded->SetBinContent(i, weightedSum_unf);//this is closest to filling on the fly, not the avg
        //            h1->SetBinContent(i, weightedSum);//this is closest to filling on the fly, not the avg
        //        }
        
        
    }
    
    h1->Scale(1.,"width");
    h1->SetLineColor(kGreen);
    
    h1_unfolded->Scale(1.,"width");
    h1_unfolded->SetLineColor(kRed);
    
    TH2D* h2_nom;
    h2_nom = (TH2D*)f->Get(histName2D.c_str());
    
    TH1D* h1_nom = (TH1D*)h2_nom->ProjectionX(Form("h1nom"),h2_nom->GetYaxis()->FindBin(pt1),h2_nom->GetYaxis()->FindBin(pt2));
    h1_nom->Scale(1.,"width");
    //        h1->Divide(h1_nom);
//    h1_unfolded->Divide(h1_nom);
          
//    h1_unfolded->Divide(h1);
//     h1->Divide(h1_nom);

    h1->Draw();
    h1_nom->Draw("same");
    h1_unfolded->Draw("same,HIST");
//
//   h2_unfolded->Divide(h2);
//   h2_unfolded->Draw();
//    h2->Draw();
//    cout<<h1_nom->GetNbinsX()<<" and "<<nBinsX<<endl;
//    cout<<h3->GetNbinsZ()<<" and "<<h2_unfolded->GetNbinsY()<<endl;
//    cout<<h3->GetNbinsZ()<<" and "<<h3_unfolded->GetNbinsZ()<<endl;
    
    
     cout<<h2->GetNbinsZ()<<" and "<<h2_unfolded->GetNbinsZ()<<endl;
     cout<<h2->GetNbinsY()<<" and "<<h2_unfolded->GetNbinsY()<<endl;
     cout<<h2->GetNbinsX()<<" and "<<h2_unfolded->GetNbinsX()<<endl;
//    
//     cout<<h3->GetNbinsZ()<<" and "<<h3_unfolded->GetNbinsZ()<<endl;
//     cout<<h3->GetNbinsY()<<" and "<<h3_unfolded->GetNbinsY()<<endl;
//     cout<<h3->GetNbinsX()<<" and "<<h3_unfolded->GetNbinsX()<<endl;
    
}
