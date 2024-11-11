#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/Vectors.hh"
#include "Rivet/Math/Matrices.hh"
#include "Rivet/Math/Units.hh"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"


#include "cmath"
#include "iostream"
#include "fstream"
#include "random"
#include "vector"
#include "stdio.h"
#include "stdlib.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"

using namespace std;

/// @Energy correlators in jets at ALICE 13TeV fast simulation
namespace Rivet {

class ALICE_fastSim : public Analysis {
public:
    
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_ENC_fastSim);
    
    /// @name Analysis methods
    ///@{
    
    /// Book histograms and initialise projections before the run
    void init() {
        // Initialise and register projections
        
        // The basic final-state projection:
        // all charged final-state particles within
        // the given eta acceptance and pTmin
        const ChargedFinalState cfs;
        declare(cfs, "CFS");
        
        //particle tracking files to read
        tr_eff_f = new TFile("/home/ar2545/palmer_scratch/tr_eff.root", "READ");
        track_eff_hist = (TH1D*)tr_eff_f->Get("tr_eff");
        tr_eff_gr = new TGraph(track_eff_hist);
        
        //Defining functions to use
        double lg_pt_eff(double p, double m, double b);
        int GetClosestXIndex(const TGraph* graph, double xValue);
        
        //Booking histograms
        book (_h_track_pt_tru, "true_track_pt", 150, 0, 150);
        book (_h_track_pt_tru_norm, "true_track_pt_norm", 150, 0, 150);
        book (_h_track_pt_tru_scale, "true_track_pt_scale", 150, 0, 150);
        
        book (_h_track_pt_det, "det_track_pt", 150, 0, 150);
        book (_h_track_pt_det_norm, "det_track_pt_norm", 150, 0, 150);
        book (_h_track_pt_det_scale, "det_track_pt_scale", 150, 0, 150);
        
        // book (_h_jet_pt_tru, "jet_pt_tru", 21, 15, 120);
        // book(_h_jet_pt_tru_norm,"jet_pt_tru_norm",21, 15, 120);
        // book( _h_jet_pt_tru_scale,"jet_pt_tru_scale",21, 15, 120);
        
        // book (_h_jet_pt_det, "jet_pt_det", 21, 15, 120);
        // book(_h_jet_pt_det_norm,"jet_pt_det_norm",21, 15, 120);
        // book( _h_jet_pt_det_scale,"jet_pt_det_scale",21, 15, 120);
        
        book (_h_jet_pt_tru, "jet_pt_tru", 25, 15, 150);
        book(_h_jet_pt_tru_norm,"jet_pt_tru_norm",25, 15, 150);
        book( _h_jet_pt_tru_scale,"jet_pt_tru_scale",25, 15, 150);
        
        book (_h_jet_pt_det, "jet_pt_det", 25, 15, 150);
        book(_h_jet_pt_det_norm,"jet_pt_det_norm",25, 15, 150);
        book( _h_jet_pt_det_scale,"jet_pt_det_scale",25, 15, 150);
        
        vector<double> xbins;
        
        Double_t from = -2;
        Double_t to = -0.4;
        Double_t width = (to-from)/20;
        for (size_t i = 0; i < 20 + 1; i++) xbins.push_back(1*pow(10, from + i*width));
        for (size_t i = 0; i < 20 + 1; i++) cout<<xbins[i]<<endl;
        //xbins[] = [0.01,0.0120226,0.0144544,0.017378,0.020893,0.0251189,0.0301995,0.0363078,0.0436516,0.0524807,0.0630957,0.0758578,0.0912011,0.109648,0.131826,0.158489,0.190546,0.229087,0.275423,0.331131,0.398107];
        
        
        //for converting to full jet pT in 20-80 range, extend jet pt range to 140 linspace(25,15,150). Usual range is linspace(21,15,120)
        book(_h_eec_pt_tru,"eec_pt_hist_tru",xbins,linspace(21,15,120));
        book(_h_eec_pt_tru_norm,"eec_pt_hist_tru_norm",xbins,linspace(21,15,120));
        book(_h_eec_pt_tru_scale,"eec_pt_hist_tru_scale",xbins,linspace(21,15,120));
        
        book(_h_eec_pt_det,"eec_pt_hist_det",xbins,linspace(21,15,120));
        book(_h_eec_pt_det_norm,"eec_pt_hist_det_norm",xbins,linspace(21,15,120));
        book(_h_eec_pt_det_scale,"eec_pt_hist_det_scale",xbins,linspace(21,15,120));
        
        book(_h_e3c_pt_tru,"e3c_pt_hist_tru",xbins,linspace(21,15,120));
        book(_h_e3c_pt_tru_norm,"e3c_pt_hist_tru_norm",xbins,linspace(21,15,120));
        book(_h_e3c_pt_tru_scale,"e3c_pt_hist_tru_scale",xbins,linspace(21,15,120));
        
        book(_h_e3c_pt_det,"e3c_pt_hist_det",xbins,linspace(21,15,120));
        book(_h_e3c_pt_det_norm,"e3c_pt_hist_det_norm",xbins,linspace(21,15,120));
        book(_h_e3c_pt_det_scale,"e3c_pt_hist_det_scale",xbins,linspace(21,15,120));
        
        //These are to get the pair efficiency of particles which we will weight the particles by when filling the ENC
        book(_h_n2_qpt_det,"npairs_qpt_hist_det",xbins,linspace(20,0,2));
        book(_h_n2_qpt,"npairs_qpt_hist",xbins,linspace(20,0,2));
        book(_h_n3_qpt_det,"npairs3_qpt_hist_det",xbins,linspace(20,0,2));
        book(_h_n3_qpt,"npairs3_qpt_hist",xbins,linspace(20,0,2));
        book(_h_n2_paireff,"n2_paireff",xbins,linspace(20,0,2));
        
    }
    
    double lg_pt_eff(double pt, double m, double b)
    {
        double lgf = m*pt + b;
        if (lgf > 0.0){return lgf;}
        else {return 0;}
    }
    
    Double_t delR(fastjet::PseudoJet ps1, fastjet::PseudoJet ps2)
    {
        double dphi = abs(ps1.phi()-ps2.phi());
        if(dphi>TMath::Pi()){dphi = (2.*TMath::Pi() - dphi);}
        
        double deta = std::fabs(ps1.eta() - ps2.eta());
        Double_t delR = std::sqrt(dphi*dphi + deta*deta);
        
        return delR;
    }
    
    int GetClosestXIndex(const TGraph* graph, double xValue)
    {
        double minDist = DBL_MAX;
        int closestIndex = -1;
        
        // Iterate over all points in the graph
        for (int i = 0; i < graph->GetN(); ++i)
        {
            double x, y;
            graph->GetPoint(i, x, y);
            
            // Calculate the distance between the current x-value and the desired x-value
            double dist = TMath::Abs(x - xValue);
            // Update the closest index if a closer point is found
            if (dist < minDist) {
                minDist = dist;
                closestIndex = i;
            }
        }
        
        return closestIndex;
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& event)
    {
        
        //Defining some functions to use
        //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        // Better linear extrapolation than the ROOT TGraph provides for pT > 150 GeV/c
        double x0 = 36; //tgraph has 82 points counting from 0. 36th point is where for you hit 10 (9.75) GeV
        double y0 = tr_eff_gr->GetPointY(x0);
        double x1 = 81; //tgraph has 82 points counting from 0. 81st point is where for you hit 150 (147.5) GeV
        double y1 = tr_eff_gr->GetPointY(x1);
        double  m = (y1 - y0) / (x1 - x0);
        double  b = y1 - m*x1;
        
        //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        Particles fsParticles = applyProjection<FinalState>(event, "CFS").particles(Cuts::abseta < 0.9 && Cuts::pt > 0.150);//Apply ALICE tracking cut
        Particles fsParticles_det = applyProjection<FinalState>(event, "CFS").particles(Cuts::abseta < 0.9 && Cuts::pt > 0.150); //apply eta cut and ALICE kinematic cut on det level particles
        PseudoJets parts, parts_det;
        
        random_device generator;
        double sigma_pt, pts, ratio;
        
        for (const Particle &p : fsParticles)
        {
            PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
            parts.push_back(pseudojet); //Truth level particles
            _h_track_pt_tru->fill(p.pt());
            _h_track_pt_tru_norm->fill(p.pt());
            _h_track_pt_tru_scale->fill(p.pt());
        }
        
        for (const Particle &p : fsParticles_det) //Det level particles
        {
            //Apply particle eff cut
            int closestIndex = GetClosestXIndex(tr_eff_gr, p.pt()); //find the index of the point with the pt value
            double lg_eff = lg_pt_eff(p.pt(),m,b);
            r = ((double) rand() / (RAND_MAX));
            
            if (p.pt()>150 &&   r > lg_eff) continue;
            else if (r > tr_eff_gr->GetPointY(closestIndex))  continue;
            else
            {
                //Apply pT smearing
                if (p.pt() < 1)
                {
                    sigma_pt = p.pt() * (-0.035 * p.pt() + 0.04);
                }
                else if (p.pt() < 60)
                {
                    sigma_pt = p.pt() * (0.00085 * p.pt() + 0.00415);
                }
                else // approximately valid for at least pt < 90 GeV
                {
                    sigma_pt = p.pt() * ((0.0015 * p.pt()) - 0.035);
                }
                normal_distribution<double> smearTrackPt(p.pt(),sigma_pt);
                pts = smearTrackPt(generator);
                ratio = pts/p.pt();
                PseudoJet p_smr_tr(p.px()*ratio,p.py()*ratio,p.pz()*ratio,p.E()*ratio);
                p_smr_tr.reset_momentum(p_smr_tr.px(), p_smr_tr.py(), p_smr_tr.pz(), sqrt(pow(p_smr_tr.px(), 2) + pow(p_smr_tr.py(), 2) + pow(p_smr_tr.pz(), 2) + pow(p.mass(), 2)));
                parts_det.push_back(p_smr_tr); //Det level particles
                _h_track_pt_det->fill(p_smr_tr.pt());
                _h_track_pt_det_norm->fill(p_smr_tr.pt());
                _h_track_pt_det_scale->fill(p_smr_tr.pt());
            }
        }
        //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        //Implementing pair efficiency: filling histograms with pairs for the EEC first as a function of charge/pt of tracks.
        double dist, dist_det, qpt, qpt_det;
        for (unsigned int i = 0; i < fsParticles; i++)
        {
            for (unsigned int j = i+1; j < fsParticles; j++)
            {
                
                dist = delR(fsParticles[i],fsParticles[j]);
                qpt = abs(fsParticles[i].charge()/fsParticles[i].pt() - fsParticles[j].charge()/fsParticles[j].pt());
                _h_n2_qpt->fill(dist,qpt);
                if(i<=fsParticles_det.size() && j<=fsParticles_det.size())
                {
                    dist_det = delR(fsParticles_det[i],fsParticles_det[j]);
                    qpt_det = abs(fsParticles_det[i].charge()/fsParticles_det[i].pt() - fsParticles_det[j].charge()/fsParticles_det[j].pt());
                    _h_n2_qpt_det->fill(dist_det,qpt_det);
                }
                //                    double ratio = dist_det/dist;
                //                    _h_n2_paireff->fill(ratio, qpt_det, qpt);
                
                
            }
        }
        
        for (unsigned int i = 0; i < fsParticles_det; i++)
        {
            for (unsigned int j = i+1; j < fsParticles_det; j++)
            {
                
                double dist_det = delR(fsParticles_det[i],fsParticles_det[j]);
                double qpt_det = abs(fsParticles_det[i].charge()/fsParticles_det[i].pt() - fsParticles_det[j].charge()/fsParticles_det[j].pt());
                _h_n2_qpt_det->fill(dist_det,qpt_det);
                
            }
        }
        
        _h_n2_paireff = (TH2D*)_h_n2_qpt_det->Clone("n2_paireff");
        _h_n2_paireff->Divide(_h_n2_qpt);
        
        
        //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        //Jet selection settings
        fastjet::Selector selector = fastjet::SelectorPtMin(15.0) * fastjet::SelectorEtaRange(-0.5, 0.5);
        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
        
        //Truth level jets
        fastjet::ClusterSequence cs(parts, jet_def);
        vector<fastjet::PseudoJet> jets = sorted_by_pt(selector(cs.inclusive_jets()));
        
        //Det level jets
        fastjet::ClusterSequence csdet(parts_det, jet_def);
        vector<fastjet::PseudoJet> jets_det = sorted_by_pt(selector(csdet.inclusive_jets()));
        
        
        
        //Jet Matching: Geometric matching
        
        //FOR MATCHED JETS EEC
        vector<fastjet::PseudoJet> constit, constit_det;
        for (unsigned int i = 0; i < jets.size(); i++)
        {
            int i_best_match = -1;
            int k_best_match = -1;
            double minDeltaR = 0.4;
            bool found_match = false;
            
            for (unsigned int k = 0; k < jets_det.size(); k++)
            {
                double dR = jet.delta_R(jet_MC);
                if (dR < minDeltaR) {
                    minDeltaR = dR;
                    i_best_match = i;
                    k_best_match = k;
                    found_match = true;
                }
                if(found_match){
                    double pTjet = jets[i_best_match].pt();
                    // double pTjet = (1.65)*jets[i].pt();
                    _h_jet_pt_tru->fill(pTjet);
                    _h_jet_pt_tru_norm->fill(pTjet);
                    _h_jet_pt_tru_scale->fill(pTjet);
                    for (unsigned int j = 0; j < jets[i_best_match].constituents().size(); j++)
                    {
                        // cout<<"entered j in true loop"<<endl;
                        for (unsigned int s = j+1; s < jets[i_best_match].constituents().size() ; s++)
                        {
                            // cout<<"entered s in true loop"<<endl;
                            PseudoJet constit1 = jets[i_best_match].constituents().at(j);
                            PseudoJet constit2 = jets[i_best_match].constituents().at(s);
                            if (constit1.pt() < 1 ||  constit2.pt() < 1) continue;
                            
                            double ee = 2*(constit1.pt()*constit2.pt())/(pTjet*pTjet);
                            double dR = delR(constit1,constit2);
                            
                            double ejjs = (3*constit1.pt()*constit1.pt()*constit2.pt())/(pow(pTjet,3));
                            double ejss = (3*constit1.pt()*constit2.pt()*constit2.pt())/(pow(pTjet,3));
                            
                            _h_eec_pt_tru->fill(dR,pTjet,ee);
                            _h_eec_pt_tru_norm->fill(dR,pTjet,ee);
                            _h_eec_pt_tru_scale->fill(dR,pTjet,ee);
                            
                            _h_e3c_pt_tru->fill(dR,pTjet,ejss);
                            _h_e3c_pt_tru_norm->fill(dR,pTjet,ejss);
                            _h_e3c_pt_tru_scale->fill(dR,pTjet,ejss);
                            
                            _h_e3c_pt_tru->fill(dR,pTjet,ejjs);
                            _h_e3c_pt_tru_norm->fill(dR,pTjet,ejjs);
                            _h_e3c_pt_tru_scale->fill(dR,pTjet,ejjs);
                            
                            // cout<<"filled"<<endl;
                            
                            for (unsigned int m = s+1; m < jets[i_best_match].constituents().size() ; m++)
                            {
                                PseudoJet constit3 = jets[i_best_match].constituents().at(m);
                                if (constit1.pt() < 1 ||  constit2.pt() < 1 || constit3.pt() <1) continue;
                                
                                double e3c_jsm = (6*constit1.pt()*constit2.pt()*constit3.pt())/(pow(pTjet,3));
                                double dRjs = delR(constit1,constit2);
                                double dRsm = delR(constit2,constit3);
                                double dRmj = delR(constit3,constit1);
                                
                                if(dRjs > dRsm && dRjs > dRmj){R_L = dRjs; }
                                else if(dRsm > dRjs && dRsm > dRmj){R_L = dRsm;}
                                else{R_L = dRmj;}
                                
                                
                                _h_e3c_pt_tru->fill(R_L,pTjet,e3c_jsm);
                                _h_e3c_pt_tru_norm->fill(R_L,pTjet,e3c_jsm);
                                _h_e3c_pt_tru_scale->fill(R_L,pTjet,e3c_jsm);
                            }
                        }
                    }
                    
                    double pTjet_det = jets_det[k_best_match].pt();
                    // double pTjet_det = (1.65)*jets_det[k].pt();
                    _h_jet_pt_det->fill(pTjet_det);
                    _h_jet_pt_det_norm->fill(pTjet);
                    _h_jet_pt_det_scale->fill(pTjet);
                    for (unsigned int j = 0; j < jets_det[k_best_match].constituents().size(); j++)
                    {
                        // cout<<"entered j in det loop"<<endl;
                        for (unsigned int s = j+1; s < jets_det[k_best_match].constituents().size(); s++)
                        {
                            //  cout<<"entered s in det loop"<<endl;
                            PseudoJet constit1 = jets_det[k_best_match].constituents().at(j);
                            PseudoJet constit2 = jets_det[k_best_match].constituents().at(s);
                            if (constit1.pt() < 1 ||  constit2.pt() < 1) continue;
                            double ee = 2*(constit1.pt()*constit2.pt())/(pTjet_det*pTjet_det);
                            double dR = delR(constit1,constit2);
                            
                            double ejjs = (3*constit1.pt()*constit1.pt()*constit2.pt())/(pow(pTjet_det,3));
                            double ejss = (3*constit1.pt()*constit2.pt()*constit2.pt())/(pow(pTjet_det,3));
                            
                            
                            _h_eec_pt_det->fill(dR,pTjet_det,ee);
                            _h_eec_pt_det_norm->fill(dR,pTjet_det,ee);
                            _h_eec_pt_det_scale->fill(dR,pTjet_det,ee);
                            
                            _h_e3c_pt_det->fill(dR,pTjet_det,ejss);
                            _h_e3c_pt_det_norm->fill(dR,pTjet_det,ejss);
                            _h_e3c_pt_det_scale->fill(dR,pTjet_det,ejss);
                            
                            _h_e3c_pt_det->fill(dR,pTjet_det,ejjs);
                            _h_e3c_pt_det_norm->fill(dR,pTjet_det,ejjs);
                            _h_e3c_pt_det_scale->fill(dR,pTjet_det,ejjs);
                            //Figuring out the charge/pt difference between two particles
                            double chargept = abs(constit1.charge()/constit1.pt() - constit2.charge()/constit2.pt());
                            //Now we want to find which bin does this belong to in our pair eff histogram so we can extract the weight that we will put in
                            //                                if bin1<chargept<bin2
                            double biny = _h_n2_paireff->GetYaxis()->FindBin(chargept);
                            hpaireff_x = _h_n2_paireff->ProjectionX(biny,biny+1); // check projection
                            
//                            hpaireff_x = _h_n2_paireff->ProjectionX(biny,biny+1); // check projection
                            
                            double binx = h->GetXaxis()->FindBin(dR); //should be dR
                            double weight_paireff = hpaireff_x->GetBinContent(binx);
                            _h_eec_pt_det_new->fill(dR,pTjet_det,ee*weight_paireff);
                            
                            for (unsigned int m = s+1; m < jets_det[k_best_match].constituents().size() ; m++)
                            {
                                PseudoJet constit3 = jets_det[k_best_match].constituents().at(m);
                                if (constit1.pt() < 1 ||  constit2.pt() < 1 || constit3.pt() <1) continue;
                                
                                double e3c_jsm = (6*constit1.pt()*constit2.pt()*constit3.pt())/(pow(pTjet,3));
                                double dRjs = delR(constit1,constit2);
                                double dRsm = delR(constit2,constit3);
                                double dRmj = delR(constit3,constit1);//Found that this was 32 on aug 30th instead of 31 and changed it
                                
                                if(dRjs > dRsm && dRjs > dRmj){R_L = dRjs; }
                                else if(dRsm > dRjs && dRsm > dRmj){R_L = dRsm;}
                                else{R_L = dRmj;}
                                
                                
                                _h_e3c_pt_det->fill(R_L,pTjet,e3c_jsm);
                                _h_e3c_pt_det_norm->fill(R_L,pTjet,e3c_jsm);
                                _h_e3c_pt_det_scale->fill(R_L,pTjet,e3c_jsm);
                            }
                        }
                    }
                }
            }
        }
        
        
        //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        
    }
    
    
    void finalize()
    {
        normalize(_h_track_pt_tru_norm, crossSection());
        normalize(_h_track_pt_det_norm, crossSection());
        
        scale(_h_track_pt_tru_scale, crossSection()/sumW());
        scale(_h_track_pt_det_scale, crossSection()/sumW());
        
        normalize(_h_jet_pt_tru_norm, crossSection());
        normalize(_h_jet_pt_det_norm, crossSection());
        
        scale(_h_jet_pt_tru_scale, crossSection()/sumW());
        scale(_h_jet_pt_det_scale, crossSection()/sumW());
        
        normalize(_h_eec_pt_tru_norm, crossSection());
        normalize(_h_eec_pt_det_norm, crossSection());
        
        scale(_h_eec_pt_tru_scale, crossSection()/sumW());
        scale(_h_eec_pt_det_scale, crossSection()/sumW());
        
        
        normalize(_h_e3c_pt_tru_norm, crossSection());
        normalize(_h_e3c_pt_det_norm, crossSection());
        
        scale(_h_e3c_pt_tru_scale, crossSection()/sumW());
        scale(_h_e3c_pt_det_scale, crossSection()/sumW());
        
        cout<< "success"<<endl;
        
        
    }
    
    //When using yoda
    Histo1DPtr _h_jet_pt_tru;
    Histo1DPtr _h_jet_pt_tru_norm;
    Histo1DPtr _h_jet_pt_tru_scale;
    
    Histo1DPtr _h_jet_pt_det;
    Histo1DPtr _h_jet_pt_det_norm;
    Histo1DPtr _h_jet_pt_det_scale;
    
    Histo1DPtr _h_track_pt_tru;
    Histo1DPtr _h_track_pt_tru_norm;
    Histo1DPtr _h_track_pt_tru_scale;
    
    Histo1DPtr _h_track_pt_det;
    Histo1DPtr _h_track_pt_det_norm;
    Histo1DPtr _h_track_pt_det_scale;
    
    Histo2DPtr _h_eec_pt_tru, _h_e3c_pt_tru;
    Histo2DPtr _h_eec_pt_tru_norm, _h_e3c_pt_tru_norm;
    Histo2DPtr _h_eec_pt_tru_scale, _h_e3c_pt_tru_scale;
    
    Histo2DPtr _h_eec_pt_det, _h_e3c_pt_det;
    Histo2DPtr _h_eec_pt_det_norm, _h_e3c_pt_det_norm;
    Histo2DPtr _h_eec_pt_det_scale, _h_e3c_pt_det_scale;
    
    Histo2DPtr _h_n2_qpt, _h_n2_qpt_det, _h_n2_paireff;
    double r, b, m;
    // double lg_pt_eff;
    // double GetClosestXIndex;
    
private:
    TH1D* track_eff_hist;
    TFile* tr_eff_f;
    TGraph* tr_eff_gr;
    
};

RIVET_DECLARE_PLUGIN(ALICE_fastSim);
}


