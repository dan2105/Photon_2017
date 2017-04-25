  
#include "TInterpreter.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TMath.h" 
#include "TFrame.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include <math.h>
#include <algorithm>
#include <list> 
#include "/home/daniel/workspace/fastjet_build/include/fastjet/ClusterSequence.hh"
#include "/home/daniel/workspace/fastjet_build/include/fastjet/PseudoJet.hh"
#include "/home/daniel/workspace/fastjet_build/include/fastjet/JetDefinition.hh"
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector> 
#include <cstdio>
#include "TRandom.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include <iostream>
#include <TROOT.h>
#include "TStopwatch.h"

using namespace std;
using namespace fastjet;




void analysis_aahbbar_signal()
{


  TFile *_file0 = new TFile("pPb_aahbB_signal_hepmc.root");  //Arquivo formato root

  TTree *T = (TTree *)_file0->Get("T");   //tree

  //=============================================================================
  //Declaração das variaveis dos eletron, muons, neutrinos e particulas carregadas
  //=============================================================================

  TStopwatch m_timer;
  double rtime;
  double ctime;

  int const NPARTMAX=10000;
  int n_part;
  int part_id[NPARTMAX];
  double part_pt[NPARTMAX];
  double part_px[NPARTMAX];
  double part_py[NPARTMAX];
  double part_pz[NPARTMAX];
  double part_eta[NPARTMAX];
  double part_phi[NPARTMAX];
  double part_energy[NPARTMAX];
  double part_mass[NPARTMAX];
  double part_charge[NPARTMAX];
  
  int const PROTONMAX=50;
  int n_proton=0;
  double proton_pt[PROTONMAX];
  double proton_px[PROTONMAX];
  double proton_py[PROTONMAX];
  double proton_pz[PROTONMAX];
  double proton_energy[PROTONMAX];



  //===============================================================
  //Quark's Branches 
  //===============================================================

  T->SetBranchAddress("n_part", &n_part);
  T->SetBranchAddress("part_id", &part_id);
  T->SetBranchAddress("part_pt", &part_pt);
  T->SetBranchAddress("part_px", &part_px);
  T->SetBranchAddress("part_py", &part_py); 
  T->SetBranchAddress("part_pz", &part_pz);
  T->SetBranchAddress("part_eta", &part_eta);
  T->SetBranchAddress("part_phi", &part_phi);
  T->SetBranchAddress("part_energy", &part_energy);
  T->SetBranchAddress("part_mass", &part_mass);
  T->SetBranchAddress("part_charge", &part_charge);

  //===============================================================
  // Protons's Branch
  //===============================================================
  T->SetBranchAddress("n_proton",&n_proton);
  T->SetBranchAddress("proton_px",&proton_px);
  T->SetBranchAddress("proton_py",&proton_py);
  T->SetBranchAddress("proton_pz",&proton_pz);
  T->SetBranchAddress("proton_pt",&proton_pt);
  T->SetBranchAddress("proton_energy",&proton_energy);


  // Declarar os histogramas no ROOT 
  TH1F* h_npart = new TH1F("npart","npart",1000,0.,1000.);
  TH1F* h_part_id = new TH1F("part_id","part_id",6000,-3000.,3000.);
  TH1F* h_part_pt = new TH1F("part_pt","part_pt",100,0.,10.);
  TH1F* h_part_px = new TH1F("part_px","part_px",100,0.,10.);
  TH1F* h_part_py = new TH1F("part_py","part_py",100,0.,10.);
  TH1F* h_part_pz = new TH1F("part_pz","part_pz",100,0.,10.);
  TH1F* h_part_eta = new TH1F("part_eta","part_eta",1000,-5.0,5.0);
  TH1F* h_part_phi = new TH1F("part_phi","part_phi",1000,-M_PI,M_PI);
  TH1F* h_part_cos = new TH1F("part_cos","part_cos",1000,-M_PI,M_PI);
  TH1F* h_part_energy = new TH1F("part_energy","part_energy",1000,0.,500.);
  TH1F* h_part_mass = new TH1F("part_mass","part_mass",100,0.,100.);
  TH1F* h_part_charge = new TH1F("part_charge","part_charge",10,-5.,5.);


  TH1F *h_jet1pt = new TH1F("jet1pt","jet1pt",400,0.,400);
  TH1F *h_jet2pt = new TH1F("jet2pt","jet2pt",400,0.,400);
  TH1F *h_invmassjet = new TH1F("invmassjet","invmassjet",2000,0.,2000);
  TH1F *h_proton1_pz = new TH1F("proton1_pz","proton1_pz",100000,-0.,10);


  //===============================================================
  //Variables contructed with the combinations of protons and quarks
  //===============================================================
  TH1F* h_proton1_xi  = new TH1F("proton1_xi","proton1_xi",1000,0.,1.);
  TH1F* h_proton2_xi  = new TH1F("proton2_xi","proton2_xi",1000,0.,1.);
  TH1F* h_mx_PbPb  = new TH1F("mx_PbPb","mx_PbPb",1000,0.,1000.);
  //===============================================================
  //t, pt information
  //===============================================================
  TH1F* h_proton1_t  = new TH1F("proton1_t","proton1_t",1000,0.,10.);
  TH1F* h_t_pt_proton1  = new TH1F("t_pt_proton1","t_pt_proton1",1000,0.,10.);
  TH1F* h_proton2_t  = new TH1F("proton2_t","proton2_t",1000,0.,10.);
  TH1F* h_t_pt_proton2  = new TH1F("t_pt_proton2","t_pt_proton2",1000,0.,10.);
  //===============================================================
  //Ratio DMF/ECM
  //===============================================================
  TH1F* h_ratio_DMF = new TH1F("ratio_DMF","ratio_DMF",500,0.,2.);
  //===============================================================

  //=================================================================
  TLorentzVector vec_part; //general Lorentz vector
  TLorentzVector proton2p4,  proton1p4, protonp4; //proton's Lorentz vector
  //==================================================================


  //==================================================================
  //FastJet
  //==================================================================
  double Rparam = 0.4;
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam, recombScheme, strategy);
  std::vector <fastjet::PseudoJet> fjInputs;
  //======================================================================


  //======================================================================
  //Tree's Loop 
  //======================================================================
  Int_t nentries = (Int_t)T->GetEntries();
  for (Int_t i=0; i<nentries; i++) {
    T->GetEntry(i);
    h_npart->Fill(n_part);   
    fjInputs.resize(0);
    //======================================================================
    
    bool found_proton1 = false;
    bool found_proton2 = false;

    for( int ipart=0; ipart < n_part; ipart++){
      vec_part.SetPxPyPzE(part_px[ipart],part_py[ipart],part_pz[ipart],part_energy[ipart]);
      
      // fjInputs.push_back(PseudoJet(part_px[ipart],part_py[ipart],part_pz[ipart],part_energy[ipart]));
          
      h_part_id->Fill(part_id[ipart]);
      h_part_pt->Fill(vec_part.Pt());
      h_part_px->Fill(vec_part.Px());
      h_part_py->Fill(vec_part.Py());
      h_part_pz->Fill(vec_part.Pz());
      h_part_energy->Fill(vec_part.E());
      h_part_eta->Fill(vec_part.Eta());
      h_part_charge->Fill(part_charge[ipart]);
      h_part_phi->Fill(vec_part.Phi());
      h_part_cos->Fill(vec_part.CosTheta());

      /*
      if(part_id[ipart]==2212 && part_pz[ipart] > 17000){
	//        double E_from_p = sqrt(pow(px[ipart],2) + pow(py[ipart],2) + pow(pz[ipart],2) + pow(0.938272046,2));
	proton1p4.SetPxPyPzE(part_px[ipart],part_py[ipart],part_pz[ipart],part_energy[ipart]);
	found_proton1 = true;
      } //+
      if(part_id[ipart]==2212 && part_pz[ipart] < -17000){
	proton2p4.SetPxPyPzE(part_px[ipart],part_py[ipart],part_pz[ipart],part_energy[ipart]);
	found_proton2 = true;
      } //-
      */     
     fjInputs.push_back(PseudoJet(part_px[ipart],part_py[ipart],part_pz[ipart],part_energy[ipart]));
    }


    //=================================================================  //pPb
    double ib1=31500;
    double ib2=31500;
    double xi_min_sel = 0.00000015;
    double xi_max_sel = 1.;	
    double ROOTs= 63000; 
    //=================================================================   


    /*
    //=================================================  //PbPb
    double ib1=19500;
    double ib2=19500;
    double xi_min_sel = 0.00000015;
    double xi_max_sel = 1.;	
    double ROOTs= 39000; 
    //=====================================================================
    */

    for(int iproton = 0; iproton < n_proton; iproton++){

      double proton1_px = proton_px[iproton];
      double proton1_py = proton_py[iproton];
      double proton1_pz = proton_pz[iproton];
      double proton1_pt = proton_pt[iproton];
      double proton1_energy = proton_energy[iproton];
      protonp4.SetPxPyPzE(proton1_px,proton1_py,proton1_pz,proton1_energy);
      
      if(protonp4.Pz()>0){
	double xi_tmp = 1. - fabs(protonp4.Pz())/ib1;
	if( (xi_tmp >= xi_min_sel) && (xi_tmp <= xi_max_sel) ){
	  proton1p4.SetPxPyPzE(proton1_px,proton1_py,proton1_pz,proton1_energy);
	  found_proton1=true;
	}
      }
      if(protonp4.Pz()<0){
	double xi_tmp =  1. - fabs(protonp4.Pz())/ib2;
	if( (xi_tmp >= xi_min_sel) && (xi_tmp <= xi_max_sel) ){
	  proton2p4.SetPxPyPzE(proton1_px,proton1_py,proton1_pz,proton1_energy);
	  found_proton2=true;
	}
      }// 4 vector Pz<0
    }// end proton loop  



   
    double xi1_proton_pb = 1. - fabs(proton1p4.Pz())/ib1;
    double xi2_proton_pb = 1. - fabs(proton2p4.Pz())/ib2;

    xi1_proton_pb*=gRandom->Gaus(1,0.02); //smearing 2%
    xi2_proton_pb*=gRandom->Gaus(1,0.02); //smearing 2%
    
    h_proton1_xi->Fill(xi1_proton_pb);
    h_proton2_xi->Fill(xi2_proton_pb);
    
    double mx = sqrt(xi1_proton_pb*xi2_proton_pb)*ROOTs;

    h_mx_PbPb->Fill(mx);

    //==================================================================================
    // t, pt information
    //=================================================================================
    double const m_p = 0.938272046;  //GeV  
    double const m_n = 0.9315;  //GeV  
    TLorentzVector E_Beam1;
    TLorentzVector E_Beam2;
  
    // E_Beam1.SetPxPyPzE(0.,0.,sqrt(pow(19500,2)-pow(m_p,2)),19500.) ;   //positive
    //E_Beam2.SetPxPyPzE(0.,0.,-sqrt(pow(19500,2)-pow(m_p,2)),19500.) ;   //negative
    
    E_Beam1.SetPxPyPzE(0.,0.,sqrt(pow(31500,2)-pow(m_p,2)),31500.) ;   //positive
    E_Beam2.SetPxPyPzE(0.,0.,-sqrt(pow(31500,2)-pow(m_n,2)),31500.) ;   //negative



    TLorentzVector t_proton1,t_proton2;
    t_proton1 = ( proton1p4 - E_Beam1 );
    t_proton2 = ( proton2p4 - E_Beam2 );

    h_proton1_t->Fill(-t_proton1.M2()); 
    h_t_pt_proton1->Fill(proton1p4.Pt(),-t_proton1.M2());   
    h_proton2_t->Fill(-t_proton2.M2()); 
    h_t_pt_proton2->Fill(proton2p4.Pt(),-t_proton2.M2());   

   
    //============================================================================
    //FastJet Analysis
    //===========================================================================
    double ptmin= 0.0;
    vector <fastjet::PseudoJet> inclusiveJets;
    vector <fastjet::PseudoJet> sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);
    inclusiveJets = clustSeq.inclusive_jets(ptmin);
    sortedJets    = sorted_by_pt(inclusiveJets);
    
    double Mjj=(sortedJets[0]+sortedJets[1]).m();
    Mjj*=gRandom->Gaus(1,0.02); //smearing 2%

    double jet1_pt = sortedJets[0].pt();
    double jet2_pt = sortedJets[1].pt();
    
    h_jet1pt->Fill(jet1_pt);
    h_jet2pt->Fill(jet2_pt);
    h_invmassjet->Fill(Mjj);

    //==============================================================
    //Dijets mass fraction (DMF)
    //==============================================================
    //Exclusive events show peak at 1. Could this be one form to diminish the backgrounds? Discussion needed!
    double ratio_DMF = Mjj/mx;
    //==============================================================
    h_ratio_DMF->Fill(ratio_DMF);
        
  
  } //Fim do loop na tree

  //===================================================================
  //Normalização
  //====================================================================
   
  //Double_t  L_int = 33000;  //pb^-1  PbPb
  Double_t  L_int = 8;  //pb^-1  pPb
  //Double_t sigma_wwinc =0.06988536;  //pb  PbPb  
  //Double_t sigma_wwinc =0.2153924e-07;  //pb  pp  
  //Double_t sigma_wwinc =0.4954e-04;  //pb  pPb   39 TeV
  Double_t sigma_wwinc =0.759e-04;  //pb  pPb   63 TeV
  
  Double_t N_ww = nentries;    
  Double_t scale1 = (sigma_wwinc * L_int)/(N_ww);  


  h_npart->Scale(scale1);     
  h_part_id->Scale(scale1);
  h_part_pt->Scale(scale1);
  h_part_px->Scale(scale1);
  h_part_py->Scale(scale1);
  h_part_pz->Scale(scale1);
  h_part_energy->Scale(scale1);
  h_part_eta->Scale(scale1);
  h_part_charge->Scale(scale1);
  h_part_phi->Scale(scale1);
  h_part_cos->Scale(scale1);

  h_invmassjet->Scale(scale1);
  h_jet1pt->Scale(scale1);
  h_jet2pt->Scale(scale1);

  h_proton1_xi->Scale(scale1);  
  h_proton2_xi->Scale(scale1);
  h_mx_PbPb->Scale(scale1);


  h_proton1_t->Scale(scale1); 
  h_t_pt_proton1->Scale(scale1);   
  h_proton2_t->Scale(scale1); 
  h_t_pt_proton2->Scale(scale1);   

  h_ratio_DMF->Scale(scale1); 






  //========================================================================
  //arquivo de saida
  //TFile* output = new TFile("pp_aahbB_signal_analysis_39TeV.root","RECREATE");
  TFile* output = new TFile("pPb_aahbB_signal_analysis_63TeV.root","RECREATE");
  //========================================================================
  h_npart->Write();     
  h_part_id->Write();
  h_part_pt->Write();
  h_part_px->Write();
  h_part_py->Write();
  h_part_pz->Write();
  h_part_energy->Write();
  h_part_eta->Write();
  h_part_charge->Write();
  h_part_phi->Write();


  //Proton information
  h_proton1_xi->Write();
  h_proton2_xi->Write();
  h_mx_PbPb->Write();
  h_proton1_t->Write(); 
  h_t_pt_proton1->Write();   
  h_proton2_t->Write(); 
  h_t_pt_proton2->Write();   


  //Jet Information
  h_invmassjet->Write();
  h_jet1pt->Write();
  h_jet2pt->Write();
  h_part_cos->Write();
  h_ratio_DMF->Write();


  output->Close();
  m_timer.Stop();
  rtime = m_timer.RealTime();
  ctime = m_timer.CpuTime();

  cout << "********************************" << endl;
  cout << "INFORMATION :: Performance : " << endl;
  cout << "RealTime= " << rtime << "seconds, CpuTime= " << ctime << "seconds" << endl;
  cout << "********************************" << endl;  




} // The end! main loop















































/*    
//===================================================================
//Normalização
//====================================================================
   
Double_t  L_int = 1000.0;
Double_t sigma_wwinc =16;  
Double_t N_ww = 100000;    
Double_t scale1 = (sigma_wwinc * L_int)/(N_ww);  
    
h_part_pt->Scale(scale1);
h_part_charge->Scale(scale1);
h_part_px->Scale(scale1);
h_part_py->Scale(scale1);
h_part_pz->Scale(scale1);
h_nb->Scale(scale1);
h_part_eta->Scale(scale1);
h_part_phi->Scale(scale1);
h_part_energy->Scale(scale1);
*/
