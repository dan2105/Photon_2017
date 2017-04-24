//.Dobbs@Cern.CH, Feb 2000
// This example shows low to use the particle and vertex iterators
//////////////////////////////////////////////////////////////////////////
// To Compile: go to the HepMC directory and type:
// gmake examples/example_UsingIterators.exe
//

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/IO_GenEvent.h"
#include "HepPDT/TableBuilder.hh"
#include "HepPDT/ParticleDataTable.hh"

#include "TROOT.h"
#include "TSystem.h"
#include "TDataType.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
//#include "/home/dan-bia/workspace/Build_clhep/include/CLHEP/Vector/LorentzVector.h"
//#include "/home/dan-bia/workspace/x86_64-slc5-gcc47-opt/include/CLHEP/Vector/LorentzVector.h"
//#include "CLHEP/Vector/LorentzVector.h"
#include <math.h>
#include <algorithm>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector> 
#include <cstdio>
#include "/home/daniel/workspace/fastjet_build/include/fastjet/ClusterSequence.hh"
#include "/home/daniel/workspace/fastjet_build/include/fastjet/PseudoJet.hh"
#include "/home/daniel/workspace/fastjet_build/include/fastjet/JetDefinition.hh"
using namespace std;
using namespace fastjet;

/*
//==============================================================
PDG-ID taken from ROOT data table
//==============================================================
11      -11     e-         e+
12      -12    nu_e       nu_ebar
13      -13     mu-        mu+
14      -14    nu_mu      nu_mubar
15      -15    tau-       tau+
16      -16    nu_tau     nu_taubar
17      -17    tau'-      tau'+
18      -18    nu'_tau    nu'_taubar
24      -24     W+         W-
37      -37     H+         H-
23       22     Z0         gamma
//===========================================================
*/


//===================================
// Declaração de classes auxiliares
//===================================
class IsW_Boson {
public:
  /// returns true if the GenParticle is a W
  bool operator()( const HepMC::GenParticle* p ) { 
    if ( p->pdg_id() == 24 ) return 1;
    return 0;
  }
};


class IsStateFinal {
public:
  /// returns true if the GenParticle does not decay
  bool operator()( const HepMC::GenParticle* p ) { 
        if ( !p->end_vertex() && p->status()==1 ) return 1;
    //if ( p->status()==1 ) return 1;
    return 0;
  }
};

//===================================
//Event selection
//===================================

class IsEventGood {
public:
  /// check this event for goodness
  bool operator()( const HepMC::GenEvent* evt ) { 
    for( HepMC::GenEvent::particle_const_iterator p 
	   = evt->particles_begin(); p != evt->particles_end(); ++p ){
      if ( !(*p)->end_vertex()  && (*p)->status()==1 ){
	//  if ((*p)->status()==1 ){
	std::cout << "Event " << evt->event_number()
		  << " is a good event." << std::endl;
	(*p)->print();
	return 1;
      }
    }
    return 0;
  }
};
//========================================
//Main analysis
///=========================================


int main(){ 
  
  //  HepMC::IO_GenEvent ascii_in("hepmc_aabB_background.dat",std::ios::in);
  HepMC::IO_GenEvent ascii_in("aacc_background_240GeV.dat",std::ios::in);    

  // declare another IO_GenEvent for writing out the good events
    HepMC::IO_GenEvent ascii_out("output.dat",std::ios::out);
  // Build HepPDT particle table
  const char infile[] ="/home/daniel/workspace/heppdt/data/particle.tbl";   
  std::ifstream pdfile( infile );
  if( !pdfile ) { 
    std::cerr << ">>> Cannot open " << infile << std::endl;
    exit(-1);
  }
  HepPDT::ParticleDataTable pdt( "Particle Table" );
  {
    // Construct table builder
    HepPDT::TableBuilder tb(pdt);
    if( !addParticleTable( pdfile, tb, true ) ) { 
      std::cout << ">> Error reading PDG pdt file " << std::endl; 
    }
  } // the tb destructor fills datacol
    // Loop over Particle Data Table
  std::ostringstream oss;
  oss << std::setw(15) << "Particle Id"
      << std::setw(22) << "Particle Name"
      << std::setw(15) << "Three-charge" << std::endl; 
  for( HepPDT::ParticleDataTable::const_iterator p = pdt.begin(); p != pdt.end(); ++p ) {
    const HepPDT::ParticleID & id = p->first;
    int pdgId = id.pid();
    int q3 = id.threeCharge();
    //double q = id.charge();
    const std::string& name = id.PDTname();
    oss << std::setw(15) << pdgId
	<< std::setw(22) << name
	<< std::setw(15) << q3 << std::endl;
  }
  std::cout << oss.str();
    
    
   
  // Declare ROOT TTree
    
  // final quarks, b and bbar in the same branch, further analysis will be employed
  int const NPARTMAX = 1000;
  int n_part=0;
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


  int const PROTONMAX = 1000;
  int n_proton =0;
  double proton_pt[PROTONMAX];
  double proton_px[PROTONMAX];
  double proton_py[PROTONMAX];
  double proton_pz[PROTONMAX];
  double proton_energy[PROTONMAX];


    

  TTree* T = new TTree("T","Tree");
  T->Branch("n_part", &n_part,"n_part/I");
  T->Branch("part_id", &part_id,"part_id/I");
  T->Branch("part_pt", &part_pt,"part_pt[n_part]/D");
  T->Branch("part_px", &part_px,"part_px[n_part]/D");
  T->Branch("part_py", &part_py,"part_py[n_part]/D"); 
  T->Branch("part_pz", &part_pz,"part_pz[n_part]/D");
  T->Branch("part_eta", &part_eta,"part_eta[n_part]/D");
  T->Branch("part_phi", &part_phi,"part_phi[n_part]/D");
  T->Branch("part_energy", &part_energy,"part_energy[n_part]/D");
  T->Branch("part_mass", &part_mass,"part_mass[n_part]/D");
  T->Branch("part_charge", &part_charge,"part_charge[n_part]/D");


 //Proton
  T->Branch("n_proton", &n_proton,"n_proton/I");
  T->Branch("proton_px", &proton_px,"proton_px[n_proton]/D");
  T->Branch("proton_py", &proton_py,"proton_py[n_proton]/D");
  T->Branch("proton_pz", &proton_pz,"proton_pz[n_proton]/D");
  T->Branch("proton_pt", &proton_pt,"proton_pt[n_proton]/D");
  T->Branch("proton_energy", &proton_energy,"proton_energy[n_proton]/D");

    


  // Declarar os histogramas no ROOT 
  TH1F* h_npart = new TH1F("npart","npart",2000,0.,2000.);
  TH1F* h_part_id = new TH1F("part_id","part_id",6000,-3000.,3000.);
  TH1F* h_part_pt = new TH1F("part_pt","part_pt",1000,0.,10.);
  TH1F* h_part_px = new TH1F("part_px","part_px",1000,0.,10.);
  TH1F* h_part_py = new TH1F("part_py","part_py",1000,0.,10.);
  TH1F* h_part_pz = new TH1F("part_pz","part_pz",1000,0.,10.);
  TH1F* h_part_eta = new TH1F("part_eta","part_eta",1000,-5.0,5.0);
  TH1F* h_part_phi = new TH1F("part_phi","part_phi",1000,-M_PI,M_PI);
  TH1F* h_part_energy = new TH1F("part_energy","part_energy",1000,0.,1000.);
  TH1F* h_part_mass = new TH1F("part_mass","part_mass",1000,0.,5.);
  TH1F* h_part_charge = new TH1F("part_charge","part_charge",20,-10.,10.);

  TH1F* h_nproton = new TH1F("nproton","nproton",100,0.,100.);
  TH1F* h_proton_pt = new TH1F("proton_pt","proton_pt",1000,0.,10.);
  TH1F* h_proton_px = new TH1F("proton_px","proton_px",1000,0.,10.);
  TH1F* h_proton_py = new TH1F("proton_py","proton_py",1000,0.,10.);
  TH1F* h_proton_pz = new TH1F("proton_pz","proton_pz",20000,-10000.,10000.);
  TH1F* h_proton_energy = new TH1F("proton_energy","proton_energy",60000,-30000.,30000.);




 // EVENT LOOP
  // declare an instance of the event selection predicate
  IsEventGood is_good_event;
  IsStateFinal isfinal;
  int icount=0;
  int num_good_events=0;
  HepMC::GenEvent* evt = ascii_in.read_next_event();
    
   
  while ( evt ) {
    icount++;
    if ( icount%50==1 ) std::cout << "Processing Event Number " << icount
				  << " its # " << evt->event_number() 
				  << std::endl;
    // Reset Tree variables per event
    
      
    n_part = 0;
    for(int ipart = 0; ipart < NPARTMAX; ++ipart) {
      part_id[ipart] =  -999.;     
      part_pt[ipart] =  -999.;
      part_px[ipart] =  -999.;
      part_py[ipart] =  -999.; 
      part_pz[ipart] =  -999.;
      part_eta[ipart] = -999.;
      part_phi[ipart] = -999.;
      part_energy[ipart] = -999.;
      part_mass[ipart] = -999.;
      part_charge[ipart] = -999.;}
      

    n_proton=0; 
    for(int iproton = 0; iproton < PROTONMAX; ++iproton) {
      proton_pt[iproton] = -999.;
      proton_px[iproton] = -999.;
      proton_py[iproton] = -999.;
      proton_pz[iproton] = -999.;
      proton_energy[iproton] = -999.;
    }

  
      
    if ( is_good_event(evt) ) {

      for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p ){

		if ( isfinal(*p) && (*p)->status() == 1 && (*p)->momentum().perp()> 0) {
		  //  	if ( isfinal(*p)) {
	  int pdg_id = (*p)->pdg_id();
	  // Could potentially be slow.
	  // Instead build map with charge per PDG id before event loop.
	  HepPDT::ParticleData * pdb;
	  pdb = pdt.particle( HepPDT::ParticleID( pdg_id ) );
	  double bcharge = pdb->charge();
	  part_id[n_part] = pdg_id; 
	  part_charge[n_part] = bcharge;
	  part_px[n_part] = (*p)->momentum().px();
	  part_pt[n_part] = (*p)->momentum().perp();
	  part_py[n_part] = (*p)->momentum().py();
	  part_pz[n_part] = (*p)->momentum().pz();
	  part_eta[n_part] = (*p)->momentum().eta(); 
	  part_phi[n_part] = (*p)->momentum().phi(); 
	  part_energy[n_part] = (*p)->momentum().e();
	  part_mass[n_part] = (*p)->momentum().m();
	    
	  h_part_id->Fill(part_id[n_part]);
	  h_part_px->Fill(part_px[n_part]);
	  h_part_py->Fill(part_py[n_part]);
	  h_part_pz->Fill(part_pz[n_part]);
	  h_part_pt->Fill(part_pt[n_part]);
	  h_part_eta->Fill(part_eta[n_part] );
	  h_part_phi->Fill(part_phi[n_part] );
	  h_part_energy->Fill(part_energy[n_part] );
	  h_part_mass->Fill(part_mass[n_part] );
	  h_part_charge->Fill(part_charge[n_part]);
	  ++n_part;
	}

	if (abs((*p)->pdg_id() == 2212) && (*p)->status() == 1 ){
	  proton_px[n_proton] = (*p)->momentum().px();
	  proton_py[n_proton] = (*p)->momentum().py();
	  proton_pz[n_proton] = (*p)->momentum().pz();
	  proton_pt[n_proton] = (*p)->momentum().perp();
	  proton_energy[n_proton] =(*p)->momentum().e();
	    
	  h_proton_px->Fill( proton_px[n_proton]);
	  h_proton_py->Fill( proton_py[n_proton]);
	  h_proton_pz->Fill(proton_pz[n_proton]);
	  h_proton_pt->Fill( proton_pt[n_proton]);
	  h_proton_energy->Fill(proton_energy[n_proton]);
	  ++n_proton;
	}
 

	 
      }
      h_npart->Fill(n_part);


      T->Fill();
      ascii_out << evt;
      ++num_good_events;
    }
    // Here it will fill TTree for all events 
    //       T->Fill();
    delete evt;
    ascii_in >> evt;
  }
  //........................................PRINT RESULT
  std::cout << num_good_events << " out of " << icount 
	    << " processed events passed the cuts. Finished." << std::endl;
  T->Print();
  // Output file
  TFile* output = new TFile("aacc_background_240GeV_hepmc.root","RECREATE");
  output->cd();
  // Write TTree and histograms to file
  T->Write();
  h_part_id->Write();
  h_npart->Write();
  h_part_pt->Write();
  h_part_px->Write();
  h_part_py->Write();
  h_part_pz->Write();
  h_part_eta->Write();
  h_part_phi->Write();
  h_part_energy->Write();
  h_part_mass->Write();
  h_part_charge->Write();


  h_nproton->Write();
  h_proton_px->Write();
  h_proton_py->Write();
  h_proton_pz->Write();
  h_proton_pt->Write();
  h_proton_energy->Write();

    

   
  output->Close();
    
  return 0;
} //end of int main



