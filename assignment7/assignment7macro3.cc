//comparing eta plots from different pdfs 
#include "TFile.h"
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;

TCanvas *c3 = new TCanvas(); //canvas for  Wminus eta distribution  
TCanvas *c4 = new TCanvas(); //canvas for  Wplus eta distribution  
TCanvas *c5 = new TCanvas(); //canvas for  Z eta distribution  
void assignment7macro2(){
	//pthat bins

//defining bins to seperate soft and hard qcd using pthat
static const int nbins =1;
static const double binedges[nbins+1] = {0., 150. };

	//initialise histograms

	// W minus histograms for different pdfs

	TH1F* sub_Wm_14000_muon_cross_section = new TH1F("sub_Wm_14000_muon_cross_section","", 50,  -10,10);
	TH1F* Wm_14000_muon_cross_section = new TH1F("Wm_14000_muon_cross_section","", 50,  -10,10);



	// W plus histograms for different pdfs

	TH1F* sub_Wp_14000_muon_cross_section = new TH1F("sub_Wp_14000_muon_cross_section","", 50,  -10,10);
	TH1F* Wp_14000_muon_cross_section = new TH1F("Wp_14000_muon_cross_section","", 50,  -10,10);



	// Z histograms for different pdfs

	TH1F* sub_Z_14000_muon_cross_section = new TH1F("sub_Z_14000_muon_cross_section","", 50,  -10,10);
	TH1F* Z_14000_muon_cross_section = new TH1F("Z_14000_muon_cross_section","", 50,  -10,10);




//*******************************************************************************
//Now looking at W-  
//*******************************************************************************
	
	

	//14000
	{
	TFile* infile = TFile::Open("./output/wm_pwgevents_14000.root", "READ");
	
	vector<double> *Luminosity;
	infile->GetObject("luminosity",Luminosity);


	//resetting sub bins (bins from etahat) to fill main histogram
	//reading files from pythia calculation
	TNtuple *muontuples = (TNtuple*)infile->Get("muons");
	Float_t type, eventNo,index, status, mother1, mother2, pAbs, pt, y, eta, id;	
	int particle_count = muontuples->GetEntries();	//number of muons
	muontuples->SetBranchAddress("eventNo", &eventNo);
	muontuples->SetBranchAddress("index", &index);
	muontuples->SetBranchAddress("mother1", &mother1);
	muontuples->SetBranchAddress("mother2",&mother2);
	muontuples->SetBranchAddress("pAbs", &pAbs);
	muontuples->SetBranchAddress("pt", &pt);
	muontuples->SetBranchAddress("y", &y);
	muontuples->SetBranchAddress("eta", &eta);
	muontuples->SetBranchAddress("id", &id);


	//filling histograms
	for (int event_counter = 0; event_counter < particle_count; event_counter++)
	{
		
		muontuples->GetEntry(event_counter);
	//if ((eta > -4.0) && (eta < -2.5)){
		
		if ((abs(mother1) == 24) || ((abs(mother1) == 13) && (abs(mother2) == 90)))
		{
			sub_Wm_14000_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wm_14000_muon_cross_section->Scale(scalebin, "width");
	Wm_14000_muon_cross_section->Add(sub_Wm_14000_muon_cross_section);
	}


	


//*******************************************************************************
//Now looking at W+  
//*******************************************************************************
	
	//14000
	{
	TFile* infile = TFile::Open("./output/wp_pwgevents_14000.root", "READ");
	
	vector<double> *Luminosity;
	infile->GetObject("luminosity",Luminosity);


	//resetting sub bins (bins from etahat) to fill main histogram
	//reading files from pythia calculation
	TNtuple *muontuples = (TNtuple*)infile->Get("muons");
	Float_t type, eventNo,index, status, mother1, mother2, pAbs, pt, y, eta, id;	
	int particle_count = muontuples->GetEntries();	//number of muons
	muontuples->SetBranchAddress("eventNo", &eventNo);
	muontuples->SetBranchAddress("index", &index);
	muontuples->SetBranchAddress("mother1", &mother1);
	muontuples->SetBranchAddress("mother2",&mother2);
	muontuples->SetBranchAddress("pAbs", &pAbs);
	muontuples->SetBranchAddress("pt", &pt);
	muontuples->SetBranchAddress("y", &y);
	muontuples->SetBranchAddress("eta", &eta);
	muontuples->SetBranchAddress("id", &id);


	//filling histograms
	for (int event_counter = 0; event_counter < particle_count; event_counter++)
	{
		
		muontuples->GetEntry(event_counter);
	//if ((eta > -4.0) && (eta < -2.5)){
		
		if ((abs(mother1) == 24) || ((abs(mother1) == 13) && (abs(mother2) == 90)))
		{
			sub_Wp_14000_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wp_14000_muon_cross_section->Scale(scalebin, "width");
	Wp_14000_muon_cross_section->Add(sub_Wp_14000_muon_cross_section);
	}


	//*******************************************************************************
//Now looking at Z  
//*******************************************************************************
	
	//14000
	{
	TFile* infile = TFile::Open("./output/z_pwgevents_14000.root", "READ");
	
	vector<double> *Luminosity;
	infile->GetObject("luminosity",Luminosity);


	//resetting sub bins (bins from etahat) to fill main histogram
	//reading files from pythia calculation
	TNtuple *muontuples = (TNtuple*)infile->Get("muons");
	Float_t type, eventNo,index, status, mother1, mother2, pAbs, pt, y, eta, id;	
	int particle_count = muontuples->GetEntries();	//number of muons
	muontuples->SetBranchAddress("eventNo", &eventNo);
	muontuples->SetBranchAddress("index", &index);
	muontuples->SetBranchAddress("mother1", &mother1);
	muontuples->SetBranchAddress("mother2",&mother2);
	muontuples->SetBranchAddress("pAbs", &pAbs);
	muontuples->SetBranchAddress("pt", &pt);
	muontuples->SetBranchAddress("y", &y);
	muontuples->SetBranchAddress("eta", &eta);
	muontuples->SetBranchAddress("id", &id);


	//filling histograms
	for (int event_counter = 0; event_counter < particle_count; event_counter++)
	{
		
		muontuples->GetEntry(event_counter);
	//if ((eta > -4.0) && (eta < -2.5)){
		
		if ((abs(mother1) == 24) || ((abs(mother1) == 13) && (abs(mother2) == 90)))
		{
			sub_Z_14000_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Z_14000_muon_cross_section->Scale(scalebin, "width");
	Z_14000_muon_cross_section->Add(sub_Z_14000_muon_cross_section);
	}



	TFile* outFile =new TFile("NLOmuonyieldmacro.root", "RECREATE");

	//plotting muons from W minus
	c3->cd();
	c3->SetGridy();

	Wm_14000_muon_cross_section->GetYaxis()->SetTitle("d#sigma/d#eta (pb/GeV/c)");
	Wm_14000_muon_cross_section->GetXaxis()->SetTitle("#eta");

	Wm_14000_muon_cross_section->Draw("SAME");

	Wm_14000_muon_cross_section->SetLineColor(kGreen);
	Wm_14000_muon_cross_section->SetMarkerStyle(27);
	Wm_14000_muon_cross_section->SetMarkerColor(kGreen);


	TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
	leg2->AddEntry(Wm_14000_muon_cross_section, "CT18NNLO", "lep");
	leg2->Draw("SAME");

	Wm_14000_muon_cross_section->SetMinimum(1);
	Wm_14000_muon_cross_section->SetMaximum(1750);
	c3->Write();


	//plotting muons from W plus
	c4->cd();
	c4->SetGridy();

	Wp_14000_muon_cross_section->GetYaxis()->SetTitle("d#sigma/d#eta (pb/GeV/c)");
	Wp_14000_muon_cross_section->GetXaxis()->SetTitle("#eta");

	Wp_14000_muon_cross_section->Draw("SAME");

	Wp_14000_muon_cross_section->SetLineColor(kGreen);
	Wp_14000_muon_cross_section->SetMarkerStyle(27);
	Wp_14000_muon_cross_section->SetMarkerColor(kGreen);

	TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
	leg->AddEntry(Wp_14000_muon_cross_section, "CT18NNLO", "lep");
	leg->Draw("SAME");

	Wp_14000_muon_cross_section->SetMinimum(1);
	Wp_14000_muon_cross_section->SetMaximum(1800);
	c4->Write();


	//plotting muons from Z
	
	c5->cd();
	c5->SetGridy();

	Z_14000_muon_cross_section->GetYaxis()->SetTitle("d#sigma/d#eta (pb/GeV/c)");
	Z_14000_muon_cross_section->GetXaxis()->SetTitle("#eta");

	Z_14000_muon_cross_section->Draw("SAME");

	Z_14000_muon_cross_section->SetLineColor(kGreen);
	Z_14000_muon_cross_section->SetMarkerStyle(27);
	Z_14000_muon_cross_section->SetMarkerColor(kGreen);

	TLegend *leg3 = new TLegend(0.6, 0.7, 0.9, 0.9);
	leg3->AddEntry(Z_14000_muon_cross_section, "CT18NNLO", "lep");
	leg3->Draw("SAME");

	Z_14000_muon_cross_section->SetMinimum(1);
	Z_14000_muon_cross_section->SetMaximum(400);
	c5->Write();

	delete outFile;
}
