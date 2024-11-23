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
	TH1F* sub_Wm_10042_muon_cross_section = new TH1F("sub_Wm_10042_muon_cross_section","", 50,  -10,10);
	TH1F* Wm_10042_muon_cross_section = new TH1F("Wm_10042_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wm_10555_muon_cross_section = new TH1F("sub_Wm_10555_muon_cross_section","", 50,  -10,10);
	TH1F* Wm_10555_muon_cross_section = new TH1F("Wm_10555_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wm_14000_muon_cross_section = new TH1F("sub_Wm_14000_muon_cross_section","", 50,  -10,10);
	TH1F* Wm_14000_muon_cross_section = new TH1F("Wm_14000_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wm_27000_muon_cross_section = new TH1F("sub_Wm_27000_muon_cross_section","", 50,  -10,10);
	TH1F* Wm_27000_muon_cross_section = new TH1F("Wm_27000_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wm_315000_muon_cross_section = new TH1F("sub_Wm_315000_muon_cross_section","", 50,  -10,10);
	TH1F* Wm_315000_muon_cross_section = new TH1F("Wm_315000_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wm_315200_muon_cross_section = new TH1F("sub_Wm_315200_muon_cross_section","", 50,  -10,10);
	TH1F* Wm_315200_muon_cross_section = new TH1F("Wm_315200_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wm_331900_muon_cross_section = new TH1F("sub_Wm_331900_muon_cross_section","", 50,  -10,10);
	TH1F* Wm_331900_muon_cross_section = new TH1F("Wm_331900_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wm_10050_muon_cross_section = new TH1F("sub_Wm_10050_muon_cross_section","", 50,  -10,10);
	TH1F* Wm_10050_muon_cross_section = new TH1F("Wm_10050_muon_cross_section","", 50,  -10,10);


	// W plus histograms for different pdfs
	TH1F* sub_Wp_10042_muon_cross_section = new TH1F("sub_Wp_10042_muon_cross_section","", 50,  -10,10);
	TH1F* Wp_10042_muon_cross_section = new TH1F("Wp_10042_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wp_10555_muon_cross_section = new TH1F("sub_Wp_10555_muon_cross_section","", 50,  -10,10);
	TH1F* Wp_10555_muon_cross_section = new TH1F("Wp_10555_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wp_14000_muon_cross_section = new TH1F("sub_Wp_14000_muon_cross_section","", 50,  -10,10);
	TH1F* Wp_14000_muon_cross_section = new TH1F("Wp_14000_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wp_27000_muon_cross_section = new TH1F("sub_Wp_27000_muon_cross_section","", 50,  -10,10);
	TH1F* Wp_27000_muon_cross_section = new TH1F("Wp_27000_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wp_315000_muon_cross_section = new TH1F("sub_Wp_315000_muon_cross_section","", 50,  -10,10);
	TH1F* Wp_315000_muon_cross_section = new TH1F("Wp_315000_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wp_315200_muon_cross_section = new TH1F("sub_Wp_315200_muon_cross_section","", 50,  -10,10);
	TH1F* Wp_315200_muon_cross_section = new TH1F("Wp_315200_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wp_331900_muon_cross_section = new TH1F("sub_Wp_331900_muon_cross_section","", 50,  -10,10);
	TH1F* Wp_331900_muon_cross_section = new TH1F("Wp_331900_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wp_10050_muon_cross_section = new TH1F("sub_Wp_10050_muon_cross_section","", 50,  -10,10);
	TH1F* Wp_10050_muon_cross_section = new TH1F("Wp_10050_muon_cross_section","", 50,  -10,10);


	// Z histograms for different pdfs
	TH1F* sub_Z_10042_muon_cross_section = new TH1F("sub_Z_10042_muon_cross_section","", 50,  -10,10);
	TH1F* Z_10042_muon_cross_section = new TH1F("Z_10042_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Z_10555_muon_cross_section = new TH1F("sub_Z_10555_muon_cross_section","", 50,  -10,10);
	TH1F* Z_10555_muon_cross_section = new TH1F("Z_10555_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Z_14000_muon_cross_section = new TH1F("sub_Z_14000_muon_cross_section","", 50,  -10,10);
	TH1F* Z_14000_muon_cross_section = new TH1F("Z_14000_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Z_27000_muon_cross_section = new TH1F("sub_Z_27000_muon_cross_section","", 50,  -10,10);
	TH1F* Z_27000_muon_cross_section = new TH1F("Z_27000_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Z_315000_muon_cross_section = new TH1F("sub_Z_315000_muon_cross_section","", 50,  -10,10);
	TH1F* Z_315000_muon_cross_section = new TH1F("Z_315000_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Z_315200_muon_cross_section = new TH1F("sub_Z_315200_muon_cross_section","", 50,  -10,10);
	TH1F* Z_315200_muon_cross_section = new TH1F("Z_315200_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Z_331900_muon_cross_section = new TH1F("sub_Z_331900_muon_cross_section","", 50,  -10,10);
	TH1F* Z_331900_muon_cross_section = new TH1F("Z_331900_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Z_10050_muon_cross_section = new TH1F("sub_Z_10050_muon_cross_section","", 50,  -10,10);
	TH1F* Z_10050_muon_cross_section = new TH1F("Z_10050_muon_cross_section","", 50,  -10,10);









//*******************************************************************************
//Now looking at W-  
//*******************************************************************************
	
	//10042
	{
	TFile* infile = TFile::Open("./output/wm_pwgevents_10042.root", "READ");
	
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
			sub_Wm_10042_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wm_10042_muon_cross_section->Scale(scalebin, "width");
	Wm_10042_muon_cross_section->Add(sub_Wm_10042_muon_cross_section);
	}
	
	//10555
	{
	TFile* infile = TFile::Open("./output/wm_pwgevents_10555.root", "READ");
	
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
			sub_Wm_10555_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wm_10555_muon_cross_section->Scale(scalebin, "width");
	Wm_10555_muon_cross_section->Add(sub_Wm_10555_muon_cross_section);
	}

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


	//27000
	{
	TFile* infile = TFile::Open("./output/wm_pwgevents_27000.root", "READ");
	
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
			sub_Wm_27000_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wm_27000_muon_cross_section->Scale(scalebin, "width");
	Wm_27000_muon_cross_section->Add(sub_Wm_27000_muon_cross_section);
	}

	//315000
	{
	TFile* infile = TFile::Open("./output/wm_pwgevents_315000.root", "READ");
	
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
			sub_Wm_315000_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wm_315000_muon_cross_section->Scale(scalebin, "width");
	Wm_315000_muon_cross_section->Add(sub_Wm_315000_muon_cross_section);
	}//315000
	{
	TFile* infile = TFile::Open("./output/wm_pwgevents_315000.root", "READ");
	
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
			sub_Wm_315000_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wm_315000_muon_cross_section->Scale(scalebin, "width");
	Wm_315000_muon_cross_section->Add(sub_Wm_315000_muon_cross_section);
	}

	//315200
	{
	TFile* infile = TFile::Open("./output/wm_pwgevents_315200.root", "READ");
	
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
			sub_Wm_315200_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wm_315200_muon_cross_section->Scale(scalebin, "width");
	Wm_315200_muon_cross_section->Add(sub_Wm_315200_muon_cross_section);
	}


	//331900
	{
	TFile* infile = TFile::Open("./output/wm_pwgevents_331900.root", "READ");
	
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
			sub_Wm_331900_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wm_331900_muon_cross_section->Scale(scalebin, "width");
	Wm_331900_muon_cross_section->Add(sub_Wm_331900_muon_cross_section);
	}

	//10050
	{
	TFile* infile = TFile::Open("../assignment6/wm_NLOmuonyield.root", "READ");
	
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
			sub_Wm_10050_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wm_10050_muon_cross_section->Scale(scalebin, "width");
	Wm_10050_muon_cross_section->Add(sub_Wm_10050_muon_cross_section);
	}









//*******************************************************************************
//Now looking at W+  
//*******************************************************************************
	
	//10042
	{
	TFile* infile = TFile::Open("./output/wp_pwgevents_10042.root", "READ");
	
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
			sub_Wp_10042_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wp_10042_muon_cross_section->Scale(scalebin, "width");
	Wp_10042_muon_cross_section->Add(sub_Wp_10042_muon_cross_section);
	}
	
	//10555
	{
	TFile* infile = TFile::Open("./output/wp_pwgevents_10555.root", "READ");
	
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
			sub_Wp_10555_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wp_10555_muon_cross_section->Scale(scalebin, "width");
	Wp_10555_muon_cross_section->Add(sub_Wp_10555_muon_cross_section);
	}

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


	//27000
	{
	TFile* infile = TFile::Open("./output/wp_pwgevents_27000.root", "READ");
	
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
			sub_Wp_27000_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wp_27000_muon_cross_section->Scale(scalebin, "width");
	Wp_27000_muon_cross_section->Add(sub_Wp_27000_muon_cross_section);
	}

	//315000
	{
	TFile* infile = TFile::Open("./output/wp_pwgevents_315000.root", "READ");
	
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
			sub_Wp_315000_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wp_315000_muon_cross_section->Scale(scalebin, "width");
	Wp_315000_muon_cross_section->Add(sub_Wp_315000_muon_cross_section);
	}//315000
	{
	TFile* infile = TFile::Open("./output/wp_pwgevents_315000.root", "READ");
	
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
			sub_Wp_315000_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wp_315000_muon_cross_section->Scale(scalebin, "width");
	Wp_315000_muon_cross_section->Add(sub_Wp_315000_muon_cross_section);
	}

	//315200
	{
	TFile* infile = TFile::Open("./output/wp_pwgevents_315200.root", "READ");
	
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
			sub_Wp_315200_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wp_315200_muon_cross_section->Scale(scalebin, "width");
	Wp_315200_muon_cross_section->Add(sub_Wp_315200_muon_cross_section);
	}


	//331900
	{
	TFile* infile = TFile::Open("./output/wp_pwgevents_331900.root", "READ");
	
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
			sub_Wp_331900_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wp_331900_muon_cross_section->Scale(scalebin, "width");
	Wp_331900_muon_cross_section->Add(sub_Wp_331900_muon_cross_section);
	}

	//10050
	{
	TFile* infile = TFile::Open("../assignment6/wp_NLOmuonyield.root", "READ");
	
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
			sub_Wp_10050_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Wp_10050_muon_cross_section->Scale(scalebin, "width");
	Wp_10050_muon_cross_section->Add(sub_Wp_10050_muon_cross_section);
	}



//*******************************************************************************
//Now looking at Z  
//*******************************************************************************
	
	//10042
	{
	TFile* infile = TFile::Open("./output/z_pwgevents_10042.root", "READ");
	
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
			sub_Z_10042_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Z_10042_muon_cross_section->Scale(scalebin, "width");
	Z_10042_muon_cross_section->Add(sub_Z_10042_muon_cross_section);
	}
	
	//10555
	{
	TFile* infile = TFile::Open("./output/z_pwgevents_10555.root", "READ");
	
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
			sub_Z_10555_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Z_10555_muon_cross_section->Scale(scalebin, "width");
	Z_10555_muon_cross_section->Add(sub_Z_10555_muon_cross_section);
	}

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


	//27000
	{
	TFile* infile = TFile::Open("./output/z_pwgevents_27000.root", "READ");
	
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
			sub_Z_27000_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Z_27000_muon_cross_section->Scale(scalebin, "width");
	Z_27000_muon_cross_section->Add(sub_Z_27000_muon_cross_section);
	}

	//315000
	{
	TFile* infile = TFile::Open("./output/z_pwgevents_315000.root", "READ");
	
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
			sub_Z_315000_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Z_315000_muon_cross_section->Scale(scalebin, "width");
	Z_315000_muon_cross_section->Add(sub_Z_315000_muon_cross_section);
	}//315000
	{
	TFile* infile = TFile::Open("./output/z_pwgevents_315000.root", "READ");
	
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
			sub_Z_315000_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Z_315000_muon_cross_section->Scale(scalebin, "width");
	Z_315000_muon_cross_section->Add(sub_Z_315000_muon_cross_section);
	}

	//315200
	{
	TFile* infile = TFile::Open("./output/z_pwgevents_315200.root", "READ");
	
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
			sub_Z_315200_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Z_315200_muon_cross_section->Scale(scalebin, "width");
	Z_315200_muon_cross_section->Add(sub_Z_315200_muon_cross_section);
	}


	//331900
	{
	TFile* infile = TFile::Open("./output/z_pwgevents_331900.root", "READ");
	
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
			sub_Z_331900_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Z_331900_muon_cross_section->Scale(scalebin, "width");
	Z_331900_muon_cross_section->Add(sub_Z_331900_muon_cross_section);
	}

	//10050
	{
	TFile* infile = TFile::Open("../assignment6/Z_NLOmuonyield.root", "READ");
	
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
			sub_Z_10050_muon_cross_section->Fill(eta);
		}


	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_Z_10050_muon_cross_section->Scale(scalebin, "width");
	Z_10050_muon_cross_section->Add(sub_Z_10050_muon_cross_section);
	}








	TFile* outFile =new TFile("NLOmuonyieldmacro.root", "RECREATE");

	//plotting muons from W minus
	c3->cd();
	c3->SetGridy();

	Wm_10042_muon_cross_section->GetYaxis()->SetTitle("d#sigma/d#eta (pb/GeV/c)");
	Wm_10042_muon_cross_section->GetXaxis()->SetTitle("#eta");

	Wm_10042_muon_cross_section->Draw("SAME");
	Wm_10555_muon_cross_section->Draw("SAME");
	Wm_14000_muon_cross_section->Draw("SAME");
	Wm_27000_muon_cross_section->Draw("SAME");
	Wm_315000_muon_cross_section->Draw("SAME");
	Wm_315200_muon_cross_section->Draw("SAME");
	Wm_331900_muon_cross_section->Draw("SAME");
	Wm_10050_muon_cross_section->Draw("SAME");

	Wm_10042_muon_cross_section->SetLineColor(kBlack);
	Wm_10042_muon_cross_section->SetMarkerStyle(28);
	Wm_10042_muon_cross_section->SetMarkerColor(kBlack);

	Wm_10555_muon_cross_section->SetLineColor(kRed);
	Wm_10555_muon_cross_section->SetMarkerStyle(25);
	Wm_10555_muon_cross_section->SetMarkerColor(kRed);

	Wm_14000_muon_cross_section->SetLineColor(kGreen);
	Wm_14000_muon_cross_section->SetMarkerStyle(27);
	Wm_14000_muon_cross_section->SetMarkerColor(kGreen);

	Wm_27000_muon_cross_section->SetLineColor(kBlue);
	Wm_27000_muon_cross_section->SetMarkerStyle(29);
	Wm_27000_muon_cross_section->SetMarkerColor(kBlue);

	Wm_315000_muon_cross_section->SetLineColor(kMagenta);
	Wm_315000_muon_cross_section->SetMarkerStyle(20);
	Wm_315000_muon_cross_section->SetMarkerColor(kMagenta);

	Wm_315200_muon_cross_section->SetLineColor(8);
	Wm_315200_muon_cross_section->SetMarkerStyle(21);
	Wm_315200_muon_cross_section->SetMarkerColor(8);

	Wm_331900_muon_cross_section->SetLineColor(kYellow);
	Wm_331900_muon_cross_section->SetMarkerStyle(23);
	Wm_331900_muon_cross_section->SetMarkerColor(kYellow);

	Wm_10050_muon_cross_section->SetLineColor(9);
	Wm_10050_muon_cross_section->SetMarkerStyle(23);
	Wm_10050_muon_cross_section->SetMarkerColor(9);

	TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
	leg2->AddEntry(Wm_10042_muon_cross_section, "CTEQ6l1", "lep");
	leg2->AddEntry(Wm_10555_muon_cross_section, "CTEQ66", "lep");
	leg2->AddEntry(Wm_14000_muon_cross_section, "CT18NNLO", "lep");
	leg2->AddEntry(Wm_27000_muon_cross_section, "MSHT20LO", "lep");
	leg2->AddEntry(Wm_315000_muon_cross_section, "NNPDF3.1LO #alpha (M_Z) = 0.118", "lep");
	leg2->AddEntry(Wm_315200_muon_cross_section, "NNPDF3.1LO #alpha (M_Z) = 0.130", "lep");
	leg2->AddEntry(Wm_331900_muon_cross_section, "NNPDF4.0LO", "lep");
	leg2->AddEntry(Wm_10050_muon_cross_section, "CTEQ6mE?", "lep");
	leg2->Draw("SAME");

	Wm_10042_muon_cross_section->SetMinimum(1);
	Wm_10042_muon_cross_section->SetMaximum(5000);
	c3->Write();


	//plotting muons from W plus
	c4->cd();
	c4->SetGridy();

	Wp_10042_muon_cross_section->GetYaxis()->SetTitle("d#sigma/d#eta (pb/GeV/c)");
	Wp_10042_muon_cross_section->GetXaxis()->SetTitle("#eta");

	Wp_10042_muon_cross_section->Draw("SAME");
	Wp_10555_muon_cross_section->Draw("SAME");
	Wp_14000_muon_cross_section->Draw("SAME");
	Wp_27000_muon_cross_section->Draw("SAME");
	Wp_315000_muon_cross_section->Draw("SAME");
	Wp_315200_muon_cross_section->Draw("SAME");
	Wp_331900_muon_cross_section->Draw("SAME");
	Wp_10050_muon_cross_section->Draw("SAME");

	Wp_10042_muon_cross_section->SetLineColor(kBlack);
	Wp_10042_muon_cross_section->SetMarkerStyle(28);
	Wp_10042_muon_cross_section->SetMarkerColor(kBlack);

	Wp_10555_muon_cross_section->SetLineColor(kRed);
	Wp_10555_muon_cross_section->SetMarkerStyle(25);
	Wp_10555_muon_cross_section->SetMarkerColor(kRed);

	Wp_14000_muon_cross_section->SetLineColor(kGreen);
	Wp_14000_muon_cross_section->SetMarkerStyle(27);
	Wp_14000_muon_cross_section->SetMarkerColor(kGreen);

	Wp_27000_muon_cross_section->SetLineColor(kBlue);
	Wp_27000_muon_cross_section->SetMarkerStyle(29);
	Wp_27000_muon_cross_section->SetMarkerColor(kBlue);

	Wp_315000_muon_cross_section->SetLineColor(kMagenta);
	Wp_315000_muon_cross_section->SetMarkerStyle(20);
	Wp_315000_muon_cross_section->SetMarkerColor(kMagenta);

	Wp_315200_muon_cross_section->SetLineColor(8);
	Wp_315200_muon_cross_section->SetMarkerStyle(21);
	Wp_315200_muon_cross_section->SetMarkerColor(8);

	Wp_331900_muon_cross_section->SetLineColor(kYellow);
	Wp_331900_muon_cross_section->SetMarkerStyle(23);
	Wp_331900_muon_cross_section->SetMarkerColor(kYellow);

	Wp_10050_muon_cross_section->SetLineColor(9);
	Wp_10050_muon_cross_section->SetMarkerStyle(23);
	Wp_10050_muon_cross_section->SetMarkerColor(9);

	TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
	leg->AddEntry(Wp_10042_muon_cross_section, "CTEQ6l1", "lep");
	leg->AddEntry(Wp_10555_muon_cross_section, "CTEQ66", "lep");
	leg->AddEntry(Wp_14000_muon_cross_section, "CT18NNLO", "lep");
	leg->AddEntry(Wp_27000_muon_cross_section, "MSHT20LO", "lep");
	leg->AddEntry(Wp_315000_muon_cross_section, "NNPDF3.1LO #alpha (M_Z) = 0.118", "lep");
	leg->AddEntry(Wp_315200_muon_cross_section, "NNPDF3.1LO #alpha (M_Z) = 0.130", "lep");
	leg->AddEntry(Wp_331900_muon_cross_section, "NNPDF4.0LO", "lep");
	leg->AddEntry(Wp_10050_muon_cross_section, "CTEQ6mE?", "lep");
	leg->Draw("SAME");

	Wp_10042_muon_cross_section->SetMinimum(1);
	Wp_10042_muon_cross_section->SetMaximum(5000);
	c4->Write();





	//plotting muons from Z
	
	c5->cd();
	c5->SetGridy();

	Z_10042_muon_cross_section->GetYaxis()->SetTitle("d#sigma/d#eta (pb/GeV/c)");
	Z_10042_muon_cross_section->GetXaxis()->SetTitle("#eta");

	Z_10042_muon_cross_section->Draw("SAME");
	Z_10555_muon_cross_section->Draw("SAME");
	Z_14000_muon_cross_section->Draw("SAME");
	Z_27000_muon_cross_section->Draw("SAME");
	Z_315000_muon_cross_section->Draw("SAME");
	Z_315200_muon_cross_section->Draw("SAME");
	Z_331900_muon_cross_section->Draw("SAME");
	Z_10050_muon_cross_section->Draw("SAME");

	Z_10042_muon_cross_section->SetLineColor(kBlack);
	Z_10042_muon_cross_section->SetMarkerStyle(28);
	Z_10042_muon_cross_section->SetMarkerColor(kBlack);

	Z_10555_muon_cross_section->SetLineColor(kRed);
	Z_10555_muon_cross_section->SetMarkerStyle(25);
	Z_10555_muon_cross_section->SetMarkerColor(kRed);

	Z_14000_muon_cross_section->SetLineColor(kGreen);
	Z_14000_muon_cross_section->SetMarkerStyle(27);
	Z_14000_muon_cross_section->SetMarkerColor(kGreen);

	Z_27000_muon_cross_section->SetLineColor(kBlue);
	Z_27000_muon_cross_section->SetMarkerStyle(29);
	Z_27000_muon_cross_section->SetMarkerColor(kBlue);

	Z_315000_muon_cross_section->SetLineColor(kMagenta);
	Z_315000_muon_cross_section->SetMarkerStyle(20);
	Z_315000_muon_cross_section->SetMarkerColor(kMagenta);

	Z_315200_muon_cross_section->SetLineColor(8);
	Z_315200_muon_cross_section->SetMarkerStyle(21);
	Z_315200_muon_cross_section->SetMarkerColor(8);

	Z_331900_muon_cross_section->SetLineColor(kYellow);
	Z_331900_muon_cross_section->SetMarkerStyle(23);
	Z_331900_muon_cross_section->SetMarkerColor(kYellow);

	Z_10050_muon_cross_section->SetLineColor(9);
	Z_10050_muon_cross_section->SetMarkerStyle(23);
	Z_10050_muon_cross_section->SetMarkerColor(9);

	TLegend *leg3 = new TLegend(0.6, 0.7, 0.9, 0.9);
	leg3->AddEntry(Z_10042_muon_cross_section, "CTEQ6l1", "lep");
	leg3->AddEntry(Z_10555_muon_cross_section, "CTEQ66", "lep");
	leg3->AddEntry(Z_14000_muon_cross_section, "CT18NNLO", "lep");
	leg3->AddEntry(Z_27000_muon_cross_section, "MSHT20LO", "lep");
	leg3->AddEntry(Z_315000_muon_cross_section, "NNPDF3.1LO #alpha (M_Z) = 0.118", "lep");
	leg3->AddEntry(Z_315200_muon_cross_section, "NNPDF3.1LO #alpha (M_Z) = 0.130", "lep");
	leg3->AddEntry(Z_331900_muon_cross_section, "NNPDF4.0LO", "lep");
	leg3->AddEntry(Z_10050_muon_cross_section, "CTEQ6mE?", "lep");
	leg3->Draw("SAME");

	Z_10042_muon_cross_section->SetMinimum(1);
	Z_10042_muon_cross_section->SetMaximum(5000);
	c5->Write();

	delete outFile;
}
