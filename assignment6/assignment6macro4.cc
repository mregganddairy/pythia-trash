//Combining all results from background and NLO into one eta plot

#include "TFile.h"
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;

TCanvas *c3 = new TCanvas(); //canvas for rapidity eta distribution  
void assignment6macro4(){
	//pthat bins

//defining bins to seperate soft and hard qcd using pthat
static const int nbins =1;
static const double binedges[nbins+1] = {0., 150. };

//creating lists of IDs of mother particles of interest
static const int b_moms[] = {511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 20523, 515, 525, 531, 533, 535, 541, 543,545, //charmed mesons
							551, 553, 555, 557,	//bbbar
							5122, 5112, 5212, 5222, 5114, 5214, 5224, 5132, 5232, 5312, 5322, 5314, 5324, 5332, 5334 //charmed baryons
							};
static const int c_moms[] = {411, 421, 10411, 10421, 413, 423, 10413, 10423, 20413, 20423, 415, 425, 431, 10431, 433, 10433, 20433, 435, //mesons
							441, 443, 445,	 //ccbar mesons
							4122, 4222, 4212, 4112, 4224, 4214, 4114, 4232, 4132, 4322, 4312, 4324, 4314, 4332, 4334, 4422 //mesons
							};
	
	//initialise histograms
	
	//total muon cross section with no cuts
	TH1F* sub_TOTAL_NLO_muon_cross_section = new TH1F("sub_TOTAL_NLO_muon_cross_section", "", 50,-10,10);
	TH1F* TOTAL_NLO_muon_cross_section = new TH1F("TOTAL_NLO_muon_cross_section", "muon differential cross section", 50,-10,10);

	//for pt distribution in central barrel
	TH1F* sub_NLO_muon_cross_section_cb = new TH1F("sub_NLO_muon_cross_section_cb","", 100, -10,10);
	TH1F* total_NLO_muon_cross_section_cb = new TH1F("total_NLO_muon_cross_section_cb","", 100, -10,10);

	//for pt distribution in forward region 
	TH1F* sub_NLO_muon_cross_section_fr = new TH1F("sub_NLO_muon_cross_section_fr","", 100, -10,10);
	TH1F* total_NLO_muon_cross_section_fr = new TH1F("total_NLO_muon_cross_section_fr","", 100, -10,10);


	//muon cross section produced from bottom or charm
	TH1F* sub_b_NLO_muon_cross_section = new TH1F("sub_b_NLO_muon_cross_section","", 50,  -10,10);
	TH1F* b_NLO_muon_cross_section = new TH1F("b_NLO_muon_cross_section","", 50,  -10,10);
	TH1F* sub_c_NLO_muon_cross_section = new TH1F("sub_c_NLO_muon_cross_section","", 50,  -10,10);
	TH1F* c_NLO_muon_cross_section = new TH1F("c_NLO_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wm_NLO_muon_cross_section = new TH1F("sub_Wm_NLO_muon_cross_section","", 50,  -10,10);
	TH1F* Wm_NLO_muon_cross_section = new TH1F("Wm_NLO_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Wp_NLO_muon_cross_section = new TH1F("sub_Wp_NLO_muon_cross_section","", 50,  -10,10);
	TH1F* Wp_NLO_muon_cross_section = new TH1F("Wp_NLO_muon_cross_section","", 50,  -10,10);

	TH1F* sub_W_NLO_muon_cross_section = new TH1F("sub_W_NLO_muon_cross_section","", 50,  -10,10);
	TH1F* W_NLO_muon_cross_section = new TH1F("W_NLO_muon_cross_section","", 50,  -10,10);

	TH1F* sub_Z_NLO_muon_cross_section = new TH1F("sub_Z_NLO_muon_cross_section","", 50,  -10,10);
	TH1F* Z_NLO_muon_cross_section = new TH1F("Z_NLO_muon_cross_section","", 50,  -10,10);

	//muons coming from muons
	TH1F* sub_muon_NLO_muon_cross_section = new TH1F("sub_muon_NLO_muon_cross_section","", 50,  -10,10);
	TH1F* muon_NLO_muon_cross_section = new TH1F("muon_NLO_muon_cross_section","", 50,  -10,10);

	//muon cross sections which don't come from bottom or charm
	TH1F* sub_other_NLO_muon_cross_section = new TH1F("sub_other_NLO_muon_cross_section","", 50,  -10,10);
	TH1F* other_NLO_muon_cross_section = new TH1F("other_NLO_muon_cross_section","", 50,  -10,10);








//*******************************************************************************
//background stuff
//*******************************************************************************
	//total muon cross section with no cuts
	TH1F* sub_TOTAL_muon_cross_section = new TH1F("sub_TOTAL_muon_cross_section", "", 50,-10,10);
	TH1F* TOTAL_muon_cross_section = new TH1F("TOTAL_muon_cross_section", "muon cross section from different mothers in 4#pi", 50,-10,10);


	//muon cross section produced from bottom or charm
	TH1F* sub_b_muon_cross_section = new TH1F("sub_b_muon_cross_section","", 50,  -10,10);
	TH1F* b_muon_cross_section = new TH1F("b_muon_cross_section","", 50,  -10,10);
	TH1F* sub_c_muon_cross_section = new TH1F("sub_c_muon_cross_section","", 50,  -10,10);
	TH1F* c_muon_cross_section = new TH1F("c_muon_cross_section","", 50,  -10,10);


	//muon cross sections which don't come from bottom or charm
	TH1F* sub_other_muon_cross_section = new TH1F("sub_other_muon_cross_section","", 50,  -10,10);
	TH1F* other_muon_cross_section = new TH1F("other_muon_cross_section","", 50,  -10,10);

	{
	//defining bins to seperate soft and hard qcd using pthat
	static const int nbins =6;
	static const double binedges[nbins+1] = {0., 14., 30., 50., 75., 100., 150. };
	
	//fetching goods
	TFile* infile = TFile::Open("/home/Josh/physics/pythia/assignment4/muonyield.root", "READ");
	
	vector<double> *binLuminosity;
	infile->GetObject("luminosity",binLuminosity);



	for (int ibin = 0; ibin < nbins; ibin++)
	{

		//resetting sub bins (bins from pthat) to fill main histogram
		sub_TOTAL_muon_cross_section->Reset();
		sub_b_muon_cross_section->Reset();
		sub_c_muon_cross_section->Reset();
		sub_other_muon_cross_section->Reset();

		//reading files from pythia calculation
		TNtuple *muontuples = (TNtuple*)infile->Get(Form("muon%d", ibin));
		Float_t type, eventNo,index, status, mother1, mother2, pAbs, pt, y, eta, id;	
		int particle_count = muontuples->GetEntries();	//number of muons
		muontuples->SetBranchAddress("binNo", &type);
		muontuples->SetBranchAddress("eventNo", &eventNo);
		muontuples->SetBranchAddress("index", &index);
		muontuples->SetBranchAddress("mother1", &mother1);
		muontuples->SetBranchAddress("mother2",&mother2);
		muontuples->SetBranchAddress("pAbs", &pAbs);
		muontuples->SetBranchAddress("pt", &pt);
		muontuples->SetBranchAddress("y", &y);
		muontuples->SetBranchAddress("eta", &eta);
		muontuples->SetBranchAddress("id", &id);

		//DEBUGGGING: check if entries have been correctly called
		//for (int entry = 0; entry < muontuples->GetEntries(); ++entry) 
		//{
		//	muontuples->GetEntry(entry);
		//	cout << "Entry " << entry << ": ";
//
//			cout << "type=" << type << " eventNo=" << eventNo << " mother1ID=" 
//				 << mother1 << " mother2ID=" << mother2 << " pAbs=" << pAbs
//				 << " pt=" << pt << " eta=" << eta << " id=" << id << endl;
//		}

		//filling histograms
		for (int event_counter = 0; event_counter < particle_count; event_counter++)
		{
			muontuples->GetEntry(event_counter);

			sub_TOTAL_muon_cross_section->Fill(eta);

			//check if mother is b or c meson/baryon.
			bool b_found = false;
			bool c_found = false;

			for (int i = 0; i < sizeof(b_moms)/sizeof(b_moms[0]); ++i)
			{
				if (abs(b_moms[i]) == abs(mother1))
				{
					b_found = true;
					break;
				}
			}

			for (int i = 0; i < sizeof(c_moms)/sizeof(c_moms[0]); ++i)
			{
				if (abs(c_moms[i]) == abs(mother1))
				{
					c_found = true;
					break;
				}
			}
			
		//checking the mother of the muon and filling the appropriate histogram
		if (b_found){sub_b_muon_cross_section->Fill(eta);}

		else if (c_found){sub_c_muon_cross_section->Fill(eta);}

	/*	else  
		{
			sub_other_muon_cross_section->Fill(pt);
			cout << " mom1 id: "<< mother1<< endl;
		}
		*/


		}
		
		//normalising bins for calculating the cross section
		double_t scalebin = 1./(*binLuminosity)[ibin];

		sub_TOTAL_muon_cross_section->Scale(scalebin, "width");
		sub_b_muon_cross_section->Scale(scalebin, "width");
		sub_c_muon_cross_section->Scale(scalebin, "width");
		sub_other_muon_cross_section->Scale(scalebin, "width");

		TOTAL_muon_cross_section->Add(sub_TOTAL_muon_cross_section);
		b_muon_cross_section->Add(sub_b_muon_cross_section);
		c_muon_cross_section->Add(sub_c_muon_cross_section);
		other_muon_cross_section->Add(sub_other_muon_cross_section);
	}	


	}





//*******************************************************************************
//Now looking at W+  
//*******************************************************************************
	
	{
	TFile* infile = TFile::Open("wp_NLOmuonyield.root", "READ");
	
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

	//DEBUGGGING: check if entries have been correctly called
	//for (int entry = 0; entry < muontuples->GetEntries(); ++entry) 
	//{
	//	muontuples->GetEntry(entry);
	//	cout << "Entry " << entry << ": ";
//
//			cout << "type=" << type << " eventNo=" << eventNo << " mother1ID=" 
//				 << mother1 << " mother2ID=" << mother2 << " pAbs=" << pAbs
//				 << " pt=" << pt << " eta=" << eta << " id=" << id << endl;
//		}

	//filling histograms
	for (int event_counter = 0; event_counter < particle_count; event_counter++)
	{
		
		muontuples->GetEntry(event_counter);
	//if ((eta > -4.0) && (eta < -2.5)){
		
		sub_TOTAL_NLO_muon_cross_section->Fill(eta);

		//check if muon is in a particular region for region plots
		if (abs(eta) < 0.9)
		{
			sub_NLO_muon_cross_section_cb->Fill(eta);
		}

		if (eta>-4.0 && eta<-2.5)
		{
			sub_NLO_muon_cross_section_fr->Fill(eta);
		}


		//check if mother is b or c meson/baryon.
		bool b_found = false;
		bool c_found = false;

		for (int i = 0; i < sizeof(b_moms)/sizeof(b_moms[0]); ++i)
		{
			if (abs(b_moms[i]) == abs(mother1))
			{
				b_found = true;
				break;
			}
		}

		for (int i = 0; i < sizeof(c_moms)/sizeof(c_moms[0]); ++i)
		{
			if (abs(c_moms[i]) == abs(mother1))
			{
				c_found = true;
				break;
			}
		}
		
		//checking the mother of the muon and filling the appropriate histogram
		if (b_found){sub_b_NLO_muon_cross_section->Fill(eta);}

		else if (c_found){sub_c_NLO_muon_cross_section->Fill(eta);}

		//checing if mother particle is a weak vector boson 

		else if ((abs(mother1) == 24) || ((abs(mother1) == 13) && (abs(mother2) == 90)))
		{
			sub_W_NLO_muon_cross_section->Fill(eta);

			sub_Wp_NLO_muon_cross_section->Fill(eta);
		}

		else if (!(abs(mother1) == 13) && !(abs(mother2) == 13))
		{
		sub_other_NLO_muon_cross_section->Fill(eta);
		//cout << " mom1 id: "<< mother1<< endl;
		}
		else if (((abs(mother1) == 13)) && !((abs(mother1) == 13) && (abs(mother2) == 90)))
		{
			sub_muon_NLO_muon_cross_section->Fill(eta);
		}
	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_TOTAL_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_NLO_muon_cross_section_cb->Scale(scalebin, "width");
	sub_NLO_muon_cross_section_fr->Scale(scalebin, "width");
	sub_b_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_c_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_other_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_W_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_Wm_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_Wp_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_Z_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_muon_NLO_muon_cross_section->Scale(scalebin, "width");
	

	TOTAL_NLO_muon_cross_section->Add(sub_TOTAL_NLO_muon_cross_section);
	total_NLO_muon_cross_section_cb->Add(sub_NLO_muon_cross_section_cb);
	total_NLO_muon_cross_section_fr->Add(sub_NLO_muon_cross_section_fr);
	b_NLO_muon_cross_section->Add(sub_b_NLO_muon_cross_section);
	c_NLO_muon_cross_section->Add(sub_c_NLO_muon_cross_section);
	other_NLO_muon_cross_section->Add(sub_other_NLO_muon_cross_section);
	W_NLO_muon_cross_section->Add(sub_W_NLO_muon_cross_section);
	Wp_NLO_muon_cross_section->Add(sub_Wp_NLO_muon_cross_section);
	Wm_NLO_muon_cross_section->Add(sub_Wm_NLO_muon_cross_section);
	Z_NLO_muon_cross_section->Add(sub_Z_NLO_muon_cross_section);
	muon_NLO_muon_cross_section->Add(sub_muon_NLO_muon_cross_section);
	}
	
	




//*******************************************************************************
//Now looking at W-  
//*******************************************************************************
	{
	sub_TOTAL_NLO_muon_cross_section->Reset();
	sub_b_NLO_muon_cross_section->Reset();
	sub_c_NLO_muon_cross_section->Reset();
	sub_other_NLO_muon_cross_section->Reset();

	sub_Z_NLO_muon_cross_section->Reset();
	sub_W_NLO_muon_cross_section->Reset();

	sub_Wp_NLO_muon_cross_section->Reset();
	sub_Wm_NLO_muon_cross_section->Reset();

	sub_NLO_muon_cross_section_cb->Reset();
	sub_NLO_muon_cross_section_fr->Reset();
	TFile* infile = TFile::Open("wm_NLOmuonyield.root", "READ");
	
	vector<double> *Luminosity;
	infile->GetObject("luminosity",Luminosity);

	




	//resetting sub bins (bins from pthat) to fill main histogram
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

	//DEBUGGGING: check if entries have been correctly called
	//for (int entry = 0; entry < muontuples->GetEntries(); ++entry) 
	//{
	//	muontuples->GetEntry(entry);
	//	cout << "Entry " << entry << ": ";
//
//			cout << "type=" << type << " eventNo=" << eventNo << " mother1ID=" 
//				 << mother1 << " mother2ID=" << mother2 << " pAbs=" << pAbs
//				 << " pt=" << pt << " eta=" << eta << " id=" << id << endl;
//		}

	//filling histograms
	for (int event_counter = 0; event_counter < particle_count; event_counter++)
	{

		muontuples->GetEntry(event_counter);

	//if ((eta > -4.0) && (eta < -2.5)){

		sub_TOTAL_NLO_muon_cross_section->Fill(eta);

		//check if muon is in a particular region for region plots
		if (abs(eta) < 0.9)
		{
			sub_NLO_muon_cross_section_cb->Fill(eta);
		}

		if (eta>-4.0 && eta<-2.5)
		{
			sub_NLO_muon_cross_section_fr->Fill(eta);
		}


		//check if mother is b or c meson/baryon.
		bool b_found = false;
		bool c_found = false;

		for (int i = 0; i < sizeof(b_moms)/sizeof(b_moms[0]); ++i)
		{
			if (abs(b_moms[i]) == abs(mother1))
			{
				b_found = true;
				break;
			}
		}

		for (int i = 0; i < sizeof(c_moms)/sizeof(c_moms[0]); ++i)
		{
			if (abs(c_moms[i]) == abs(mother1))
			{
				c_found = true;
				break;
			}
		}
		
		//checking the mother of the muon and filling the appropriate histogram
		if (b_found){sub_b_NLO_muon_cross_section->Fill(eta);}

		else if (c_found){sub_c_NLO_muon_cross_section->Fill(eta);}

		//checing if mother particle is a weak vector boson 

		else if ((abs(mother1) == 24) || ((abs(mother1) == 13) && (abs(mother2) == 90)))

		{
			sub_W_NLO_muon_cross_section->Fill(eta);
			
			sub_Wm_NLO_muon_cross_section->Fill(eta);
		}

		else if (!(abs(mother1) == 13) && !(abs(mother2) == 13))
		{
		sub_other_NLO_muon_cross_section->Fill(eta);
		//cout << " mom1 id: "<< mother1<< endl;
		}
		else if (((abs(mother1) == 13)) && !((abs(mother1) == 13) && (abs(mother2) == 90))  && !((abs(mother1) == 13) && (abs(mother2) == 13)))

		{
			cout << mother1 << "  " << mother2 << endl;
			sub_muon_NLO_muon_cross_section->Fill(eta);
		}
	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_TOTAL_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_NLO_muon_cross_section_cb->Scale(scalebin, "width");
	sub_NLO_muon_cross_section_fr->Scale(scalebin, "width");
	sub_b_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_c_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_other_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_W_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_Wm_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_Wp_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_Z_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_muon_NLO_muon_cross_section->Scale(scalebin, "width");
	

	TOTAL_NLO_muon_cross_section->Add(sub_TOTAL_NLO_muon_cross_section);
	total_NLO_muon_cross_section_cb->Add(sub_NLO_muon_cross_section_cb);
	total_NLO_muon_cross_section_fr->Add(sub_NLO_muon_cross_section_fr);
	b_NLO_muon_cross_section->Add(sub_b_NLO_muon_cross_section);
	c_NLO_muon_cross_section->Add(sub_c_NLO_muon_cross_section);
	other_NLO_muon_cross_section->Add(sub_other_NLO_muon_cross_section);
	W_NLO_muon_cross_section->Add(sub_W_NLO_muon_cross_section);
	Wp_NLO_muon_cross_section->Add(sub_Wp_NLO_muon_cross_section);
	Wm_NLO_muon_cross_section->Add(sub_Wm_NLO_muon_cross_section);
	Z_NLO_muon_cross_section->Add(sub_Z_NLO_muon_cross_section);
	muon_NLO_muon_cross_section->Add(sub_muon_NLO_muon_cross_section);
	}

//*******************************************************************************
//Now looking at Z  
//*******************************************************************************
	{
	sub_TOTAL_NLO_muon_cross_section->Reset();
	sub_b_NLO_muon_cross_section->Reset();
	sub_c_NLO_muon_cross_section->Reset();
	sub_other_NLO_muon_cross_section->Reset();

	sub_Z_NLO_muon_cross_section->Reset();
	sub_W_NLO_muon_cross_section->Reset();

	sub_Wp_NLO_muon_cross_section->Reset();
	sub_Wm_NLO_muon_cross_section->Reset();

	sub_NLO_muon_cross_section_cb->Reset();
	sub_NLO_muon_cross_section_fr->Reset();
	TFile* infile = TFile::Open("Z_NLOmuonyield.root", "READ");
	
	vector<double> *Luminosity;
	infile->GetObject("luminosity",Luminosity);

	




	//resetting sub bins (bins from pthat) to fill main histogram
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

	//DEBUGGGING: check if entries have been correctly called
	//for (int entry = 0; entry < muontuples->GetEntries(); ++entry) 
	//{
	//	muontuples->GetEntry(entry);
	//	cout << "Entry " << entry << ": ";
//
//			cout << "type=" << type << " eventNo=" << eventNo << " mother1ID=" 
//				 << mother1 << " mother2ID=" << mother2 << " pAbs=" << pAbs
//				 << " pt=" << pt << " eta=" << eta << " id=" << id << endl;
//		}

	//filling histograms
	for (int event_counter = 0; event_counter < particle_count; event_counter++)
	{

		muontuples->GetEntry(event_counter);
		
	//if ((eta > -4.0) && (eta < -2.5)){
		
		sub_TOTAL_NLO_muon_cross_section->Fill(eta);

		//check if muon is in a particular region for region plots
		if (abs(eta) < 0.9)
		{
			sub_NLO_muon_cross_section_cb->Fill(eta);
		}

		if (eta>-4.0 && eta<-2.5)
		{
			sub_NLO_muon_cross_section_fr->Fill(eta);
		}


		//check if mother is b or c meson/baryon.
		bool b_found = false;
		bool c_found = false;

		for (int i = 0; i < sizeof(b_moms)/sizeof(b_moms[0]); ++i)
		{
			if (abs(b_moms[i]) == abs(mother1))
			{
				b_found = true;
				break;
			}
		}

		for (int i = 0; i < sizeof(c_moms)/sizeof(c_moms[0]); ++i)
		{
			if (abs(c_moms[i]) == abs(mother1))
			{
				c_found = true;
				break;
			}
		}
		
		//checking the mother of the muon and filling the appropriate histogram
		if (b_found){sub_b_NLO_muon_cross_section->Fill(eta);}

		else if (c_found){sub_c_NLO_muon_cross_section->Fill(eta);}

		//checing if mother particle is a weak vector boson 
		else if ((abs(mother1) == 23) || ((abs(mother1) == 13) && (abs(mother2) == 90)) || ((abs(mother1) == 13) && (abs(mother2) == 13)))
		{
			sub_Z_NLO_muon_cross_section->Fill(eta);
		}

		else if (!(abs(mother1) == 13) && !(abs(mother2) == 13))
		{
		sub_other_NLO_muon_cross_section->Fill(eta);
		//cout << " mom1 id: "<< mother1<< endl;
		}
		else if (((abs(mother1) == 13)) && !((abs(mother1) == 13) && (abs(mother2) == 90))  && !((abs(mother1) == 13) && (abs(mother2) == 13)))
		{
			sub_muon_NLO_muon_cross_section->Fill(eta);
		}
	//}
	}
	
	//normalising bins for calculating the cross section
	double_t scalebin = 1./(*Luminosity)[0];

	sub_TOTAL_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_NLO_muon_cross_section_cb->Scale(scalebin, "width");
	sub_NLO_muon_cross_section_fr->Scale(scalebin, "width");
	sub_b_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_c_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_other_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_W_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_Wm_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_Wp_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_Z_NLO_muon_cross_section->Scale(scalebin, "width");
	sub_muon_NLO_muon_cross_section->Scale(scalebin, "width");
	

	TOTAL_NLO_muon_cross_section->Add(sub_TOTAL_NLO_muon_cross_section);
	total_NLO_muon_cross_section_cb->Add(sub_NLO_muon_cross_section_cb);
	total_NLO_muon_cross_section_fr->Add(sub_NLO_muon_cross_section_fr);
	b_NLO_muon_cross_section->Add(sub_b_NLO_muon_cross_section);
	c_NLO_muon_cross_section->Add(sub_c_NLO_muon_cross_section);
	other_NLO_muon_cross_section->Add(sub_other_NLO_muon_cross_section);
	W_NLO_muon_cross_section->Add(sub_W_NLO_muon_cross_section);
	Wp_NLO_muon_cross_section->Add(sub_Wp_NLO_muon_cross_section);
	Wm_NLO_muon_cross_section->Add(sub_Wm_NLO_muon_cross_section);
	Z_NLO_muon_cross_section->Add(sub_Z_NLO_muon_cross_section);
	muon_NLO_muon_cross_section->Add(sub_muon_NLO_muon_cross_section);

	TOTAL_NLO_muon_cross_section->Add(TOTAL_muon_cross_section);
	}





	TFile* outFile =new TFile("NLOmuonyieldmacro.root", "RECREATE");

	//plotting muons from different sources
	c3->cd();
	c3->SetLogy();
	c3->SetGridy();

	TOTAL_NLO_muon_cross_section->GetYaxis()->SetTitle("d#sigma/d#eta (pb/GeV/c)");
	TOTAL_NLO_muon_cross_section->GetXaxis()->SetTitle("#eta (GeV/c)");

	TOTAL_NLO_muon_cross_section->Draw("SAME");
//	c_NLO_muon_cross_section->Draw("SAME");
	//c_muon_cross_section ->Draw("SAME");
	Z_NLO_muon_cross_section->Draw("SAME");
	W_NLO_muon_cross_section->Draw("SAME");
	Wm_NLO_muon_cross_section->Draw("SAME");
	Wp_NLO_muon_cross_section->Draw("SAME");
//	b_NLO_muon_cross_section->Draw("SAME");
	//b_muon_cross_section ->Draw("SAME");
//	other_NLO_muon_cross_section->Draw("SAME");
//	other_muon_cross_section ->Draw("SAME");
//	muon_NLO_muon_cross_section->Draw("SAME");

	TOTAL_NLO_muon_cross_section->SetLineColor(kBlack);
//	b_NLO_muon_cross_section->SetLineColor(kRed);
	//b_muon_cross_section->SetLineColor(kRed);
//	c_NLO_muon_cross_section->SetLineColor(kGreen);
	//c_muon_cross_section->SetLineColor(kGreen);
	W_NLO_muon_cross_section->SetLineColor(kBlue);
	Wm_NLO_muon_cross_section->SetLineColor(7);
	Wp_NLO_muon_cross_section->SetLineColor(8);
	Z_NLO_muon_cross_section->SetLineColor(kYellow);
//	other_NLO_muon_cross_section->SetLineColor(kMagenta);
//	muon_NLO_muon_cross_section->SetLineColor(20);

//	b_NLO_muon_cross_section->SetMarkerStyle(25);
//	b_NLO_muon_cross_section->SetMarkerColor(kRed);
	//b_muon_cross_section->SetMarkerStyle(25);
	//b_muon_cross_section->SetMarkerColor(kRed);
//	c_NLO_muon_cross_section->SetMarkerStyle(26);
//	c_NLO_muon_cross_section->SetMarkerColor(kGreen);
	//c_muon_cross_section->SetMarkerStyle(26);
	//c_muon_cross_section->SetMarkerColor(kGreen);
	W_NLO_muon_cross_section->SetMarkerStyle(28);
	W_NLO_muon_cross_section->SetMarkerColor(kBlue);
	Z_NLO_muon_cross_section->SetMarkerStyle(23);
	Z_NLO_muon_cross_section->SetMarkerColor(kYellow);
	TOTAL_NLO_muon_cross_section->SetMarkerStyle(29);
	TOTAL_NLO_muon_cross_section->SetMarkerColor(kBlack);
	Wm_NLO_muon_cross_section->SetMarkerStyle(20);
	Wm_NLO_muon_cross_section->SetMarkerColor(7);
	Wp_NLO_muon_cross_section->SetMarkerStyle(24);
	Wp_NLO_muon_cross_section->SetMarkerColor(8);
//	muon_NLO_muon_cross_section->SetMarkerStyle(39);
//	muon_NLO_muon_cross_section->SetMarkerColor(20);



	TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
	leg2->AddEntry(TOTAL_NLO_muon_cross_section,"Total", "lep");
//	leg2->AddEntry(b_NLO_muon_cross_section,"W/Z/#gamma*(probably) #rightarrow b #rightarrow #mu", "lep");
//	leg2->AddEntry(c_NLO_muon_cross_section,"W/Z/#gamma*(probably) #rightarrow c #rightarrow #mu", "lep");
//	leg2->AddEntry(b_muon_cross_section,"b #rightarrow #mu", "lep");
//	leg2->AddEntry(c_muon_cross_section,"c #rightarrow #mu", "lep");
	leg2->AddEntry(Z_NLO_muon_cross_section, "Z/#gamma* #rightarrow #mu", "lep");
	leg2->AddEntry(W_NLO_muon_cross_section, "W #rightarrow #mu", "lep");
	leg2->AddEntry(Wp_NLO_muon_cross_section, "W+ #rightarrow #mu", "lep");
	leg2->AddEntry(Wm_NLO_muon_cross_section, "W- #rightarrow #mu", "lep");
//	leg2->AddEntry(other_NLO_muon_cross_section,"other #rightarrow #mu", "lep");
//	leg2->AddEntry(muon_NLO_muon_cross_section,"#mu #rightarrow #mu + X", "lep");
	leg2->Draw("SAME");

	TOTAL_NLO_muon_cross_section->SetMinimum(1);
	TOTAL_NLO_muon_cross_section->SetMaximum(5000000000);
	c3->Write();


	delete outFile;
}
