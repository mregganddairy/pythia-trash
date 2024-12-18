//Comparison of PYTHIA and POWHEG calculations in electroweak processes
//Could very much be cleaned up.

#include "TFile.h"
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;


TCanvas *c4 = new TCanvas(); //canvas for Wp comparison plots
TCanvas *c5 = new TCanvas(); //canvas for Wm comparison plots
TCanvas *c6 = new TCanvas(); //canvas for Z comparison plots

							 
void assignment6macro2(){
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
	TH1F* sub_TOTAL_NLO_muon_cross_section = new TH1F("sub_TOTAL_NLO_muon_cross_section", "", 50,0,100);
	TH1F* TOTAL_NLO_muon_cross_section = new TH1F("TOTAL_NLO_muon_cross_section", "NLO muon differential cross section different mothers (mainly W processes)", 50,0,100);

	//for pt distribution in central barrel
	TH1F* sub_NLO_muon_cross_section_cb = new TH1F("sub_NLO_muon_cross_section_cb","", 100, 0, 100);
	TH1F* total_NLO_muon_cross_section_cb = new TH1F("total_NLO_muon_cross_section_cb","", 100, 0, 100);

	//for pt distribution in forward region 
	TH1F* sub_NLO_muon_cross_section_fr = new TH1F("sub_NLO_muon_cross_section_fr","", 100, 0, 100);
	TH1F* total_NLO_muon_cross_section_fr = new TH1F("total_NLO_muon_cross_section_fr","", 100, 0, 100);


	//muon cross section produced from bottom or charm
	TH1F* sub_b_NLO_muon_cross_section = new TH1F("sub_b_NLO_muon_cross_section","", 50,  0, 100);
	TH1F* b_NLO_muon_cross_section = new TH1F("b_NLO_muon_cross_section","", 50,  0, 100);
	TH1F* sub_c_NLO_muon_cross_section = new TH1F("sub_c_NLO_muon_cross_section","", 50,  0, 100);
	TH1F* c_NLO_muon_cross_section = new TH1F("c_NLO_muon_cross_section","", 50,  0, 100);

	TH1F* sub_Wm_NLO_muon_cross_section = new TH1F("sub_Wm_NLO_muon_cross_section","", 50,  0, 100);
	TH1F* Wm_NLO_muon_cross_section = new TH1F("Wm_NLO_muon_cross_section","", 50,  0, 100);

	TH1F* sub_Wp_NLO_muon_cross_section = new TH1F("sub_Wp_NLO_muon_cross_section","", 50,  0, 100);
	TH1F* Wp_NLO_muon_cross_section = new TH1F("Wp_NLO_muon_cross_section","", 50,  0, 100);

	TH1F* sub_W_NLO_muon_cross_section = new TH1F("sub_W_NLO_muon_cross_section","", 50,  0, 100);
	TH1F* W_NLO_muon_cross_section = new TH1F("W_NLO_muon_cross_section","", 50,  0, 100);

	TH1F* sub_Z_NLO_muon_cross_section = new TH1F("sub_Z_NLO_muon_cross_section","", 50,  0, 100);
	TH1F* Z_NLO_muon_cross_section = new TH1F("Z_NLO_muon_cross_section","", 50,  0, 100);

	//muons coming from muons
	TH1F* sub_muon_NLO_muon_cross_section = new TH1F("sub_muon_NLO_muon_cross_section","", 50,  0, 100);
	TH1F* muon_NLO_muon_cross_section = new TH1F("muon_NLO_muon_cross_section","", 50,  0, 100);

	//muon cross sections which don't come from bottom or charm
	TH1F* sub_other_NLO_muon_cross_section = new TH1F("sub_other_NLO_muon_cross_section","", 50,  0, 100);
	TH1F* other_NLO_muon_cross_section = new TH1F("other_NLO_muon_cross_section","", 50,  0, 100);


//*******************************************************************************
//LO stuff  
//*******************************************************************************
	TH1F* sub_TOTAL_muon_cross_section = new TH1F("sub_TOTAL_muon_cross_section", "", 50,0,100);
	TH1F* TOTAL_muon_cross_section = new TH1F("TOTAL_muon_cross_section", "LO differential cross section from different mother particles (weak processes)", 50,0,100);

	//for pt distribution in central barrel
	TH1F* sub_muon_cross_section_cb = new TH1F("sub_muon_cross_section_cb","", 100, 0, 100);
	TH1F* total_muon_cross_section_cb = new TH1F("total_muon_cross_section_cb","", 100, 0, 100);

	//for pt distribution in forward region 
	TH1F* sub_muon_cross_section_fr = new TH1F("sub_muon_cross_section_fr","", 100, 0, 100);
	TH1F* total_muon_cross_section_fr = new TH1F("total_muon_cross_section_fr","", 100, 0, 100);


	//muon cross section produced from bottom or charm
	TH1F* sub_b_muon_cross_section = new TH1F("sub_b_muon_cross_section","", 50,  0, 100);
	TH1F* b_muon_cross_section = new TH1F("b_muon_cross_section","", 50,  0, 100);
	TH1F* sub_c_muon_cross_section = new TH1F("sub_c_muon_cross_section","", 50,  0, 100);
	TH1F* c_muon_cross_section = new TH1F("c_muon_cross_section","", 50,  0, 100);

	TH1F* sub_Wm_muon_cross_section = new TH1F("sub_Wm_muon_cross_section","", 50,  0, 100);
	TH1F* Wm_muon_cross_section = new TH1F("Wm_muon_cross_section","", 50,  0, 100);

	TH1F* sub_Wp_muon_cross_section = new TH1F("sub_Wp_muon_cross_section","", 50,  0, 100);
	TH1F* Wp_muon_cross_section = new TH1F("Wp_muon_cross_section","", 50,  0, 100);

	TH1F* sub_W_muon_cross_section = new TH1F("sub_W_muon_cross_section","", 50,  0, 100);
	TH1F* W_muon_cross_section = new TH1F("W_muon_cross_section","", 50,  0, 100);

	TH1F* sub_Z_muon_cross_section = new TH1F("sub_Z_muon_cross_section","", 50,  0, 100);
	TH1F* Z_muon_cross_section = new TH1F("Z_muon_cross_section","", 50,  0, 100);

	//muons coming from muons
	TH1F* sub_muon_muon_cross_section = new TH1F("sub_muon_muon_cross_section","", 50,  0, 100);
	TH1F* muon_muon_cross_section = new TH1F("muon_muon_cross_section","", 50,  0, 100);

	TH1F* sub_muon_annihl_cross_section = new TH1F("sub_muon_annihl_cross_section","", 50,  0, 100);
	TH1F* muon_annihl_cross_section = new TH1F("muon_annihl_cross_section","", 50,  0, 100);

	TH1F* sub_muon_kak_cross_section = new TH1F("sub_muon_kak_cross_section","", 50,  0, 100);
	TH1F* muon_kak_cross_section = new TH1F("muon_kak_cross_section","", 50,  0, 100);

	//muon cross sections which don't come from bottom or charm
	TH1F* sub_other_muon_cross_section = new TH1F("sub_other_muon_cross_section","", 50,  0, 100);
	TH1F* other_muon_cross_section = new TH1F("other_muon_cross_section","", 50,  0, 100);


	{
	static const int nbins =7;
	static const double binedges[nbins+1] = {0., 14.,20., 30., 50., 75., 100., 150. };

	TFile* infile = TFile::Open("/home/Josh/physics/pythia/assignment5/muonyield.root", "READ");
	
	vector<double> *binLuminosity;
	infile->GetObject("luminosity",binLuminosity);


	for (int ibin = 0; ibin < nbins; ibin++)
	{

		//resetting sub bins (bins from pthat) to fill main histogram
		sub_TOTAL_muon_cross_section->Reset();
		sub_b_muon_cross_section->Reset();
		sub_c_muon_cross_section->Reset();
		sub_other_muon_cross_section->Reset();

		sub_Z_muon_cross_section->Reset();
		sub_W_muon_cross_section->Reset();

		sub_Wp_muon_cross_section->Reset();
		sub_Wm_muon_cross_section->Reset();

		sub_muon_cross_section_cb->Reset();
		sub_muon_cross_section_fr->Reset();

		sub_muon_muon_cross_section->Reset();
		sub_muon_annihl_cross_section->Reset();
		sub_muon_kak_cross_section->Reset();

		//reading files from pythia calculation
		TNtuple *muontuples = (TNtuple*)infile->Get(Form("muon%d", ibin));
		Float_t type, eventNo,index, status, mother1, mother2, mother11, mother111, pAbs, pt, y, eta, id;	
		int particle_count = muontuples->GetEntries();	//number of muons
		muontuples->SetBranchAddress("binNo", &type);
		muontuples->SetBranchAddress("eventNo", &eventNo);
		muontuples->SetBranchAddress("index", &index);
		muontuples->SetBranchAddress("status", &status);
		muontuples->SetBranchAddress("mother1", &mother1);
		muontuples->SetBranchAddress("mother2",&mother2);
		muontuples->SetBranchAddress("mother11", &mother11);
		muontuples->SetBranchAddress("mother111",&mother111);
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

			sub_TOTAL_muon_cross_section->Fill(pt);

			//check if muon is in a particular region for region plots
			if (abs(eta) < 0.9)
			{
				sub_muon_cross_section_cb->Fill(pt);
			}

			if (eta>-4.0 && eta<-2.5)
			{
				sub_muon_cross_section_fr->Fill(pt);
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
			if (b_found){sub_b_muon_cross_section->Fill(pt);
			}

			else if (c_found){sub_c_muon_cross_section->Fill(pt);}

			//checing if mother particle is a weak vector boson 

			else if (abs(mother1) == 24 || (abs(mother11) ==24) || (abs(mother111) ==24))
			{
				sub_W_muon_cross_section->Fill(pt);

				if (mother1 == 24 || (mother11 ==24) || (mother111 ==24))
				{
				sub_Wp_muon_cross_section->Fill(pt);
				}

				else if (mother1 == -24 || (mother11 ==-24) || (mother111 ==-24))
				{
				sub_Wm_muon_cross_section->Fill(pt);
				}


			}
			else if(abs(mother1) == 23 || (abs(mother11) ==23) || (abs(mother111) ==23))
			{
				sub_Z_muon_cross_section->Fill(pt);
			}

			else if (!(abs(mother1) == 13) && !(abs(mother2) == 13))
			{
			sub_other_muon_cross_section->Fill(pt);
	//		cout << " mom1 id: "<< mother1<< endl;
			}
			else if ((abs(mother1) == 13) && !(abs(mother11) ==23 || abs(mother11)==24) && !(abs(mother111) ==23 || abs(mother111)==24))
			{
				sub_muon_muon_cross_section->Fill(pt);
				//cout << "muon->muon moms:" << mother1 << "and" << mother2 << endl;
				if ((abs(mother2) == 13))
				{
					sub_muon_annihl_cross_section->Fill(pt);
				}
			
				else if (abs(mother2) == 90)
				{
					sub_muon_kak_cross_section->Fill(pt);
					//cout << "status: "<< status<<endl;
				}
			}

		}
		
		//normalising bins for calculating the cross section
		double_t scalebin = 1./(*binLuminosity)[ibin];

		sub_TOTAL_muon_cross_section->Scale(scalebin, "width");
		sub_muon_cross_section_cb->Scale(scalebin, "width");
		sub_muon_cross_section_fr->Scale(scalebin, "width");
		sub_b_muon_cross_section->Scale(scalebin, "width");
		sub_c_muon_cross_section->Scale(scalebin, "width");
		sub_other_muon_cross_section->Scale(scalebin, "width");
		sub_W_muon_cross_section->Scale(scalebin, "width");
		sub_Wm_muon_cross_section->Scale(scalebin, "width");
		sub_Wp_muon_cross_section->Scale(scalebin, "width");
		sub_Z_muon_cross_section->Scale(scalebin, "width");
		sub_muon_muon_cross_section->Scale(scalebin, "width");

		sub_muon_kak_cross_section->Scale(scalebin, "width");
		sub_muon_annihl_cross_section->Scale(scalebin, "width");
		

		TOTAL_muon_cross_section->Add(sub_TOTAL_muon_cross_section);
		total_muon_cross_section_cb->Add(sub_muon_cross_section_cb);
		total_muon_cross_section_fr->Add(sub_muon_cross_section_fr);
		b_muon_cross_section->Add(sub_b_muon_cross_section);
		c_muon_cross_section->Add(sub_c_muon_cross_section);
		other_muon_cross_section->Add(sub_other_muon_cross_section);
		W_muon_cross_section->Add(sub_W_muon_cross_section);
		Wp_muon_cross_section->Add(sub_Wp_muon_cross_section);
		Wm_muon_cross_section->Add(sub_Wm_muon_cross_section);
		Z_muon_cross_section->Add(sub_Z_muon_cross_section);
		muon_muon_cross_section->Add(sub_muon_muon_cross_section);

		muon_kak_cross_section->Add(sub_muon_kak_cross_section);
		muon_annihl_cross_section->Add(sub_muon_annihl_cross_section);
	}	

	}

//*******************************************************************************
//Now looking at W+  
//*******************************************************************************
	
	{
	TFile* infile = TFile::Open("wp_NLOmuonyield.root", "READ");
	
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
		
		sub_TOTAL_NLO_muon_cross_section->Fill(pt);

		//check if muon is in a particular region for region plots
		if (abs(eta) < 0.9)
		{
			sub_NLO_muon_cross_section_cb->Fill(pt);
		}

		if (eta>-4.0 && eta<-2.5)
		{
			sub_NLO_muon_cross_section_fr->Fill(pt);
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
		if (b_found){sub_b_NLO_muon_cross_section->Fill(pt);}

		else if (c_found){sub_c_NLO_muon_cross_section->Fill(pt);}

		//checing if mother particle is a weak vector boson 

		else if ((abs(mother1) == 24) || ((abs(mother1) == 13) && (abs(mother2) == 90)))
		{
			sub_W_NLO_muon_cross_section->Fill(pt);

			sub_Wp_NLO_muon_cross_section->Fill(pt);
		}

		else if (!(abs(mother1) == 13) && !(abs(mother2) == 13))
		{
		sub_other_NLO_muon_cross_section->Fill(pt);
		//cout << " mom1 id: "<< mother1<< endl;
		}
		else if (((abs(mother1) == 13)) && !((abs(mother1) == 13) && (abs(mother2) == 90)))
		{
			sub_muon_NLO_muon_cross_section->Fill(pt);
		}
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

		sub_TOTAL_NLO_muon_cross_section->Fill(pt);

		//check if muon is in a particular region for region plots
		if (abs(eta) < 0.9)
		{
			sub_NLO_muon_cross_section_cb->Fill(pt);
		}

		if (eta>-4.0 && eta<-2.5)
		{
			sub_NLO_muon_cross_section_fr->Fill(pt);
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
		if (b_found){sub_b_NLO_muon_cross_section->Fill(pt);}

		else if (c_found){sub_c_NLO_muon_cross_section->Fill(pt);}

		//checing if mother particle is a weak vector boson 

		else if ((abs(mother1) == 24) || ((abs(mother1) == 13) && (abs(mother2) == 90)))

		{
			sub_W_NLO_muon_cross_section->Fill(pt);
			
			sub_Wm_NLO_muon_cross_section->Fill(pt);
		}

		else if (!(abs(mother1) == 13) && !(abs(mother2) == 13))
		{
		sub_other_NLO_muon_cross_section->Fill(pt);
		//cout << " mom1 id: "<< mother1<< endl;
		}
		else if (((abs(mother1) == 13)) && !((abs(mother1) == 13) && (abs(mother2) == 90))  && !((abs(mother1) == 13) && (abs(mother2) == 13)))

		{
			cout << mother1 << "  " << mother2 << endl;
			sub_muon_NLO_muon_cross_section->Fill(pt);
		}
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
		
		
		sub_TOTAL_NLO_muon_cross_section->Fill(pt);

		//check if muon is in a particular region for region plots
		if (abs(eta) < 0.9)
		{
			sub_NLO_muon_cross_section_cb->Fill(pt);
		}

		if (eta>-4.0 && eta<-2.5)
		{
			sub_NLO_muon_cross_section_fr->Fill(pt);
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
		if (b_found){sub_b_NLO_muon_cross_section->Fill(pt);}

		else if (c_found){sub_c_NLO_muon_cross_section->Fill(pt);}

		//checing if mother particle is a weak vector boson 
		else if ((abs(mother1) == 23) || ((abs(mother1) == 13) && (abs(mother2) == 90)) || ((abs(mother1) == 13) && (abs(mother2) == 13)))
		{
			sub_Z_NLO_muon_cross_section->Fill(pt);
		}

		else if (!(abs(mother1) == 13) && !(abs(mother2) == 13))
		{
		sub_other_NLO_muon_cross_section->Fill(pt);
		//cout << " mom1 id: "<< mother1<< endl;
		}
		else if (((abs(mother1) == 13)) && !((abs(mother1) == 13) && (abs(mother2) == 90))  && !((abs(mother1) == 13) && (abs(mother2) == 13)))
		{
			sub_muon_NLO_muon_cross_section->Fill(pt);
		}
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





	TFile* outFile =new TFile("NLOmuonyieldmacro.root", "RECREATE");

	c4->cd();
	
	auto Rm = new TRatioPlot(Wm_muon_cross_section, Wm_NLO_muon_cross_section);

	Rm->Draw();  

	// Customize the TRatioPlot appearance (optional)
	Rm->GetUpperRefYaxis()->SetTitle("d#sigma/dpt (pb/GeV/c)");
	Rm->GetLowerRefYaxis()->SetTitle("Ratio LO/NLO");
	Rm->GetLowerRefXaxis()->SetTitle("pT (GeV)");
	Wm_muon_cross_section->SetTitle("W- Cross Section Comparison");


// Set marker styles 
	Rm->GetUpperPad()->cd();  // Access upper pad
	Wm_NLO_muon_cross_section->SetLineColor(kRed);    
	Wm_NLO_muon_cross_section->SetLineStyle(1);       
//	Wm_NLO_muon_cross_section->SetMarkerStyle(25);    
//	Wm_NLO_muon_cross_section->SetMarkerColor(kRed);    

	Wm_muon_cross_section->SetLineColor(kBlue);       
	Wm_muon_cross_section->SetLineStyle(1);           
//	Wm_muon_cross_section->SetMarkerStyle(24);        
//	Wm_muon_cross_section->SetMarkerColor(kBlue);    





	// Draw the legend on the upper pad
	TLegend *leg3 = new TLegend(0.6, 0.7, 0.9, 0.9);  // Position it as needed
	leg3->AddEntry(Wm_NLO_muon_cross_section, "NLO", "lep");
	leg3->AddEntry(Wm_muon_cross_section, "LO", "lep");
	Rm->GetUpperPad()->cd();  // Make sure you're in the upper pad
	leg3->Draw();

	// Update the canvas to reflect changes
	c4->Update();







c5->cd();
	
	auto Rp = new TRatioPlot(Wp_muon_cross_section, Wp_NLO_muon_cross_section);

	Rp->Draw();  

	// Customize the TRatioPlot appearance (optional)
	Rp->GetUpperRefYaxis()->SetTitle("d#sigma/dpt (pb/GeV/c)");
	Rp->GetLowerRefYaxis()->SetTitle("Ratio LO/NLO");
	Rp->GetLowerRefXaxis()->SetTitle("pT (GeV)");
	Wp_muon_cross_section->SetTitle("W+ Cross Section Comparison");


// Set marker styles 
	Rp->GetUpperPad()->cd();  // Access upper pad
	Wp_NLO_muon_cross_section->SetLineColor(kRed);    
	Wp_NLO_muon_cross_section->SetLineStyle(1);       
//	Wp_NLO_muon_cross_section->SetMarkerStyle(25);    
//	Wp_NLO_muon_cross_section->SetMarkerColor(kRed);    

	Wp_muon_cross_section->SetLineColor(kBlue);       
	Wp_muon_cross_section->SetLineStyle(1);           
//	Wp_muon_cross_section->SetMarkerStyle(24);        
//	Wp_muon_cross_section->SetMarkerColor(kBlue);    





	// Draw the legend on the upper pad
	TLegend *leg4 = new TLegend(0.6, 0.7, 0.9, 0.9);  // Position it as needed
	leg4->AddEntry(Wp_NLO_muon_cross_section, "NLO", "lep");
	leg4->AddEntry(Wp_muon_cross_section, "LO", "lep");
	Rp->GetUpperPad()->cd();  // Make sure you're in the upper pad
	leg4->Draw();

	// Update the canvas to reflect changes
	c5->Update();







c6->cd();
	
	auto Rz = new TRatioPlot(Z_muon_cross_section, Z_NLO_muon_cross_section);

	Rz->Draw();  

	// Customize the TRatioPlot appearance (optional)
	Rz->GetUpperRefYaxis()->SetTitle("d#sigma/dpt (pb/GeV/c)");
	Rz->GetLowerRefYaxis()->SetTitle("Ratio LO/NLO");
	Rz->GetLowerRefXaxis()->SetTitle("pT (GeV)");
	Z_muon_cross_section->SetTitle("Z Cross Section Comparison");


// Set marker styles 
	Rz->GetUpperPad()->cd();  // Access upper pad
	Z_NLO_muon_cross_section->SetLineColor(kRed);    
	Z_NLO_muon_cross_section->SetLineStyle(1);       
//	Z_NLO_muon_cross_section->SetMarkerStyle(25);    
//	Z_NLO_muon_cross_section->SetMarkerColor(kRed);    

	Z_muon_cross_section->SetLineColor(kBlue);       
	Z_muon_cross_section->SetLineStyle(1);           
//	Z_muon_cross_section->SetMarkerStyle(24);        
//	Z_muon_cross_section->SetMarkerColor(kBlue);    





	// Draw the legend on the upper pad
	TLegend *leg5 = new TLegend(0.6, 0.7, 0.9, 0.9);  // Position it as needed
	leg5->AddEntry(Wp_NLO_muon_cross_section, "NLO", "lep");
	leg5->AddEntry(Wp_muon_cross_section, "LO", "lep");
	Rz->GetUpperPad()->cd();  // Make sure you're in the upper pad
	leg5->Draw();

	// Update the canvas to reflect changes
	c6->Update();

	delete outFile;
}
