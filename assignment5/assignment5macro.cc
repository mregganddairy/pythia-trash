//total eta distribution of muons and also distribution from different mothers 

#include "TFile.h"
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;

TCanvas *c2 = new TCanvas(); //canvas for pt distribution of muons from electroweak particles 
TCanvas *c3 = new TCanvas(); //canvas for rapidity pt distribution  
void assignment4macro(){
	//pthat bins

//defining bins to seperate soft and hard qcd using pthat
static const int nbins =6;
static const double binedges[nbins+1] = {0., 14., 30., 50., 75., 100., 150. };

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
	TH1F* sub_TOTAL_muon_cross_section = new TH1F("sub_TOTAL_muon_cross_section", "", 50,0,100);
	TH1F* TOTAL_muon_cross_section = new TH1F("TOTAL_muon_cross_section", "cross section from different mothers", 50,0,100);

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

	//muon cross sections which don't come from bottom or charm
	TH1F* sub_other_muon_cross_section = new TH1F("sub_other_muon_cross_section","", 50,  0, 100);
	TH1F* other_muon_cross_section = new TH1F("other_muon_cross_section","", 50,  0, 100);

	TFile* infile = TFile::Open("muonyield.root", "READ");
	
	vector<double> *binLuminosity;
	infile->GetObject("luminosity",binLuminosity);

	for (int ibin = 0; ibin < nbins; ibin++)
	{

		//resetting sub bins (bins from pthat) to fill main histogram
		sub_TOTAL_muon_cross_section->Reset();
		sub_b_muon_cross_section->Reset();
		sub_c_muon_cross_section->Reset();
		sub_other_muon_cross_section->Reset();


		sub_muon_cross_section_cb->Reset();
		sub_muon_cross_section_fr->Reset();

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
		if (b_found){sub_b_muon_cross_section->Fill(pt);}

		else if (c_found){sub_c_muon_cross_section->Fill(pt);}

		else
		{
			sub_other_muon_cross_section->Fill(pt);
			cout << " mom1 id: "<< mother1<< endl;
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

		TOTAL_muon_cross_section->Add(sub_TOTAL_muon_cross_section);
		total_muon_cross_section_cb->Add(sub_muon_cross_section_cb);
		total_muon_cross_section_fr->Add(sub_muon_cross_section_fr);
		b_muon_cross_section->Add(sub_b_muon_cross_section);
		c_muon_cross_section->Add(sub_c_muon_cross_section);
		other_muon_cross_section->Add(sub_other_muon_cross_section);
	}	
	TFile* outFile =new TFile("muonyieldmacro.root", "RECREATE");



	//checking contributions in different regions
	c2->cd();
//	total_muon_cross_section_cb->SetMinimum(0.);
	total_muon_cross_section_fr->GetXaxis()->SetTitle("pt (GeV)");
	total_muon_cross_section_fr->GetYaxis()->SetTitle("d#sigma/dpt (pb/GeV/c)");
	total_muon_cross_section_fr->SetLineColor(kRed);
	total_muon_cross_section_fr->Draw("SAME");

//	total_muon_cross_section_cb->SetMinimum(0.);
	total_muon_cross_section_cb->GetXaxis()->SetTitle("pt (GeV)");
	total_muon_cross_section_cb->SetLineColor(kBlue);
	total_muon_cross_section_cb->GetYaxis()->SetTitle("d#sigma/dpt (pb/GeV/c)");
	total_muon_cross_section_cb->Draw("SAME");


	TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
	leg->AddEntry(total_muon_cross_section_cb,"total_muon_cross_section_cb", "lep");
	leg->AddEntry(total_muon_cross_section_fr,"total_muon_cross_section_fr", "lep");
	leg->Draw("SAME");

	c2->Write();
	


	//plotting muons from different sources
	c3->cd();

	TOTAL_muon_cross_section->GetYaxis()->SetTitle("d#sigma/dpt (pb/GeV/c)");
	TOTAL_muon_cross_section->GetXaxis()->SetTitle("pt (GeV)");

	TOTAL_muon_cross_section->Draw("SAME");
	b_muon_cross_section->Draw("SAME");
	c_muon_cross_section->Draw("SAME");
	other_muon_cross_section->Draw("SAME");

	TOTAL_muon_cross_section->SetLineColor(kBlack);
	b_muon_cross_section->SetLineColor(kBlue);
	c_muon_cross_section->SetLineColor(kRed);
	other_muon_cross_section->SetLineColor(kGreen);

	TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
	leg2->AddEntry(TOTAL_muon_cross_section,"Total", "lep");
	leg2->AddEntry(b_muon_cross_section,"bottom hadron moms", "lep");
	leg2->AddEntry(c_muon_cross_section,"charm hadron moms", "lep");
	leg2->AddEntry(other_muon_cross_section,"other moms", "lep");
	leg2->Draw("SAME");

	c3->Write();


	delete outFile;
}
