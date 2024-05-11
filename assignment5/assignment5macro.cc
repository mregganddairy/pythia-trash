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

TCanvas *c1 = new TCanvas(); //canvas for eta distribution
TCanvas *c2 = new TCanvas(); // canvas for pseudorapidity pt distribution 
TCanvas *c3 = new TCanvas(); // canvas for rapidity pt distribution  
void assignment4macro(){
	//pthat bins

//defining bins to seperate soft and hard qcd using pthat
static const int nbins =6;
static const double binedges[nbins+1] = {0., 14., 30., 50., 75., 100., 150. };

//creating lists of IDs of mother particles of interest
//static const int b_moms[] = {511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 20523, 515, 525, 531, 533, 535, 541, 543,545, //mesons
//							5122, 5112, 5212, 5222, 5114, 5214, 5224, 5132, 5232, 5312, 5322, 5314, 5324, 5332, 5334}; //baryons
//static const int c_moms[] = {411, 421, 10411, 10421, 413, 423, 10413, 10423, 20413, 20423, 415, 425, 431, 10431, 433, 10433, 20433, 435, //mesons
//							4122, 4222, 4212, 4112, 4224, 4214, 4114, 4232, 4132, 4322, 4312, 4324, 4314, 4332, 4334, 4422}; \mesons
	
	//initialise histograms
	

	//for eta distribution
	TH1F* sub_muon_cross_section = new TH1F("sub_muon_cross_section","", 60, -15, 15);
	TH1F* total_muon_cross_section = new TH1F("total_muon_cross_section","pseudorapidity distribution", 60, -15, 15);


	//for y distribution
	TH1F* sub_muon_cross_section_y = new TH1F("sub_muon_cross_section_y","", 60, -15, 15);
	TH1F* total_muon_cross_section_y = new TH1F("total_muon_cross_section_y","rapidity distribution", 60, -15, 15);
	
	//for pt distribution in central barrel
	TH1F* sub_muon_cross_section_cb = new TH1F("sub_muon_cross_section_cb","", 100, 0, 100);
	TH1F* total_muon_cross_section_cb = new TH1F("total_muon_cross_section_cb","", 100, 0, 100);


	//for pt distribution in forward region 
	TH1F* sub_muon_cross_section_fr = new TH1F("sub_muon_cross_section_fr","", 100, 0, 100);
	TH1F* total_muon_cross_section_fr = new TH1F("total_muon_cross_section_fr","", 100, 0, 100);


	//muon cross section produced from bottom and charm
//	TH1F* b_muon_cross_section = new TH1F("b_muon_cross_section","", 50,  0, 100);
//	TH1F* c_muon_cross_section = new TH1F("c_muon_cross_section","", 50,  0, 100);



	TFile* infile = TFile::Open("muonyield.root", "READ");
	
	vector<double> *binLuminosity;
	infile->GetObject("luminosity",binLuminosity);

	for (int ibin = 0; ibin < nbins; ibin++)
	{

		//resetting sub bins (bins from pthat) to fill main histogram
		sub_muon_cross_section->Reset();
		sub_muon_cross_section_y->Reset();
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

			sub_muon_cross_section->Fill(eta);

			sub_muon_cross_section_y->Fill(y);
			
			if (abs(eta) < 0.9)
			{
				sub_muon_cross_section_cb->Fill(pt);
			}

			if (eta>-4.0 && eta<-2.5)
			{
				sub_muon_cross_section_fr->Fill(pt);
			}
		}
		
		//normalising bins for calculating the cross section
		double_t scalebin = 1./(*binLuminosity)[ibin];
		sub_muon_cross_section->Scale(scalebin, "width");
		sub_muon_cross_section_y->Scale(scalebin, "width");
		sub_muon_cross_section_cb->Scale(scalebin, "width");
		sub_muon_cross_section_fr->Scale(scalebin, "width");

		total_muon_cross_section->Add(sub_muon_cross_section);
		total_muon_cross_section_y->Add(sub_muon_cross_section_y);
		total_muon_cross_section_cb->Add(sub_muon_cross_section_cb);
		total_muon_cross_section_fr->Add(sub_muon_cross_section_fr);
	}	
	TFile* outFile =new TFile("muonyieldmacro.root", "RECREATE");


	c1->cd();
//	total_muon_cross_section->SetMinimum(0.);

	total_muon_cross_section->GetXaxis()->SetTitle("#eta");
	total_muon_cross_section->GetYaxis()->SetTitle("d#sigma/d#eta (pb)");

	total_muon_cross_section->SetLineColor(kBlue);
	total_muon_cross_section->Draw("SAME");

	total_muon_cross_section_y->SetLineColor(kRed);
	total_muon_cross_section_y->Draw("SAME");	

	TLegend *leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
	leg1->AddEntry(total_muon_cross_section,"pseudorapidity ", "l");
	leg1->AddEntry(total_muon_cross_section_y, "rapidity", "l");
	leg1->Draw("SAME");

	c1->Write();

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
	leg->AddEntry(total_muon_cross_section_cb,"total_muon_cross_section_cb", "l");
	leg->AddEntry(total_muon_cross_section_fr,"total_muon_cross_section_fr", "l");
	leg->Draw("SAME");

	c2->Write();
	delete outFile;
}
