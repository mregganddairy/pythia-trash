#include "TFile.h"
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;

TCanvas *c1 = new TCanvas();
void assignment3macro(){
	//pthat bins
	static const int nbins = 2;
	static const double binedges[nbins+1] = {0., 14., 150.};

	//initialise histograms
	TH1F* muon_yield = new TH1F("muon_yield","muon yield; #hat{p}_{T} (GeV/c; #frac{dN}{dp_{T}}", 150, 0, 15);
	TH1F* hard_muon_yield = new TH1F("hard_muon_yield","", 150, 0, 15);
	TH1F* soft_muon_yield = new TH1F("soft_muon_yield","", 150, 0, 15);


	TFile* infile = TFile::Open("soft_and_hard_qcd.root", "READ");

	for (int ibin = 0; ibin < nbins; ibin++)
	{

		//reading files from pythia calculation
		//

		TNtuple *muontuples = (TNtuple*)infile->Get(Form("muon%d", ibin));
		Float_t type, eventNo, pAbs, pt, eta, id;	
		int particle_count = muontuples->GetEntries();	//number of muons
		muontuples->SetBranchAddress("type", &type);
		muontuples->SetBranchAddress("eventNo", &eventNo);
		muontuples->SetBranchAddress("pAbs", &pAbs);
		muontuples->SetBranchAddress("pt", &pt);
		muontuples->SetBranchAddress("eta", &eta);
		muontuples->SetBranchAddress("id", &id);

		//DEBUGGGING: check if entries have been correctly called
		for (int entry = 0; entry < muontuples->GetEntries(); ++entry) 
		{
			muontuples->GetEntry(entry);
			cout << "Entry " << entry << ": ";

			cout << "type=" << type << " eventNo=" << eventNo << " pAbs=" << pAbs
				 << " pt=" << pt << " eta=" << eta << " id=" << id << endl;
		}

		//filling histograms
		for (int event_counter = 0; event_counter < particle_count; event_counter++)
		{
			muontuples->GetEntry(event_counter);
			if (ibin == 0) {
				soft_muon_yield->Fill(pt);
			}
			else{
				hard_muon_yield->Fill(pt);
			}
		}
	}	
	TFile* outFile =new TFile("soft_and_hard_qcdMacro.root", "RECREATE");


	hard_muon_yield->Draw("SAME");
	
	soft_muon_yield->Draw("SAME");
	c1->Write();

	delete outFile;
}
