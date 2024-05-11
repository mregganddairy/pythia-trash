#include "Pythia8/Pythia.h"
#include <cmath>
#include "TCanvas.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include <TNtuple.h>

using namespace Pythia8;

//attempt at implimenting factorisation theorem for soft and hard qcd

int main(){
	//Turn softQCD on/off
	bool softQCD= false;

	Pythia pythia;

	//creating ROOT file for histograms
	TFile* outFile = new TFile("soft_and_hard_qcd.root", "RECREATE");
	TH1F* hard_muon_yield = new TH1F("hard_muon_yield","", 50, 0, 50);
	TH1F* soft_muon_yield = new TH1F("soft_muon_yield","", 50, 0, 50);

	TH1F* total_muon_yield = new TH1F("total_muon_yield","muon yield;pT;dN/pT", 50, 0, 50);

	pythia.readString("Beams:eCM = 13600.");
	
	//defining bins to seperate soft and hard qcd using pthat
	static const int nbins =2;
	static const double binedges[nbins+1] = {0., 14.,150.};

	//tuple with relevant information on muons
	vector<TNtuple*> muontuples(nbins);

	for (int i=0; i < nbins; ++i){
		muontuples[i] = new TNtuple("mu_stuff", "mu_stuff", "type:eventNo:pAbs:pt:eta:id");
	//muon number
		}
	
	int nevents = 1000;

	for (int ibin = 0; ibin < nbins; ++ibin){
		if (ibin == 0){
			pythia.readString("HardQCD:all = off");
			pythia.readString("SoftQCD:all = on"); //find out why only non diffractive processes are used.
		}
		else {
			pythia.readString("HardQCD:all = on");
			pythia.readString("SoftQCD:all = off");
		}
	
		//setting limits on pthat for soft and hard qcd.
		pythia.settings.parm("PhaseSpace:pTHatMin", binedges[ibin]);
		pythia.settings.parm("PhaseSpace:pTHatMax",binedges[ibin+1]);

		pythia.init();
		

		//begining event loop
		for (int iEvent = 0; iEvent < nevents; ++iEvent)
		{
			if (!pythia.next()) continue;

			double pThat = pythia.info.pTHat();

			if (ibin == 0 && pythia.info.isNonDiffractive() && pThat > binedges[ibin+1]) continue;

			//begin particle loop
			for (int i=0; i < pythia.event.size();++i)
			{
				if (abs(pythia.event[i].id()) == 13 && pythia.event[i].isFinal())
				{
					double particlePAbs = pythia.event[i].pAbs();
					double particlePt = pythia.event[i].pT();
					double particlePseudoRapidity = pythia.event[i].eta();
					double particleID = pythia.event[i].id();
					
					//filling tuple bin entries
					muontuples[ibin]->Fill(ibin, iEvent, particlePAbs, particlePt, particlePseudoRapidity, particleID);
					if (ibin == 0)
					{
						soft_muon_yield->Fill(particlePt);
					}
					else
					{
						hard_muon_yield->Fill(particlePt);
					}
					total_muon_yield->Fill(particlePt);
				}
			}
			

		}

	}
	//things that haven't been accounte for:
	//luminosity and normalisation of cross sections
	//different types of soft processes

	TCanvas *c1 = new TCanvas();

	hard_muon_yield->Draw("SAME");
//	hard_muon_yield->Write();
	
	soft_muon_yield->Draw("SAME");
	soft_muon_yield->Write();


	total_muon_yield->Draw("SAME");
	total_muon_yield->Write();

	c1->Write();

	//Read Tuples to output file
	for (int ibin = 0; ibin < nbins; ++ibin)
	{
		muontuples[ibin]->Write(Form("muon%d", ibin));
	}

//check to see if tuples were filled correctly 
for (int i = 0; i < muontuples.size(); ++i) 
{
    TNtuple* tuple = muontuples[i];
    cout << "Entries in muontuples[" << i << "]:" << endl;

	Float_t type, eventNo, pAbs, pt, eta, id;	
	tuple->SetBranchAddress("type", &type);
	tuple->SetBranchAddress("eventNo", &eventNo);
	tuple->SetBranchAddress("pAbs", &pAbs);
	tuple->SetBranchAddress("pt", &pt);
	tuple->SetBranchAddress("eta", &eta);
	tuple->SetBranchAddress("id", &id);




	
    for (int entry = 0; entry < tuple->GetEntries(); ++entry) 
	{
        tuple->GetEntry(entry);
        cout << "Entry " << entry << ": ";
        //printing variables in tuple. 
        cout << "type=" << type << " eventNo=" << eventNo << " pAbs=" << pAbs
             << " pt=" << pt << " eta=" << eta << " id=" << id << endl;
    }
}



	delete outFile;
	
	return 0;
	
}
