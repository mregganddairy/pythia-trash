#include "Pythia8/Pythia.h"
#include "Pythia8/PythiaParallel.h"
#include <cmath>
#include "TCanvas.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1D.h"
#include <TNtuple.h>

using namespace Pythia8;

//attempt at implimenting factorisation theorem for soft and hard qcd

int main()
{

	PythiaParallel pythia;

	//processes which will be on for all bins
//	pythia.readString("WeakDoubleBoson:all = on");
	pythia.readString("WeakSingleBoson:ffbar2W =on");
	pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
//	pythia.readString("WeakZ0:gmZmode =2");
	pythia.readString("Parallelism:numThreads = 4");
	
	//force Z and W to decay to muons to get better stats
/*	pythia.readString("23:onMode = off");
	pythia.readString("23:onIfAny = 13 -13"); 
	
	pythia.readString("24:onMode = off");
	pythia.readString("24:onIfAny =  13 14"); 
	*/

	//creating ROOT file for histograms
	TFile* outFile = new TFile("muonyield.root", "RECREATE");
	TH1F* hard_muon_yield = new TH1F("hard_muon_yield","", 150, 0, 150);
	TH1F* soft_muon_yield = new TH1F("soft_muon_yield","", 150, 0, 150);
	TH1F* sum_hard_muon_yield = new TH1F("sum_hard_muon_yield","", 150, 0, 150);

	TH1F* total_muon_yield = new TH1F("total_muon_yield","muon yield;pT;dN/pT", 150, 0, 150);

	pythia.readString("Beams:eCM = 13600.");
	
	//defining bins to seperate soft and hard qcd using pthat
	static const int nbins =7;
	static const double binedges[nbins+1] = {0., 14., 20., 30., 50., 75., 100., 150. };

	//tuple with relevant information on muons
	vector<TNtuple*> muontuples(nbins);


	//defining branches within the tree

	vector<double> binLuminosity(nbins); //luminosity from generated process sigma to calculate cross sections

	for (int i=0; i < nbins; ++i){
		muontuples[i] = new TNtuple("mu_stuff", "mu_stuff", "binNo:eventNo:index:status:mother1:mother2:mother11:mother111:pAbs:pt:y:eta:id");
	//muon number

		}
	
	int nevents = 10000;

	for (int ibin = 0; ibin < nbins; ++ibin)
	{
	/*	if (ibin == 0)
		{
		//	pythia.readString("HardQCD:all = off");
		//	pythia.readString("SoftQCD:nonDiffractive = on"); //find out why only non diffractive processes are used.

			pythia.readString("WeakSingleBoson:all = on"); //switching on low energy EW boson production
			pythia.readString("WeakBosonAndParton:all = off"); //switching off high energy EW boson production
		}
		
		//switching off softqcd but leaving single boson production on.
		else if (ibin == 1)
		{
		//	pythia.readString("HardQCD:all = on");
		//	pythia.readString("SoftQCD:nonDiffractive = off"); //find out why only non diffractive processes are used.
			pythia.readString("WeakSingleBoson:all = on"); //switching on low energy EW boson production
			pythia.readString("WeakBosonAndParton:all = off"); //switching off high energy EW boson production
		}

		//switching on boson and parton production
		else 
		{
		//	pythia.readString("HardQCD:all = on");
		//	pythia.readString("SoftQCD:nonDiffractive = off");

			pythia.readString("WeakSingleBoson:all = off");
			pythia.readString("WeakBosonAndParton:all = on"); //switching on high energy EW boson production
		}
		*/

	
		//setting limits on pthat for soft and hard qcd.
		pythia.settings.parm("PhaseSpace:pTHatMin", binedges[ibin]);
		pythia.settings.parm("PhaseSpace:pTHatMax",binedges[ibin+1]);

		pythia.init();
		
		hard_muon_yield->Reset(); //restart hard process binning for new bin 

		int event_count = 0; // to account for softqcd being dumb, we count the number of events in each bin seperately.
							 

		pythia.run(nevents, [&](Pythia* pythiaPtr)
		{
			//giving reference to the instance that generated the event.
			Event& event = pythiaPtr->event;
			const Info& info = pythiaPtr->info;

			double pThat = info.pTHat();

			if (ibin == 0 && info.isNonDiffractive() && pThat > binedges[ibin+1]) return;//apparently softqcd is stupid or something and doesn't have an upper limit on its pThat or something. So contribution above pThat max need to be manually removed.

			event_count++;

			//begin particle loop
			for (int i=0; i < event.size();++i)
			{
				if (abs(event[i].id()) == 13 && event[i].isFinal())
				{
					double particlemother1 = event[event[i].mother1()].id();
					double particlemother2 =event[event[i].mother2()].id();
					double particlemother11 =event[event[event[i].mother1()].mother1()].id();
					double particlemother111 =event[event[event[event[i].mother1()].mother1()].mother1()].id();
					double particlePAbs = event[i].pAbs();
					double particleStatus = event[i].status();
					double particlePt = event[i].pT();
					double particleRapidity = event[i].y();
					double particlePseudoRapidity = event[i].eta();
					double particleID = event[i].id();
					double eventNo = event_count;


					//cout<< "testing moms: OG particle: " <<event[i].id() <<"("<<event[i].index()<< ")" << ". mothers: " << event[event[i].mother1()].id()<<"("<< event[i].mother1()<<")" << " and "<< event[event[i].mother2()].id()<<"("<< event[i].mother2()<< ")"<< ". mothers of mothers: " << event[event[event[i].mother1()].mother1()].id()<<"("<< event[event[i].mother1()].mother1()<< ")" <<" and " << event[event[event[i].mother1()].mother2()].id()<<"("<< event[event[i].mother1()].mother2()<<")"<< endl;

					//cout <<"testing daughters: OG particle: " <<event[i].id() <<". daughter particles: "<< event[i].daughter1() <<" and "<< event[i].daughter2() << ". daughter particles: "<<event[event[i].daughter1()].daughter1() << " and " << event[event[i].daughter1()].daughter2() << endl;
					
					////man fuck this noise//////
					/*if ((abs(particlemother1) == 13) && (abs(particlemother2 == 13)))
					{
					cout << "newdaughterList: [";
					std::cout << event[event[event[i].mother1()].daughter1()].id() << ", " << event[event[event[i].mother1()].daughter2()].id(); 
					cout << "]" << endl;
					cout << "momentum of daughter and mother:" << event[event[event[i].mother1()].daughter1()].pT() << " and pT_m1=" << event[event[i].mother1()].pT() << " and pT_d2=" << event[event[event[i].mother1()].daughter2()].pT() << endl;
					}
					*/

					
					//filling tuple bin entries
					muontuples[ibin]->Fill(ibin,eventNo,i, particleStatus, particlemother1, particlemother2, particlemother11, 
							particlemother111, 
							particlePAbs, particlePt, particleRapidity, particlePseudoRapidity, particleID);
					
//////////////////////////////////////////////////////////////////////////////
	/*				//test filling bins
					if (ibin == 0)
					{
						soft_muon_yield->Fill(particlePt);
					}
					else
					{
						hard_muon_yield->Fill(particlePt);
					}
					*/
//////////////////////////////////////////////////////////////////////////////
				}
			}

		});
		binLuminosity[ibin] = event_count/(pythia.sigmaGen()*pow(10,9));//integrated luminosity used for normalisation/calculation of the cross sections

		cout <<"bin number"<< ibin<<" cross section: " << pythia.sigmaGen()*pow(10,9)<<endl;
		if (ibin ==0)
		{
			soft_muon_yield->Scale(1/binLuminosity[ibin] , "width");
		}
		else
		{
			hard_muon_yield->Scale(1/binLuminosity[ibin] , "width");
			hard_muon_yield->Draw();
			hard_muon_yield->Write();
			sum_hard_muon_yield->Add(hard_muon_yield); //adding contributions from each bin to get the total muon cross section
		}

	}

//plotting stuff now	

	TCanvas *c1 = new TCanvas();

	sum_hard_muon_yield->SetMarkerStyle(2);
	sum_hard_muon_yield->SetLineColor(kRed);
	sum_hard_muon_yield->Draw("SAME");

	soft_muon_yield->SetMarkerStyle(3);
	soft_muon_yield->SetLineColor(kGreen);
	soft_muon_yield->Draw("SAME");

	total_muon_yield->Add(sum_hard_muon_yield);
	total_muon_yield->SetLineColor(kBlue);
	total_muon_yield->Add(soft_muon_yield);

	total_muon_yield->SetMarkerStyle(1);
	total_muon_yield->Draw("SAME");

	soft_muon_yield->Write();
	sum_hard_muon_yield->Write();

	TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
	leg->AddEntry(sum_hard_muon_yield,"sum_hard_muon_yield", "l");
	leg->AddEntry(soft_muon_yield,"soft_muon_yield", "l");
	leg->AddEntry(total_muon_yield,"total_muon_yield", "l");
	leg->Draw("SAME");

	//total_muon_yield->Write();

	c1->Write();

	//Read Tuples to output file
	for (int ibin = 0; ibin < nbins; ++ibin)
	{
		muontuples[ibin]->Write(Form("muon%d", ibin));
	}


/*//////////////////////////////////////////////////////////////////////////////////////////
	//DEBUGGING: checkings stuffs
	for (int i = 0; i < muontuples.size(); ++i)
	 {
		TNtuple* tuples = muontuples[i];
		Float_t binNo, eventNo,index, mother1, mother2, pAbs, pt, eta, id;	
		int particle_count = tuples->GetEntries();	//number of muons
		tuples->SetBranchAddress("binNo", &binNo);
		tuples->SetBranchAddress("eventNo", &eventNo);
		tuples->SetBranchAddress("mother1", &mother1);
		tuples->SetBranchAddress("mother2",&mother2);
		tuples->SetBranchAddress("pAbs", &pAbs);
		tuples->SetBranchAddress("pt", &pt);
		tuples->SetBranchAddress("eta", &eta);
		tuples->SetBranchAddress("id", &id);	 
	 
	  
	 
		 for (int entry = 0; entry < tuples->GetEntries(); ++entry)
		 {
			 tuples->GetEntry(entry);
			 cout << "Entry " << entry << ": ";
			 //printing variables in tuple. 
			 cout << "binNo=" << binNo << " eventNo=" << eventNo << " mother1ID=" 
				 << mother1 << " mother2ID=" << mother2 << " pAbs=" << pAbs
				  << " pt=" << pt << " eta=" << eta << " id=" << id << endl;
		 }
	 }
*////////////////////////////////////////////////////////////////////////////

	outFile->WriteObject(&binLuminosity, "luminosity");

	delete outFile;
	
	return 0;
	
}
