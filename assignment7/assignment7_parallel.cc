//Showering powheg files that use various different PDFs
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
#include "Pythia8Plugins/PowhegHooks.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Generator
  PythiaParallel pythia;

  // Load configuration file
  pythia.readFile("assignment7.cmnd");
  pythia.readString("Parallelism:numThreads = 4");
  //pythia.readString("Main:debug = on");

  // Read in main settings.
  int nEvent      = pythia.settings.mode("Main:numberOfEvents");
  int nError      = pythia.settings.mode("Main:timesAllowErrors");
  // Read in key POWHEG matching settings.
  int vetoMode    = pythia.settings.mode("POWHEG:veto");
  int MPIvetoMode = pythia.settings.mode("POWHEG:MPIveto");
  bool loadHooks  = (vetoMode > 0 || MPIvetoMode > 0);
  // Read in shower settings.
  int showerModel = pythia.settings.mode("PartonShowers:model");


  
	TNtuple* muontuples;
	vector<double> Luminosity(1); //luminosity from generated process sigma to calculate cross sections

	TFile* outFile = new TFile("wm_pwgevents_14000.root", "RECREATE");
	muontuples = new TNtuple("mu_stuff", "mu_stuff", "eventNo:index:status:mother1:mother2:daughter1:daughter2:pAbs:pt:y:eta:id");


  // Add in user hooks for shower vetoing.
	shared_ptr<PowhegHooks> powhegHooks;
	powhegHooks = make_shared<PowhegHooks>();
    //pythia.setUserHooksPtr((UserHooksPtr)powhegHooks);

  if (loadHooks) {

    // Set showers to start at the kinematical limit.
    if (vetoMode > 0) {
      if (showerModel == 1 || showerModel == 3) {
        // Pythia and Dire have separate settings for ISR and FSR.
        pythia.readString("SpaceShower:pTmaxMatch = 2");
        pythia.readString("TimeShower:pTmaxMatch = 2");
      } else if (showerModel == 2) {
        // Vincia only has one common setting for both ISR and FSR.
        pythia.readString("Vincia:pTmaxMatch = 2");
      }
    }

    // Set MPI to start at the kinematical limit.
    if (MPIvetoMode > 0) {
      pythia.readString("MultipartonInteractions:pTmaxMatch = 2");
    }

    
  }

  // Initialise and list settings
  if (!pythia.init([&](Pythia* pythiaInstance)
	  {
	  pythiaInstance->setUserHooksPtr(powhegHooks);
	  return true;
	  })) return 1;

  // Counters for number of ISR/FSR emissions vetoed
  unsigned long int nISRveto = 0, nFSRveto = 0;

  // Begin event loop; generate until nEvent events are processed
  // or end of LHEF file
  int event_count = 0;// to account for potentially vetoed events
  int iEvent = 0, iError = 0;

    // Generate the next event
      pythia.run(nEvent, [&](Pythia* pythiaPtr) 
	  {

	  Event& event =pythiaPtr->event;
	  const Info& info = pythiaPtr->info;


      // If failure because reached end of file then exit event loop

      // Otherwise count event failure and continue/exit as necessary
	  if(!pythiaPtr->next())
	  {
		  if (info.atEndOfFile()) return;
		  cout << "Warning: event " << iEvent << " failed" << endl;
		  if (++iError == nError) {
			cout << "Error: too many event failures... exiting" << endl;
			}
	  return;
	  }


    /*
     * Process dependent checks and analysis may be inserted here
     */ 
	 event_count++;
     for (int i=0; i < event.size();++i)
			{
				if (abs(event[i].id()) == 13 && event[i].isFinal())
				{
					double particlemother1 = event[event[i].mother1()].id();
					double particlemother2 =event[event[i].mother2()].id();
					double particledaughter1 =event[i].daughter1();
					double particledaughter2 =event[i].daughter2();
					double particlePAbs = event[i].pAbs();
					double particleStatus = event[i].status();
					double particlePt = event[i].pT();
					double particleRapidity = event[i].y();
					double particlePseudoRapidity = event[i].eta();
					double particleID = event[i].id();
					
					//filling tuple entries
					muontuples->Fill( iEvent,i, particleStatus, particlemother1, particlemother2, particledaughter1, particledaughter2, 
							particlePAbs, particlePt, particleRapidity, particlePseudoRapidity, particleID);
				}
				if (abs(event[event[i].mother1()].id()) == 13 && event[event[i].mother1()].isFinal())
				{
					cout <<"pdg id: " << event[i].id() << endl;
				}
			}	
    // Update ISR/FSR veto counters
    if (loadHooks) {
      nISRveto += powhegHooks->getNISRveto();
      nFSRveto += powhegHooks->getNFSRveto();
    }

	
    // If nEvent is set, check and exit loop if necessary
    ++iEvent;
    if (nEvent != 0 && iEvent == nEvent) return;

  }); // End of event loop.

	  // Statistics, histograms and veto information
  pythia.stat();
  cout << "Number of ISR emissions vetoed: " << nISRveto << endl;
  cout << "Number of FSR emissions vetoed: " << nFSRveto << endl;
  cout << endl;


  
	Luminosity[0] = event_count/(pythia.sigmaGen()*pow(10,9));//integrated luminosity used for normalisation/calculation of the cross sections
	muontuples->Write("muons");
	


	outFile->WriteObject(&Luminosity, "luminosity"); 
	delete outFile;


  // Done.
  return 0;
}
