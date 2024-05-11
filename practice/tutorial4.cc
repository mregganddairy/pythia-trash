//Attempt at plotting the dN/dpT distribution of muons produced in a pp collision
//with CM energy = 13.6TeV, pT>4GeV and various cuts on the psuedorapidity.
#include "Pythia/Pythia.h"

using namespace Pythia8;

int main()
	{
	Pythia pythia;	//Switch on process
	pythia.readString("Beam:eCM = 13600.")
	pythia.readString("PhaseSpace:pTHatMin = 4.")

	Hist pT ("muon transverse momentum", 100, 4., 200.);	//booking histograms

	//Begin event loop.
	for (int iEvent = 0; iEvent < 1000; ++iEvent)
		{
		pythia.next();	//Generate events

//		int imuon = 0;	//initiliaze value to search for muon

		//begin particle loop.
		for (int i = 0; i < pythia.event.size(); ++i)
			{
			if (pythia.event[i].id() == 13) pT.fill( pythia.event[i].pT() );
			}
		}

	cout << pT << eta;

	return 0;
	}
