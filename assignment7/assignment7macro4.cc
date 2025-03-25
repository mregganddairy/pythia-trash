#include <TFile.h>
#include <TNtuple.h>
#include <TH1D.h>
#include <TString.h>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>

//this function processes one file: it opens the ROOT file, fills the provided histogram, and scales it.
double processFile(const std::string &filename, TH1D* hist) 
{
    TFile* infile = TFile::Open(filename.c_str(), "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return -1;
    }

    std::vector<double>* Luminosity = nullptr;
    infile->GetObject("luminosity", Luminosity);
    
    //getting the TNtuple with muon information
    TNtuple *muontuples = (TNtuple*)infile->Get("muons");
    if (!muontuples) {
        std::cerr << "No TNtuple 'muons' found in " << filename << std::endl;
        infile->Close();
        delete infile;
        return -1;
    }
    
    Float_t eventNo, index, mother1, mother2, pAbs, pt, y, eta, id;
    int particle_count = muontuples->GetEntries();
    muontuples->SetBranchAddress("eventNo", &eventNo);
    muontuples->SetBranchAddress("index", &index);
    muontuples->SetBranchAddress("mother1", &mother1);
    muontuples->SetBranchAddress("mother2", &mother2);
    muontuples->SetBranchAddress("pAbs", &pAbs);
    muontuples->SetBranchAddress("pt", &pt);
    muontuples->SetBranchAddress("y", &y);
    muontuples->SetBranchAddress("eta", &eta);
    muontuples->SetBranchAddress("id", &id);
    
    for (int i = 0; i < particle_count; ++i) {
        muontuples->GetEntry(i);
        if ((std::abs(mother1) == 24) || ((std::abs(mother1)==13) && (std::abs(mother2)==90))) {
            hist->Fill(eta);
        }
    }
    
    //normalize using luminosity
    if (Luminosity && !Luminosity->empty()) {
        double scalebin = 1.0 / (*Luminosity)[0];
        hist->Scale(scalebin, "width");
    }
    
    double integral = hist->Integral();
    
    infile->Close();
    delete infile;
    return integral;
}

int assignment7macro4() 
{
    //define histogram binning parameters
    const int nbins = 5;
    const double xlow = -5.0;
    const double xhigh = 5.0;

	std::string PDFID= "274";
	const int numVariations = 58;
    
	TFile* outFile =new TFile("NLOMuonEtaDistMacro.root", "RECREATE");
	TCanvas *c1 = new TCanvas(); //canvas for eta distribution  

	TGraphAsymmErrors* WpUnc = new TGraphAsymmErrors(nbins);
	TGraphAsymmErrors* WmUnc = new TGraphAsymmErrors(nbins);
	TGraphAsymmErrors* ZUnc = new TGraphAsymmErrors(nbins);
	TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
	//loop over each boson
	for (int k = 0; k < 3; ++k)
	{
		std::string boson;
		if (k == 0)
		{
			boson = "wm";
		}
		else if (k == 1)
		{
			boson = "wp";
		}
		else if (k == 2)
		{
			boson = "z";
		}

		//processing central file
		TString centralFilename = PDFID+"00_output/"+PDFID+"00_pythia_output/"+boson+"_pwgevents_"+PDFID+"00.root";
		TH1D* centralHist = new TH1D("centralHist", "Central Cross Section; #eta; d#sigma/d#eta", nbins, xlow, xhigh);
		if (processFile(centralFilename.Data(), centralHist) < 0) {
			std::cerr << "Error processing central file." << std::endl;
			return 1;
		}
		
		//processing uncertainty files
		std::vector<TH1D*> variationHists;
		for (int i = 1; i <= numVariations; ++i) {
			//getting names of files for each pdf 
			std::ostringstream oss;
			oss << PDFID+"00_output/" << PDFID << "00_pythia_output/" << boson << "_pwgevents_" << PDFID << (i < 10 ? "0" : "") << i << ".root";
			std::string varFilename = oss.str();
			
			//create a histogram for this variation (clone the binning of centralHist)
			TH1D* hVar = new TH1D(Form("variation_%d", i), "Variation", nbins, xlow, xhigh);
			if (processFile(varFilename, hVar) < 0) {
				std::cerr << "Error processing file: " << varFilename << std::endl;
				continue;
			}
			variationHists.push_back(hVar);
		}
		
		const int numPairs = variationHists.size() / 2;  
		TH1D* uncertaintyHist = (TH1D*)centralHist->Clone("uncertaintyHist");
		TH1D* upperUncertaintyHist = (TH1D*)centralHist->Clone("uncertaintyHist");
		TH1D* lowerUncertaintyHist = (TH1D*)centralHist->Clone("uncertaintyHist");
		uncertaintyHist->Reset();
		
		//compute the uncertainty in each bin
		for (int bin = 1; bin <= nbins; bin++) {
			double sumUpperUncert = 0.0;
			double sumLowerUncert = 0.0;
			double zero = centralHist->GetBinContent(bin);
			for (int i = 0; i < numPairs; ++i) {
				//getting plus and minus variations from paired files
				double plus  = variationHists[2*i]->GetBinContent(bin);
				double minus = variationHists[2*i + 1]->GetBinContent(bin);
				//compute the difference according to the asymmetric Hessian prescription
				sumUpperUncert += std::pow(std::max(std::max(plus - zero, minus - zero), 0.0),2);
				sumLowerUncert += std::pow(std::max(std::max(zero - plus, zero - minus), 0.0),2);
			}
			double upperUncert = std::sqrt(sumUpperUncert);
			double lowerUncert = std::sqrt(sumLowerUncert);
			upperUncertaintyHist->SetBinContent(bin, upperUncert);
			lowerUncertaintyHist->SetBinContent(bin, lowerUncert);
		}
		
		c1->cd();
		c1->SetGridy();

		WmUnc->GetYaxis()->SetTitle("d#sigma/d#eta (pb/GeV/c)");
		WmUnc->GetXaxis()->SetTitle("#eta");


		//plotting graphs
		for (int bin = 1; bin <= nbins; ++bin)
		{
			double x = centralHist->GetBinCenter(bin);
			double y = centralHist->GetBinContent(bin);
			double upperErr = upperUncertaintyHist->GetBinContent(bin);
			double lowerErr = lowerUncertaintyHist->GetBinContent(bin);

			if (k==0)
			{
			WmUnc->SetPoint(bin-1, x, y);
			WmUnc->SetPointError(bin-1, 0, 0, upperErr, lowerErr);
			}
			else if (k==1)
			{
			WpUnc->SetPoint(bin-1, x, y);
			WpUnc->SetPointError(bin-1, 0, 0, upperErr, lowerErr);
			}
			else if (k==2)
			{
			ZUnc->SetPoint(bin-1, x, y);
			ZUnc->SetPointError(bin-1, 0, 0, upperErr, lowerErr);
			}
		}

		if (k==0)
		{
		WmUnc->SetTitle("POGGERS");
		WmUnc->SetMarkerStyle(27);
		WmUnc->SetMarkerColor(kRed);
		WmUnc->SetMinimum(0);
		WmUnc->SetMaximum(2000);
		WmUnc->Draw("AP");
		leg->AddEntry(WmUnc, "W^{-}", "lep");
		}
		else if (k==1)
		{
		WpUnc->SetMarkerStyle(25);
		WpUnc->SetMarkerColor(kBlue);
		WpUnc->Draw("P SAME");
		leg->AddEntry(WpUnc, "W^{+}", "lep");
		}
		else if (k==2)
		{
		ZUnc->SetMarkerStyle(26);
		ZUnc->SetMarkerColor(kGreen);
		ZUnc->Draw("P SAME");
		leg->AddEntry(ZUnc, "Z", "lep");
		leg->Draw("SAME");
		}
		
		
		// Clean up histograms
		delete centralHist;
		delete uncertaintyHist;
		for (auto h : variationHists) {
			delete h;
		}
	}
	outFile->cd();
	c1->Write();
	outFile->Close();
	delete outFile;

	return 1;
    
}

