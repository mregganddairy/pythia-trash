//Calculating rapidity cross secion WITH UNCERTAINTY.
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
    const int nbins = 50;
    const double xlow = -5.0;
    const double xhigh = 5.0;

	std::string PDFID= "274";
	const int numVariations = 64; 
    
	TCanvas *c1 = new TCanvas(); //canvas for eta distribution  
	TCanvas *c2 = new TCanvas(); //canvas for eta uncertainty ratios 
	TCanvas *c3 = new TCanvas(); //canvas for ratio of W+/W-
	TCanvas *c4 = new TCanvas(); //Canvas for ratio of W/Z
	TCanvas *c5 = new TCanvas(); //Canvas for A asymmetry
	TCanvas *c6 = new TCanvas(); //Canvas for all ratios

	//plots with uncertainty
	TGraphAsymmErrors* WpUnc = new TGraphAsymmErrors(nbins);
	TGraphAsymmErrors* WmUnc = new TGraphAsymmErrors(nbins);
	TGraphAsymmErrors* ZUnc = new TGraphAsymmErrors(nbins);


	//Uncerainty ratio plots
	TGraph* UpperWpUncRat = new TGraph(nbins); 
	TGraph* LowerWpUncRat = new TGraph(nbins);
	TGraph* SymmWpUncRat = new TGraph(nbins);

	TGraph* UpperWmUncRat = new TGraph(nbins); 
	TGraph* LowerWmUncRat = new TGraph(nbins); 
	TGraph* SymmWmUncRat = new TGraph(nbins);

	TGraph* UpperZUncRat = new TGraph(nbins); 
	TGraph* LowerZUncRat = new TGraph(nbins); 
	TGraph* SymmZUncRat = new TGraph(nbins);
	
	//Ratios of calculated values and their uncertainties
	TGraphErrors* R_pm = new TGraphErrors(nbins);
	vector<double> W_pUnc(nbins, 0);

	TGraphErrors* R_ZW = new TGraphErrors(nbins);
	vector<double> W_mUnc(nbins, 0);


	TGraphErrors* A_pm = new TGraphErrors(nbins);
	vector<double> Z_Unc(nbins, 0);


	TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
	TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
	TLegend *leg3 = new TLegend(0.6, 0.7, 0.9, 0.9);
	TLegend *leg4 = new TLegend(0.6, 0.7, 0.9, 0.9);
	TLegend *leg5 = new TLegend(0.6, 0.7, 0.9, 0.9);
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
		TH1D* symmUncertaintyHist = (TH1D*)centralHist->Clone("uncertaintyHist");
		uncertaintyHist->Reset();
		
		//compute the uncertainty in each bin
		for (int bin = 1; bin <= nbins; bin++) {
			double sumUpperUncert = 0.0;
			double sumLowerUncert = 0.0;
			double sumSymmUncert = 0.0;
			double zero = centralHist->GetBinContent(bin);
			for (int i = 0; i < numPairs; ++i) {
				//getting plus and minus variations from paired files
				double plus  = variationHists[2*i]->GetBinContent(bin);
				double minus = variationHists[2*i + 1]->GetBinContent(bin);
				//compute the difference according to the asymmetric Hessian prescription
				sumUpperUncert += std::pow(std::max(std::max(plus - zero, minus - zero), 0.0),2);
				sumLowerUncert += std::pow(std::max(std::max(zero - plus, zero - minus), 0.0),2);
				sumSymmUncert += std::pow(plus - minus,2);

			}
			double upperUncert = std::sqrt(sumUpperUncert);
			double lowerUncert = std::sqrt(sumLowerUncert);

			double SymmUncert = std::sqrt(sumSymmUncert)*0.5;

			upperUncertaintyHist->SetBinContent(bin, upperUncert);
			lowerUncertaintyHist->SetBinContent(bin, lowerUncert);
			symmUncertaintyHist->SetBinContent(bin, SymmUncert);
		}
		

		WmUnc->GetYaxis()->SetTitle("d#sigma/dy (pb)");
		WmUnc->GetXaxis()->SetTitle("y");


		//plotting graphs
		for (int bin = 1; bin <= nbins; ++bin)
		{
			double x = centralHist->GetBinCenter(bin);
			double y = centralHist->GetBinContent(bin);
			double upperErr = upperUncertaintyHist->GetBinContent(bin);
			double lowerErr = lowerUncertaintyHist->GetBinContent(bin);
			double symmErr = symmUncertaintyHist->GetBinContent(bin);

			if (k==0)
			{
				WmUnc->SetPoint(bin-1, x, y);
				WmUnc->SetPointError(bin-1, 0, 0, upperErr, lowerErr);
				W_mUnc[bin-1] = symmErr;
				

				UpperWmUncRat->SetPoint(bin-1, x, 100*(1-std::abs(y-upperErr)/y));
				LowerWmUncRat->SetPoint(bin-1, x, 100*(1-std::abs(y-lowerErr)/y));
				SymmWmUncRat->SetPoint(bin-1, x, 100*(1-std::abs(y-symmErr)/y));
			}
			else if (k==1)
			{
				WpUnc->SetPoint(bin-1, x, y);
				WpUnc->SetPointError(bin-1, 0, 0, upperErr, lowerErr);
				SymmWpUncRat->SetPoint(bin-1, x, 100*(1-std::abs(y-symmErr)/y));
				W_pUnc[bin-1] = symmErr;


			}
			else if (k==2)
			{
				ZUnc->SetPoint(bin-1, x, y);
				ZUnc->SetPointError(bin-1, 0, 0, upperErr, lowerErr);
				SymmZUncRat->SetPoint(bin-1, x, 100*(1-std::abs(y-symmErr)/y));
				Z_Unc[bin-1] = symmErr;
			}
		}


		c1->cd();
		c1->SetGridy();

		if (k==0)
		{
			
			WmUnc->SetTitle(" ");
			WmUnc->SetMarkerStyle(45);
			WmUnc->SetMarkerColor(kRed);
			WmUnc->SetLineColor(kRed);
			WmUnc->SetMarkerSize(0.8);
			WmUnc->SetMinimum(0);
			WmUnc->SetMaximum(2000);
			WmUnc->Draw("AP");
			leg->AddEntry(WmUnc, "#mu #leftarrow W^{-}", "lep");
			
		}
		else if (k==1)
		{

			WpUnc->SetMarkerStyle(44);
			WpUnc->SetMarkerColor(kBlue);
			WpUnc->SetLineColor(kBlue);
			WpUnc->SetMarkerSize(0.8);
			WpUnc->Draw("P SAME");
			leg->AddEntry(WpUnc, "#mu #leftarrow W^{+}", "lep");
		}
		else if (k==2)
		{

			ZUnc->SetMarkerStyle(43);
			ZUnc->SetMarkerColor(kMagenta);
			ZUnc->SetLineColor(kMagenta);
			ZUnc->SetMarkerSize(0.7);
			ZUnc->Draw("P SAME");
			leg->AddEntry(ZUnc, "#mu #leftarrow Z", "lep");
			leg->Draw("SAME");
		}
		
		//Plotting uncertainty percentage
		c2->cd();
		c2->SetGridy();
		
		if (k==0)
		{
			SymmWmUncRat->SetTitle(" ");
			SymmWmUncRat->SetMarkerStyle(45);
			SymmWmUncRat->SetMarkerColor(kRed);
			SymmWmUncRat->SetLineColor(kRed);
			SymmWmUncRat->Draw("AP");
			SymmWmUncRat->GetYaxis()->SetTitle("Percentage Uncertainty (%)");
			SymmWmUncRat->GetXaxis()->SetTitle("y");
			leg2->AddEntry(SymmWmUncRat, "#mu #leftarrow W^{-}", "lep");

		}
		else if (k==1)
		{
			SymmWpUncRat->SetMarkerStyle(44);
			SymmWpUncRat->SetMarkerColor(kBlue);
			SymmWpUncRat->SetLineColor(kBlue);
			SymmWpUncRat->Draw("P SAME");
			leg2->AddEntry(SymmWpUncRat, "#mu #leftarrow W^{+}", "lep");

		}
		else if (k==2)
		{
			SymmZUncRat->SetMarkerStyle(43);
			SymmZUncRat->SetMarkerColor(kMagenta);
			SymmZUncRat->SetLineColor(kMagenta);
			SymmZUncRat->Draw("P SAME");
			leg2->AddEntry(SymmZUncRat, "#mu #leftarrow Z", "lep");

			leg2->Draw();
		}
		//Clean up histograms (include or shit breaks)
		delete centralHist;
		delete uncertaintyHist;
		for (auto h : variationHists){delete h;}
	}

	//Now calculating and setting values for ratios
	
	for (int bin = 1; bin <= nbins; bin++) 
	{

	double Wp_x = WpUnc->GetPointX(bin-1);
	double Wp_y = WpUnc->GetPointY(bin-1);

	double Wm_x = WmUnc->GetPointX(bin-1);
	double Wm_y = WmUnc->GetPointY(bin-1);

	double Z_x = ZUnc->GetPointX(bin-1);
	double Z_y = ZUnc->GetPointY(bin-1);

	R_pm->SetPoint(bin-1, Wm_x, Wp_y/Wm_y);
	R_pm->SetPointError(bin-1, 0, std::sqrt(std::pow(Wp_y/Wm_y, 2)*(std::pow(W_mUnc[bin-1]/Wm_y,2) + std::pow(W_pUnc[bin-1]/Wp_y,2))));

	A_pm->SetPoint(bin-1, Wm_x, ((Wp_y-Wm_y)/(Wp_y + Wm_x)));
	A_pm->SetPointError(bin-1, 0, std::sqrt(std::pow((Wp_y-Wm_y)/(Wm_y + Wp_y), 2)*(std::pow((std::pow(W_mUnc[bin-1],2)+std::pow(W_pUnc[bin-1],2))/(Wp_y-Wm_y),2) + (std::pow((std::pow(W_mUnc[bin-1],2)+std::pow(W_pUnc[bin-1],2))/(Wp_y+Wm_y),2)))));

	
	R_ZW->SetPoint(bin-1, Wm_x, ((Z_y)/(Wp_y + Wm_x)));
	R_ZW->SetPointError(bin-1, 0, std::sqrt(std::pow((Z_y)/(Wm_y+Wp_y), 2)*((std::pow(Z_Unc[bin-1]/Z_y,2)) + (std::pow((std::pow(W_mUnc[bin-1],2)+std::pow(W_pUnc[bin-1],2))/(Wp_y+Wm_y),2)))));

	}
	
	c3->cd(); //plotting R_pm
	R_pm->SetTitle(" ");
	R_pm->SetMarkerStyle(45);
	R_pm->SetMarkerColor(kRed);
	R_pm->SetLineColor(kRed);
	R_pm->Draw("AP");
	R_pm->GetYaxis()->SetTitle("R_{#pm}");
	R_pm->GetXaxis()->SetTitle("y");

	c4->cd(); //plotting A_pm
	A_pm->SetTitle(" ");
	A_pm->SetMarkerStyle(44);
	A_pm->SetMarkerColor(kBlue);
	A_pm->SetLineColor(kBlue);
	A_pm->Draw("AP");
	A_pm->GetYaxis()->SetTitle("A_{#pm}");
	A_pm->GetXaxis()->SetTitle("y");

	c5->cd(); //plotting R_Z/W
	R_ZW->SetTitle(" ");
	R_ZW->SetMarkerStyle(43);
	R_ZW->SetMarkerColor(kMagenta);
	R_ZW->SetLineColor(kMagenta);
	R_ZW->Draw("AP");
	R_ZW->GetYaxis()->SetTitle("R_{Z/W}");
	R_ZW->GetXaxis()->SetTitle("y");

	c3->Write();
	c4->Write();
	c5->Write();


	TGraphAsymmErrors* pain = new TGraphAsymmErrors(1);
	c6->cd(); //Plotting all ratios
			  
	pain->SetMaximum(1.7);
	pain->SetMinimum(-1.7);
	pain->GetYaxis()->SetTitle("Ratio");
	pain->GetXaxis()->SetTitle("y");
	pain->Draw("AP");

	R_ZW->Draw("P SAME");
	leg4->AddEntry(R_ZW, "R_{Z/W}", "lep");
	R_pm->Draw("P SAME");
	leg4->AddEntry(R_pm, "R_{#pm}", "lep");
	A_pm->Draw("P SAME");
	leg4->AddEntry(A_pm, "A_{#pm}", "lep");
	leg4->Draw();
	

	TFile* outFile =new TFile("NLOMuonEtaDistMacro.root", "RECREATE");
	c1->Write();
	c2->Write();
	c6->Write();
	outFile->Close();
	delete outFile;

	return 1;
    
}

