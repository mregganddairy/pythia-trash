//Calculating rapidity cross secions and ratios WITH UNCERTAINTY using the reweighting method, CORRECTLY.
#include <TFile.h>
#include <TNtuple.h>
#include <TH1D.h>
#include <TString.h>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>

//this function processes one file: it opens the ROOT file, fills the provided histogram, and scales it.
double processFile(const std::string &filename, TH1D* hist, int VarNo) 
{
    TFile* infile = TFile::Open(filename.c_str(), "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return -1;
    }

    std::vector<double>* Luminosity = nullptr;
    infile->GetObject("luminosity", Luminosity);

	std::vector<std::vector<double>>* eventweights = nullptr;
    infile->GetObject("Eventweights", eventweights);
    
    //getting the TNtuple with muon information
    TNtuple *muontuples = (TNtuple*)infile->Get("muons");
    if (!muontuples) {
        std::cerr << "No TNtuple 'muons' found in " << filename << std::endl;
        infile->Close();
        delete infile;
        return -1;
    }

    Float_t  eventNo, index, mother1, mother2, pAbs, pt, y, eta, id;
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
			int   EventNo = (int)eventNo;
			double w0 = ((*eventweights)[EventNo][0]);
			double wi = ((*eventweights)[EventNo][VarNo]);
			double w = wi/w0;
            hist->Fill(eta, w);
        }
    }
    
    //normalize uing luminosity
    if (Luminosity && !Luminosity->empty()) {
        double scalebin = 1.0 / (*Luminosity)[0];
        hist->Scale(scalebin, "width");
    }
    
    double integral = hist->Integral();
    
    infile->Close();
    delete infile;
	delete eventweights;
    return integral;
}

int assignment8macro() 
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
	TCanvas *c7 = new TCanvas(); //Canvas for uncertaity ratios of all ratios

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
	TGraph* R_pmUnc = new TGraph(nbins);

	TGraphErrors* R_ZW = new TGraphErrors(nbins);
	TGraph* R_ZWUnc = new TGraph(nbins);

	TGraphErrors* A_pm = new TGraphErrors(nbins);
	TGraph* A_pmUnc = new TGraph(nbins);

	//Uncertainty vectors for boson cross sections
	vector<double> W_pUnc(nbins, 0);
	vector<double> W_mUnc(nbins, 0);
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
			boson = "wp";
		}
		else if (k == 1)
		{
			boson = "wm";
		}
		else if (k == 2)
		{
			boson = "z";
		}

		//processing central file
	//	TString centralFilename = PDFID+"00_output/"+PDFID+"00_pythia_output/"+boson+"_pwgevents_"+PDFID+"00.root";
		TString centralFilename = PDFID+"00_pythia_reweighted_output/"+PDFID+"00_pythia_output/"+boson+"_pwgevents_"+PDFID+"00.root";
		TH1D* centralHist = new TH1D("centralHist", "Central Cross Section; #eta; d#sigma/d#eta", nbins, xlow, xhigh);
		if (processFile(centralFilename.Data(), centralHist, 0) < 0) {
			std::cerr << "Error processing central file." << std::endl;
			return 1;
		}
		
		//processing uncertainty files
		std::vector<TH1D*> variationHists;
		for (int i = 1; i <= numVariations; ++i) {
			//getting names of files for each pdf 
			//std::ostringstream oss;
			//oss << PDFID+"00_output/" << PDFID << "00_pythia_output/" << boson << "_pwgevents_" << PDFID << (i < 10 ? "0" : "") << i << ".root";
			//oss << PDFID+"00_pythia_reweighted_output/" << PDFID << "00_pythia_output/" << boson << "_pwgevents_" << PDFID << "00.root";
			//std::string varFilename = oss.str();
			
			//create a histogram for this variation (clone the binning of centralHist)
			TH1D* hVar = new TH1D(Form("variation_%d", i), "Variation", nbins, xlow, xhigh);
			if (processFile(centralFilename.Data(), hVar, i) < 0) 
			{
				std::cerr << "Error processing file: " << centralFilename << std::endl;
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
		for (int bin = 1; bin <= nbins; bin++)
		{
			double sumUpperUncert = 0.0;
			double sumLowerUncert = 0.0;
			double sumSymmUncert = 0.0;
			double zero = centralHist->GetBinContent(bin);
			for (int i = 0; i < numPairs; ++i)
			{
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
			if (bin==20)
			{
				//cout <<"Symmetric uncertainty: "<< SymmUncert << endl;
			}

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

			if (bin==20)
			{
			//cout << "Uncertainty ratio:	" << symmErr/y << endl;
			}

			if (k==0)
			{
				WmUnc->SetPoint(bin-1, x, y);
				WmUnc->SetPointError(bin-1, 0, 0, upperErr, lowerErr);
				W_mUnc[bin-1] = symmErr;

				UpperWmUncRat->SetPoint(bin-1, x, 100*(1-(y-upperErr)/y));
				LowerWmUncRat->SetPoint(bin-1, x, -100*(1-(y-lowerErr)/y));
				SymmWmUncRat->SetPoint(bin-1, x, 100*(1-std::abs(y-symmErr)/y));
				std::cout << "bin: "<<bin << "	central:  " << y << endl;
				std::cout << "bin: "<<bin << "	upperErr: " << upperErr << endl;
				std::cout << "bin: "<<bin << "	lowerErr: " << lowerErr << endl;
			}
			else if (k==1)
			{
				WpUnc->SetPoint(bin-1, x, y);
				WpUnc->SetPointError(bin-1, 0, 0, upperErr, lowerErr);
				UpperWpUncRat->SetPoint(bin-1, x, 100*(1-(y-upperErr)/y));
				LowerWpUncRat->SetPoint(bin-1, x, -100*(1-(y-lowerErr)/y));
				SymmWpUncRat->SetPoint(bin-1, x, 100*(1-std::abs(y-symmErr)/y));
				W_pUnc[bin-1] = symmErr;


			}
			else if (k==2)
			{
				ZUnc->SetPoint(bin-1, x, y);
				ZUnc->SetPointError(bin-1, 0, 0, upperErr, lowerErr);
				UpperZUncRat->SetPoint(bin-1, x, 100*(1-(y-upperErr)/y));
				LowerZUncRat->SetPoint(bin-1, x, -100*(1-(y-lowerErr)/y));
				SymmZUncRat->SetPoint(bin-1, x, 100*(1-std::abs(y-symmErr)/y));
				Z_Unc[bin-1] = symmErr;
			}
		}


		c1->cd();
		c1->SetGridy();

		if (k==0)
		{
			
			WmUnc->SetTitle(" ");
			WmUnc->SetMarkerStyle(73);
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

			WpUnc->SetMarkerStyle(72);
			WpUnc->SetMarkerColor(kBlue);
			WpUnc->SetLineColor(kBlue);
			WpUnc->SetMarkerSize(0.8);
			WpUnc->Draw("P SAME");
			leg->AddEntry(WpUnc, "#mu #leftarrow W^{+}", "lep");
		}
		else if (k==2)
		{

			ZUnc->SetMarkerStyle(71);
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
			UpperWmUncRat->SetTitle(" ");
			UpperWmUncRat->SetMarkerStyle(73);
			UpperWmUncRat->SetMarkerColor(kRed);
			UpperWmUncRat->SetLineColor(kRed);
			UpperWmUncRat->Draw("AP");
			UpperWmUncRat->GetYaxis()->SetTitle("Percentage Uncertainty (%)");
			UpperWmUncRat->GetXaxis()->SetTitle("y");
			UpperWmUncRat->SetMaximum(5);
			UpperWmUncRat->SetMinimum(-5);
			leg2->AddEntry(UpperWmUncRat, "#mu #leftarrow W^{-}", "lep");

			LowerWmUncRat->SetMarkerStyle(73);
			LowerWmUncRat->SetMarkerColor(kRed);
			LowerWmUncRat->SetLineColor(kRed);
			LowerWmUncRat->Draw("P SAME");

		}
		else if (k==1)
		{
			UpperWpUncRat->SetMarkerStyle(72);
			UpperWpUncRat->SetMarkerColor(kBlue);
			UpperWpUncRat->SetLineColor(kBlue);
			UpperWpUncRat->Draw("P SAME");
			leg2->AddEntry(UpperWpUncRat, "#mu #leftarrow W^{+}", "lep");

			LowerWpUncRat->SetMarkerStyle(72);
			LowerWpUncRat->SetMarkerColor(kBlue);
			LowerWpUncRat->SetLineColor(kBlue);
			LowerWpUncRat->Draw("P SAME");

		}
		else if (k==2)
		{
			UpperZUncRat->SetMarkerStyle(71);
			UpperZUncRat->SetMarkerColor(kMagenta);
			UpperZUncRat->SetLineColor(kMagenta);
			UpperZUncRat->Draw("P SAME");
			leg2->AddEntry(UpperZUncRat, "#mu #leftarrow Z", "lep");

			LowerZUncRat->SetMarkerStyle(71);
			LowerZUncRat->SetMarkerColor(kMagenta);
			LowerZUncRat->SetLineColor(kMagenta);
			LowerZUncRat->Draw("P SAME");

			leg2->Draw();
		}
		//Clean up histograms (include or shit breaks)
		delete centralHist;
		delete uncertaintyHist;
		for (auto h : variationHists){delete h;}
	}

	//Now calculating and setting values for ratios and uncertainty percentages of uncertainties
	
	for (int bin = 1; bin <= nbins; bin++) 
	{

	double Wp_x = WpUnc->GetPointX(bin-1);
	double Wp_y = WpUnc->GetPointY(bin-1);

	double Wm_x = WmUnc->GetPointX(bin-1);
	double Wm_y = WmUnc->GetPointY(bin-1);

	double Z_x = ZUnc->GetPointX(bin-1);
	double Z_y = ZUnc->GetPointY(bin-1);

	double R_pm_value = Wp_y/Wm_y;
	double R_pm_uncertainty = std::sqrt(std::pow(Wp_y/Wm_y, 2)*(std::pow(W_mUnc[bin-1]/Wm_y,2) + std::pow(W_pUnc[bin-1]/Wp_y,2)));

	double A_pm_value = ((Wp_y-Wm_y)/(Wp_y + Wm_y));
	double A_pm_uncertainty = std::sqrt(std::pow(((1.0-A_pm_value)/(Wp_y+Wm_y))*W_pUnc[bin-1],2) + std::pow(((1.0+A_pm_value)/(Wp_y+Wm_y))*W_mUnc[bin-1],2));

	double R_ZW_value = ((Z_y)/(Wp_y + Wm_y));
	double R_ZW_uncertainty = std::sqrt(std::pow(R_ZW_value,2)*((std::pow(W_pUnc[bin-1],2) +std::pow(W_mUnc[bin-1],2))/std::pow(Wp_y+Wm_y,2) + std::pow(Z_Unc[bin-1]/Z_y, 2)));

	R_pm->SetPoint(bin-1, Wm_x, R_pm_value);
	R_pm->SetPointError(bin-1, 0, R_pm_uncertainty);

	A_pm->SetPoint(bin-1, Wm_x, A_pm_value);
	A_pm->SetPointError(bin-1, 0, A_pm_uncertainty);

	
	R_ZW->SetPoint(bin-1, Wm_x, R_ZW_value);
	R_ZW->SetPointError(bin-1, 0, R_ZW_uncertainty);

	//Now calculating and setting values for ratios and uncertainty percentages of uncertainties
	
	R_pmUnc->SetPoint(bin-1, Wm_x, 100*((std::abs(R_pm_uncertainty)/std::abs(R_pm_value))));
	A_pmUnc->SetPoint(bin-1, Wm_x, 100*((std::abs(A_pm_uncertainty)/std::abs(A_pm_value))));
	R_ZWUnc->SetPoint(bin-1, Wm_x, 100*((std::abs(R_ZW_uncertainty)/std::abs(R_ZW_value))));
	
	//std::cout << "Wplus:	 " << bin <<  WpUnc->GetPointY(bin-1) << std::endl;
	//std::cout << "Wminus:	 " << bin << WmUnc->GetPointY(bin-1) << std::endl;
	//std::cout << "Z:		 " << bin << ZUnc->GetPointY(bin-1) << std::endl;

	}
	
	//Plotting ratios on separate graphs
	
	c3->cd(); //plotting R_pm
	R_pm->SetTitle(" ");
	R_pm->SetMarkerStyle(23);
	R_pm->SetMarkerColor(kOrange);
	R_pm->SetLineColor(kOrange);
	R_pm->Draw("AP");
	R_pm->GetYaxis()->SetTitle("R_{#pm}");
	R_pm->GetXaxis()->SetTitle("y");

	c4->cd(); //plotting A_pm
	A_pm->SetTitle(" ");
	A_pm->SetMarkerStyle(28);
	A_pm->SetMarkerColor(kAzure);
	A_pm->SetLineColor(kAzure);
	A_pm->Draw("AP");
	A_pm->GetYaxis()->SetTitle("A_{#pm}");
	A_pm->GetXaxis()->SetTitle("y");

	c5->cd(); //plotting R_Z/W
	R_ZW->SetTitle(" ");
	R_ZW->SetMarkerStyle(47);
	R_ZW->SetMarkerColor(kViolet);
	R_ZW->SetLineColor(kViolet);
	R_ZW->Draw("AP");
	R_ZW->GetYaxis()->SetTitle("R_{Z/W}");
	R_ZW->GetXaxis()->SetTitle("y");


	//Plotting all ratios
	c6->cd(); 
			  
	TGraphErrors* R_ZW_2 = (TGraphErrors*)R_ZW->Clone("R_ZW_2");

	R_ZW_2->SetMarkerStyle(47);
	R_ZW_2->SetMarkerColor(kViolet);
	R_ZW_2->SetLineColor(kViolet);
	R_ZW_2->SetMaximum(1.7);
	R_ZW_2->SetMinimum(-1.7);
	R_ZW_2->GetYaxis()->SetTitle("Ratio");
	R_ZW_2->GetXaxis()->SetTitle("y");
	R_ZW_2->Draw("AP");

	leg4->AddEntry(R_ZW_2, "R_{Z/W}", "lep");
	R_pm->Draw("P SAME");
	leg4->AddEntry(R_pm, "R_{#pm}", "lep");
	A_pm->Draw("P SAME");
	leg4->AddEntry(A_pm, "A_{#pm}", "lep");
	leg4->Draw();
	
	//plotting uncertaintyi ratios
	c7->cd();

	R_ZWUnc->SetTitle("  ");
	R_ZWUnc->SetMarkerStyle(47);
	R_ZWUnc->SetMarkerColor(kViolet);
	R_ZWUnc->SetMinimum(0);
	R_ZWUnc->SetMaximum(40);
	R_ZWUnc->SetLineColor(kViolet);
	R_ZWUnc->GetYaxis()->SetTitle("Uncerainty Percentage (%)");
	R_ZWUnc->GetXaxis()->SetTitle("#eta");
	R_ZWUnc->Draw("AP");
	leg5->AddEntry(R_ZWUnc, "R_{Z/W}", "lep");

	R_pmUnc->SetMarkerStyle(23);
	R_pmUnc->SetMarkerColor(kOrange);
	R_pmUnc->SetLineColor(kOrange);
	R_pmUnc->Draw("P SAME");
	leg5->AddEntry(R_pmUnc, "R_{#pm}", "lep");


	A_pmUnc->SetMarkerStyle(28);
	A_pmUnc->SetMarkerColor(kAzure);
	A_pmUnc->SetLineColor(kAzure);
	A_pmUnc->Draw("P SAME");
	leg5->AddEntry(A_pmUnc, "A_{#pm}", "lep");
	leg5->Draw();

	TFile* outFile =new TFile("NLOMuonEtaDistMacro.root", "RECREATE");
	c1->Write();
	c2->Write();
	c3->Write();
	c4->Write();
	c5->Write();
	c6->Write();
	c7->Write();
	outFile->Close();
	//delete outFile;

	return 1;
    
}

