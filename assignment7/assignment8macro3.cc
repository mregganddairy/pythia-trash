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

int assignment8macro3() 
{
    //define histogram binning parameters
    const int nbins = 50;
	const int nbosons =3;
    const double xlow = -5.0;
    const double xhigh = 5.0;
	const int nPDFs = 3;
	double_t binwidth = (xhigh-xlow)/(2*nbins);
	cout<<binwidth;
	//

	std::vector<std::string> PDFIDs= {"274", "140", "3311"};
	const std::vector<int> AllnumVariations = {64, 58, 99}; 
	//const std::vector<int> AllnumVariations = {8, 8, 4}; 

	
	//vectors with the values and uncertainty of each of each boson
	double_t x_values[nbins];

	double_t y_values[nPDFs][nbosons][nbins];
	double_t y_upper_uncertainty_values[nPDFs][nbosons][nbins];
	double_t y_lower_uncertainty_values[nPDFs][nbosons][nbins];

	double_t y_ratio_values[nPDFs][nbosons][nbins];
	double_t y_ratio_upper_uncertainty_values[nPDFs][nbosons][nbins];
	double_t y_ratio_lower_uncertainty_values[nPDFs][nbosons][nbins];

	double_t binwidth_vector[nbins];
	double_t stat_uncertainty[nbins];

			
	TCanvas *c1 = new TCanvas(); //canvas for eta of W+
	TCanvas *c2 = new TCanvas(); //canvas for eta of W-
	TCanvas *c3 = new TCanvas(); //canvas for eta of Z
	TCanvas *c4 = new TCanvas(); //Canvas for pt of W+
	TCanvas *c5 = new TCanvas(); //Canvas for pt of W-
	TCanvas *c6 = new TCanvas(); //Canvas for pt of Z

	TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
	TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
	TLegend *leg3 = new TLegend(0.6, 0.7, 0.9, 0.9);
	TLegend *leg4 = new TLegend(0.6, 0.7, 0.9, 0.9);
	TLegend *leg5 = new TLegend(0.6, 0.7, 0.9, 0.9);
	TLegend *leg6 = new TLegend(0.6, 0.7, 0.9, 0.9);

	//reference values for calculating ratios
	TGraph* RefUpperWmUnc = new TGraph(nbins);
	TGraph* RefLowerWmUnc =  new TGraph(nbins);
	TGraph* RefCentreWmUnc = new TGraph(nbins);

	TGraph* RefUpperWpUnc = new TGraph(nbins);
	TGraph* RefLowerWpUnc =  new TGraph(nbins);
	TGraph* RefCentreWpUnc = new TGraph(nbins);

	TGraph* RefUpperZUnc = new TGraph(nbins);
	TGraph* RefLowerZUnc =  new TGraph(nbins);
	TGraph* RefCentreZUnc = new TGraph(nbins);


	//plots with uncertainty
	int j = 0;
	for (std::string PDFID : PDFIDs)
	{
		int numVariations = AllnumVariations[j];
		TGraphAsymmErrors* WpUnc = new TGraphAsymmErrors(nbins);
		TGraphAsymmErrors* WmUnc = new TGraphAsymmErrors(nbins);
		TGraphAsymmErrors* ZUnc = new TGraphAsymmErrors(nbins);


		TGraph* WpCentre = new TGraph(nbins);
		TGraph* WmCentre = new TGraph(nbins);
		TGraph* ZCentre = new TGraph(nbins);


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
		

		//Uncertainty vectors for boson cross sections
		vector<double> W_pUnc(nbins, 0);
		vector<double> W_mUnc(nbins, 0);
		vector<double> Z_Unc(nbins, 0);

		//loop over each boson
		for (int k = 0; k < nbosons; ++k)
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
				double upperUncert, lowerUncert, SymmUncert;

				double zero = centralHist->GetBinContent(bin);
				if (PDFID == "3311")
				{
					double noOfVars = numVariations;

					for (int i = 0; i < numVariations; ++i)
					{
						//getting plus and minus variations from paired files
						double var  = variationHists[i]->GetBinContent(bin);

						sumUpperUncert += std::pow((var -zero), 2);
					}
					upperUncert = std::sqrt((1/(noOfVars-1))*sumUpperUncert);
					lowerUncert =  upperUncert;
					SymmUncert = upperUncert;
				}
				
				//Hessian pdf calculation
				else
				{
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

					upperUncert = std::sqrt(sumUpperUncert);
					lowerUncert = std::sqrt(sumLowerUncert);
					SymmUncert = std::sqrt(sumSymmUncert)*0.5;
				}


				y_lower_uncertainty_values[j][k][bin-1] = lowerUncert;
				y_upper_uncertainty_values[j][k][bin-1] = upperUncert;

				upperUncertaintyHist->SetBinContent(bin, upperUncert);
				lowerUncertaintyHist->SetBinContent(bin, lowerUncert);
				symmUncertaintyHist->SetBinContent(bin, SymmUncert);
			}
			

			WmUnc->GetYaxis()->SetTitle("d#sigma/d#eta (pb)");
			WmUnc->GetXaxis()->SetTitle("#eta");
			WpUnc->GetYaxis()->SetTitle("d#sigma/d#eta (pb)");
			WpUnc->GetXaxis()->SetTitle("#eta");
			ZUnc->GetYaxis()->SetTitle("d#sigma/d#eta (pb)");
			ZUnc->GetXaxis()->SetTitle("#eta");


			//plotting graphs
			for (int bin = 1; bin <= nbins; ++bin)
			{
				y_values[j][k][bin-1] = centralHist->GetBinContent(bin);
				x_values[bin-1] = centralHist->GetBinCenter(bin);

				y_ratio_values[j][k][bin-1] = y_values[j][k][bin-1]/y_values[0][k][bin-1];

				y_ratio_lower_uncertainty_values[j][k][bin-1] = y_ratio_values[j][k][bin-1]*((y_lower_uncertainty_values[j][k][bin-1])/y_values[j][k][bin-1]);
				y_ratio_upper_uncertainty_values[j][k][bin-1] = y_ratio_values[j][k][bin-1]*((y_upper_uncertainty_values[j][k][bin-1])/y_values[j][k][bin-1]);

				binwidth_vector[bin-1] = binwidth;
				stat_uncertainty[bin-1] = 0.;

				double x = centralHist->GetBinCenter(bin);
				double y = centralHist->GetBinContent(bin);
				double upperErr = upperUncertaintyHist->GetBinContent(bin);
				double lowerErr = lowerUncertaintyHist->GetBinContent(bin);
				double symmErr = symmUncertaintyHist->GetBinContent(bin);
				
				double refy =1.; //reference pdf value for ratio plot


				if (k==0)
				{

					WmUnc->SetPoint(bin-1, x, y);
					WmUnc->SetPointError(bin-1, 0, 0,lowerErr, upperErr);
					W_mUnc[bin-1] = symmErr;

					if (j != 0)
					{
						refy = RefCentreWmUnc->GetPointY(bin-1);
						WmCentre->SetPoint(bin-1, x, y/refy);
						UpperWmUncRat->SetPoint(bin-1, x, y/refy*(1+(upperErr)/y));
						LowerWmUncRat->SetPoint(bin-1, x, y/refy*(1-(lowerErr)/y));

						cout << y/refy << "*" <<"(1-(" << lowerErr << ")/"<< y <<"))= "<< y/refy*(1-(lowerErr)/y)<<endl;

						cout << y_ratio_values[j][k][bin-1] << "*" <<"(1-("
							<< y_lower_uncertainty_values[j][k][bin-1] << ")/"<< y_values[j][k][bin-1] <<"))= "
							<<y_ratio_values[j][k][bin-1]*(1-(y_lower_uncertainty_values[j][k][bin-1])/y_values[j][k][bin-1])<< endl;

						SymmWmUncRat->SetPoint(bin-1, x, y/refy+(std::abs(symmErr)/refy));
					}

					else if (j==0)
					{
						WmCentre->SetPoint(bin-1, x, y/refy);
						UpperWmUncRat->SetPoint(bin-1, x, y/y*(1+(upperErr)/y));
						cout << y << "/" << y << "*" <<"(1+(" << upperErr << ")/"<< y <<"))= "<< y/y*(1+(upperErr)/y)<< endl;
						LowerWmUncRat->SetPoint(bin-1, x, y/y*(1-(lowerErr)/y));
						SymmWmUncRat->SetPoint(bin-1, x, y/y+(std::abs(symmErr)/y));
					}


					std::cout << "bin: "<<bin << "	central:  " << y << endl;
					std::cout << "bin: "<<bin << "	upperErr: " << upperErr << endl;
					std::cout << "bin: "<<bin << "	lowerErr: " << lowerErr << endl;
				}

				else if (k==1)
				{

					WpUnc->SetPoint(bin-1, x, y);
					WpUnc->SetPointError(bin-1, 0, 0,lowerErr, upperErr);
					W_pUnc[bin-1] = symmErr;

					if (j != 0)
					{
						refy = RefCentreWpUnc->GetPointY(bin);
						WpCentre->SetPoint(bin-1, x, y/refy);
						UpperWpUncRat->SetPoint(bin-1, x, y/refy*(1+(upperErr)/y));
						LowerWpUncRat->SetPoint(bin-1, x, y/refy*(1-(lowerErr)/y));
						SymmWpUncRat->SetPoint(bin-1, x, y/refy+(std::abs(symmErr)/y));
					}

					else if (j==0)
					{
						UpperWpUncRat->SetPoint(bin-1, x, y/refy*(1+(upperErr)/y));
						LowerWpUncRat->SetPoint(bin-1, x, y/refy*(1-(lowerErr)/y));
						SymmWpUncRat->SetPoint(bin-1, x, y/y+(std::abs(symmErr)/y));
					}


				}

				else if (k==2)
				{
					ZUnc->SetPoint(bin-1, x, y);
					ZUnc->SetPointError(bin-1, 0, 0,lowerErr, upperErr);
					Z_Unc[bin-1] = symmErr;

					if (j != 0)
					{
						refy = RefCentreWmUnc->GetPointY(bin);
						ZCentre->SetPoint(bin-1, x, y/refy);
						UpperZUncRat->SetPoint(bin-1, x, y/refy*(1+(upperErr)/y));
						LowerZUncRat->SetPoint(bin-1, x, y/refy*(1-(lowerErr)/y));
						SymmZUncRat->SetPoint(bin-1, x, y/refy+(std::abs(symmErr)/refy));
					}

					else if (j==0)
					{
						UpperZUncRat->SetPoint(bin-1, x, y/refy*(1+(upperErr)/y));
						LowerZUncRat->SetPoint(bin-1, x, y/refy*(1-(lowerErr)/y));
						SymmZUncRat->SetPoint(bin-1, x, y/y+(std::abs(symmErr)/y));
					}
				}

			}


			if (j==0)
			{
				RefCentreWmUnc = (TGraph*)WmCentre->Clone("RefCentreWmUnc");
				RefCentreWpUnc = (TGraph*)WpCentre->Clone("RefCentreWpUnc");
				RefCentreZUnc = (TGraph*)ZCentre->Clone("RefCentreZUnc");

				for (int bin = 1; bin <= nbins; ++bin)
				{
					double x = centralHist->GetBinCenter(bin);
					double y = centralHist->GetBinContent(bin);

					WmCentre->SetPoint(bin-1, x, y/y);
					WpCentre->SetPoint(bin-1, x, y/y);
					ZCentre->SetPoint(bin-1, x, y/y);

				}
			

			}
			

			std::vector<int> colours = {2, 4, 6};
			std::vector<int> markers = {73,72,71};
			std::vector<int> fills = {3002,3004,3005};

			std::vector<const char*> PDFnames = {"CT18NNLO", "MSHT20NNLO", "NNPDF40NNLO"}; 

			//plotting different boson eta distributions
			//include central values and shade uncertainty properly
			if (k==0)
			{
				c1->cd();
				c1->SetGridy();
				
				WmUnc->SetTitle("W^{-}");
				WmUnc->SetMarkerStyle(markers[j]);
				WmUnc->SetMarkerColor(colours[j]);
				WmUnc->SetLineColor(colours[j]);
				WmUnc->SetMarkerSize(0.8);
				WmUnc->SetMinimum(0);
				WmUnc->SetMaximum(2000);
				if (j==0) WmUnc->Draw("AP");
				else WmUnc->Draw("P SAME");
				leg->AddEntry(WmUnc, PDFnames[j], "lep");
				if (j==(PDFIDs.size()-1)) leg->Draw("SAME");
				
			}
			else if (k==1)
			{
				c2->cd();
				c2->SetGridy();

				WpUnc->SetTitle("W^{+}");
				WpUnc->SetMarkerStyle(markers[j]);
				WpUnc->SetMarkerColor(colours[j]);
				WpUnc->SetLineColor(colours[j]);
				WpUnc->SetMarkerSize(0.8);
				if (j==0) WpUnc->Draw("AP");
				else WpUnc->Draw("P SAME");
				leg2->AddEntry(WpUnc, PDFnames[j], "lep");
				if (j==(PDFIDs.size()-1)) leg2->Draw("SAME");
			}
			else if (k==2)
			{
				c3->cd();
				c3->SetGridy();

				ZUnc->SetTitle("Z");
				ZUnc->SetMarkerStyle(markers[j]);
				ZUnc->SetMarkerColor(colours[j]);
				ZUnc->SetLineColor(colours[j]);
				ZUnc->SetMarkerSize(0.7);
				if (j==0) ZUnc->Draw("AP");
				else ZUnc->Draw("P SAME");
				leg3->AddEntry(ZUnc, PDFnames[j], "lep");
				if (j==(PDFIDs.size()-1)) leg3->Draw("SAME");

			}
			
			//Plotting uncertainty percentage for different bosons
			
			if (k==0)
			{
				c4->cd();
				c4->SetGridy();

				TGraphMultiErrors* Wmratio = new TGraphMultiErrors("Wmratio", "W^{-}", nbins, x_values, y_ratio_values[j][k], binwidth_vector, binwidth_vector, stat_uncertainty, stat_uncertainty);
				Wmratio->AddYError(nbins, y_ratio_lower_uncertainty_values[j][k], y_ratio_upper_uncertainty_values[j][k]);

				Wmratio->SetMarkerStyle(markers[j]);
				Wmratio->GetAttLine(0)->SetLineColor(kRed);
				Wmratio->GetAttFill(0)->SetFillStyle(fills[j]);
				Wmratio->GetAttFill(0)->SetFillColor(colours[j]);
				Wmratio->SetMarkerColor(colours[j]);

				Wmratio->GetAttLine(1)->SetLineColor(colours[j]);
				Wmratio->GetAttFill(1)->SetFillStyle(fills[j]);
				Wmratio->GetAttFill(1)->SetFillColor(colours[j]);
				Wmratio->SetMaximum(1.1);
				Wmratio->SetMinimum(.9);

				if (j==0) Wmratio->Draw("APS s=0.0 ; ; 5");
				else Wmratio->Draw("PS SAME s=0.0 ; ; 5");

				leg4->AddEntry(Wmratio, PDFnames[j], "fp");
				if (j==(PDFIDs.size()-1)) leg4->Draw();

				/*
				c5->cd();
				c5->SetGridy();
				UpperWmUncRat->SetTitle("W^{-}");
				UpperWmUncRat->SetMarkerStyle(73);
				UpperWmUncRat->SetMarkerColor(colours[j]);
				UpperWmUncRat->SetLineColor(colours[j]);
				if (j==0) UpperWmUncRat->Draw("AP");
				else  UpperWmUncRat->Draw("P SAME");
				UpperWmUncRat->GetYaxis()->SetTitle("Percentage Uncertainty (%)");
				UpperWmUncRat->GetXaxis()->SetTitle("y");
				UpperWmUncRat->SetMaximum(1.1);
				UpperWmUncRat->SetMinimum(.9);


				leg4->AddEntry(UpperWmUncRat, PDFnames[j], "lep");
				WmCentre->SetMarkerSize(52);
				WmCentre->Draw("P SAME");

				LowerWmUncRat->SetMarkerStyle(73);
				LowerWmUncRat->SetMarkerColor(colours[j]);
				LowerWmUncRat->SetLineColor(colours[j]);
				LowerWmUncRat->Draw("P SAME");
				if (j==(PDFIDs.size()-1)) leg4->Draw();
				*/
				
			}

			else if (k==1)

			{
				c5->cd();
				c5->SetGridy();
				TGraphMultiErrors* Wpratio = new TGraphMultiErrors("Wpratio", "W^{+}", nbins, x_values, y_ratio_values[j][k], binwidth_vector, binwidth_vector, stat_uncertainty, stat_uncertainty);
				Wpratio->AddYError(nbins, y_ratio_lower_uncertainty_values[j][k], y_ratio_upper_uncertainty_values[j][k]);

				Wpratio->SetMarkerStyle(markers[j]);
				Wpratio->GetAttLine(0)->SetLineColor(kRed);
				Wpratio->GetAttFill(0)->SetFillStyle(fills[j]);
				Wpratio->GetAttFill(0)->SetFillColor(colours[j]);
				Wpratio->SetMarkerColor(colours[j]);

				Wpratio->GetAttLine(1)->SetLineColor(colours[j]);
				Wpratio->GetAttFill(1)->SetFillStyle(fills[j]);
				Wpratio->GetAttFill(1)->SetFillColor(colours[j]);
				Wpratio->SetMaximum(1.1);
				Wpratio->SetMinimum(.9);

				if (j==0) Wpratio->Draw("APS s=0.0 ; ; 5");
				else Wpratio->Draw("PS SAME s=0.0 ; ; 5");

				leg5->AddEntry(Wpratio, PDFnames[j], "fp");
				if (j==(PDFIDs.size()-1)) leg5->Draw();

				/*
				UpperWpUncRat->SetTitle("W^{+}");
				UpperWpUncRat->SetMarkerStyle(72);
				UpperWpUncRat->SetMarkerColor(colours[j]);
				UpperWpUncRat->SetLineColor(colours[j]);
				UpperWpUncRat->SetMaximum(1.1);
				UpperWpUncRat->SetMinimum(.9);
				if (j==0) UpperWpUncRat->Draw("AP");
				else  UpperWpUncRat->Draw("P SAME");
				WpCentre->Draw("P SAME");
				leg5->AddEntry(UpperWpUncRat, PDFnames[j], "lep");

				LowerWpUncRat->SetMarkerStyle(72);
				LowerWpUncRat->SetMarkerColor(colours[j]);
				LowerWpUncRat->SetLineColor(colours[j]);
				LowerWpUncRat->Draw("P SAME");
				if (j==(PDFIDs.size()-1)) leg5->Draw();
				*/
			}

			else if (k==2)

			{
				c6->cd();
				c6->SetGridy();
				TGraphMultiErrors* Zratio = new TGraphMultiErrors("Wpratio", "W^{+}", nbins, x_values, y_ratio_values[j][k], binwidth_vector, binwidth_vector, stat_uncertainty, stat_uncertainty);
				Zratio->AddYError(nbins, y_ratio_lower_uncertainty_values[j][k], y_ratio_upper_uncertainty_values[j][k]);

				Zratio->SetMarkerStyle(markers[j]);
				Zratio->GetAttLine(0)->SetLineColor(kRed);
				Zratio->GetAttFill(0)->SetFillStyle(fills[j]);
				Zratio->GetAttFill(0)->SetFillColor(colours[j]);
				Zratio->SetMarkerColor(colours[j]);

				Zratio->GetAttLine(1)->SetLineColor(colours[j]);
				Zratio->GetAttFill(1)->SetFillStyle(fills[j]);
				Zratio->GetAttFill(1)->SetFillColor(colours[j]);
				Zratio->SetMaximum(1.1);
				Zratio->SetMinimum(.9);

				if (j==0) Zratio->Draw("APS s=0.0 ; ; 5");
				else Zratio->Draw("PS SAME s=0.0 ; ; 5");

				leg6->AddEntry(Zratio, PDFnames[j], "fp");
				if (j==(PDFIDs.size()-1)) leg6->Draw();

				/*
				UpperZUncRat->SetTitle("Z");
				UpperZUncRat->SetMarkerStyle(71);
				UpperZUncRat->SetMarkerColor(colours[j]);
				UpperZUncRat->SetLineColor(colours[j]);
				UpperZUncRat->SetMaximum(1.1);
				UpperZUncRat->SetMinimum(.9);
				if (j==0) UpperZUncRat->Draw("AP");
				else  UpperZUncRat->Draw("P SAME");
				leg6->AddEntry(UpperZUncRat, PDFnames[j], "lep");
				ZCentre->Draw("P SAME");

				LowerZUncRat->SetMarkerStyle(71);
				LowerZUncRat->SetMarkerColor(colours[j]);
				LowerZUncRat->SetLineColor(colours[j]);
				LowerZUncRat->Draw("P SAME");

				if (j==(PDFIDs.size()-1)) leg6->Draw();
				*/
			}

			//Clean up histograms (include or shit breaks)
			delete centralHist;
			delete uncertaintyHist;
			for (auto h : variationHists){delete h;}
		}
		j += 1;
	}


	TFile* outFile =new TFile("NLOMuonEtaDistComparision.root", "RECREATE");
	c1->Write();
	c2->Write();
	c3->Write();
	c4->Write();
	c5->Write();
	c6->Write();
	outFile->Close();
	//delete outFile;

	return 1;
    
}

