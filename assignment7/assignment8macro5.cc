//Calculating invariant mass cross secions and ratios WITH UNCERTAINTY using the reweighting method, CORRECTLY.
//Specifically looking at the percentage uncertainty of each pdf and comparing the bosons
#include <TFile.h>
#include <TNtuple.h>
#include <TH1D.h>
#include <TString.h>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>

#include <TPad.h>
#include <TLatex.h>

using namespace std;

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

    Float_t  eventNo, index, mother1, mother2, pAbs, pt, y, eta, id, particleCharge,mothermass1, mothermass2;
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
	muontuples->SetBranchAddress("charge", &particleCharge);
	muontuples->SetBranchAddress("mothermass1", &mothermass1);
	muontuples->SetBranchAddress("mothermass2", &mothermass2);


    
    for (int i = 0; i < particle_count; ++i)
	{
        muontuples->GetEntry(i);
        if ((((mother1 == 24) && (particleCharge==1)) || ((mother1==13) && (std::abs(mother2)==90))) //check if particle has W mom
        || (((mother1 == -24) && (particleCharge==-1)) || ((mother1==-13) && (std::abs(mother2)==90))) //check if particle has W mom
        || (((std::abs(mother1) == 23) && (particleCharge==1))  || ((std::abs(mother1)==13) && (std::abs(mother2)==90) && (particleCharge==1)))) //check if particle has a Z/gamma mom
		{
			//if ((std::abs(eta) <2.4) && (pt > 25)){
				int   EventNo = (int)eventNo;
				double w0 = ((*eventweights)[EventNo][0]);
				double wi = ((*eventweights)[EventNo][VarNo]);
				double w = wi/w0;
				hist->Fill(mothermass1, w);
			//}
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

int assignment8macro5() 
{
    //define histogram binning parameters
    const int nbins = 25;
	const int nbosons =3;
    const double xlow = 30.;
    const double xhigh = 120;
	const int nPDFs = 3;
	double_t binwidth = (xhigh-xlow)/(2*nbins);
	cout<<binwidth;
	

	std::vector<std::string> PDFIDs= {"140", "274", "3311"};
	//std::vector<std::string> PDFIDs= {"140"};//, "274", "3311"};
	const std::vector<int> AllnumVariations = {58, 64, 99}; 
	//const std::vector<int> AllnumVariations = {2};//, 2, 2}; 
	//const std::vector<int> AllTrueNumVariations = {58};//, 64, 99}; 
	const std::vector<int> AllTrueNumVariations = {58, 64, 99}; 
	int scaleVars = 6;
	//int scaleVars = 2;

	
	//vectors with the values and uncertainty of each of each boson
	double_t x_values[nbins];

	double_t y_values[nPDFs][nbosons][nbins];
	double_t y_upper_uncertainty_values[nPDFs][nbosons][nbins];
	double_t y_lower_uncertainty_values[nPDFs][nbosons][nbins];

	double_t y_upper_total_uncertainty_values[nPDFs][nbosons][nbins];
	double_t y_lower_total_uncertainty_values[nPDFs][nbosons][nbins];

	double_t y_lower_scale_uncertainty_values[nPDFs][nbosons][nbins];
	double_t y_upper_scale_uncertainty_values[nPDFs][nbosons][nbins];

	double_t y_ratio_values[nPDFs][nbosons][nbins];

	double_t y_ratio_upper_uncertainty_values[nPDFs][nbosons][nbins];
	double_t y_ratio_lower_uncertainty_values[nPDFs][nbosons][nbins];

	double_t y_ratio_upper_scale_uncertainty_values[nPDFs][nbosons][nbins];
	double_t y_ratio_lower_scale_uncertainty_values[nPDFs][nbosons][nbins];
	
	double_t y_ratio_upper_total_uncertainty_values[nPDFs][nbosons][nbins];
	double_t y_ratio_lower_total_uncertainty_values[nPDFs][nbosons][nbins];

	double_t binwidth_vector[nbins];
	double_t stat_uncertainty[nbins];

			
	TCanvas *c1 = new TCanvas(); //canvas for eta of W+
	TCanvas *c2 = new TCanvas(); //canvas for eta of W-
	TCanvas *c3 = new TCanvas(); //canvas for eta of Z
	TCanvas *c4 = new TCanvas(); //canvas for eta of W+
	TCanvas *c5 = new TCanvas(); //canvas for eta of W-
	TCanvas *c6 = new TCanvas(); //canvas for eta of Z

								 

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
		int TruenumVariations = AllTrueNumVariations[j]; //for scale uncertainty calculation
													   
		TGraphAsymmErrors* WpUnc = new TGraphAsymmErrors(nbins);
		TGraphAsymmErrors* WmUnc = new TGraphAsymmErrors(nbins);
		TGraphAsymmErrors* ZUnc = new TGraphAsymmErrors(nbins);


		//Uncertainty vectors for boson cross sections
		vector<double> W_pUnc(nbins, 0);
		vector<double> W_mUnc(nbins, 0);
		vector<double> Z_Unc(nbins, 0);

		//loop over each boson
		for (int k = 2; k < nbosons; ++k)
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

			TString centralFilename = PDFID+"00_pythia_reweighted_output/"+boson+"_pwgevents_"+PDFID+"00.root";
			TH1D* centralHist = new TH1D("centralHist", "Central Cross Section; #eta; d#sigma/d#eta", nbins, xlow, xhigh);
			if (processFile(centralFilename.Data(), centralHist, 0) < 0) {
				std::cerr << "Error processing central file." << std::endl;
				return 1;
			}
			
			//processing uncertainty files
			std::vector<TH1D*> variationHists;
			std::vector<TH1D*> scaleVarHists;

			for (int i = 1; i <= numVariations; ++i)
			{
				cout << i<< endl;
				//create a histogram for this variation (clone the binning of centralHist)
				TH1D* hVar = new TH1D(Form("variation_%d", i), "Variation", nbins, xlow, xhigh);
				if (processFile(centralFilename.Data(), hVar, i) < 0) 
				{
					std::cerr << "Error processing file: " << centralFilename << std::endl;
					continue;
				}
				variationHists.push_back(hVar);
				if (i==0) scaleVarHists.push_back(hVar);
			}

			//calculating variation distributions for scales variations
			for (int i = TruenumVariations+1; i <= TruenumVariations+scaleVars; ++i)
			{
				cout << i<< endl;
				TH1D* hVar = new TH1D(Form("variation_%d", i), "Variation", nbins, xlow, xhigh);
				if (processFile(centralFilename.Data(), hVar, i) < 0) 
				{
					std::cerr << "Error processing file: " << centralFilename << std::endl;
					continue;
				}
				if (i > numVariations || i==0) scaleVarHists.push_back(hVar);
				
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
				std::vector<double> ScaleVars; //values for scales in one bin

				double ScaleUpperUncert = 0.0;
				double ScaleLowerUncert = 0.0;

				double sumUpperUncert = 0.0;
				double sumLowerUncert = 0.0;
				double sumSymmUncert = 0.0;
				double upperUncert, lowerUncert, SymmUncert;

				double zero = centralHist->GetBinContent(bin);


				//Scale uncertainty
				for (int i = 0; i < scaleVars; ++i)
				{
					double scaleVar = scaleVarHists[i]->GetBinContent(bin);
					ScaleVars.push_back(scaleVar);
				}
				//Take the largest and smallest values to be the upper and lower uncertainty respectively.
				ScaleUpperUncert = *max_element(ScaleVars.begin(),ScaleVars.end()) - zero;

				//cout << "central value " << bin << ": " << zero << endl;
				//cout << "scale variation max value " << ": " << ScaleVars[4] << endl;
				ScaleLowerUncert = *min_element(ScaleVars.begin(),ScaleVars.end()) - zero;



				
				//NNPDF/Gaussian pdf calculation 
				if ((PDFID == "3311") ||(PDFID == "3317"))
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

				
				y_lower_scale_uncertainty_values[j][k][bin-1] = std::abs(ScaleLowerUncert);
				y_upper_scale_uncertainty_values[j][k][bin-1] = ScaleUpperUncert;

				if (j==0) 
				{
				y_lower_uncertainty_values[j][k][bin-1] = lowerUncert/1.645;
				y_upper_uncertainty_values[j][k][bin-1] = upperUncert/1.645;

				y_lower_total_uncertainty_values[j][k][bin-1] = std::sqrt(std::pow(lowerUncert/1.645,2) +std::pow(ScaleLowerUncert,2));
				y_upper_total_uncertainty_values[j][k][bin-1] = std::sqrt(std::pow(upperUncert/1.645,2) +std::pow(ScaleUpperUncert,2));

				}

				else
				{
				y_lower_uncertainty_values[j][k][bin-1] = lowerUncert;
				y_upper_uncertainty_values[j][k][bin-1] = upperUncert;

				y_lower_total_uncertainty_values[j][k][bin-1] = std::sqrt(std::pow(lowerUncert,2) +std::pow(ScaleLowerUncert,2));
				y_upper_total_uncertainty_values[j][k][bin-1] = std::sqrt(std::pow(upperUncert,2) +std::pow(ScaleUpperUncert,2));
				}
				
				upperUncertaintyHist->SetBinContent(bin, y_upper_total_uncertainty_values[j][k][bin-1]);
				lowerUncertaintyHist->SetBinContent(bin, y_lower_total_uncertainty_values[j][k][bin-1]);
				symmUncertaintyHist->SetBinContent(bin, SymmUncert);
			}
			

			WmUnc->GetYaxis()->SetTitle("d#sigma/dy (pb)");
			WmUnc->GetXaxis()->SetTitle("rapidity (y)");
			WpUnc->GetYaxis()->SetTitle("d#sigma/d#eta (pb)");
			ZUnc->GetYaxis()->SetTitle("d#sigma/dM_{Z} (pb/Gev/c^{2})");
			ZUnc->GetXaxis()->SetTitle("M_{Z} (GeV/c^{2})");


			//plotting graphs
			for (int bin = 1; bin <= nbins; ++bin)
			{
				y_values[j][k][bin-1] = centralHist->GetBinContent(bin);
				x_values[bin-1] = centralHist->GetBinCenter(bin);
				
				//slightly shifting points so that the uncertainty lines don't completely overlap
				if (j==1){x_values[bin-1] += 0.02;}
				if (j==2){x_values[bin-1] += -0.02;}
				

				y_ratio_values[j][k][bin-1] = y_values[j][k][bin-1]/y_values[0][k][bin-1];

				y_ratio_lower_uncertainty_values[j][k][bin-1] = y_ratio_values[j][k][bin-1]*((y_lower_uncertainty_values[j][k][bin-1])/y_values[j][k][bin-1]);
				y_ratio_upper_uncertainty_values[j][k][bin-1] = y_ratio_values[j][k][bin-1]*((y_upper_uncertainty_values[j][k][bin-1])/y_values[j][k][bin-1]);


				y_ratio_lower_scale_uncertainty_values[j][k][bin-1] = y_ratio_values[j][k][bin-1]*((y_lower_scale_uncertainty_values[j][k][bin-1])/y_values[j][k][bin-1]);
				y_ratio_upper_scale_uncertainty_values[j][k][bin-1] = y_ratio_values[j][k][bin-1]*((y_upper_scale_uncertainty_values[j][k][bin-1])/y_values[j][k][bin-1]);
				cout << "lower_scale unc_rat: " << y_ratio_lower_scale_uncertainty_values[j][k][bin-1] << endl;
				cout << "upper_scale unc_rat: " << y_ratio_upper_scale_uncertainty_values[j][k][bin-1] << endl;

				y_ratio_lower_total_uncertainty_values[j][k][bin-1] = y_ratio_values[j][k][bin-1]*((y_lower_total_uncertainty_values[j][k][bin-1])/y_values[j][k][bin-1]);
				y_ratio_upper_total_uncertainty_values[j][k][bin-1] = y_ratio_values[j][k][bin-1]*((y_upper_total_uncertainty_values[j][k][bin-1])/y_values[j][k][bin-1]);

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
				}

				else if (k==1)
				{
					WpUnc->SetPoint(bin-1, x, y);
					WpUnc->SetPointError(bin-1, 0, 0,lowerErr, upperErr);
					W_pUnc[bin-1] = symmErr;
				}

				else if (k==2)
				{
					ZUnc->SetPoint(bin-1, x, y);
					ZUnc->SetPointError(bin-1, 0, 0,lowerErr, upperErr);
					Z_Unc[bin-1] = symmErr;
				}

			}

			

			std::vector<int> colours = {2, 6, 4};
			std::vector<int> markercolours = {2, 6, 4};
			std::vector<int> markers = {73,72,71};
			std::vector<int> fills = {3002,3004,3005};

			std::vector<const char*> PDFnames = {"CT18NNLO", "MSHT20NNLO", "NNPDF40NNLO"}; 
			//std::vector<const char*> PDFnames = {"MSHT20NLO", "CT18NLO", "NNPDF40NLO"}; 

			//plotting different boson eta distributions
			//include central values and shade uncertainty properly
			if (j==0)
			{
					ZUnc->SetTitle("");
					ZUnc->SetMarkerStyle(markers[k]);
					ZUnc->SetMarkerColor(markercolours[k]);
					ZUnc->SetFillStyle(fills[k]);
					ZUnc->SetLineColor(markercolours[k]);
					ZUnc->SetMarkerSize(0.8);
					
					
					ZUnc->GetYaxis()->SetRangeUser(0., 250.);
					ZUnc->GetXaxis()->SetRangeUser(xlow, xhigh);
					/*
				if (k==0)
				{
					c1->cd();
					//upperPad1->Draw();
					//upperPad1->cd();
					
					//WmUnc->SetTitle("W^{-}");
					if (j==0)
					{
						//random text stuff to go into plot
						TLatex randomtexts;
						TLatex texts;
						texts.SetTextSize(0.035);
						randomtexts.SetTextSize(0.025);
						randomtexts.DrawLatex(-3.8, 2000, "Muon Spectrometer");
						randomtexts.DrawLatex(-0.9, 2000, "Central Barrel");
						texts.DrawLatex(-4.5, 2350, "p-p @ #sqrt{s} = 13.6 TeV");
						texts.DrawLatex(-4.5, 2250, "POWHEG+PYTHIA8");

						WmUnc->Draw("ASP");
					}
					else WmUnc->Draw("SP SAME");
					leg->AddEntry(WmUnc, "#mu #leftarrow W^{-}", "lep");
					if (j==(PDFIDs.size()-1)) leg->Draw("SAME");
					
				}
				else if (k==1)
				{
					c1->cd();
					//upperPad2->Draw();
					//upperPad2->cd();

					//WpUnc->SetTitle("W^{+}");
					WpUnc->SetMarkerStyle(markers[k]);
					WpUnc->SetMarkerColor(markercolours[k]);
					WpUnc->SetLineColor(markercolours[k]);
					WpUnc->SetMarkerSize(0.8);
					WpUnc->GetXaxis()->SetRangeUser(xlow, xhigh);
					if (j==0) WpUnc->Draw("SP SAME");
					//else WpUnc->Draw("SP SAME");
					leg->AddEntry(WpUnc, "#mu #leftarrow W^{+}", "lep");
					if (j==(PDFIDs.size()-1)) leg2->Draw("SAME");
				}
				*/
				if (k==2)
				{
					c1->cd();
					c1->SetGridy();
					//upperPad3->Draw();
					//upperPad3->cd();

					//random text stuff to go into plot

					ZUnc->SetMarkerStyle(markers[k]);
					ZUnc->SetMarkerColor(markercolours[k]);
					ZUnc->SetLineColor(markercolours[k]);
					ZUnc->SetMarkerSize(0.8);
					ZUnc->GetXaxis()->SetRangeUser(xlow, xhigh);

					ZUnc->GetXaxis()->SetTitleSize(0.045);
					ZUnc->GetYaxis()->SetTitleSize(0.05);
					if (j==0) ZUnc->Draw("ASP");

					TLatex texts;
					texts.SetTextAlign(12);
					texts.SetTextSize(0.035);
					texts.DrawLatex(100, 100, "p-p @ #sqrt{s} = 13.6 TeV");
					texts.DrawLatex(100, 70, "POWHEG+PYTHIA8");
					texts.DrawLatex(100, 50, "CT18NNLO");

					c1->Write();

					//leg->AddEntry(ZUnc, "#mu #leftarrow Z", "lep");

				}
			}
			
			//Plotting uncertainty percentage for different bosons
			
			if (k==0)
			{
				c4->cd();
				//lowerPad1->Draw();
				//lowerPad1->cd();

				TGraphMultiErrors* Wmratio = new TGraphMultiErrors("Wmratio", "W^{-}", nbins, x_values,
						y_ratio_values[j][k], binwidth_vector, binwidth_vector, y_ratio_lower_scale_uncertainty_values[j][k], y_ratio_upper_scale_uncertainty_values[j][k]);
				Wmratio->AddYError(nbins, y_ratio_lower_uncertainty_values[j][k], y_ratio_upper_uncertainty_values[j][k]);

				Wmratio->SetMarkerStyle(markers[j]);
				Wmratio->GetAttLine(0)->SetLineColor(colours[j]);
				Wmratio->GetAttFill(0)->SetFillStyle(fills[j]);
				Wmratio->GetAttFill(0)->SetFillColor(colours[j]);
				Wmratio->SetMarkerColor(colours[j]);
				Wmratio->GetXaxis()->SetTitle("#eta");
				Wmratio->SetLineColor(colours[j]);
				Wmratio->SetFillStyle(fills[j]);
				Wmratio->SetFillColor(colours[j]);

				Wmratio->GetAttLine(1)->SetLineColor(colours[j]);
				Wmratio->GetAttFill(1)->SetFillStyle(fills[j]);
				Wmratio->GetAttFill(1)->SetFillColor(colours[j]);
				Wmratio->SetMaximum(1.15);
				Wmratio->SetMinimum(.85);
				Wmratio->SetMarkerSize(0.8);
				Wmratio->GetXaxis()->SetRangeUser(xlow, xhigh);
				Wmratio->GetYaxis()->SetTitle("ratio to CT18NNLO");

				if (j==0) Wmratio->Draw("APSE s=0.0 ; ; 5");
				else Wmratio->Draw("PSE SAME s=0.0 ; ; 5");

				leg4->AddEntry(Wmratio, PDFnames[j], "PEFL");
				if (j==(PDFIDs.size()-1)) leg4->Draw();
				
			}

			else if (k==1)

			{
				c5->cd();
				//lowerPad2->Draw();
				//lowerPad2->cd();

				TGraphMultiErrors* Wpratio = new TGraphMultiErrors("Wpratio", "W^{+}", nbins, x_values,
						y_ratio_values[j][k], binwidth_vector, binwidth_vector, y_ratio_lower_scale_uncertainty_values[j][k], y_ratio_upper_scale_uncertainty_values[j][k]);



				Wpratio->AddYError(nbins, y_ratio_lower_uncertainty_values[j][k], y_ratio_upper_uncertainty_values[j][k]);

				Wpratio->SetMarkerStyle(markers[j]);
				Wpratio->GetAttLine(0)->SetLineColor(colours[j]);
				Wpratio->GetAttFill(0)->SetFillStyle(fills[j]);
				Wpratio->GetAttFill(0)->SetFillColor(colours[j]);
				Wpratio->SetMarkerColor(colours[j]);
				Wpratio->SetLineColor(colours[j]);
				Wpratio->SetFillStyle(fills[j]);
				Wpratio->SetFillColor(colours[j]);

				Wpratio->GetAttLine(1)->SetLineColor(colours[j]);
				Wpratio->GetAttFill(1)->SetFillStyle(fills[j]);
				Wpratio->GetAttFill(1)->SetFillColor(colours[j]);
				Wpratio->SetMaximum(1.15);
				Wpratio->SetMinimum(.85);
				Wpratio->SetMarkerSize(0.8);
				Wpratio->GetXaxis()->SetRangeUser(xlow, xhigh);
				Wpratio->GetXaxis()->SetTitle("#eta");
				Wpratio->GetYaxis()->SetTitle("#frac{d#sigma_{PDF_{i}}}{dM}/#frac{d#sigma_{MSHT}}{dM}");

				if (j==0) Wpratio->Draw("APSE s=0.0 ; ; 5");
				else Wpratio->Draw("PSE SAME s=0.0 ; ; 5");

				leg5->AddEntry(Wpratio, PDFnames[j], "PEFL");
				if (j==(PDFIDs.size()-1)) leg5->Draw();

			}

			else if (k==2)

			{
				c6->cd();
				//lowerPad3->Draw();
				//lowerPad3->cd();

				TGraphMultiErrors* Zratio = new TGraphMultiErrors("Zratio", "Z", nbins, x_values,
						y_ratio_values[j][k], binwidth_vector, binwidth_vector, y_ratio_lower_scale_uncertainty_values[j][k],  y_ratio_upper_scale_uncertainty_values[j][k]);

				Zratio->AddYError(nbins, y_ratio_lower_uncertainty_values[j][k], y_ratio_upper_uncertainty_values[j][k]);

				Zratio->SetMarkerStyle(markers[j]);
				Zratio->GetAttLine(0)->SetLineColor(colours[j]);
				Zratio->GetAttFill(0)->SetFillStyle(fills[j]);
				Zratio->GetAttFill(0)->SetFillColor(colours[j]);
				Zratio->SetMarkerColor(colours[j]);
				Zratio->SetLineColor(colours[j]);
				Zratio->SetFillStyle(fills[j]);
				Zratio->SetFillColor(colours[j]);

				Zratio->GetAttLine(1)->SetLineColor(colours[j]);
				Zratio->GetAttFill(1)->SetFillStyle(fills[j]);
				Zratio->GetAttFill(1)->SetFillColor(colours[j]);
				Zratio->SetMaximum(1.2);
				Zratio->SetMinimum(.8);
				Zratio->SetMarkerSize(0.8);
				Zratio->GetXaxis()->SetRangeUser(xlow, xhigh);
				Zratio->GetXaxis()->SetTitle("M_{Z} (GeV/c^{2})");
				Zratio->GetXaxis()->SetTitleSize(0.035);
				Zratio->GetYaxis()->SetTitle("#frac{d#sigma_{PDF_{i}}}{dM}/#frac{d#sigma_{MSHT}}{dM}");
				Zratio->GetYaxis()->SetTitleSize(0.045);

				if (j==0) Zratio->Draw("APSE s=0.0 ; ; 5");
				else Zratio->Draw("PSE SAME s=0.0 ; ; 5");

				leg6->AddEntry(Zratio, PDFnames[j], "PEFL");
				if (j==(PDFIDs.size()-1)) leg6->Draw();

			}

			//Clean up histograms (getting rid of any potential memory leaks)
			delete centralHist;
			delete uncertaintyHist;
			for (auto h : variationHists){delete h;}
			TFile* outFileScales = new TFile("Scalevariations.root", "RECREATE");
			for (auto h : scaleVarHists){h->Write();}
			outFileScales->Close();
			for (auto h : scaleVarHists){delete h;}
		}
		j += 1;
	}


	TFile* outFile =new TFile("NLOMuonInvMassDistComparisionWithScales.root", "RECREATE");
	c2->Write();
	c3->Write();
	c4->Write();
	c5->Write();
	c6->Write();
	outFile->Close();
	//delete outFile;

	return 1;
    
}

