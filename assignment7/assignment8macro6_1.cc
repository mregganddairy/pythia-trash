//Calculating transverse momentum cross secions WITH UNCERTAINTY using the reweighting method, CORRECTLY.
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

    Float_t  eventNo, index, mother1, mother2, pAbs, pt, y, eta, id, particleCharge;
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


    
    for (int i = 0; i < particle_count; ++i)
	{
        muontuples->GetEntry(i);
        if ((((mother1 == 24) && (particleCharge==1)) || ((mother1==13) && (std::abs(mother2)==90))) //check if particle has W mom
        || (((mother1 == -24) && (particleCharge==-1)) || ((mother1==-13) && (std::abs(mother2)==90))) //check if particle has W mom
        || (((std::abs(mother1) == 23) && (particleCharge==1))  || ((std::abs(mother1)==13) && (std::abs(mother2)==90) && (particleCharge==1)))) //check if particle has a Z/gamma mom
		{
			if ((std::abs(eta) <4.0) && (std::abs(eta) >2.5)){
				int   EventNo = (int)eventNo;
				double w0 = ((*eventweights)[EventNo][0]);
				double wi = ((*eventweights)[EventNo][VarNo]);
				double w = wi/w0;
				hist->Fill(pt, w);
			}
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

//Ratio functions
double ratiofunction(double Wp, double Wm, double Z, int ratFuncNo)
{
	if (ratFuncNo == 0)
	{
		//Apm
		double output = (Wp-Wm)/(Wp+Wm);
		return output;
	}
	if (ratFuncNo == 1)
	{
		//Rpm
		double output = Wm/Wp;
		return output;
	}
	if (ratFuncNo == 2)
	{
		//RWZ
		double output = (Z)/(Wp+Wm);
		return output;
	}
	else
	{
		cout << "shits the bed";
		exit(0);
	}
}


int assignment8macro6_1() 
{
    //define histogram binning parameters
    const int nbins = 50;
	const int nbosons =3;
    const double xlow = 0.;
    const double xhigh = 100.;
	const int nPDFs = 3;
	double_t binwidth = (xhigh-xlow)/(2*nbins);
	cout<<binwidth;
	

	std::vector<std::string> PDFIDs= {"140", "274", "3311"};
	const std::vector<int> AllnumVariations = {58, 64, 99}; 
	//const std::vector<int> AllnumVariations = {58, 2, 2}; 
	const std::vector<int> AllTrueNumVariations = {58, 64, 99}; 
	int scaleVars = 6;
	//int scaleVars = 1;

	
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



	//plots with uncertainty
	int j = 0;
	for (std::string PDFID : PDFIDs)
	{
		int numVariations = AllnumVariations[j];
		int TruenumVariations = AllTrueNumVariations[j]; //for scale uncertainty calculation
													   
		
		/*
		TGraphAsymmErrors* ZUnc = new TGraphAsymmErrors(nbins);
		*/


		TGraphAsymmErrors* RpmUnc = new TGraphAsymmErrors(nbins);
		TGraphAsymmErrors* ApmUnc = new TGraphAsymmErrors(nbins);
		TGraphAsymmErrors* RZWUnc = new TGraphAsymmErrors(nbins);


		vector<TString> centralFilenames;
		vector<TH1D*> centralHists;

		vector<vector<TH1D*>> BIGvariationHists;
		vector<vector<TH1D*>> BIGscaleVarHists;

		//loop over each boson
		for (int k = 0; k < nbosons; ++k)
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

			//processing all bosons at once
			//

			TString centralFilename = PDFID+"00_pythia_reweighted_output/"+boson+"_pwgevents_"+PDFID+"00.root";
			TH1D* centralHist = new TH1D("centralHist", "Central Cross Section; #eta; d#sigma/d#eta", nbins, xlow, xhigh);
			if (processFile(centralFilename.Data(), centralHist, 0) < 0) {
				std::cerr << "Error processing central file." << std::endl;
				return 1;
			}

			centralFilenames.push_back(centralFilename);
			centralHists.push_back(centralHist);
			
			//processing uncertainty files
			vector<TH1D*> variationHists;
			vector<TH1D*> scaleVarHists;

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
				cout << j << boson << ": " << i <<endl;
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
				cout << j << boson << ": " << i <<endl;
				
			}
			BIGvariationHists.push_back(variationHists);
			BIGscaleVarHists.push_back(scaleVarHists);
			
			/*
			//Clean up histograms (getting rid of any potential memory leaks)
			delete centralHist;
			//delete uncertaintyHist;
			for (auto h : variationHists){delete h;}
			for (auto h : scaleVarHists){delete h;}
			*/
			
		}
			

		for (int ratFuncNo = 0; ratFuncNo <= 2; ratFuncNo++)
		{

			const int numPairs = BIGvariationHists.at(ratFuncNo).size() / 2;  
			cout << numPairs <<endl;
			TH1D* uncertaintyHist = (TH1D*)centralHists.at(ratFuncNo)->Clone("uncertaintyHist");
			TH1D* upperUncertaintyHist = (TH1D*)centralHists.at(ratFuncNo)->Clone("uncertaintyHist");
			TH1D* lowerUncertaintyHist = (TH1D*)centralHists.at(ratFuncNo)->Clone("uncertaintyHist");
			TH1D* symmUncertaintyHist = (TH1D*)centralHists.at(ratFuncNo)->Clone("uncertaintyHist");
			uncertaintyHist->Reset();


			//compute the uncertainty in each bin
			for (int bin = 1; bin <= nbins; bin++)
			{
				vector<double> ScaleVars; //values for scales in one bin

				double ScaleUpperUncert = 0.0;
				double ScaleLowerUncert = 0.0;

				double sumUpperUncert = 0.0;
				double sumLowerUncert = 0.0;
				double sumSymmUncert = 0.0;
				double upperUncert, lowerUncert, SymmUncert;

				double zerowm = centralHists.at(0)->GetBinContent(bin);
				double zerowp = centralHists.at(1)->GetBinContent(bin);
				double zeroz = centralHists.at(2)->GetBinContent(bin);

				double zero = ratiofunction(zerowp, zerowm, zeroz, ratFuncNo);



				//Scale uncertainty
				for (int i = 0; i < scaleVars; ++i)
				{
					double scaleVarwm = BIGscaleVarHists[0][i]->GetBinContent(bin);
					double scaleVarwp = BIGscaleVarHists[1][i]->GetBinContent(bin);
					double scaleVarz = BIGscaleVarHists[2][i]->GetBinContent(bin);

					double scale_var = ratiofunction(scaleVarwp,scaleVarwm, scaleVarz, ratFuncNo);

					ScaleVars.push_back(scale_var);
				}
				
				//Take the largest and smallest values to be the upper and lower uncertainty respectively.
				ScaleUpperUncert = *max_element(ScaleVars.begin(),ScaleVars.end()) - zero;
				ScaleLowerUncert = *min_element(ScaleVars.begin(),ScaleVars.end()) - zero;


				
				//NNPDF/Monte Carlo pdf calculation 
				if ((PDFID == "3311") ||(PDFID == "3317"))
				{
					double noOfVars = numVariations;

					for (int i = 0; i < numVariations; ++i)
					{
						//getting plus and minus variations from paired files
						double varwm  = BIGvariationHists[0][i]->GetBinContent(bin);
						double varwp  = BIGvariationHists[1][i]->GetBinContent(bin);
						double varz  = BIGvariationHists[2][i]->GetBinContent(bin);

						double var = ratiofunction(varwp, varwm, varz, ratFuncNo);

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
						double pluswm  = BIGvariationHists.at(0)[2*i]->GetBinContent(bin);
						double pluswp  = BIGvariationHists.at(1)[2*i]->GetBinContent(bin);
						double plusz  = BIGvariationHists.at(2)[2*i]->GetBinContent(bin);

						double minuswm = BIGvariationHists.at(0)[2*i + 1]->GetBinContent(bin);
						double minuswp = BIGvariationHists.at(1)[2*i + 1]->GetBinContent(bin);
						double minusz = BIGvariationHists.at(2)[2*i + 1]->GetBinContent(bin);

						double plus = ratiofunction(pluswp, pluswm, plusz, ratFuncNo);
						double minus = ratiofunction(minuswp, minuswm, minusz, ratFuncNo);


						//compute the difference according to the asymmetric Hessian prescription
						sumUpperUncert += std::pow(std::max(std::max(plus - zero, minus - zero), 0.0),2);
						sumLowerUncert += std::pow(std::max(std::max(zero - plus, zero - minus), 0.0),2);
						sumSymmUncert += std::pow(plus - minus,2);

					}

					upperUncert = std::sqrt(sumUpperUncert);
					lowerUncert = std::sqrt(sumLowerUncert);
					SymmUncert = std::sqrt(sumSymmUncert)*0.5;
				}

				
				y_lower_scale_uncertainty_values[j][ratFuncNo][bin-1] = std::abs(ScaleLowerUncert);
				y_upper_scale_uncertainty_values[j][ratFuncNo][bin-1] = ScaleUpperUncert;

				if (j==0) 
				{
				y_lower_uncertainty_values[j][ratFuncNo][bin-1] = lowerUncert/1.645;
				y_upper_uncertainty_values[j][ratFuncNo][bin-1] = upperUncert/1.645;

				y_lower_total_uncertainty_values[j][ratFuncNo][bin-1] = std::sqrt(std::pow(lowerUncert/1.645,2) +std::pow(ScaleLowerUncert,2));
				y_upper_total_uncertainty_values[j][ratFuncNo][bin-1] = std::sqrt(std::pow(upperUncert/1.645,2) +std::pow(ScaleUpperUncert,2));

				}

				else
				{
				y_lower_uncertainty_values[j][ratFuncNo][bin-1] = lowerUncert;
				y_upper_uncertainty_values[j][ratFuncNo][bin-1] = upperUncert;

				y_lower_total_uncertainty_values[j][ratFuncNo][bin-1] = std::sqrt(std::pow(lowerUncert,2) +std::pow(ScaleLowerUncert,2));
				y_upper_total_uncertainty_values[j][ratFuncNo][bin-1] = std::sqrt(std::pow(upperUncert,2) +std::pow(ScaleUpperUncert,2));
				}
				
				upperUncertaintyHist->SetBinContent(bin, y_upper_total_uncertainty_values[j][ratFuncNo][bin-1]);
				lowerUncertaintyHist->SetBinContent(bin, y_lower_total_uncertainty_values[j][ratFuncNo][bin-1]);
				symmUncertaintyHist->SetBinContent(bin, SymmUncert);
			}

			

		TGraphAsymmErrors* ApmUnc = new TGraphAsymmErrors(nbins);
		TGraphAsymmErrors* RpmUnc = new TGraphAsymmErrors(nbins);
		TGraphAsymmErrors* RZWUnc = new TGraphAsymmErrors(nbins);

			ApmUnc->GetYaxis()->SetTitle("A_{#pm}(p_{T}) ");
			ApmUnc->GetXaxis()->SetTitle("p_{T} (GeV/c)");

			RpmUnc->GetYaxis()->SetTitle("R_{#pm}(p_{T}) ");
			RpmUnc->GetXaxis()->SetTitle("p_{T} (GeV/c)");

			RZWUnc->GetYaxis()->SetTitle("R_{Z/W}(p_{T}) ");
			RZWUnc->GetXaxis()->SetTitle("p_{T} (GeV/c)");



			//plotting graphs
			for (int bin = 1; bin <= nbins; ++bin)
			{
				

				double wm = centralHists.at(0)->GetBinContent(bin);
				double wp = centralHists.at(1)->GetBinContent(bin);
				double z = centralHists.at(2)->GetBinContent(bin);

				y_values[j][ratFuncNo][bin-1] = ratiofunction(wp, wm, z, ratFuncNo);
				x_values[bin-1] = centralHists.at(ratFuncNo)->GetBinCenter(bin);
				
				//slightly shifting points so that the uncertainty lines don't completely overlap
				if (j==1){x_values[bin-1] += 0.02;}
				if (j==2){x_values[bin-1] += -0.02;}
				

				y_ratio_values[j][ratFuncNo][bin-1] = y_values[j][ratFuncNo][bin-1]/y_values[0][ratFuncNo][bin-1];

				y_ratio_lower_uncertainty_values[j][ratFuncNo][bin-1] = y_ratio_values[j][ratFuncNo][bin-1]*((y_lower_uncertainty_values[j][ratFuncNo][bin-1])/y_values[j][ratFuncNo][bin-1]);
				y_ratio_upper_uncertainty_values[j][ratFuncNo][bin-1] = y_ratio_values[j][ratFuncNo][bin-1]*((y_upper_uncertainty_values[j][ratFuncNo][bin-1])/y_values[j][ratFuncNo][bin-1]);


				y_ratio_lower_scale_uncertainty_values[j][ratFuncNo][bin-1] = y_ratio_values[j][ratFuncNo][bin-1]*((y_lower_scale_uncertainty_values[j][ratFuncNo][bin-1])/y_values[j][ratFuncNo][bin-1]);
				y_ratio_upper_scale_uncertainty_values[j][ratFuncNo][bin-1] = y_ratio_values[j][ratFuncNo][bin-1]*((y_upper_scale_uncertainty_values[j][ratFuncNo][bin-1])/y_values[j][ratFuncNo][bin-1]);
				cout << "lower_scale unc_rat: " << y_ratio_lower_scale_uncertainty_values[j][ratFuncNo][bin-1] << endl;
				cout << "upper_scale unc_rat: " << y_ratio_upper_scale_uncertainty_values[j][ratFuncNo][bin-1] << endl;

				y_ratio_lower_total_uncertainty_values[j][ratFuncNo][bin-1] = y_ratio_values[j][ratFuncNo][bin-1]*((y_lower_total_uncertainty_values[j][ratFuncNo][bin-1])/y_values[j][ratFuncNo][bin-1]);
				y_ratio_upper_total_uncertainty_values[j][ratFuncNo][bin-1] = y_ratio_values[j][ratFuncNo][bin-1]*((y_upper_total_uncertainty_values[j][ratFuncNo][bin-1])/y_values[j][ratFuncNo][bin-1]);

				binwidth_vector[bin-1] = binwidth;
				stat_uncertainty[bin-1] = 0.;

				double x = centralHists.at(ratFuncNo)->GetBinCenter(bin);
				double y = ratiofunction(wp, wm, z, ratFuncNo);
				double upperErr = upperUncertaintyHist->GetBinContent(bin);
				double lowerErr = lowerUncertaintyHist->GetBinContent(bin);
				double symmErr = symmUncertaintyHist->GetBinContent(bin);
				
				double refy =1.; //reference pdf value for ratio plot
								 


				if (ratFuncNo==0)
				{
					ApmUnc->SetPoint(bin-1, x, y);
					ApmUnc->SetPointError(bin-1, 0, 0,lowerErr, upperErr);
				}

				else if (ratFuncNo==1)
				{
					RpmUnc->SetPoint(bin-1, x, y);
					RpmUnc->SetPointError(bin-1, 0, 0,lowerErr, upperErr);
				}

				else if (ratFuncNo==2)
				{
					RZWUnc->SetPoint(bin-1, x, y);
					RZWUnc->SetPointError(bin-1, 0, 0,lowerErr, upperErr);
				}

			}


			std::vector<int> colours = {3, 2, 4};
			std::vector<int> markercolours = {3, 2, 4};
			std::vector<int> markers = {73,72,71};
			std::vector<int> fills = {3002,3004,3005};

			std::vector<const char*> PDFnames = {"CT18NNLO", "MSHT20NNLO", "NNPDF40NNLO"}; 
			//std::vector<const char*> PDFnames = {"MSHT20NLO", "CT18NLO", "NNPDF40NLO"}; 

			//plotting different boson eta distributions
			//include central values and shade uncertainty properly
			if (j==0)
			{
				if (ratFuncNo==0)
				{
					c1->cd();
					c1->SetGridy();
					//upperPad1->Draw();
					//upperPad1->cd();
					
					ApmUnc->SetTitle("");
					ApmUnc->SetMarkerStyle(markers[ratFuncNo]);
					ApmUnc->SetMarkerColor(markercolours[ratFuncNo]);
					ApmUnc->SetFillStyle(fills[ratFuncNo]);
					ApmUnc->SetLineColor(markercolours[ratFuncNo]);
					ApmUnc->SetMarkerSize(0.8);
					ApmUnc->GetXaxis()->SetTitleSize(0.05);
					ApmUnc->GetYaxis()->SetTitleSize(0.04);
					ApmUnc->GetXaxis()->SetRangeUser(xlow, xhigh);
					if (j==0)
					{
						//random text stuff to go into plot
						TLatex texts;
						texts.SetTextSize(0.035);
						ApmUnc->Draw("ASP");

						texts.DrawLatex(65, 0.17, "p-p @ #sqrt{s} = 13.6 TeV");
						texts.DrawLatex(65, 0.15, "POWHEG+PYTHIA8");
						texts.DrawLatex(65, 0.13, "CT18NNLO");
						texts.DrawLatex(65, 0.11, "-2.5<y<-4.0");
					}

					else ApmUnc->Draw("SP SAME");
//					leg->AddEntry(WmUnc, "#mu #leftarrow W^{-}", "lep");
//					if (j==(PDFIDs.size()-1)) leg->Draw("SAME");
					
				}
				else if (ratFuncNo==1)
				{
					c2->cd();
					c2->SetGridy();
					//upperPad2->Draw();
					//upperPad2->cd();

					RpmUnc->SetTitle("");
					RpmUnc->SetMarkerStyle(markers[ratFuncNo]);
					RpmUnc->SetMarkerColor(markercolours[ratFuncNo]);
					RpmUnc->SetLineColor(markercolours[ratFuncNo]);
					RpmUnc->SetMarkerSize(0.8);
					RpmUnc->GetXaxis()->SetRangeUser(xlow, xhigh);
					if (j==0)
					{	RpmUnc->Draw("ASP");
						TLatex texts;
						texts.SetTextSize(0.035);

						texts.DrawLatex(10, 1.15, "p-p @ #sqrt{s} = 13.6 TeV");
						texts.DrawLatex(10, 1.1, "POWHEG+PYTHIA8");
						texts.DrawLatex(10, 0.95, "CT18NNLO");
						texts.DrawLatex(10, 0.9, "-2.5<y<-4.0");
					}
					//else RpmUnc->Draw("SP SAME");
					//leg->AddEntry(WpUnc, "#mu #leftarrow W^{+}", "lep");
					//if (j==(PDFIDs.size()-1)) leg2->Draw("SAME");
				}
				else if (ratFuncNo==2)
				{
					c3->cd();
					c3->SetGridy();
					//upperPad3->Draw();
					//upperPad3->cd();

				//	ZUnc->SetTitle("Z");
					RZWUnc->SetTitle("");
					RZWUnc->SetMarkerStyle(markers[ratFuncNo]);
					RZWUnc->SetMarkerColor(markercolours[ratFuncNo]);
					RZWUnc->SetLineColor(markercolours[ratFuncNo]);
					RZWUnc->SetMarkerSize(0.8);
					RZWUnc->GetXaxis()->SetRangeUser(xlow, xhigh);
					if (j==0) 
					{	
						RZWUnc->Draw("ASP");
						TLatex texts;
						texts.SetTextSize(0.035);

						texts.DrawLatex(10, .25, "p-p @ #sqrt{s} = 13.6 TeV");
						texts.DrawLatex(10, .23, "POWHEG+PYTHIA8");
						texts.DrawLatex(10, .21, "CT18NNLO");
						texts.DrawLatex(10, .19, "-4.0<y<-2.5");
					}
					//else RZWUnc->Draw("SP SAME");
//					leg->AddEntry(RZWUnc, "#mu #leftarrow Z", "lep");
//					leg->Draw("SAME");

				}
			}
			
			//Plotting uncertainty percentage for different bosons
			
			if (ratFuncNo==0)
			{
				c4->cd();
				//lowerPad1->Draw();
				//lowerPad1->cd();

				TGraphMultiErrors* Apmratio = new TGraphMultiErrors("Apmratio", "A_{#pm}", nbins, x_values,
						y_ratio_values[j][ratFuncNo], binwidth_vector, binwidth_vector, y_ratio_lower_scale_uncertainty_values[j][ratFuncNo], y_ratio_upper_scale_uncertainty_values[j][ratFuncNo]);
				Apmratio->AddYError(nbins, y_ratio_lower_uncertainty_values[j][ratFuncNo], y_ratio_upper_uncertainty_values[j][ratFuncNo]);

				Apmratio->SetMarkerStyle(markers[j]);
				Apmratio->GetAttLine(0)->SetLineColor(colours[j]);
				Apmratio->GetAttFill(0)->SetFillStyle(fills[j]);
				Apmratio->GetAttFill(0)->SetFillColor(colours[j]);
				Apmratio->SetMarkerColor(colours[j]);
				Apmratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
				Apmratio->GetXaxis()->SetTitleSize(0.05);
				Apmratio->GetYaxis()->SetTitleSize(0.04);
				Apmratio->SetLineColor(colours[j]);
				Apmratio->SetFillStyle(fills[j]);
				Apmratio->SetFillColor(colours[j]);

				Apmratio->GetAttLine(1)->SetLineColor(colours[j]);
				Apmratio->GetAttFill(1)->SetFillStyle(fills[j]);
				Apmratio->GetAttFill(1)->SetFillColor(colours[j]);
				Apmratio->SetMaximum(1.15);
				Apmratio->SetMinimum(.85);
				Apmratio->SetMarkerSize(0.8);
				Apmratio->GetXaxis()->SetRangeUser(xlow, xhigh);
				Apmratio->GetYaxis()->SetTitle("#frac{dA^{PDF_{i}}_{#pm}}{dp_{T}}/#frac{dA_{#pm}^{CT18}}{dp_{T}}");

				if (j==0) Apmratio->Draw("APSE s=0.0 ; ; 5");
				else Apmratio->Draw("PSE SAME s=0.0 ; ; 5");

				leg4->AddEntry(Apmratio, PDFnames[j], "PEFL");
				if (j==(PDFIDs.size()-1)) leg4->Draw();
				
			}

			else if (ratFuncNo==1)

			{
				c5->cd();
				//lowerPad2->Draw();
				//lowerPad2->cd();

				TGraphMultiErrors* Rpmratio = new TGraphMultiErrors("Rpmratio", "R_{#pm}", nbins, x_values,
						y_ratio_values[j][ratFuncNo], binwidth_vector, binwidth_vector, y_ratio_lower_scale_uncertainty_values[j][ratFuncNo], y_ratio_upper_scale_uncertainty_values[j][ratFuncNo]);



				Rpmratio->AddYError(nbins, y_ratio_lower_uncertainty_values[j][ratFuncNo], y_ratio_upper_uncertainty_values[j][ratFuncNo]);

				Rpmratio->SetMarkerStyle(markers[j]);
				Rpmratio->GetAttLine(0)->SetLineColor(colours[j]);
				Rpmratio->GetAttFill(0)->SetFillStyle(fills[j]);
				Rpmratio->GetAttFill(0)->SetFillColor(colours[j]);
				Rpmratio->SetMarkerColor(colours[j]);
				Rpmratio->SetLineColor(colours[j]);
				Rpmratio->SetFillStyle(fills[j]);
				Rpmratio->SetFillColor(colours[j]);

				Rpmratio->GetAttLine(1)->SetLineColor(colours[j]);
				Rpmratio->GetAttFill(1)->SetFillStyle(fills[j]);
				Rpmratio->GetAttFill(1)->SetFillColor(colours[j]);
				Rpmratio->SetMaximum(1.15);
				Rpmratio->SetMinimum(.85);
				Rpmratio->SetMarkerSize(0.8);
				Rpmratio->GetXaxis()->SetRangeUser(xlow, xhigh);
				Rpmratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
				Rpmratio->GetXaxis()->SetTitleSize(0.05);
				Rpmratio->GetYaxis()->SetTitle("#frac{dR_{#pm}^{PDF_{i}}}{dp_{T}}/#frac{dR_{#pm}^{CT18}}{dp_{T}}");
				Rpmratio->GetYaxis()->SetTitleSize(0.04);

				if (j==0) Rpmratio->Draw("APSE s=0.0 ; ; 5");
				else Rpmratio->Draw("PSE SAME s=0.0 ; ; 5");

				leg5->AddEntry(Rpmratio, PDFnames[j], "PEFL");
				if (j==(PDFIDs.size()-1)) leg5->Draw();

			}

			else if (ratFuncNo==2)

			{
				c6->cd();
				//lowerPad3->Draw();
				//lowerPad3->cd();

				TGraphMultiErrors* RWZratio = new TGraphMultiErrors("RWZratio", "R_{Z/W}", nbins, x_values,
						y_ratio_values[j][ratFuncNo], binwidth_vector, binwidth_vector, y_ratio_lower_scale_uncertainty_values[j][ratFuncNo],  y_ratio_upper_scale_uncertainty_values[j][ratFuncNo]);

				RWZratio->AddYError(nbins, y_ratio_lower_uncertainty_values[j][ratFuncNo], y_ratio_upper_uncertainty_values[j][ratFuncNo]);

				RWZratio->SetMarkerStyle(markers[j]);
				RWZratio->GetAttLine(0)->SetLineColor(colours[j]);
				RWZratio->GetAttFill(0)->SetFillStyle(fills[j]);
				RWZratio->GetAttFill(0)->SetFillColor(colours[j]);
				RWZratio->SetMarkerColor(colours[j]);
				RWZratio->SetLineColor(colours[j]);
				RWZratio->SetFillStyle(fills[j]);
				RWZratio->SetFillColor(colours[j]);

				RWZratio->GetAttLine(1)->SetLineColor(colours[j]);
				RWZratio->GetAttFill(1)->SetFillStyle(fills[j]);
				RWZratio->GetAttFill(1)->SetFillColor(colours[j]);
				RWZratio->SetMaximum(1.15);
				RWZratio->SetMinimum(.85);
				RWZratio->SetMarkerSize(0.8);
				RWZratio->GetXaxis()->SetRangeUser(xlow, xhigh);
				RWZratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
				RWZratio->GetYaxis()->SetTitle("#frac{dR_{Z/W}^{PDF_{i}}}{dp_{T}}/#frac{dR_{Z/W}^{CT18}}{dp_{T}}");
				RWZratio->GetXaxis()->SetTitleSize(0.05);
				RWZratio->GetYaxis()->SetTitleSize(0.04);

				if (j==0) RWZratio->Draw("APSE s=0.0 ; ; 5");
				else RWZratio->Draw("PSE SAME s=0.0 ; ; 5");

				leg6->AddEntry(RWZratio, PDFnames[j], "PEFL");
				if (j==(PDFIDs.size()-1)) leg6->Draw();

			}

		}
		j += 1;

	}


	TFile* outFile =new TFile("NLORatiosMuonPTDistComparisionWithScalesForward.root", "RECREATE");
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

