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

int assignment8macro2() 
{
    //define histogram binning parameters
    const int nbins = 50;
    const double xlow = -5.0;
    const double xhigh = 5.0;

	std::vector<std::string> PDFIDs= {"274", "140", "3311"};
	const std::vector<int> AllnumVariations = {64, 58, 99}; 
    
			
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
		for (int k = 0; k < 1; ++k)
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



			std::vector<int> colours = {2, 4, 6};
			std::vector<int> markers = {73,72,71};
			std::vector<const char*> PDFnames = {"CT18NNLO", "MSHT20NNLO", "NNPDF40NNLO"}; 

			//plotting different boson eta distributions
			
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

				WmUnc->SetTitle("W^{-}");
				UpperWmUncRat->SetTitle(" ");
				UpperWmUncRat->SetMarkerStyle(73);
				UpperWmUncRat->SetMarkerColor(colours[j]);
				UpperWmUncRat->SetLineColor(colours[j]);
				if (j==0) UpperWmUncRat->Draw("AP");
				else  UpperWmUncRat->Draw("P SAME");
				UpperWmUncRat->GetYaxis()->SetTitle("Percentage Uncertainty (%)");
				UpperWmUncRat->GetXaxis()->SetTitle("y");
				UpperWmUncRat->SetMaximum(10);
				UpperWmUncRat->SetMinimum(-10);
				leg4->AddEntry(UpperWmUncRat, PDFnames[j], "lep");

				LowerWmUncRat->SetMarkerStyle(73);
				LowerWmUncRat->SetMarkerColor(colours[j]);
				LowerWmUncRat->SetLineColor(colours[j]);
				LowerWmUncRat->Draw("P SAME");
				if (j==(PDFIDs.size()-1)) leg4->Draw();
			}
			else if (k==1)
			{
				c5->cd();
				c5->SetGridy();
				WpUnc->SetTitle("W^{+}");
				UpperWpUncRat->SetMarkerStyle(72);
				UpperWpUncRat->SetMarkerColor(colours[j]);
				UpperWpUncRat->SetLineColor(colours[j]);
				UpperWpUncRat->SetMaximum(10);
				UpperWpUncRat->SetMinimum(-10);
				if (j==0) UpperWpUncRat->Draw("AP");
				else  UpperWpUncRat->Draw("P SAME");
				leg5->AddEntry(UpperWpUncRat, PDFnames[j], "lep");

				LowerWpUncRat->SetMarkerStyle(72);
				LowerWpUncRat->SetMarkerColor(colours[j]);
				LowerWpUncRat->SetLineColor(colours[j]);
				LowerWpUncRat->Draw("P SAME");
				if (j==(PDFIDs.size()-1)) leg5->Draw();
			}
			else if (k==2)
			{
				c6->cd();
				c6->SetGridy();
				ZUnc->SetTitle("Z");
				UpperZUncRat->SetMarkerStyle(71);
				UpperZUncRat->SetMarkerColor(colours[j]);
				UpperZUncRat->SetLineColor(colours[j]);
				UpperZUncRat->SetMaximum(10);
				UpperZUncRat->SetMinimum(-10);
				if (j==0) UpperZUncRat->Draw("AP");
				else  UpperZUncRat->Draw("P SAME");
				leg6->AddEntry(UpperZUncRat, PDFnames[j], "lep");

				LowerZUncRat->SetMarkerStyle(71);
				LowerZUncRat->SetMarkerColor(colours[j]);
				LowerZUncRat->SetLineColor(colours[j]);
				LowerZUncRat->Draw("P SAME");

				if (j==(PDFIDs.size()-1)) leg6->Draw();
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

