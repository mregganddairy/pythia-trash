#include <LHAPDF/LHAPDF.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <iostream>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <cmath>


//calculates the central value and uncertainty of a pdf set
int PDFxf(const std::string pdfSetName, int nPoints, double x[500], double PDFValues[11][500], double PDFUpperuncertainty[11][500], double PDFLoweruncertainty[11][500], bool NNPDF) 
{
	
	//pdf sets
	std::vector<LHAPDF::PDF*> pdfs = LHAPDF::mkPDFs(pdfSetName); //0th member of the set
	//central values
	LHAPDF::PDF* pdf = LHAPDF::mkPDF(pdfSetName,0);
	
	
	//const double Q =80.4;//units of GeV
	//const double Q =std::sqrt(10.);//units of GeV
	const double Q =std::sqrt(80.4);//units of GeV
	double Q2 =std::pow(Q, 2.0);
	const double xMin = 1e-6, xMax = 0.9;

	int IDs[11] = {21, 1, 2, 3, 4, 5, -1, -2, -3, -4, -5};
	
	//upper bound of parameter variations
	double  PDFUppervariations[11][nPoints];

	//lower bound of parameter  variations
	double  PDFLowervariations[11][nPoints];

	for (int i = 0; i < nPoints; ++i)
	{
		x[i] = xMin * std::pow(xMax / xMin, double(i) / (nPoints - 1));

		for ( int k = 0; k < 11; k++)
		{
			PDFValues[k][i]=pdf->xfxQ2(IDs[k],x[i],Q2);

			//Setting initial uncertainty to zero because c++ is fucking stupid and sometimes sets values to NaN
			PDFLoweruncertainty[k][i]=0.0;
			PDFUpperuncertainty[k][i]=0.0;

			//uncertainty loop
			if (NNPDF==true)
			{
				double pdfsize =pdfs.size();
				/*
				for (size_t j = 1; j < pdfs.size(); j+=2)
				{
					PDFLowervariations[k][i]=pdfs[j+1]->xfxQ2(IDs[k], x[i], Q2);
					PDFUppervariations[k][i]=pdfs[j]->xfxQ2(IDs[k], x[i], Q2);

					//calculating asymmetric uncertainty
					PDFUpperuncertainty[k][i] += std::pow(std::max(std::max(PDFUppervariations[k][i] - PDFValues[k][i], PDFLowervariations[k][i] - PDFValues[k][i]), 0.0),2);
					PDFLoweruncertainty[k][i] += std::pow(std::max(std::max(PDFValues[k][i] - PDFUppervariations[k][i], PDFValues[k][i] - PDFLowervariations[k][i]), 0.0),2);
				}
				*/

				//member set loop
				for (size_t j = 1; j < pdfs.size(); j++)
				{
					PDFUppervariations[k][i]=pdfs[j]->xfxQ2(IDs[k], x[i], Q2);
					PDFLowervariations[k][i]=pdfs[j]->xfxQ2(IDs[k], x[i], Q2);

					//calculating asymmetric uncertainty
					PDFUpperuncertainty[k][i] += std::pow(PDFUppervariations[k][i] - PDFValues[k][i], 2);
					PDFLoweruncertainty[k][i] += std::pow(PDFUppervariations[k][i] - PDFValues[k][i], 2);
				}

				//PDFUpperuncertainty[k][i] = std::sqrt((1/(pdfs.size()-2))*PDFUpperuncertainty[k][i]);
				//PDFLoweruncertainty[k][i] = std::sqrt((1/(pdfs.size()-2))*PDFLoweruncertainty[k][i]);
				PDFUpperuncertainty[k][i] = (std::sqrt(1/(pdfsize-2)))*(std::sqrt(PDFUpperuncertainty[k][i]));
				PDFLoweruncertainty[k][i] = (std::sqrt(1/(pdfsize-2)))*(std::sqrt(PDFLoweruncertainty[k][i]));
			}
			else
			{
				for (size_t j = 1; j < pdfs.size(); j+=2)
				{
					PDFLowervariations[k][i]=pdfs[j+1]->xfxQ2(IDs[k], x[i], Q2);
					PDFUppervariations[k][i]=pdfs[j]->xfxQ2(IDs[k], x[i], Q2);

					//calculating asymmetric uncertainty
					PDFUpperuncertainty[k][i] += std::pow(std::max(std::max(PDFUppervariations[k][i] - PDFValues[k][i], PDFLowervariations[k][i] - PDFValues[k][i]), 0.0),2);
					PDFLoweruncertainty[k][i] += std::pow(std::max(std::max(PDFValues[k][i] - PDFUppervariations[k][i], PDFValues[k][i] - PDFLowervariations[k][i]), 0.0),2);
				}

				PDFUpperuncertainty[k][i] = std::sqrt(PDFUpperuncertainty[k][i]);
				PDFLoweruncertainty[k][i] = std::sqrt(PDFLoweruncertainty[k][i]);
			}
		}
	}
	
    return 0;
}

int main()
{
	//TGraph stuff
	const int nPoints = 500;
	double PDFUpperuncertainty[3][11][nPoints];
	double PDFLoweruncertainty[3][11][nPoints];
	double PDFValues[3][11][nPoints];
	const std::vector<std::string> pdfSet = {"MSHT20nnlo_as118","CT18NNLO","NNPDF40_nnlo_as_01180"};
	std::vector<const char*> charpdfSet = {"MSHT20nnlo_as118","CT18NNLO","NNPDF40_nnlo_as_01180"};
	std::vector<const char*> charpdfTitles = {"MSHT20NNLO","CT18NNLO","NNPDF40NNLO"};
	std::vector<const char*> partons = {"gluon", "down", "up", "strange", "charm", "bottom", "anti-down", "anti-up", "anti-strange", "anti-charm", "anti-bottom"};
	std::vector<std::string> partonstring = {"gluon", "down", "up", "strange", "charm", "bottom", "anti-down", "anti-up", "anti-strange", "anti-charm", "anti-bottom"};

	int noOfPartons = partons.size();
	int noOfPDFs = pdfSet.size();

	double x[nPoints];

	std::vector<std::vector<TGraph*>> Pgraphs(noOfPDFs, std::vector<TGraph*>(noOfPartons)); 
	std::vector<std::vector<TGraph*>> Mgraphs(noOfPDFs, std::vector<TGraph*>(noOfPartons)); 
	std::vector<std::vector<TGraph*>> Shade(noOfPDFs, std::vector<TGraph*>(noOfPartons)); 
	
	std::vector<std::string> colours = {"kGreen", "kCyan", "kOrange", "kMagenta", "kBlue", "kRed"};

	//Making Canvases
	std::vector<TCanvas*> c(3);
	std::vector<TCanvas*> cparton(noOfPartons);


	
	for (int j=0; j<3; j++)
	{
		c[j] = new TCanvas(charpdfSet[j], charpdfSet[j], 800, 600);
		c[j]->cd();

		if (j==2) PDFxf(pdfSet[j], nPoints, x, PDFValues[j], PDFUpperuncertainty[j], PDFLoweruncertainty[j], true);
		else PDFxf(pdfSet[j], nPoints, x, PDFValues[j], PDFUpperuncertainty[j], PDFLoweruncertainty[j], false);

		for (int k=0; k<noOfPartons; k++)
		{
			Pgraphs.at(j).at(k) = new TGraph(nPoints, x, PDFUpperuncertainty[j][k]);
			Mgraphs.at(j).at(k) = new TGraph(nPoints, x, PDFLoweruncertainty[j][k]);
			Shade.at(j).at(k) = new TGraph(2*nPoints); 

			//Setting position of lower and upperbounds on plots
			for (int i=0;i<nPoints;i++)
			{
				if (j==0)
				{
				Pgraphs[j][k]->SetPoint(i, x[i], PDFUpperuncertainty[j][k][i]/1.645 +PDFValues[j][k][i]);
				Mgraphs[j][k]->SetPoint(i, x[i], -PDFLoweruncertainty[j][k][i]/1.645 +PDFValues[j][k][i]);
				}
				else
				{
				Pgraphs[j][k]->SetPoint(i, x[i], PDFUpperuncertainty[j][k][i] +PDFValues[j][k][i]);
				Mgraphs[j][k]->SetPoint(i, x[i], -PDFLoweruncertainty[j][k][i] +PDFValues[j][k][i]);
				}

				//Stuff for creating shading between lower and upper bounds
				Shade.at(j).at(k)->SetPoint(i,x[i], PDFUpperuncertainty[j][k][i]+PDFValues[j][k][i]);
				Shade.at(j).at(k)->SetPoint(nPoints+i,x[nPoints-i-1], -PDFLoweruncertainty[j][k][nPoints-i-1]+PDFValues[j][k][nPoints-i-1]);
			}
		}

	//Plotting stuff
    c[j]->SetLogx();
    //c->SetLogy();

	Pgraphs[j][3]->Draw("ALP");
	Mgraphs[j][3]->Draw("L SAME");
	Shade[j][3]->Draw("f");

	for (int k=0; k<8; k++)
	{
		if (k==3) continue;
		Pgraphs[j][k]->Draw("L SAME");
		Mgraphs[j][k]->Draw("L SAME");
		//Pgraphs[j][k]->SetLineStyle(9);
		//Mgraphs[j][k]->SetLineStyle(9);
		Shade[j][k]->Draw("f");
	}


	//int IDs[11] = {21, 1, 2, 3, 4, 5, -1, -2, -3, -4, -5};
	Pgraphs[j][0]->SetLineColor(kMagenta);
	Mgraphs[j][0]->SetLineColor(kMagenta);
	Shade[j][0]->SetFillColor(kMagenta);
	Shade[j][0]->SetFillStyle(3001);
	Pgraphs[j][1]->SetLineColor(kBlue);
	Mgraphs[j][1]->SetLineColor(kBlue);
	Shade[j][1]->SetFillColor(kBlue);
	Shade[j][1]->SetFillStyle(3365);
	Pgraphs[j][2]->SetLineColor(kRed);
	Mgraphs[j][2]->SetLineColor(kRed);
	Shade[j][2]->SetFillColor(kRed);
	Shade[j][2]->SetFillStyle(3325);
	Pgraphs[j][3]->SetLineColor(kGreen);
	Mgraphs[j][3]->SetLineColor(kGreen);
	Shade[j][3]->SetFillColor(kGreen);
	Shade[j][3]->SetFillStyle(3001);
	Pgraphs[j][4]->SetLineColor(kCyan);
	Mgraphs[j][4]->SetLineColor(kCyan);
	Shade[j][4]->SetFillColor(kCyan);
	Shade[j][4]->SetFillStyle(3001);
	Pgraphs[j][5]->SetLineColor(kOrange);
	Mgraphs[j][5]->SetLineColor(kOrange);
	Shade[j][5]->SetFillColor(kOrange);
	Shade[j][5]->SetFillStyle(3001);
	Pgraphs[j][6]->SetLineColor(kRed);
	Mgraphs[j][6]->SetLineColor(kRed);
	Shade[j][6]->SetFillColor(kRed);
	Shade[j][6]->SetFillStyle(3356);
	Pgraphs[j][7]->SetLineColor(kBlue);
	Mgraphs[j][7]->SetLineColor(kBlue);
	Shade[j][7]->SetFillColor(kBlue);
	Shade[j][7]->SetFillStyle(3352);

	Pgraphs[j][3]->GetXaxis()->SetTitle("x");
    Pgraphs[j][3]->GetYaxis()->SetTitle("xf(x, Q^{2})");
    Pgraphs[j][3]->SetTitle((pdfSet[j]+" at Q = M_{W} ").c_str());
    TLegend* legend = new TLegend(0.9, 0.55, 0.7, 0.9);

    legend->AddEntry(Shade[j][0], "g", "f");
    legend->AddEntry(Shade[j][2], "u", "f");
    legend->AddEntry(Shade[j][1], "d", "f");
    legend->AddEntry(Shade[j][3], "s", "f");
    legend->AddEntry(Shade[j][4], "c", "f");
    legend->AddEntry(Shade[j][5], "b", "f");

    legend->AddEntry(Shade[j][6], "#bar{u}", "f");
    legend->AddEntry(Shade[j][7], "#bar{d}", "f");
	legend->Draw();
    c[j]->SaveAs((pdfSet[j] + "PDFs.pdf").c_str());
	}
	
	//plotting ratios for uncertainty comparisons
	std::vector<std::vector<TGraph*>> centralratios(3, std::vector<TGraph*>(11));
	std::vector<std::vector<TGraph*>> Pratios(3, std::vector<TGraph*>(11));
	std::vector<std::vector<TGraph*>> Mratios(3, std::vector<TGraph*>(11));
	std::vector<std::vector<TGraph*>> shaderatios(3, std::vector<TGraph*>(11));

	std::vector<TLegend*> partonlegend(11);
	

	TFile* outFile = new TFile("PDFRatios.root","RECREATE");
	for (int k = 0; k < 11; ++k)
	{
		cparton.at(k) = new TCanvas(partons[k], partons[k], 800, 600);
		partonlegend[k]= new TLegend(0.75, 0.7, 0.35, 0.9);
		for (int j= 0; j<3; ++j)
		{
			centralratios.at(j).at(k) = new TGraph(nPoints, x, PDFValues[j][k]);
			Pratios.at(j).at(k) = new TGraph(nPoints, x, PDFValues[j][k]);
			Mratios.at(j).at(k) = new TGraph(nPoints, x, PDFValues[j][k]);
			shaderatios.at(j).at(k) = new TGraph(2*nPoints);

			centralratios[j][k]->SetLineColor(j+2);
			Pratios[j][k]->SetLineColor(j+2);
			Pratios[j][k]->SetLineStyle(9);
			Mratios[j][k]->SetLineColor(j+2);
			Mratios[j][k]->SetLineStyle(9);
			shaderatios[j][k]->SetFillColor(j+2);

			partonlegend[k]->AddEntry(shaderatios[j][k], charpdfTitles[j] , "f");

			for (int i =0; i <nPoints; ++i)
			{
				centralratios.at(j).at(k)->SetPoint(i, x[i], PDFValues[j][k][i]/PDFValues[2][k][i]);
				Pratios.at(j).at(k)->SetPoint(i, x[i], (PDFUpperuncertainty[j][k][i] +PDFValues[j][k][i])/PDFValues[2][k][i]);
				Mratios.at(j).at(k)->SetPoint(i, x[i], (-PDFLoweruncertainty[j][k][i] +PDFValues[j][k][i])/PDFValues[2][k][i]);

				shaderatios.at(j).at(k)->SetPoint(i, x[i], (PDFUpperuncertainty[j][k][i]+PDFValues[j][k][i])/PDFValues[2][k][i]);
				shaderatios.at(j).at(k)->SetPoint(nPoints+i, x[nPoints-i-1], (-PDFLoweruncertainty[j][k][nPoints-i-1]+PDFValues[j][k][nPoints-i-1])/PDFValues[2][k][nPoints-i-1]);
			}

			cparton[k]->cd();
			cparton[k]->SetLogx();

			if (j==0) 
			{
				centralratios[j][k]->Draw("ALP");
				centralratios[j][k]->SetTitle((partonstring[k]+" at Q^{2} = M_{W}^{2} GeV^{2}").c_str());
				centralratios[j][k]->GetXaxis()->SetTitle("x");
				centralratios[j][k]->GetYaxis()->SetTitle("xf(x, Q^{2}) ratio to NNPDF4.0NNLO");
				centralratios[j][k]->GetYaxis()->CenterTitle(true);
				centralratios[j][k]->GetYaxis()->SetTitleSize(0.045);
				centralratios[j][k]->GetXaxis()->SetTitleSize(0.05);
				shaderatios[j][k]->SetFillStyle(3001);
			}
			else centralratios[j][k]->Draw("L SAME");
			if (j==1) shaderatios[j][k]->SetFillStyle(3001);
			if (j==2) shaderatios[j][k]->SetFillStyle(3021);
			if (j==0) shaderatios[j][k]->SetFillStyle(3022);
			//if (j==2) shaderatios[j][k]->SetFillStyle(3352);
			
			Pratios[j][k]->Draw("L SAME");
			Mratios[j][k]->Draw("L SAME");
			shaderatios[j][k]->Draw("f");

			centralratios[0][k]->SetMaximum(1.5);
			centralratios[0][k]->SetMinimum(0.5);
			partonlegend[k]->Draw();

		}
		cparton.at(k)->SaveAs((partonstring.at(k)+ "Q=M_ratioPDFs.pdf").c_str());
		cparton.at(k)->Write();
	}
	outFile->Close();

    return 0;
}
