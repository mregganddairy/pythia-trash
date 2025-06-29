#include <LHAPDF/LHAPDF.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <iostream>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <cmath>

int main() {
	const std::string pdfSet = "MSHT20nnlo_as118";
	std::vector<LHAPDF::PDF*> pdfs = LHAPDF::mkPDFs(pdfSet); //0th member of the set
	LHAPDF::PDF* pdf = LHAPDF::mkPDF(pdfSet,0);
	
	const double Q =80.4;//units of GeV
	double Q2 =std::pow(Q, 2.0);
	const int nPoints = 500;
	const double xMin = 1e-4, xMax = 1.0;

	//arrays for x and PDF values
	double x[nPoints], gluon[nPoints], up[nPoints], down[nPoints], strange[nPoints], charm[nPoints], bottom[nPoints], aup[nPoints], adown[nPoints], astrange[nPoints], acharm[nPoints], abottom[nPoints];

	//upper bound of parameter variations
	double  Pgluon[nPoints], Pup[nPoints], Pdown[nPoints], Pstrange[nPoints], Pcharm[nPoints], Pbottom[nPoints], Paup[nPoints], Padown[nPoints], Pastrange[nPoints], Pacharm[nPoints], Pabottom[nPoints];

	//upper bound of uncertainty
	double  Pdeltag[nPoints], Pdeltau[nPoints], Pdeltad[nPoints], Pdeltas[nPoints], Pdeltac[nPoints], Pdeltab[nPoints], Pdeltaau[nPoints], Pdeltaad[nPoints], Pdeltaas[nPoints], Pdeltaac[nPoints], Pdeltaab[nPoints];

	//lower bound of parameter  variations
	double  Mgluon[nPoints], Mup[nPoints], Mdown[nPoints], Mstrange[nPoints], Mcharm[nPoints], Mbottom[nPoints], Maup[nPoints], Madown[nPoints], Mastrange[nPoints], Macharm[nPoints], Mabottom[nPoints];

	//lower bound of uncertainty
	double  Mdeltag[nPoints], Mdeltau[nPoints], Mdeltad[nPoints], Mdeltas[nPoints], Mdeltac[nPoints], Mdeltab[nPoints], Mdeltaau[nPoints], Mdeltaad[nPoints], Mdeltaas[nPoints], Mdeltaac[nPoints], Mdeltaab[nPoints];

	for (int i = 0; i < nPoints; ++i)
	{

		x[i] = xMin * std::pow(xMax / xMin, double(i) / (nPoints - 1));
		gluon[i] = pdf->xfxQ2(21, x[i], Q2); //Gluon (ID 21)
		down[i] = pdf->xfxQ2(1, x[i], Q2); //Down (ID 1)
		up[i] = pdf->xfxQ2(2, x[i], Q2); //Up (ID 2)
		strange[i] = pdf->xfxQ2(3, x[i], Q2); //Strange (ID 3)
		charm[i] = pdf->xfxQ2(4, x[i], Q2); //Charm (ID 4)
		bottom[i] = pdf->xfxQ2(5, x[i], Q2); //Bottom (ID 5)

		adown[i] = pdf->xfxQ2(-1, x[i], Q2); //anti-Down (ID 1)
		aup[i] = pdf->xfxQ2(-2, x[i], Q2); //anti-Up (ID 2)
		astrange[i] = pdf->xfxQ2(-3, x[i], Q2); //anti-Strange (ID 3)
		acharm[i] = pdf->xfxQ2(-4, x[i], Q2); //anti-Charm (ID 4)
		abottom[i] = pdf->xfxQ2(-5, x[i], Q2); //anti-Bottom (ID 5)

		//Setting initial uncertainty to zero because c++ is fucking stupid and sometimes sets values to NaN
		Mdeltag[i] = 0.0;
		Pdeltag[i] = 0.0;

		Mdeltau[i] = 0.0;
		Pdeltau[i] = 0.0;

		Mdeltad[i] = 0.0;
		Pdeltad[i] = 0.0;

		Mdeltas[i] = 0.0;
		Pdeltas[i] = 0.0;

		Mdeltac[i] = 0.0;
		Pdeltac[i] = 0.0;

		Mdeltab[i] = 0.0;
		Pdeltab[i] = 0.0;

		Mdeltaau[i] = 0.0;
		Pdeltaau[i] = 0.0;

		Mdeltaad[i] = 0.0;
		Pdeltaad[i] = 0.0;

		Mdeltaas[i] = 0.0;
		Pdeltaas[i] = 0.0;

		Mdeltaac[i] = 0.0;
		Pdeltaac[i] = 0.0;

		Mdeltaab[i] = 0.0;
		Pdeltaab[i] = 0.0;

		//uncertainty loop
		for (size_t j = 1; j < pdfs.size(); j+=2)
		{

			Pgluon[i] = pdfs[j]->xfxQ2(21, x[i], Q2); //Gluon (ID 21)
			Pdown[i] = pdfs[j]->xfxQ2(1, x[i], Q2); //Down (ID 1)
			Pup[i] = pdfs[j]->xfxQ2(2, x[i], Q2); //Up (ID 2)
			Pstrange[i] = pdfs[j]->xfxQ2(3, x[i], Q2); //Strange (ID 3)
			Pcharm[i] = pdfs[j]->xfxQ2(4, x[i], Q2); //Charm (ID 4)
			Pbottom[i] = pdfs[j]->xfxQ2(5, x[i], Q2); //Bottom (ID 5)

			Padown[i] = pdfs[j]->xfxQ2(-1, x[i], Q2); //anti-Down (ID 1)
			Paup[i] = pdfs[j]->xfxQ2(-2, x[i], Q2); //anti-Up (ID 2)
			Pastrange[i] = pdfs[j]->xfxQ2(-3, x[i], Q2); //anti-Strange (ID 3)
			Pacharm[i] = pdfs[j]->xfxQ2(-4, x[i], Q2); //anti-Charm (ID 4)
			Pabottom[i] = pdfs[j]->xfxQ2(-5, x[i], Q2); //anti-Bottom (ID 5)

			Mgluon[i] = pdfs[j+1]->xfxQ2(21, x[i], Q2); //Gluon (ID 21)
			Mdown[i] = pdfs[j+1]->xfxQ2(1, x[i], Q2); //Down (ID 1)
			Mup[i] = pdfs[j+1]->xfxQ2(2, x[i], Q2); //Up (ID 2)
			Mstrange[i] = pdfs[j+1]->xfxQ2(3, x[i], Q2); //Strange (ID 3)
			Mcharm[i] = pdfs[j+1]->xfxQ2(4, x[i], Q2); //Charm (ID 4)
			Mbottom[i] = pdfs[j+1]->xfxQ2(5, x[i], Q2); //Bottom (ID 5)

			Madown[i] = pdfs[j+1]->xfxQ2(-1, x[i], Q2); //anti-Down (ID 1)
			Maup[i] = pdfs[j+1]->xfxQ2(-2, x[i], Q2); //anti-Up (ID 2)
			Mastrange[i] = pdfs[j+1]->xfxQ2(-3, x[i], Q2); //anti-Strange (ID 3)
			Macharm[i] = pdfs[j+1]->xfxQ2(-4, x[i], Q2); //anti-Charm (ID 4)
			Mabottom[i] = pdfs[j+1]->xfxQ2(-5, x[i], Q2); //anti-Bottom (ID 5)

			//calculating asymmetric uncertainty
/*
			std::cout << "Pcharm[" << i << "] = " << Pcharm[i] << std::endl;
			std::cout << "Mcharm[" << i << "] = " << Mcharm[i] << std::endl;
			std::cout << "charm[" << i << "] = " << charm[i] << std::endl;
				std::cout << "test[" << i << "] = " << std::max(gluon[i] - Pgluon[i], gluon[i] - Mgluon[i]) << std::endl;
			*/
			Pdeltag[i] += std::pow(std::max(std::max(Pgluon[i] - gluon[i], Mgluon[i] - gluon[i]), 0.0),2);
			Mdeltag[i] += std::pow(std::max(std::max(gluon[i] - Pgluon[i], gluon[i] - Mgluon[i]), 0.0),2);


			Pdeltau[i] += std::pow(std::max(std::max(Pup[i] - up[i], Mup[i] - up[i]), 0.0),2);
			Mdeltau[i] += std::pow(std::max(std::max(up[i] - Pup[i], up[i] - Mup[i]), 0.0),2);

			Pdeltad[i] += std::pow(std::max(std::max(Pdown[i] - down[i], Mdown[i] - down[i]), 0.0),2);
			Mdeltad[i] += std::pow(std::max(std::max(down[i] - Pdown[i], down[i] - Mdown[i]), 0.0),2);

			Pdeltas[i] += std::pow(std::max(std::max(Pstrange[i] - strange[i], Mstrange[i] - strange[i]), 0.0),2);
			Mdeltas[i] += std::pow(std::max(std::max(strange[i] - Pstrange[i], strange[i] - Mstrange[i]), 0.0),2);

			Pdeltac[i] += std::pow(std::max(std::max(Pcharm[i] - charm[i], Mcharm[i] - charm[i]), 0.0),2);
			Mdeltac[i] += std::pow(std::max(std::max(charm[i] - Pcharm[i], charm[i] - Mcharm[i]), 0.0),2);

			Pdeltab[i] += std::pow(std::max(std::max(Pbottom[i] - bottom[i], Mbottom[i] - bottom[i]), 0.0),2);
			Mdeltab[i] += std::pow(std::max(std::max(bottom[i] - Pbottom[i], bottom[i] - Mbottom[i]), 0.0),2);


			Pdeltaau[i] += std::pow(std::max(std::max(Paup[i] - aup[i], Maup[i] - aup[i]), 0.0),2);
			Mdeltaau[i] += std::pow(std::max(std::max(aup[i] - Paup[i], aup[i] - Maup[i]), 0.0),2);

			Pdeltaad[i] += std::pow(std::max(std::max(Padown[i] - adown[i], Madown[i] - adown[i]), 0.0),2);
			Mdeltaad[i] += std::pow(std::max(std::max(adown[i] - Padown[i], adown[i] - Madown[i]), 0.0),2);

			Pdeltaas[i] += std::pow(std::max(std::max(Pastrange[i] - astrange[i], Mastrange[i] - astrange[i]), 0.0),2);
			Mdeltaas[i] += std::pow(std::max(std::max(astrange[i] - Pastrange[i], astrange[i] - Mastrange[i]), 0.0),2);

			Pdeltaac[i] += std::pow(std::max(std::max(Pacharm[i] - acharm[i], Macharm[i] - acharm[i]), 0.0),2);
			Mdeltaac[i] += std::pow(std::max(std::max(acharm[i] - Pacharm[i], acharm[i] - Macharm[i]), 0.0),2);

			Pdeltaab[i] += std::pow(std::max(std::max(Pabottom[i] - abottom[i], Mabottom[i] - abottom[i]), 0.0),2);
			Mdeltaab[i] += std::pow(std::max(std::max(abottom[i] - Pabottom[i], abottom[i] - Mabottom[i]), 0.0),2);

		}

		Pdeltag[i] = std::sqrt(Pdeltag[i]);
		Mdeltag[i] = std::sqrt(Mdeltag[i]);

		Pdeltau[i] = std::sqrt(Pdeltau[i]);
		Mdeltau[i] = std::sqrt(Mdeltau[i]);
		
		std::cout << "Pdeltag[" << i << "] = " << Pdeltag[i] << std::endl;
		std::cout << "Mdeltag[" << i << "] = " << Mdeltag[i] << std::endl;

		Pdeltad[i] = std::sqrt(Pdeltad[i]);
		Mdeltad[i] = std::sqrt(Mdeltad[i]);

		Pdeltas[i] = std::sqrt(Pdeltas[i]);
		Mdeltas[i] = std::sqrt(Mdeltas[i]);

		Pdeltac[i] = std::sqrt(Pdeltac[i]);
		Mdeltac[i] = std::sqrt(Mdeltac[i]);

		Pdeltab[i] = std::sqrt(Pdeltab[i]);
		Mdeltab[i] = std::sqrt(Mdeltab[i]);

		Pdeltaau[i] = std::sqrt(Pdeltaau[i]);
		Mdeltaau[i] = std::sqrt(Mdeltaau[i]);

		Pdeltaad[i] = std::sqrt(Pdeltaad[i]);
		Mdeltaad[i] = std::sqrt(Mdeltaad[i]);

		Pdeltaas[i] = std::sqrt(Pdeltaas[i]);
		Mdeltaas[i] = std::sqrt(Mdeltaas[i]);

		Pdeltaac[i] = std::sqrt(Pdeltaac[i]);
		Mdeltaac[i] = std::sqrt(Mdeltaac[i]);

		Pdeltaab[i] = std::sqrt(Pdeltaab[i]);
		Mdeltaab[i] = std::sqrt(Mdeltaab[i]);

	}

	// TGraph stuff
	
	TGraph* gPGluon= new TGraph(nPoints, x, Pdeltag);
	TGraph* gMGluon= new TGraph(nPoints, x, Mdeltag);

	TGraph* gPUp= new TGraph(nPoints, x, Pdeltau);
	TGraph* gMUp= new TGraph(nPoints, x, Mdeltau);

	TGraph* gPDown= new TGraph(nPoints, x, Pdeltad);
	TGraph* gMDown= new TGraph(nPoints, x, Mdeltad);

	TGraph* gPStrange= new TGraph(nPoints, x, Pdeltas);
	TGraph* gMStrange= new TGraph(nPoints, x, Mdeltas);

	TGraph* gPCharm= new TGraph(nPoints, x, Pdeltac);
	TGraph* gMCharm= new TGraph(nPoints, x, Mdeltac);

	TGraph* gPBottom= new TGraph(nPoints, x, Pdeltab);
	TGraph* gMBottom= new TGraph(nPoints, x, Mdeltab);

	TGraph* gPAUp= new TGraph(nPoints, x, Pdeltaau);
	TGraph* gMAUp= new TGraph(nPoints, x, Mdeltaau);

	TGraph* gPADown= new TGraph(nPoints, x, Pdeltaad);
	TGraph* gMADown= new TGraph(nPoints, x, Mdeltaad);

	//Setting position of lower and upperbounds on plots
	
	for (int i=0;i<nPoints;i++)
	{
	gPGluon->SetPoint(i,x[i],Pdeltag[i]/5 +gluon[i]/5);
	gMGluon->SetPoint(i,x[i],-Mdeltag[i]/5 +gluon[i]/5);

	gPUp->SetPoint(i,x[i],Pdeltau[i] +up[i]);
	gMUp->SetPoint(i,x[i],-Mdeltau[i] +up[i]);

	gPDown->SetPoint(i,x[i],Pdeltad[i] +down[i]);
	gMDown->SetPoint(i,x[i],-Mdeltad[i] +down[i]);

	gPStrange->SetPoint(i,x[i],Pdeltas[i] +strange[i]);
	gMStrange->SetPoint(i,x[i],-Mdeltas[i] +strange[i]);

	gPCharm->SetPoint(i,x[i],Pdeltac[i] +charm[i]);
	gMCharm->SetPoint(i,x[i],-Mdeltac[i] +charm[i]);

	gPBottom->SetPoint(i,x[i],Pdeltab[i] +bottom[i]);
	gMBottom->SetPoint(i,x[i],-Mdeltab[i] +bottom[i]);

	gPAUp->SetPoint(i,x[i],Pdeltaau[i] +aup[i]);
	gMAUp->SetPoint(i,x[i],-Mdeltaau[i] +aup[i]);

	gPADown->SetPoint(i,x[i],Pdeltaad[i] +adown[i]);
	gMADown->SetPoint(i,x[i],-Mdeltaad[i] +adown[i]);

	}

	//Stuff for creating shading between lower and upper bounds
	TGraph *gluonshade = new TGraph(2*nPoints);
	TGraph *upshade = new TGraph(2*nPoints);
	TGraph *downshade = new TGraph(2*nPoints);
	TGraph *strangeshade = new TGraph(2*nPoints);
	TGraph *charmshade = new TGraph(2*nPoints);
	TGraph *bottomshade = new TGraph(2*nPoints);
	TGraph *aupshade = new TGraph(2*nPoints);
	TGraph *adownshade = new TGraph(2*nPoints);
	for (int i=0;i<nPoints;i++)
	{
	gluonshade->SetPoint(i,x[i], Pdeltag[i]/5+gluon[i]/5);
	gluonshade->SetPoint(nPoints+i,x[nPoints-i-1],-Mdeltag[nPoints-i-1]/5+gluon[nPoints-i-1]/5);

	upshade->SetPoint(i,x[i],Pdeltau[i]+up[i]);
	upshade->SetPoint(nPoints+i,x[nPoints-i-1],-Mdeltau[nPoints-i-1]+up[nPoints-i-1]);

	downshade->SetPoint(i,x[i],Pdeltad[i]+down[i]);
	downshade->SetPoint(nPoints+i,x[nPoints-i-1],-Mdeltad[nPoints-i-1]+down[nPoints-i-1]);

	strangeshade->SetPoint(i,x[i],Pdeltas[i]+strange[i]);
	strangeshade->SetPoint(nPoints+i,x[nPoints-i-1],-Mdeltas[nPoints-i-1]+strange[nPoints-i-1]);

	charmshade->SetPoint(i,x[i],Pdeltac[i]+charm[i]);
	charmshade->SetPoint(nPoints+i,x[nPoints-i-1],-Mdeltac[nPoints-i-1]+charm[nPoints-i-1]);

	bottomshade->SetPoint(i,x[i],Pdeltab[i]+bottom[i]);
	bottomshade->SetPoint(nPoints+i,x[nPoints-i-1],-Mdeltab[nPoints-i-1]+bottom[nPoints-i-1]);

	aupshade->SetPoint(i,x[i],Pdeltaau[i]+aup[i]);
	aupshade->SetPoint(nPoints+i,x[nPoints-i-1],-Mdeltaau[nPoints-i-1]+aup[nPoints-i-1]);

	adownshade->SetPoint(i,x[i],Pdeltaad[i]+adown[i]);
	adownshade->SetPoint(nPoints+i,x[nPoints-i-1],-Mdeltaad[nPoints-i-1]+adown[nPoints-i-1]);
	}

    TCanvas* c = new TCanvas("c", "PDFs", 800, 600);
    c->SetLogx();
    //c->SetLogy();

	gPStrange->Draw("ALP");
	gMStrange->Draw("L SAME");
	gMStrange->SetLineStyle(9);
	gPStrange->SetLineStyle(9);
	gMStrange->SetLineColor(kGreen);
	gPStrange->SetLineColor(kGreen);
	strangeshade->SetFillStyle(3001);
	strangeshade->SetFillColor(kGreen);
	strangeshade->Draw("f");

	gPStrange->SetMinimum(0);

	gPCharm->Draw("L SAME");
	gMCharm->Draw("L SAME");
	gMCharm->SetLineStyle(9);
	gPCharm->SetLineStyle(9);
	gMCharm->SetLineColor(kCyan);
	gPCharm->SetLineColor(kCyan);
	charmshade->SetFillStyle(3001);
	charmshade->SetFillColor(kCyan);
	charmshade->Draw("f");

	gPBottom->Draw("L SAME");
	gMBottom->Draw("L SAME");
	gMBottom->SetLineStyle(9);
	gPBottom->SetLineStyle(9);
	gMBottom->SetLineColor(kOrange);
	gPBottom->SetLineColor(kOrange);
	bottomshade->SetFillStyle(3001);
	bottomshade->SetFillColor(kOrange);
	bottomshade->Draw("f");


	gPGluon->Draw("L SAME");
	gMGluon->Draw("L SAME");
	gMGluon->SetLineStyle(9);
	gPGluon->SetLineStyle(9);
	gMGluon->SetLineColor(kMagenta);
	gPGluon->SetLineColor(kMagenta);
	gluonshade->SetFillStyle(3001);
	gluonshade->SetFillColor(kMagenta);
	gluonshade->Draw("f");

	gPUp->Draw("L SAME");
	gMUp->Draw("L SAME");
	gMUp->SetLineStyle(9);
	gPUp->SetLineStyle(9);
	gPUp->SetLineColor(kBlue);
	gMUp->SetLineColor(kBlue);
	upshade->SetFillStyle(3365);
	upshade->SetFillColor(kBlue);
	upshade->Draw("f");

	gPDown->Draw("L SAME");
	gMDown->Draw("L SAME");
	gMDown->SetLineStyle(9);
	gPDown->SetLineStyle(9);
	gMDown->SetLineColor(kRed);
	gPDown->SetLineColor(kRed);
	downshade->SetFillStyle(3325);
	downshade->SetFillColor(kRed);
	downshade->Draw("f");


	gPAUp->Draw("L SAME");
	gMAUp->Draw("L SAME");
	gMAUp->SetLineStyle(9);
	gPAUp->SetLineStyle(9);
	gPAUp->SetLineColor(kBlue);
	gMAUp->SetLineColor(kBlue);
	aupshade->SetFillStyle(3356);
	aupshade->SetFillColor(kBlue);
	aupshade->Draw("f");


	gPADown->Draw("L SAME");
	gMADown->Draw("L SAME");
	gMADown->SetLineStyle(9);
	gPADown->SetLineStyle(9);
	gMADown->SetLineColor(kRed);
	gPADown->SetLineColor(kRed);
	adownshade->SetFillStyle(3352);
	adownshade->SetFillColor(kRed);
	adownshade->Draw("f");

    gPStrange->GetXaxis()->SetTitle("x");
    gPStrange->GetYaxis()->SetTitle("xf(x, Q^{2})");
    gPStrange->SetTitle((pdfSet+" at Q = M_{W} ").c_str());
    TLegend* legend = new TLegend(0.9, 0.55, 0.7, 0.9);

    legend->AddEntry(gluonshade, "g/5", "f");
    legend->AddEntry(upshade, "u", "f");
    legend->AddEntry(downshade, "d", "f");
    legend->AddEntry(strangeshade, "s", "f");
    legend->AddEntry(charmshade, "c", "f");
    legend->AddEntry(bottomshade, "b", "f");

    legend->AddEntry(aupshade, "#bar{u}", "f");
    legend->AddEntry(adownshade, "#bar{d}", "f");
   /* legend->AddEntry(gAStrangeUnc, "#bar{s}", "f");
    legend->AddEntry(gACharmUnc, "#bar{c}", "f");
    legend->AddEntry(gABottomUnc, "#bar{b}", "f");
	*/
    legend->Draw();

	
	
	/*
	TGraphAsymmErrors* gGluonUnc = new TGraphAsymmErrors(nPoints, x, gluon, nullptr, nullptr, Mdeltag, Pdeltag);

	TGraphAsymmErrors* gUpUnc = new TGraphAsymmErrors(nPoints, x, up, nullptr, nullptr, Mdeltau, Pdeltau);
	TGraphAsymmErrors* gDownUnc = new TGraphAsymmErrors(nPoints, x, down, nullptr, nullptr, Mdeltad, Pdeltad);
	TGraphAsymmErrors* gCharmUnc = new TGraphAsymmErrors(nPoints, x, charm, nullptr, nullptr, Mdeltac, Pdeltac);
	TGraphAsymmErrors* gStrangeUnc = new TGraphAsymmErrors(nPoints, x, strange, nullptr, nullptr, Mdeltas, Pdeltas);
	TGraphAsymmErrors* gBottomUnc = new TGraphAsymmErrors(nPoints, x, bottom, nullptr, nullptr, Mdeltab, Pdeltab);

	TGraphAsymmErrors* gAUpUnc = new TGraphAsymmErrors(nPoints, x, aup, nullptr, nullptr, Mdeltaau, Pdeltaau);
	TGraphAsymmErrors* gADownUnc = new TGraphAsymmErrors(nPoints, x, adown, nullptr, nullptr, Mdeltaad, Pdeltaad);
	TGraphAsymmErrors* gACharmUnc = new TGraphAsymmErrors(nPoints, x, acharm, nullptr, nullptr, Mdeltaac, Pdeltaac);
	TGraphAsymmErrors* gAStrangeUnc = new TGraphAsymmErrors(nPoints, x, astrange, nullptr, nullptr, Mdeltaas, Pdeltaas);
	TGraphAsymmErrors* gABottomUnc = new TGraphAsymmErrors(nPoints, x, abottom, nullptr, nullptr, Mdeltaab, Pdeltaab);

	//stylling
    gGluonUnc->SetFillColorAlpha(kRed, 0.3);
    gUpUnc->SetFillColorAlpha(kBlue, 0.3);
    gDownUnc->SetFillColorAlpha(kGreen, 0.3);
    gStrangeUnc->SetFillColorAlpha(kMagenta, 0.3);
    gCharmUnc->SetFillColorAlpha(kCyan, 0.3);
    gBottomUnc->SetFillColorAlpha(kOrange, 0.3);
	gGluonUnc->SetLineColor(kRed);  
	gUpUnc->SetLineColor(kBlue);
	gDownUnc->SetLineColor(kGreen);
	gStrangeUnc->SetLineColor(kMagenta);
	gCharmUnc->SetLineColor(kCyan);
	gBottomUnc->SetLineColor(kOrange);

	gGluonUnc->SetFillStyle(1004);
	gUpUnc->SetFillStyle(1004);
	gDownUnc->SetFillStyle(1004);
	gStrangeUnc->SetFillStyle(1004);
	gCharmUnc->SetFillStyle(1004);
	gBottomUnc->SetFillStyle(1004);

    gAUpUnc->SetFillColorAlpha(kBlue, 0.3);
    gADownUnc->SetFillColorAlpha(kGreen, 0.3);
    gAStrangeUnc->SetFillColorAlpha(kMagenta, 0.3);
    gACharmUnc->SetFillColorAlpha(kCyan, 0.3);
    gABottomUnc->SetFillColorAlpha(kOrange, 0.3);
	gAUpUnc->SetLineColor(kBlue);
	gADownUnc->SetLineColor(kGreen);
	gAStrangeUnc->SetLineColor(kMagenta);
	gACharmUnc->SetLineColor(kCyan);
	gABottomUnc->SetLineColor(kOrange);

    gAUpUnc->SetLineStyle(2);
    gADownUnc->SetLineStyle(2);
    gAStrangeUnc->SetLineStyle(2);
    gACharmUnc->SetLineStyle(2);
    gABottomUnc->SetLineStyle(2);
	gAUpUnc->SetFillStyle(3005);
	gADownUnc->SetFillStyle(3005);
	gAStrangeUnc->SetFillStyle(3005);
	gACharmUnc->SetFillStyle(3005);
	gABottomUnc->SetFillStyle(3005);

    TCanvas* c = new TCanvas("c", "PDFs", 800, 600);
    c->SetLogx();
    //c->SetLogy();

	gUpUnc->Draw("ALP");
	gGluonUnc->Draw("L SAME");
	gStrangeUnc->Draw("L SAME");
	gDownUnc->Draw("L SAME");
	gCharmUnc->Draw("L SAME");
	gBottomUnc->Draw("L SAME");

	gAUpUnc->Draw("L SAME");
	gAStrangeUnc->Draw("L SAME");
	gADownUnc->Draw("L SAME");
	gACharmUnc->Draw("L SAME");
	gABottomUnc->Draw("L SAME");

    gUpUnc->GetXaxis()->SetTitle("x");
    gUpUnc->GetYaxis()->SetTitle("xf(x, Q^{2})");
    gUpUnc->SetTitle((pdfSet+" at Q = 80.3 GeV").c_str());

    
    // Create a legend
    TLegend* legend = new TLegend(0.15, 0.55, 0.4, 0.85);
    legend->AddEntry(gGluonUnc, "g", "f");
    legend->AddEntry(gUpUnc, "u", "f");
    legend->AddEntry(gDownUnc, "d", "f");
    legend->AddEntry(gStrangeUnc, "s", "f");
    legend->AddEntry(gCharmUnc, "c", "f");
    legend->AddEntry(gBottomUnc, "b", "f");

    legend->AddEntry(gAUpUnc, "#bar{u}", "f");
    legend->AddEntry(gADownUnc, "#bar{d}", "f");
    legend->AddEntry(gAStrangeUnc, "#bar{s}", "f");
    legend->AddEntry(gACharmUnc, "#bar{c}", "f");
    legend->AddEntry(gABottomUnc, "#bar{b}", "f");
    legend->AddEntry(gAup, "#bar{u}", "l");
    legend->AddEntry(gAdown, "#bar{d}", "l");
    legend->AddEntry(gAstrange, "#bar{s}", "l");
    legend->AddEntry(gAcharm, "#bar{c}", "l");
    legend->AddEntry(gAbottom, "#bar{b}", "l");
    legend->Draw();

	*/
    // Save the plot
    c->SaveAs((pdfSet + "PDFs.pdf").c_str());

    return 0;
}
