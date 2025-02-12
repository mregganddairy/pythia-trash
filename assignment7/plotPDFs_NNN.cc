#include <LHAPDF/LHAPDF.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <iostream>
#include <cmath>
#include <string>

int main() {
    // Choose your PDF set
    const std::string pdfSet = "NNPDF40_lo_as_01180";
    std::vector<LHAPDF::PDF*> pdfs = LHAPDF::mkPDFs(pdfSet);
    LHAPDF::PDF* centralPDF = LHAPDF::mkPDF(pdfSet, 0);
    
    const double Q = 10;
    double Q2 = Q * Q;
    const int nPoints = 500;
    const double xMin = 1e-4, xMax = 1.0;
    
    // Arrays for x and central PDF values (for particles and anti-particles)
    double x[nPoints];
    double gluon[nPoints], up[nPoints], down[nPoints], strange[nPoints], charm[nPoints], bottom[nPoints];
    double aup[nPoints], adown[nPoints], astrange[nPoints], acharm[nPoints], abottom[nPoints];
    
    // Arrays for uncertainties (symmetric errors)
    double errGluon[nPoints], errUp[nPoints], errDown[nPoints], errStrange[nPoints], errCharm[nPoints], errBottom[nPoints];
    double errAup[nPoints], errAdown[nPoints], errAstrange[nPoints], errAcharm[nPoints], errAbottom[nPoints];
    
    // Loop over x points
    for (int i = 0; i < nPoints; i++) {
        x[i] = xMin * std::pow(xMax / xMin, double(i) / (nPoints - 1));
        
        // Central values from member 0
        gluon[i]  = centralPDF->xfxQ2(21, x[i], Q2);
        down[i]   = centralPDF->xfxQ2(1, x[i], Q2);
        up[i]     = centralPDF->xfxQ2(2, x[i], Q2);
        strange[i]= centralPDF->xfxQ2(3, x[i], Q2);
        charm[i]  = centralPDF->xfxQ2(4, x[i], Q2);
        bottom[i] = centralPDF->xfxQ2(5, x[i], Q2);
        
        adown[i] = centralPDF->xfxQ2(-1, x[i], Q2);
        aup[i]   = centralPDF->xfxQ2(-2, x[i], Q2);
        astrange[i] = centralPDF->xfxQ2(-3, x[i], Q2);
        acharm[i]   = centralPDF->xfxQ2(-4, x[i], Q2);
        abottom[i]  = centralPDF->xfxQ2(-5, x[i], Q2);
        
        // Initialize sums for RMS calculation
        double sum2_g = 0, sum2_u = 0, sum2_d = 0, sum2_s = 0, sum2_c = 0, sum2_b = 0;
        double sum2_au = 0, sum2_ad = 0, sum2_as = 0, sum2_ac = 0, sum2_ab = 0;
        
        // Number of error members (assumes pdfs[0] is central)
        int N = pdfs.size() - 1;
        for (size_t j = 1; j < pdfs.size(); j++) {
            double val_g  = pdfs[j]->xfxQ2(21, x[i], Q2);
            double val_d  = pdfs[j]->xfxQ2(1, x[i], Q2);
            double val_u  = pdfs[j]->xfxQ2(2, x[i], Q2);
            double val_s  = pdfs[j]->xfxQ2(3, x[i], Q2);
            double val_c  = pdfs[j]->xfxQ2(4, x[i], Q2);
            double val_b  = pdfs[j]->xfxQ2(5, x[i], Q2);
            
            double val_ad = pdfs[j]->xfxQ2(-1, x[i], Q2);
            double val_au = pdfs[j]->xfxQ2(-2, x[i], Q2);
            double val_as = pdfs[j]->xfxQ2(-3, x[i], Q2);
            double val_ac = pdfs[j]->xfxQ2(-4, x[i], Q2);
            double val_ab = pdfs[j]->xfxQ2(-5, x[i], Q2);
            
            sum2_g  += std::pow(val_g - gluon[i], 2);
            sum2_d  += std::pow(val_d - down[i], 2);
            sum2_u  += std::pow(val_u - up[i], 2);
            sum2_s  += std::pow(val_s - strange[i], 2);
            sum2_c  += std::pow(val_c - charm[i], 2);
            sum2_b  += std::pow(val_b - bottom[i], 2);
            
            sum2_ad += std::pow(val_ad - adown[i], 2);
            sum2_au += std::pow(val_au - aup[i], 2);
            sum2_as += std::pow(val_as - astrange[i], 2);
            sum2_ac += std::pow(val_ac - acharm[i], 2);
            sum2_ab += std::pow(val_ab - abottom[i], 2);
        }
        // Compute RMS error (if N > 0)
        errGluon[i]  = (N > 0) ? std::sqrt(sum2_g / N)  : 0;
        errDown[i]   = (N > 0) ? std::sqrt(sum2_d / N)  : 0;
        errUp[i]     = (N > 0) ? std::sqrt(sum2_u / N)  : 0;
        errStrange[i]= (N > 0) ? std::sqrt(sum2_s / N)  : 0;
        errCharm[i]  = (N > 0) ? std::sqrt(sum2_c / N)  : 0;
        errBottom[i] = (N > 0) ? std::sqrt(sum2_b / N)  : 0;
        
        errAdown[i] = (N > 0) ? std::sqrt(sum2_ad / N) : 0;
        errAup[i]   = (N > 0) ? std::sqrt(sum2_au / N) : 0;
        errAstrange[i]= (N > 0) ? std::sqrt(sum2_as / N) : 0;
        errAcharm[i] = (N > 0) ? std::sqrt(sum2_ac / N) : 0;
        errAbottom[i]= (N > 0) ? std::sqrt(sum2_ab / N) : 0;
    }
    
    // Create TGraphErrors for each flavor (particles)
    TGraphErrors* gGluon  = new TGraphErrors(nPoints, x, gluon, nullptr, errGluon);
    TGraphErrors* gUp     = new TGraphErrors(nPoints, x, up, nullptr, errUp);
    TGraphErrors* gDown   = new TGraphErrors(nPoints, x, down, nullptr, errDown);
    TGraphErrors* gStrange= new TGraphErrors(nPoints, x, strange, nullptr, errStrange);
    TGraphErrors* gCharm  = new TGraphErrors(nPoints, x, charm, nullptr, errCharm);
    TGraphErrors* gBottom = new TGraphErrors(nPoints, x, bottom, nullptr, errBottom);
    
    // And for the anti-partons
    TGraphErrors* gAUp    = new TGraphErrors(nPoints, x, aup, nullptr, errAup);
    TGraphErrors* gADown  = new TGraphErrors(nPoints, x, adown, nullptr, errAdown);
    TGraphErrors* gAStrange= new TGraphErrors(nPoints, x, astrange, nullptr, errAstrange);
    TGraphErrors* gACharm = new TGraphErrors(nPoints, x, acharm, nullptr, errAcharm);
    TGraphErrors* gABottom= new TGraphErrors(nPoints, x, abottom, nullptr, errAbottom);
    
    // Set line and fill colors/styles for particles
    gGluon->SetLineColor(kRed);
    gUp->SetLineColor(kBlue);
    gDown->SetLineColor(kGreen);
    gStrange->SetLineColor(kMagenta);
    gCharm->SetLineColor(kCyan);
    gBottom->SetLineColor(kOrange);
    
    gGluon->SetFillColorAlpha(kRed, 0.3);
    gUp->SetFillColorAlpha(kBlue, 0.3);
    gDown->SetFillColorAlpha(kGreen, 0.3);
    gStrange->SetFillColorAlpha(kMagenta, 0.3);
    gCharm->SetFillColorAlpha(kCyan, 0.3);
    gBottom->SetFillColorAlpha(kOrange, 0.3);
    
    // For anti-particles, use the same color but a different (e.g. dashed) line and a hatched fill style
    gAUp->SetLineColor(kBlue);    gAUp->SetLineStyle(2);    gAUp->SetFillStyle(3005);    gAUp->SetFillColorAlpha(kBlue, 0.3);
    gADown->SetLineColor(kGreen);  gADown->SetLineStyle(2);  gADown->SetFillStyle(3005);  gADown->SetFillColorAlpha(kGreen, 0.3);
    gAStrange->SetLineColor(kMagenta); gAStrange->SetLineStyle(2); gAStrange->SetFillStyle(3005); gAStrange->SetFillColorAlpha(kMagenta, 0.3);
    gACharm->SetLineColor(kCyan);   gACharm->SetLineStyle(2);   gACharm->SetFillStyle(3005);   gACharm->SetFillColorAlpha(kCyan, 0.3);
    gABottom->SetLineColor(kOrange); gABottom->SetLineStyle(2); gABottom->SetFillStyle(3005); gABottom->SetFillColorAlpha(kOrange, 0.3);
    
    // Create canvas and draw
    TCanvas* c = new TCanvas("c", "PDFs", 800, 600);
    c->SetLogx();
    
    // Draw central particle curves with errors first
    gUp->Draw("ALP");
    gGluon->Draw("L SAME");
    gStrange->Draw("L SAME");
    gDown->Draw("L SAME");
    gCharm->Draw("L SAME");
    gBottom->Draw("L SAME");
    
    // Then draw the anti-particles
    gAUp->Draw("L SAME");
    gAStrange->Draw("L SAME");
    gADown->Draw("L SAME");
    gACharm->Draw("L SAME");
    gABottom->Draw("L SAME");
    
    // Set axis titles and overall title
    gUp->GetXaxis()->SetTitle("x");
    gUp->GetYaxis()->SetTitle("xf(x, Q^{2})");
    std::string title = pdfSet + " at Q = " + "0.01" + " TeV";
    gUp->SetTitle(title.c_str());
    
    // Create and draw legend
    TLegend* legend = new TLegend(0.15, 0.55, 0.4, 0.85);
    legend->AddEntry(gGluon, "g", "lf");
    legend->AddEntry(gUp, "u", "lf");
    legend->AddEntry(gDown, "d", "lf");
    legend->AddEntry(gStrange, "s", "lf");
    legend->AddEntry(gCharm, "c", "lf");
    legend->AddEntry(gBottom, "b", "lf");
    legend->AddEntry(gAUp, "#bar{u}", "lf");
    legend->AddEntry(gADown, "#bar{d}", "lf");
    legend->AddEntry(gAStrange, "#bar{s}", "lf");
    legend->AddEntry(gACharm, "#bar{c}", "lf");
    legend->AddEntry(gABottom, "#bar{b}", "lf");
    legend->Draw();
    
    // Save plot as PDF
    c->SaveAs((pdfSet + "PDFs.pdf").c_str());
    
    return 0;
}

