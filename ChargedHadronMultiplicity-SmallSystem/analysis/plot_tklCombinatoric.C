#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

// ---------------------------------------------------------------------------
// Build a grid of named TPads that share the canvas margins properly.
// ---------------------------------------------------------------------------
void CanvasPartition(TCanvas *c, int Nx, int Ny, float lMargin, float rMargin, float bMargin, float tMargin)
{
    if (!c)
        return;

    const float hStep = (1.f - lMargin - rMargin) / Nx;
    const float vStep = (1.f - bMargin - tMargin) / Ny;

    for (int i = 0; i < Nx; ++i)
    {
        float hposl = (i == 0) ? 0.f : lMargin + i * hStep;
        float hposr = (i == Nx - 1) ? 1.f : lMargin + (i + 1) * hStep;
        float hmarl = (i == 0) ? lMargin / (hposr - hposl) : 0.f;
        float hmarr = (i == Nx - 1) ? rMargin / (hposr - hposl) : 0.f;

        for (int j = 0; j < Ny; ++j)
        {
            float vposd = (j == 0) ? 0.f : bMargin + j * vStep;
            float vposu = (j == Ny - 1) ? 1.f : bMargin + (j + 1) * vStep;
            float vmard = (j == 0) ? bMargin / (vposu - vposd) : 0.f;
            float vmaru = (j == Ny - 1) ? tMargin / (vposu - vposd) : 0.f;

            c->cd(0);
            TString name = TString::Format("pad_%d_%d", i, j);
            delete c->FindObject(name.Data()); // remove stale pad if any

            auto *pad = new TPad(name, "", hposl, vposd, hposr, vposu);
            pad->SetLeftMargin(hmarl);
            pad->SetRightMargin(hmarr);
            pad->SetBottomMargin(vmard);
            pad->SetTopMargin(vmaru);
            pad->SetBorderMode(0);
            pad->SetBorderSize(0);
            pad->SetFrameBorderMode(0);
            pad->Draw();
        }
    }
}

// ---------------------------------------------------------------------------
// Convert axis-fraction coordinates to pad NDC (accounts for margins).
// ---------------------------------------------------------------------------
static double XtoPad(double x)
{
    double lm = gPad->GetLeftMargin(), rm = gPad->GetRightMargin();
    return lm + x * (1. - lm - rm);
}
static double YtoPad(double y)
{
    double bm = gPad->GetBottomMargin(), tm = gPad->GetTopMargin();
    return bm + y * (1. - bm - tm);
}

template <typename T> static std::string to_string_with_precision(const T a_value, const int n = 2)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

static std::string EtaPhiRangeNameSuffix(const TH2D *h2, int ieta, int iphi)
{
    const auto *etaAx = h2->GetXaxis();
    const auto *phiAx = h2->GetYaxis();
    const double etaLow = etaAx->GetBinLowEdge(ieta + 1);
    const double etaHigh = etaAx->GetBinUpEdge(ieta + 1);
    const double phiLow = phiAx->GetBinLowEdge(iphi + 1);
    const double phiHigh = phiAx->GetBinUpEdge(iphi + 1);

    std::ostringstream nm;
    nm << "eta" << to_string_with_precision(etaLow) << "-" << to_string_with_precision(etaHigh) << "_phi" << to_string_with_precision(phiLow) << "-" << to_string_with_precision(phiHigh);
    return nm.str();
}

// ---------------------------------------------------------------------------
// Retrieve a cell histogram; tries the requested mode then falls back to dR.
// ---------------------------------------------------------------------------
static TH1D *GetCellHist(TFile *fin, const char *dirname, const std::string &prefix, const std::string &mode, const std::string &tag, int ieta, int iphi, const TH2D *hGrid)
{
    if (!fin || !hGrid)
        return nullptr;

    std::vector<std::string> dir_candidates;
    dir_candidates.push_back(dirname);                // legacy inclusive output
    dir_candidates.push_back(prefix + "_" + dirname); // multiplicity-binned outputs

    for (const auto &dir_name : dir_candidates)
    {
        auto *dir = fin->GetDirectory(dir_name.c_str());
        if (!dir)
            continue;

        auto tryGet = [&](const std::string &m) -> TH1D *
        {
            std::string n = prefix + "_" + m + "_" + tag + "_" + EtaPhiRangeNameSuffix(hGrid, ieta, iphi);
            return dynamic_cast<TH1D *>(dir->Get(n.c_str()));
        };

        if (auto *h = tryGet(mode))
            return h; // exact match
        if (mode != "dR")
            if (auto *h = tryGet("dR"))
                return h; // legacy fallback
    }
    return nullptr;
}

static TH2D *GetPrimaryGrid(TFile *fin, const std::string &prefix)
{
    if (!fin)
        return nullptr;
    for (const char *suffix : {"_Sig", "_Bkg", "_Sub"})
    {
        const std::string hname = prefix + suffix;
        if (auto *h = dynamic_cast<TH2D *>(fin->Get(hname.c_str())))
        {
            return h;
        }
    }
    return nullptr;
}
static TH2D *GetRawCountGrid(TFile *fin, const std::string &prefix)
{
    if (!fin)
        return nullptr;
    return dynamic_cast<TH2D *>(fin->Get((prefix + "_Sub_AbsDphiCount").c_str()));
}

static constexpr double kDirectRawCountCut = 0.015;




// ---------------------------------------------------------------------------
// Integral between [xmin, xmax] using bin edges.
// ---------------------------------------------------------------------------
static double IntegralInRange(TH1D *h, double xmin, double xmax)
{
    if (!h || xmax <= xmin)
        return 0.;
    int b1 = h->GetXaxis()->FindBin(xmin + 1e-9);
    int b2 = h->GetXaxis()->FindBin(xmax - 1e-9);
    return (b2 >= b1) ? h->Integral(b1, b2) : 0.;
}

static TH1D *MakeFractionalHist(TH1D *h, TH1D *hRef, const std::string &newname)
{
    if (!h || !hRef)
        return nullptr;

    TH1D *hFrac = dynamic_cast<TH1D *>(h->Clone(newname.c_str()));
    if (!hFrac)
        return nullptr;
    hFrac->SetDirectory(nullptr);

    const int nbins = hFrac->GetNbinsX();
    for (int i = 1; i <= nbins; ++i)
    {
        const double denom = hRef->GetBinContent(i);
        if (denom > 0.0)
        {
            hFrac->SetBinContent(i, h->GetBinContent(i) / denom);
            hFrac->SetBinError(i, h->GetBinError(i) / denom);
        }
        else
        {
            hFrac->SetBinContent(i, 0.0);
            hFrac->SetBinError(i, 0.0);
        }
    }
    return hFrac;
}

static TH1D *MakeInverseHist(TH1D *h, const std::string &newname)
{
    if (!h)
        return nullptr;

    TH1D *hInv = dynamic_cast<TH1D *>(h->Clone(newname.c_str()));
    if (!hInv)
        return nullptr;
    hInv->SetDirectory(nullptr);

    const int nbins = hInv->GetNbinsX();
    for (int i = 1; i <= nbins; ++i)
    {
        const double value = h->GetBinContent(i);
        const double error = h->GetBinError(i);
        if (value != 0.0)
        {
            hInv->SetBinContent(i, 1.0 / value);
            hInv->SetBinError(i, error / (value * value));
        }
        else
        {
            hInv->SetBinContent(i, 0.0);
            hInv->SetBinError(i, 0.0);
        }
    }
    return hInv;
}

static void draw2Dhistogram_adjZaxis(TH2 *hist,                                   //
                                     bool logz,                                   //
                                     const std::string &xtitle,                   //
                                     const std::string &ytitle,                   //
                                     const std::string &ztitle,                   //
                                     std::pair<double, double> zrange,            //
                                     const std::vector<std::string> &addinfo,     //
                                     const std::string &drawoption = "colz text", //
                                     const std::string &paintTextFmt = ".3g",     //
                                     bool show_bin_value_with_error = true,       //
                                     bool show_bin_value_only = false,            //
                                     const std::string &filename = "plot2d")
{
    if (!hist)
        return;

    gStyle->SetPalette(kLightTemperature);
    // gStyle->SetPalette(kThermometer);
    // gStyle->SetPalette(kRainbow);
    // gStyle->SetPalette(kTemperatureMap);

    TCanvas *c = new TCanvas("c_2d", "c_2d", 900, 800);
    c->cd();
    gPad->SetRightMargin(0.21);
    c->SetLogz(logz);

    gStyle->SetPaintTextFormat(paintTextFmt.c_str());
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    hist->GetYaxis()->SetTitle(ytitle.c_str());
    hist->GetZaxis()->SetTitle(ztitle.c_str());
    hist->GetZaxis()->SetTitleOffset(1.5);
    if (zrange.first < zrange.second)
        hist->GetZaxis()->SetRangeUser(zrange.first, zrange.second);

    std::string drawopt = drawoption;
    if (show_bin_value_with_error)
    {
        const std::vector<std::string> tokens_to_strip = {"TEXT", "text"};
        for (const auto &tok : tokens_to_strip)
        {
            size_t pos = std::string::npos;
            while ((pos = drawopt.find(tok)) != std::string::npos)
                drawopt.erase(pos, tok.size());
        }
        if (drawopt.find("colz") == std::string::npos && drawopt.find("COLZ") == std::string::npos)
            drawopt += " colz";
    }
    else if (show_bin_value_only)
    {
        const std::vector<std::string> tokens_to_strip = {"TEXT", "text"};
        for (const auto &tok : tokens_to_strip)
        {
            size_t pos = std::string::npos;
            while ((pos = drawopt.find(tok)) != std::string::npos)
                drawopt.erase(pos, tok.size());
        }
        if (drawopt.find("colz") == std::string::npos && drawopt.find("COLZ") == std::string::npos)
            drawopt += " colz";
    }
    hist->Draw(drawopt.c_str());

    if (show_bin_value_with_error)
    {
        TLatex tbin_val;
        tbin_val.SetTextAlign(22);
        tbin_val.SetTextFont(42);
        tbin_val.SetTextSize(0.018);

        TLatex tbin_err;
        tbin_err.SetTextAlign(22);
        tbin_err.SetTextFont(42);
        tbin_err.SetTextSize(0.015);
        // tbin_err.SetTextColor(kRed + 1);

        for (int ibx = 1; ibx <= hist->GetNbinsX(); ++ibx)
        {
            for (int iby = 1; iby <= hist->GetNbinsY(); ++iby)
            {
                const double val = hist->GetBinContent(ibx, iby);
                const double err = hist->GetBinError(ibx, iby);
                const double x = hist->GetXaxis()->GetBinCenter(ibx);
                const double y = hist->GetYaxis()->GetBinCenter(iby);
                const double y_offset = 0.18 * hist->GetYaxis()->GetBinWidth(iby);

                tbin_val.DrawLatex(x, y + y_offset, Form("%.2f", val));
                tbin_err.DrawLatex(x, y - y_offset, Form("#pm%.3g", err));
            }
        }
    }
    else if (show_bin_value_only)
    {
        TLatex tbin_val;
        tbin_val.SetTextAlign(22);
        tbin_val.SetTextFont(42);
        tbin_val.SetTextSize(0.014);

        for (int ibx = 1; ibx <= hist->GetNbinsX(); ++ibx)
        {
            for (int iby = 1; iby <= hist->GetNbinsY(); ++iby)
            {
                const double val = hist->GetBinContent(ibx, iby);
                const double x = hist->GetXaxis()->GetBinCenter(ibx);
                const double y = hist->GetYaxis()->GetBinCenter(iby);
                tbin_val.DrawLatex(x, y, Form("%.3g", val));
            }
        }
    }

    // for now, the format is only good for 1 line
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    for (size_t i = 0; i < addinfo.size(); ++i)
    {
        latex.DrawLatex(gPad->GetLeftMargin(), (1 - gPad->GetTopMargin() + 0.01) - 0.045 * i, addinfo[i].c_str());
    }

    c->SaveAs((filename + ".png").c_str());
    c->SaveAs((filename + ".pdf").c_str());
    delete c;
}

static void draw1Dhistogram(TH1 *hist,                                 //
                            const bool logy,                           //
                            const std::string &xtitle,                 //
                            const std::string &ytitle,                 //
                            const std::vector<std::string> &addinfo,   //
                            const std::string &drawoption = "hist e1", //
                            const std::string &filename = "plot1d")
{
    if (!hist)
        return;

    TCanvas *c = new TCanvas("c_1d", "c_1d", 900, 800);
    c->SetLogy(logy);
    c->cd();
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.05);
    // gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.08);

    hist->SetTitle("");
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    hist->GetYaxis()->SetTitle(ytitle.c_str());
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->GetYaxis()->SetRangeUser((logy ? hist->GetMinimum(0) * 0.5 : 0.0), hist->GetMaximum() * (logy ? 10.0 : 1.4));
    hist->SetLineWidth(2);
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(kFullCircle);
    hist->SetMarkerSize(1.0);
    hist->SetMarkerColor(kBlack);
    hist->Draw(drawoption.c_str());

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    // set right- and bottom aligned for now, which is only good for 1 line
    latex.SetTextAlign(kHAlignRight + kVAlignBottom);
    for (size_t i = 0; i < addinfo.size(); ++i)
    {
        latex.DrawLatex(1 - gPad->GetRightMargin(), (1 - gPad->GetTopMargin() + 0.01) - 0.045 * i, addinfo[i].c_str());
    }

    c->SaveAs((filename + ".png").c_str());
    c->SaveAs((filename + ".pdf").c_str());
    delete c;
}

static void drawComparisonWithRatio(std::pair<TH1D *, std::string> hDenInput, //
                                    std::pair<TH1D *, std::string> hNumInput, //
                                    const bool logy,                          //
                                    const std::string &xtitle,                //
                                    const std::string &ytitle,                //
                                    const std::string &ratio_ytitle,          //
                                    const std::vector<std::string> &addinfo,  //
                                    const std::string &filename)              //
{
    if (!hDenInput.first || !hNumInput.first)
        return;

    TH1D *hDen = dynamic_cast<TH1D *>(hDenInput.first->Clone((filename + "_den").c_str()));
    TH1D *hNum = dynamic_cast<TH1D *>(hNumInput.first->Clone((filename + "_num").c_str()));
    TH1D *hRatio = dynamic_cast<TH1D *>(hNumInput.first->Clone((filename + "_ratio").c_str()));
    if (!hDen || !hNum || !hRatio)
    {
        delete hDen;
        delete hNum;
        delete hRatio;
        return;
    }

    hDen->SetDirectory(nullptr);
    hNum->SetDirectory(nullptr);
    hRatio->SetDirectory(nullptr);
    hRatio->Reset("ICES");
    hRatio->Divide(hNum, hDen, 1.0, 1.0, "B");

    hDen->SetLineColor(kBlack);
    hDen->SetMarkerColor(kBlack);
    hDen->SetMarkerStyle(kFullCircle);
    hDen->SetMarkerSize(0.8);
    hDen->SetLineWidth(2);

    hNum->SetLineColor(kRed + 1);
    hNum->SetMarkerColor(kRed + 1);
    hNum->SetMarkerStyle(kOpenCircle);
    hNum->SetMarkerSize(0.8);
    hNum->SetLineWidth(2);

    hRatio->SetLineColor(kRed + 1);
    hRatio->SetMarkerColor(kRed + 1);
    hRatio->SetMarkerStyle(kOpenCircle);
    hRatio->SetMarkerSize(0.8);
    hRatio->SetLineWidth(2);

    double ymax = std::max(hDen->GetMaximum(), hNum->GetMaximum());
    if (ymax <= 0.0)
        ymax = 1.0;

    double ymin = std::min(hDen->GetMinimum(0), hNum->GetMinimum(0));

    TCanvas *c = new TCanvas((filename + "_canvas").c_str(), (filename + "_canvas").c_str(), 900, 800);
    c->cd();
    TPad *pad_top = new TPad("pad_top", "pad_top", 0.0, 0.32, 1.0, 1.0);
    TPad *pad_bottom = new TPad("pad_bottom", "pad_bottom", 0.0, 0.0, 1.0, 0.32);
    pad_top->SetLeftMargin(0.14);
    pad_top->SetRightMargin(0.05);
    pad_top->SetTopMargin(0.05);
    pad_top->SetBottomMargin(0.025);
    pad_bottom->SetLeftMargin(0.14);
    pad_bottom->SetRightMargin(0.05);
    pad_bottom->SetTopMargin(0.025);
    pad_bottom->SetBottomMargin(0.35);
    pad_top->Draw();
    pad_bottom->Draw();

    float toptitlesize = 0.05;
    pad_top->cd();
    pad_top->SetLogy(logy);
    hDen->SetTitle("");
    hDen->GetYaxis()->SetTitle(ytitle.c_str());
    hDen->GetYaxis()->SetTitleOffset(1.35);
    hDen->GetYaxis()->SetTitleSize(toptitlesize);
    hDen->GetYaxis()->SetLabelSize(toptitlesize);
    hDen->GetXaxis()->SetLabelSize(0.0);
    hDen->GetXaxis()->SetTitleSize(0.0);
    // hDen->SetMinimum((logy ? ymin * 0.1 : 0.0));
    hDen->SetMaximum((logy ? 100.0 : 1.5) * ymax);
    hDen->Draw("e1");
    hNum->Draw("e1 same");

    TLegend *leg = new TLegend(gPad->GetLeftMargin() + 0.03,   //
                               1 - gPad->GetTopMargin() - 0.2, //
                               gPad->GetLeftMargin() + 0.15,   //
                               1 - gPad->GetTopMargin() - 0.05 //
    );
    leg->SetHeader(addinfo[0].c_str(), "L"); // assuming addinfo has only 1 entry
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hDen, hDenInput.second.c_str(), "pel");
    leg->AddEntry(hNum, hNumInput.second.c_str(), "pel");
    leg->Draw();

    // TLatex latex;
    // latex.SetNDC();
    // latex.SetTextFont(42);
    // latex.SetTextSize(0.035);
    // latex.SetTextAlign(kHAlignRight + kVAlignBottom);
    // for (size_t i = 0; i < addinfo.size(); ++i)
    // {
    //     latex.DrawLatex(1 - pad_top->GetRightMargin(), (1 - pad_top->GetTopMargin() + 0.01) - 0.045 * i, addinfo[i].c_str());
    // }

    float textsize_scalefactor = 2.1;
    pad_bottom->cd();
    hRatio->SetTitle("");
    hRatio->GetXaxis()->SetTitle(xtitle.c_str());
    hRatio->GetYaxis()->SetTitle(ratio_ytitle.c_str());
    hRatio->GetXaxis()->SetTitleSize(toptitlesize * textsize_scalefactor);
    hRatio->GetXaxis()->SetLabelSize(toptitlesize * textsize_scalefactor);
    // hRatio->GetXaxis()->SetTitleOffset(1.1);
    hRatio->GetYaxis()->SetTitleSize(toptitlesize * textsize_scalefactor);
    hRatio->GetYaxis()->SetLabelSize(toptitlesize * textsize_scalefactor);
    hRatio->GetYaxis()->SetTitleOffset(0.6);
    hRatio->GetYaxis()->SetNdivisions(505);
    // hRatio->SetMinimum(0);
    hRatio->SetMaximum(1);
    hRatio->Draw("e1");

    TLine *line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");

    c->SaveAs((filename + ".png").c_str());
    c->SaveAs((filename + ".pdf").c_str());

    delete line;
    delete leg;
    delete hRatio;
    delete hNum;
    delete hDen;
    delete c;
}

// ===========================================================================
static void draw_comparison(std::vector<TH1 *> histograms,               //
                            const std::vector<std::string> &colors,      //
                            std::vector<std::string> labels,             //
                            const std::string &xAxisTitle,               //
                            const std::string &yAxisTitle,               //
                            const std::pair<float, float> &yRange,       //
                            const bool logy,                             //
                            const bool normalize,                        //
                            const int xAxisNdivisions,                   //
                            const std::pair<bool, float> &referenceline, //
                            const std::string &outputFileName            //
)
{
    if (histograms.empty())
        return;

    if (normalize)
    {
        for (auto *hist : histograms)
        {
            if (hist && hist->Integral(-1, -1) > 0.0)
                hist->Scale(1. / hist->Integral(-1, -1));
        }
    }

    float yMin = yRange.first, yMax = yRange.second;
    if (yRange.first < 0)
    {
        float histMin = 1E6;
        for (auto *hist : histograms)
        {
            if (hist && hist->GetMinimum(0) < histMin)
                histMin = hist->GetMinimum(0);
        }
        yMin = logy ? histMin * 0.5f : 0.f;
    }

    if (yRange.second < 0)
    {
        float histMax = -1E6;
        for (auto *hist : histograms)
        {
            if (hist && hist->GetMaximum() > histMax)
                histMax = hist->GetMaximum();
        }
        yMax = logy ? histMax * 250.f : histMax * 1.5f;
    }

    TCanvas *c = new TCanvas("c_comparison", "Comparison", 800, 700);
    c->SetLogy(logy);
    c->cd();
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.07);
    for (size_t i = 0; i < histograms.size(); ++i)
    {
        if (!histograms[i])
            continue;
        histograms[i]->SetLineColor(TColor::GetColor(colors[i].c_str()));
        histograms[i]->SetMarkerColor(TColor::GetColor(colors[i].c_str()));
        histograms[i]->SetMarkerStyle(20 + static_cast<int>(i));
        histograms[i]->SetMarkerSize(1.0);
        histograms[i]->SetLineWidth(2);
        if (i == 0)
        {
            histograms[i]->SetTitle("");
            histograms[i]->GetXaxis()->SetTitle(xAxisTitle.c_str());
            histograms[i]->GetXaxis()->SetTitleOffset(1.4);
            histograms[i]->GetXaxis()->SetNdivisions(xAxisNdivisions);
            histograms[i]->GetYaxis()->SetTitle(yAxisTitle.c_str());
            histograms[i]->GetYaxis()->SetRangeUser(yMin, yMax);
            histograms[i]->GetYaxis()->SetTitleOffset(1.5);
            histograms[i]->Draw("HIST E1");
        }
        else
        {
            histograms[i]->Draw("HIST E1 SAME");
        }
    }

    // Draw reference line if requested
    if (referenceline.first)
    {
        TLine *line = new TLine(histograms[0]->GetXaxis()->GetXmin(), referenceline.second, histograms[0]->GetXaxis()->GetXmax(), referenceline.second);
        line->SetLineColor(kGray + 1);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");
    }

    TLegend *legend = new TLegend(gPad->GetLeftMargin() + 0.04,                                        //
                                  (1 - gPad->GetTopMargin()) - 0.035 - 0.04 * (histograms.size() + 1), //
                                  gPad->GetLeftMargin() + 0.14,                                        //
                                  (1 - gPad->GetTopMargin()) - 0.035                                   //
    );
    legend->SetTextAlign(kHAlignLeft + kVAlignCenter);
    legend->SetTextSize(0.035);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    for (size_t i = 0; i < labels.size() && i < histograms.size(); ++i)
    {
        if (histograms[i])
            legend->AddEntry(histograms[i], labels[i].c_str(), "pel");
    }
    legend->Draw();

    c->SaveAs(Form("%s.png", outputFileName.c_str()));
    c->SaveAs(Form("%s.pdf", outputFileName.c_str()));
    delete legend;
    delete c;
}

// ===========================================================================
// Alternative version of draw_comparison for the INTT doublet DCA
static void draw_comparison_alt(std::vector<TH1 *> histograms,               //
                                const std::vector<std::string> &colors,      //
                                const std::vector<int> &markers,             //
                                const std::vector<float> &markerSizes,       //
                                std::vector<std::string> labels,             //
                                const std::string &xAxisTitle,               //
                                const std::string &yAxisTitle,               //
                                const std::pair<float, float> &yRange,       //
                                const bool logy,                             //
                                const bool normalize,                        //
                                const int xAxisNdivisions,                   //
                                const std::pair<bool, float> &referenceline, //
                                const std::string &outputFileName            //
)
{
    // make a dummy histogram for axis setting
    TH1 *hDummy = new TH1D("hDummy", "", 1, histograms[0]->GetXaxis()->GetXmin(), histograms[0]->GetXaxis()->GetXmax());

    if (histograms.empty())
        return;

    if (normalize)
    {
        for (auto *hist : histograms)
        {
            if (hist && hist->Integral(-1, -1) > 0.0)
                hist->Scale(1. / hist->Integral(-1, -1));
        }
    }

    float yMin = yRange.first, yMax = yRange.second;
    if (yRange.first < 0)
    {
        float histMin = std::numeric_limits<float>::max();
        for (auto *hist : histograms)
        {
            std::cout << "hist min: " << (hist ? hist->GetMinimum(0) : -1) << std::endl;
            if (hist && hist->GetMinimum(0) < histMin)
                histMin = hist->GetMinimum(0);
        }
        yMin = logy ? histMin * 0.5f : 0.f;
    }

    if (yRange.second < 0)
    {
        float histMax = -std::numeric_limits<float>::max();
        for (auto *hist : histograms)
        {
            if (hist && hist->GetMaximum() > histMax)
                histMax = hist->GetMaximum();
        }
        yMax = logy ? histMax * 2.f : histMax * 1.5f;
    }

    std::cout << "yMin: " << yMin << ", yMax: " << yMax << std::endl;

    TCanvas *c = new TCanvas("c_comparison", "Comparison", 800, 700);
    c->SetLogy(logy);
    c->cd();
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.07);
    // draw dummy first to set axis
    hDummy->SetTitle("");
    hDummy->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hDummy->GetXaxis()->SetTitleOffset(1.4);
    hDummy->GetXaxis()->SetNdivisions(xAxisNdivisions);
    hDummy->GetYaxis()->SetTitle(yAxisTitle.c_str());
    hDummy->GetYaxis()->SetRangeUser(yMin, yMax);
    hDummy->GetYaxis()->SetTitleOffset(1.5);
    hDummy->Draw("AXIS");

    for (size_t i = 0; i < histograms.size(); ++i)
    {
        if (!histograms[i])
            continue;
        histograms[i]->SetLineColor(TColor::GetColor(colors[i].c_str()));
        histograms[i]->SetMarkerColor(TColor::GetColor(colors[i].c_str()));
        histograms[i]->SetMarkerStyle(markers[i]);
        histograms[i]->SetMarkerSize(markerSizes[i]);
        histograms[i]->SetLineWidth(2);
        histograms[i]->Draw("E1 SAME");
    }

    // Draw reference line if requested
    if (referenceline.first)
    {
        TLine *line = new TLine(histograms[0]->GetXaxis()->GetXmin(), referenceline.second, histograms[0]->GetXaxis()->GetXmax(), referenceline.second);
        line->SetLineColor(kGray + 1);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");
    }

    TLegend *legend = new TLegend(gPad->GetLeftMargin() + 0.25,                                        //
                                  (1 - gPad->GetTopMargin()) - 0.035 - 0.04 * (histograms.size() + 1), //
                                  gPad->GetLeftMargin() + 0.35,                                        //
                                  (1 - gPad->GetTopMargin()) - 0.035                                   //
    );
    legend->SetTextAlign(kHAlignLeft + kVAlignCenter);
    legend->SetTextSize(0.035);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    for (size_t i = 0; i < labels.size() && i < histograms.size(); ++i)
    {
        if (histograms[i])
            legend->AddEntry(histograms[i], labels[i].c_str(), "pel");
    }
    legend->Draw();

    c->SaveAs(Form("%s.png", outputFileName.c_str()));
    c->SaveAs(Form("%s.pdf", outputFileName.c_str()));
    delete legend;
    delete c;
}

struct ProjectionRatioOutputs
{
    std::string prefix;
    std::string label;
    TH1D *eta_ratio = nullptr;
    TH1D *phi_ratio = nullptr;
    TH1D *eta_ratio_good = nullptr;
    TH1D *phi_ratio_good = nullptr;
};

static bool StartsWith(const std::string &value, const std::string &prefix) { return value.rfind(prefix, 0) == 0; }

static bool EndsWith(const std::string &value, const std::string &suffix) { return value.size() >= suffix.size() && value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0; }

static int ExtractMultiplicityBinNumber(const std::string &base_prefix, const std::string &object_name)
{
    const std::string probe_prefix = base_prefix + "_multPercentileBin";
    const std::string suffix = "_Sub";
    if (!StartsWith(object_name, probe_prefix) || !EndsWith(object_name, suffix))
        return -1;

    const std::string bin_text = object_name.substr(probe_prefix.size(), object_name.size() - probe_prefix.size() - suffix.size());
    if (bin_text.empty())
        return -1;
    for (char ch : bin_text)
    {
        if (!std::isdigit(static_cast<unsigned char>(ch)))
            return -1;
    }
    return std::stoi(bin_text);
}

static std::vector<std::string> CollectAvailablePrefixes(TFile *fin, const std::string &base_prefix)
{
    std::vector<std::pair<int, std::string>> multiplicity_prefixes;
    if (!fin)
        return {base_prefix};

    TIter next(fin->GetListOfKeys());
    while (auto *key = dynamic_cast<TKey *>(next()))
    {
        const std::string name = key->GetName();
        const int bin_number = ExtractMultiplicityBinNumber(base_prefix, name);
        if (bin_number >= 0)
            multiplicity_prefixes.push_back({bin_number, base_prefix + "_multPercentileBin" + std::to_string(bin_number)});
    }

    std::sort(multiplicity_prefixes.begin(), multiplicity_prefixes.end(), [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });
    multiplicity_prefixes.erase(std::unique(multiplicity_prefixes.begin(), multiplicity_prefixes.end(), [](const auto &lhs, const auto &rhs) { return lhs.first == rhs.first; }), multiplicity_prefixes.end());

    std::vector<std::string> prefixes = {base_prefix};
    for (const auto &entry : multiplicity_prefixes)
        prefixes.push_back(entry.second);
    return prefixes;
}

static std::string PrefixSuffixForOutput(const std::string &base_prefix, const std::string &current_prefix)
{
    if (current_prefix == base_prefix)
        return "";
    if (StartsWith(current_prefix, base_prefix))
        return current_prefix.substr(base_prefix.size());
    return "_" + current_prefix;
}

static std::string PrefixDisplayLabel(const std::string &base_prefix, const std::string &current_prefix)
{
    if (current_prefix == base_prefix)
        return "Inclusive";
    if (StartsWith(current_prefix, base_prefix + "_"))
        return current_prefix.substr(base_prefix.size() + 1);
    return current_prefix;
}

static int ExtractMultiplicityBinFromPrefix(const std::string &base_prefix, const std::string &current_prefix)
{
    const std::string probe_prefix = base_prefix + "_multPercentileBin";
    if (!StartsWith(current_prefix, probe_prefix))
        return -1;

    const std::string bin_text = current_prefix.substr(probe_prefix.size());
    if (bin_text.empty())
        return -1;
    for (char ch : bin_text)
    {
        if (!std::isdigit(static_cast<unsigned char>(ch)))
            return -1;
    }
    return std::stoi(bin_text);
}

static std::string FormatPercentileValue(double value)
{
    std::ostringstream out;
    if (std::fabs(value - std::round(value)) < 1e-6)
        out << static_cast<int>(std::round(value));
    else
        out << std::fixed << std::setprecision(1) << value;
    return out.str();
}

static std::map<int, std::pair<double, double>> LoadPercentileBoundariesByBin(const bool is_simulation)
{
    std::map<int, std::pair<double, double>> boundaries;
    const std::string multiplicity_boundary_prefix = is_simulation ? "MC_" : "Data_";
    const std::string input_file = "figure/figure-NInttClusterCrossing/" + multiplicity_boundary_prefix + "NInttClusterPercentileBoundaries.root";

    TFile fin(input_file.c_str(), "READ");
    if (fin.IsZombie())
        return boundaries;

    auto *tree = dynamic_cast<TTree *>(fin.Get("trigger_percentile_intervals"));
    if (!tree)
        return boundaries;

    int percentile_bin = -1;
    double percentile_low = 0.0;
    double percentile_high = 0.0;
    tree->SetBranchAddress("percentile_bin", &percentile_bin);
    tree->SetBranchAddress("percentile_low", &percentile_low);
    tree->SetBranchAddress("percentile_high", &percentile_high);

    const Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i)
    {
        tree->GetEntry(i);
        boundaries[percentile_bin] = {percentile_low, percentile_high};
    }

    return boundaries;
}

static std::string PrefixDisplayLabel(const std::string &base_prefix, const std::string &current_prefix, const std::map<int, std::pair<double, double>> &percentile_boundaries)
{
    if (current_prefix == base_prefix)
        return "Inclusive";

    const int bin_number = ExtractMultiplicityBinFromPrefix(base_prefix, current_prefix);
    if (bin_number >= 0)
    {
        auto it = percentile_boundaries.find(bin_number);
        if (it != percentile_boundaries.end())
        {
            return "N_{INTT clusters} " + FormatPercentileValue(it->second.first) + "-" + FormatPercentileValue(it->second.second) + "%";
        }
    }

    if (StartsWith(current_prefix, base_prefix + "_"))
        return current_prefix.substr(base_prefix.size() + 1);
    return current_prefix;
}

static double DecodeTruthVtxZHistToken(const std::string &token, bool &is_infinite, bool &ok)
{
    is_infinite = false;
    ok = false;
    if (token == "mInf")
    {
        is_infinite = true;
        ok = true;
        return -std::numeric_limits<double>::infinity();
    }
    if (token == "pInf")
    {
        is_infinite = true;
        ok = true;
        return std::numeric_limits<double>::infinity();
    }

    std::string value = token;
    double sign = 1.0;
    if (!value.empty() && value[0] == 'm')
    {
        sign = -1.0;
        value.erase(0, 1);
    }
    else if (!value.empty() && value[0] == 'p')
    {
        value.erase(0, 1);
    }

    std::replace(value.begin(), value.end(), 'p', '.');
    char *endptr = nullptr;
    const double parsed = std::strtod(value.c_str(), &endptr);
    if (endptr && *endptr == 0)
    {
        ok = true;
        return sign * parsed;
    }
    return 0.0;
}

struct TruthVtxZDCAHistBin
{
    int bin_index = -1;
    double xlow = 0.0;
    double xhigh = 0.0;
    TH1D *hist = nullptr;
};

static bool ParseTruthVtxZDCAHistName(const std::string &hist_name, const std::string &prefix, TruthVtxZDCAHistBin &info)
{
    const std::string probe_prefix = prefix + "_SilSeedDCA3D_truthVtxZ_bin";
    if (!StartsWith(hist_name, probe_prefix))
        return false;

    const std::string remainder = hist_name.substr(probe_prefix.size());
    if (remainder.size() < 4 || !std::isdigit(static_cast<unsigned char>(remainder[0])) || !std::isdigit(static_cast<unsigned char>(remainder[1])) || remainder[2] != '_')
        return false;

    const int bin_index = std::stoi(remainder.substr(0, 2));
    const std::string encoded_bounds = remainder.substr(3);
    const size_t sep = encoded_bounds.find('_');
    if (sep == std::string::npos)
        return false;

    bool low_inf = false;
    bool high_inf = false;
    bool low_ok = false;
    bool high_ok = false;
    const double low_raw = DecodeTruthVtxZHistToken(encoded_bounds.substr(0, sep), low_inf, low_ok);
    const double high_raw = DecodeTruthVtxZHistToken(encoded_bounds.substr(sep + 1), high_inf, high_ok);
    if (!low_ok || !high_ok)
        return false;

    const double low = low_inf ? (high_raw - 10.0) : low_raw;
    const double high = high_inf ? (low_raw + 10.0) : high_raw;
    if (!(high > low))
        return false;

    info.bin_index = bin_index;
    info.xlow = low;
    info.xhigh = high;
    return true;
}

static std::vector<TruthVtxZDCAHistBin> CollectTruthVtxZDCAHistBins(TFile *fin, const std::string &prefix)
{
    std::vector<TruthVtxZDCAHistBin> bins;
    if (!fin)
        return bins;

    TIter next(fin->GetListOfKeys());
    while (auto *key = dynamic_cast<TKey *>(next()))
    {
        const std::string hist_name = key->GetName();
        TruthVtxZDCAHistBin info;
        if (!ParseTruthVtxZDCAHistName(hist_name, prefix, info))
            continue;

        info.hist = dynamic_cast<TH1D *>(fin->Get(hist_name.c_str()));
        if (!info.hist)
            continue;
        bins.push_back(info);
    }

    std::sort(bins.begin(), bins.end(), [](const TruthVtxZDCAHistBin &lhs, const TruthVtxZDCAHistBin &rhs) { return lhs.bin_index < rhs.bin_index; });
    return bins;
}

static TH2D *BuildTruthVtxZDCAHeatmap(TFile *fin, const std::string &prefix)
{
    std::vector<TruthVtxZDCAHistBin> bins = CollectTruthVtxZDCAHistBins(fin, prefix);
    if (bins.empty())
        return nullptr;

    TH1D *href = bins.front().hist;
    if (!href)
        return nullptr;

    std::vector<double> xedges;
    xedges.reserve(bins.size() + 1);
    xedges.push_back(bins.front().xlow);
    for (const auto &entry : bins)
        xedges.push_back(entry.xhigh);

    std::vector<double> yedges;
    yedges.reserve(href->GetNbinsX() + 1);
    for (int ibin = 1; ibin <= href->GetNbinsX(); ++ibin)
        yedges.push_back(href->GetXaxis()->GetBinLowEdge(ibin));
    yedges.push_back(href->GetXaxis()->GetBinUpEdge(href->GetNbinsX()));

    TH2D *h2 = new TH2D((prefix + "_SilSeedDCA3D_truthVtxZ_2D").c_str(), (prefix + "_SilSeedDCA3D_truthVtxZ_2D").c_str(), static_cast<int>(bins.size()), xedges.data(), href->GetNbinsX(), yedges.data());
    h2->SetDirectory(nullptr);
    h2->Sumw2();

    for (size_t ix = 0; ix < bins.size(); ++ix)
    {
        TH1D *h1 = bins[ix].hist;
        if (!h1)
            continue;

        const double norm = h1->Integral(1, h1->GetNbinsX());
        for (int iy = 1; iy <= h1->GetNbinsX(); ++iy)
        {
            const double value = (norm > 0.0) ? (h1->GetBinContent(iy) / norm) : 0.0;
            const double error = (norm > 0.0) ? (h1->GetBinError(iy) / norm) : 0.0;
            h2->SetBinContent(static_cast<int>(ix) + 1, iy, value);
            h2->SetBinError(static_cast<int>(ix) + 1, iy, error);
        }
    }

    return h2;
}

static std::vector<std::pair<TH1D *, std::string>> BuildMergedAbsTruthVtxZDCAHists(TFile *fin, const std::string &prefix)
{
    std::vector<std::pair<TH1D *, std::string>> merged_hists;
    const std::vector<TruthVtxZDCAHistBin> bins = CollectTruthVtxZDCAHistBins(fin, prefix);
    if (bins.empty())
        return merged_hists;

    auto format_abs_bound = [](double value) -> std::string
    {
        std::ostringstream out;
        out << std::fixed << std::setprecision((std::fabs(value - std::round(value)) < 1e-6) ? 0 : 1) << value;
        return out.str();
    };

    const size_t npairs = (bins.size() + 1) / 2;
    merged_hists.reserve(npairs);
    for (size_t i = 0; i < npairs; ++i)
    {
        const size_t j = bins.size() - 1 - i;
        TH1D *href = bins[i].hist;
        if (!href)
            continue;

        TH1D *hmerge = dynamic_cast<TH1D *>(href->Clone((prefix + "_SilSeedDCA3D_absTruthVtxZ_bin" + std::to_string(i)).c_str()));
        if (!hmerge)
            continue;
        hmerge->SetDirectory(nullptr);
        if (j != i && bins[j].hist)
            hmerge->Add(bins[j].hist);

        const double abs_low = (j != i) ? std::min(std::fabs(bins[i].xhigh), std::fabs(bins[j].xlow)) : std::min(std::fabs(bins[i].xlow), std::fabs(bins[i].xhigh));
        const double abs_high = (j != i) ? std::max(std::fabs(bins[i].xlow), std::fabs(bins[j].xhigh)) : std::max(std::fabs(bins[i].xlow), std::fabs(bins[i].xhigh));
        merged_hists.push_back({hmerge, "|z_{truth}|=" + format_abs_bound(abs_low) + "-" + format_abs_bound(abs_high) + " cm"});
    }

    return merged_hists;
}

static ProjectionRatioOutputs plot_tklCombinatoric_single(const std::string &inputfile = "test_OO_combinatoric.root",     //
                                                          const std::string &prefix = "tkl_Combinatoric",                 //
                                                          const std::string &mode = "absdPhi",                            //
                                                          bool is_simulation = false,                                     //
                                                          double zvtx_cut_min = -10.0,                                    //
                                                          double zvtx_cut_max = 10.0,                                     //
                                                          double raw_count_xmin = 0.0,                                    //
                                                          double raw_count_xmax = kDirectRawCountCut,                     //
                                                          double plot_xmin = 0.0,                                         //
                                                          double plot_xmax = 0.2,                                         //
                                                          const std::string &outstem = "plot_tklCombinatoric",            //
                                                          const std::string &plotdir = "./figure/figure-tklCombinatoric", //
                                                          bool zero_x_error = true,                                       //
                                                          const std::string &ratio_draw_option = "colz text",             //
                                                          const std::string &ratio_text_format = ".3f",                   //
                                                          bool use_rounded_zmax = true,                                   //
                                                          const std::string &display_label = "")
{
    ProjectionRatioOutputs ratio_outputs;
    ratio_outputs.prefix = prefix;
    ratio_outputs.label = display_label.empty() ? prefix : display_label;

    gStyle->SetOptStat(0);
    TGaxis::SetMaxDigits(3);
    if (zero_x_error)
        gStyle->SetErrorX(0.f);

    gSystem->mkdir(plotdir.c_str(), true);
    const std::string single_canvas_dir = plotdir + "/single-canvas";
    gSystem->mkdir(single_canvas_dir.c_str(), true);

    auto *fin = TFile::Open(inputfile.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "[plot_tklCombinatoric] Cannot open: " << inputfile << "\n";
        return ratio_outputs;
    }

    auto *hGrid = GetPrimaryGrid(fin, prefix);
    if (!hGrid)
    {
        std::cerr << "[plot_tklCombinatoric] Cannot load " << prefix << "_Sig/_Bkg/_Sub for grid axes.\n";
        fin->Close();
        return ratio_outputs;
    }
    const int nEta = hGrid->GetXaxis()->GetNbins();
    const int nPhi = hGrid->GetYaxis()->GetNbins();

    auto *hSub = dynamic_cast<TH2D *>(fin->Get((prefix + "_Sub").c_str()));
    if (!hSub)
    {
        std::cerr << "[plot_tklCombinatoric] Missing " << prefix << "_Sub.\n";
        fin->Close();
        return ratio_outputs;
    }

    auto *hSilSeed = dynamic_cast<TH2D *>(fin->Get((prefix + "_SilSeedEtaPhi").c_str()));
    if (!hSilSeed)
    {
        std::cerr << "[plot_tklCombinatoric] Missing " << prefix << "_SilSeedEtaPhi.\n";
        fin->Close();
        return ratio_outputs;
    }
    auto *hSilSeed_nMVTX3nINTT2 = dynamic_cast<TH2D *>(fin->Get((prefix + "_SilSeedEtaPhi_nMVTX3nINTT2").c_str()));
    auto *hSilSeed_nMVTX3nINTT1 = dynamic_cast<TH2D *>(fin->Get((prefix + "_SilSeedEtaPhi_nMVTX3nINTT1").c_str()));
    auto *hSilSeed_nMVTX2nINTT2 = dynamic_cast<TH2D *>(fin->Get((prefix + "_SilSeedEtaPhi_nMVTX2nINTT2").c_str()));
    auto *hSilSeed_nMVTX2nINTT1 = dynamic_cast<TH2D *>(fin->Get((prefix + "_SilSeedEtaPhi_nMVTX2nINTT1").c_str()));
    auto *hSilSeed_nMVTXnINTTOther = dynamic_cast<TH2D *>(fin->Get((prefix + "_SilSeedEtaPhi_nMVTXnINTTOther").c_str()));
    auto *hINTTClusterVsAssociatedSilSeed = dynamic_cast<TH2D *>(fin->Get((prefix + "_INTTClusterVsAssociatedSilSeed").c_str()));
    auto *hINTTClusterVsAllSilSeed = dynamic_cast<TH2D *>(fin->Get((prefix + "_INTTClusterVsAllSilSeed").c_str()));
    auto *hSeedCrossing_selectedCrossings = dynamic_cast<TH1D *>(fin->Get((prefix + "_SeedCrossing_selectedCrossings").c_str()));
    auto *hSeedCrossing_selectedCrossings_largeDiffNClusterSeeds = dynamic_cast<TH1D *>(fin->Get((prefix + "_SeedCrossing_selectedCrossings_largeDiffNClusterSeeds").c_str()));
    auto *hDCA3D_Sig = dynamic_cast<TH1D *>(fin->Get((prefix + "_DCA3D_Sig").c_str()));
    auto *hDCA3D_Bkg = dynamic_cast<TH1D *>(fin->Get((prefix + "_DCA3D_Bkg").c_str()));
    auto *hSilSeedDCA3D = dynamic_cast<TH1D *>(fin->Get((prefix + "_SilSeedDCA3D").c_str()));
    TH2D *hSilSeedDCA3D_truthVtxZ_2D = is_simulation ? BuildTruthVtxZDCAHeatmap(fin, prefix) : nullptr;
    std::vector<std::pair<TH1D *, std::string>> hSilSeedDCA3D_truthVtxZ_absMerged = is_simulation ? BuildMergedAbsTruthVtxZDCAHists(fin, prefix) : std::vector<std::pair<TH1D *, std::string>>{};

    if (!hSilSeed_nMVTX3nINTT2 || !hSilSeed_nMVTX3nINTT1 || !hSilSeed_nMVTX2nINTT2 || !hSilSeed_nMVTX2nINTT1 || !hSilSeed_nMVTXnINTTOther)
    {
        std::cerr << "[plot_tklCombinatoric] Missing one or more silicon-seed category histograms." << std::endl;
        fin->Close();
        return ratio_outputs;
    }

    // Sum the high-quality seed categories: (3+2) + (3+1) + (2+2).
    auto *hSilSeed_good = dynamic_cast<TH2D *>(hSilSeed_nMVTX3nINTT2->Clone((prefix + "_SilSeedEtaPhi_good").c_str()));
    if (!hSilSeed_good)
    {
        std::cerr << "[plot_tklCombinatoric] Failed to build " << prefix << "_SilSeedEtaPhi_good." << std::endl;
        fin->Close();
        return ratio_outputs;
    }
    hSilSeed_good->SetDirectory(nullptr);
    hSilSeed_good->Add(hSilSeed_nMVTX3nINTT1);
    hSilSeed_good->Add(hSilSeed_nMVTX2nINTT2);

    // this is cluster on the first layer, not "silicon seeds"
    auto *hSeed = dynamic_cast<TH2D *>(fin->Get((prefix + "_SeedEtaPhi").c_str()));
    if (!hSeed)
    {
        std::cerr << "[plot_tklCombinatoric] Missing " << prefix << "_SeedEtaPhi.\n";
        fin->Close();
        return ratio_outputs;
    }

    TH2D *hPrimaryCharged = nullptr;
    if (is_simulation)
    {
        hPrimaryCharged = dynamic_cast<TH2D *>(fin->Get((prefix + "_PrimaryChargedEtaPhi").c_str()));
        if (!hPrimaryCharged)
        {
            std::cerr << "[plot_tklCombinatoric] Missing " << prefix << "_PrimaryChargedEtaPhi for simulation mode.\n";
            fin->Close();
            return ratio_outputs;
        }
    }

    bool use_direct_raw_count_grid = false;
    TH2D *h2_nraw = nullptr;
    if (auto *hRawDirect = GetRawCountGrid(fin, prefix))
    {
        h2_nraw = dynamic_cast<TH2D *>(hRawDirect->Clone((prefix + "_NRawEtaPhi").c_str()));
        if (h2_nraw)
        {
            h2_nraw->SetDirectory(nullptr);
            use_direct_raw_count_grid = true;
        }
    }
    if (!h2_nraw)
    {
        h2_nraw = dynamic_cast<TH2D *>(hSub->Clone((prefix + "_NRawEtaPhi").c_str()));
        h2_nraw->SetDirectory(nullptr);
        h2_nraw->Reset("ICES");
        std::cerr << "[plot_tklCombinatoric] Missing " << prefix << "_Sub_AbsDphiCount. Falling back to per-cell integration for h2_nraw.\n";
    }
    if (h2_nraw->GetNbinsX() != hSilSeed->GetNbinsX() || h2_nraw->GetNbinsY() != hSilSeed->GetNbinsY())
    {
        std::cerr << "[plot_tklCombinatoric] Incompatible binning between " << prefix << "_Sub and " << prefix << "_SilSeedEtaPhi.\n";
        delete h2_nraw;
        fin->Close();
        return ratio_outputs;
    }
    if (is_simulation && (h2_nraw->GetNbinsX() != hPrimaryCharged->GetNbinsX() || h2_nraw->GetNbinsY() != hPrimaryCharged->GetNbinsY()))
    {
        std::cerr << "[plot_tklCombinatoric] Incompatible binning between " << prefix << "_NRawEtaPhi and " << prefix << "_PrimaryChargedEtaPhi.\n";
        delete h2_nraw;
        fin->Close();
        return ratio_outputs;
    }
    auto check_same_binning = [](const TH2D *ha, const TH2D *hb) -> bool { return ha && hb && ha->GetNbinsX() == hb->GetNbinsX() && ha->GetNbinsY() == hb->GetNbinsY(); };
    if (!check_same_binning(hSilSeed_nMVTX3nINTT2, hSilSeed) || !check_same_binning(hSilSeed_nMVTX3nINTT1, hSilSeed) || !check_same_binning(hSilSeed_nMVTX2nINTT2, hSilSeed) || !check_same_binning(hSilSeed_nMVTX2nINTT1, hSilSeed) || !check_same_binning(hSilSeed_nMVTXnINTTOther, hSilSeed))
    {
        std::cerr << "[plot_tklCombinatoric] Incompatible binning between silicon-seed categories and inclusive silicon-seed histogram.\n";
        delete h2_nraw;
        fin->Close();
        return ratio_outputs;
    }

    const bool etaAscending = hGrid->GetXaxis()->GetBinCenter(1) <= hGrid->GetXaxis()->GetBinCenter(nEta);
    const bool phiAscending = hGrid->GetYaxis()->GetBinCenter(1) <= hGrid->GetYaxis()->GetBinCenter(nPhi);
    const bool use_custom_xrange = std::isfinite(plot_xmin) && std::isfinite(plot_xmax) && (plot_xmax > plot_xmin);

    // Compute a shared y-axis max across all available signal/background cell histograms.
    double global_ymax = 0.0;
    double ref_bin_size = -1.0;
    for (int ieta = 0; ieta < nEta; ++ieta)
    {
        for (int iphi = 0; iphi < nPhi; ++iphi)
        {
            auto *hSig = GetCellHist(fin, "dR_Sig", prefix, mode, "sig", ieta, iphi, hGrid);
            auto *hBkg = GetCellHist(fin, "dR_Bkg", prefix, mode, "bkg", ieta, iphi, hGrid);

            if (hSig)
            {
                global_ymax = std::max(global_ymax, hSig->GetMaximum());
                if (ref_bin_size < 0.0)
                    ref_bin_size = hSig->GetXaxis()->GetBinWidth(1);
            }
            if (hBkg)
            {
                global_ymax = std::max(global_ymax, hBkg->GetMaximum());
                if (ref_bin_size < 0.0)
                    ref_bin_size = hBkg->GetXaxis()->GetBinWidth(1);
            }
        }
    }
    const double shared_ymax = (global_ymax > 0.0) ? 1.3 * global_ymax : 1.0;
    const std::string yAxisTitle = (ref_bin_size > 0.0) ? Form("Counts / (%.2g radian)", ref_bin_size) : "Counts / (bin size) [radian]";

    auto *c = new TCanvas("c_tklCombinatoric", "c_tklCombinatoric", 5500, 5000);
    const float lM = 0.03f, rM = 0.02f, bM = 0.03f, tM = 0.02f;
    CanvasPartition(c, nEta, nPhi, lM, rM, bM, tM);

    auto *pad00 = dynamic_cast<TPad *>(c->FindObject("pad_0_0"));
    if (!pad00)
    {
        std::cerr << "[plot_tklCombinatoric] pad_0_0 not found.\n";
        delete h2_nraw;
        fin->Close();
        delete c;
        return ratio_outputs;
    }

    std::vector<TH1D *> frames;
    std::vector<TLine *> cutLines;

    int nDrawn = 0;

    auto styleHist = [](TH1D *h, Style_t ms, Color_t col, float markerSize)
    {
        h->SetMarkerStyle(ms);
        h->SetMarkerSize(markerSize);
        h->SetLineColor(col);
        h->SetMarkerColor(col);
        h->Draw("E1 SAME");
    };

    auto createAndDrawCutLine = [&](double x, int linewidth) -> TLine *
    {
        auto *line = new TLine(x, 0.0, x, shared_ymax);
        line->SetLineStyle(3);
        line->SetLineWidth(linewidth);
        line->SetLineColor(kBlue + 2);
        line->Draw("SAME");
        cutLines.push_back(line);
        return line;
    };

    auto SaveSingleCellPlot = [&](bool isdata, TH1D *hSig, TH1D *hBkg, int ieta, int iphi, double etaLow, double etaHigh, double phiLow, double phiHigh, double nraw, std::string str_mode, std::pair<double, double> cutrange)
    {
        if (!hSig && !hBkg)
            return;

        const std::string range_suffix = EtaPhiRangeNameSuffix(hGrid, ieta, iphi);
        const std::string cname = "c_single_" + range_suffix;
        TCanvas csingle(cname.c_str(), cname.c_str(), 800, 700);
        csingle.cd();
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.05);
        // gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.06);
        gPad->SetTicks(1, 1);

        TH1D *frame_single = (TH1D *)(hSig ? hSig : hBkg)->Clone((cname + "_frame").c_str());
        frame_single->SetDirectory(nullptr);
        frame_single->Reset("ICES");
        // frame_single->SetMinimum(0.0);
        frame_single->SetMaximum(shared_ymax);
        if (use_custom_xrange)
            frame_single->GetXaxis()->SetRangeUser(plot_xmin, plot_xmax);
        frame_single->GetXaxis()->SetTitle("INTT doublet |#Delta#phi| [radian]");
        frame_single->GetYaxis()->SetTitle(yAxisTitle.c_str());
        frame_single->GetYaxis()->SetTitleOffset(1.5);
        frame_single->Draw("AXIS");

        if (hSig)
            styleHist(hSig, kFullCircle, kBlack, 1.0f);
        if (hBkg)
            styleHist(hBkg, kOpenCircle, kRed + 1, 1.0f);

        TLegend leg(gPad->GetLeftMargin() + 0.05,     //
                    1. - gPad->GetTopMargin() - 0.22, //
                    gPad->GetLeftMargin() + 0.2,      //
                    1. - gPad->GetTopMargin() - 0.04  //
        );
        leg.SetHeader(Form("%s. #eta=[%s,%s], #phi=[%s,%s]", isdata ? "Data" : "Simulation", to_string_with_precision(etaLow).c_str(), to_string_with_precision(etaHigh).c_str(), to_string_with_precision(phiLow).c_str(), to_string_with_precision(phiHigh).c_str()));
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.035);
        if (hSig)
            leg.AddEntry(hSig, "Unrotated", "pel");
        if (hBkg)
            leg.AddEntry(hBkg, "Rotated", "pel");
        leg.AddEntry((TObject *)0, Form("N^{Uncorr.}_{Doublets}=%.1f (%s=[%s,%s])", nraw, str_mode.c_str(), to_string_with_precision(cutrange.first).c_str(), to_string_with_precision(cutrange.second).c_str()), "");
        leg.Draw();

        createAndDrawCutLine(cutrange.first, 2);
        createAndDrawCutLine(cutrange.second, 2);

        const std::string outbase = single_canvas_dir + "/" + outstem + "_" + range_suffix;
        csingle.SaveAs((outbase + ".pdf").c_str());
        csingle.SaveAs((outbase + ".png").c_str());
        delete frame_single;
    };

    // a lambda to create a string for the mode
    auto modeStringLambda = [&]() -> std::string
    {
        if (mode == "absdPhi")
            return "|#Delta#phi|";
        else if (mode == "dR")
            return "#DeltaR";
        else
            return mode; // as-is if unrecognized, e.g. for custom modes
    };
    std::string modeString = modeStringLambda();
    std::string multiplicity_text = "inclusive";
    const std::string multiplicity_prefix_1 = "Cluster multiplicity percentile ";
    const std::string multiplicity_prefix_2 = "N_{INTT clusters} ";
    if (StartsWith(ratio_outputs.label, multiplicity_prefix_1))
    {
        multiplicity_text = ratio_outputs.label.substr(multiplicity_prefix_1.size());
    }
    else if (StartsWith(ratio_outputs.label, multiplicity_prefix_2))
    {
        multiplicity_text = ratio_outputs.label.substr(multiplicity_prefix_2.size());
    }
    else if (!ratio_outputs.label.empty() && ratio_outputs.label != "Inclusive")
    {
        multiplicity_text = ratio_outputs.label;
    }
    const std::string sample_label = Form("%s, %s#leqZ_{vtx}^{silicon}#leq%s [cm], %s", is_simulation ? "Simulation" : "Data", to_string_with_precision(zvtx_cut_min, 1).c_str(), to_string_with_precision(zvtx_cut_max, 1).c_str(), multiplicity_text.c_str());

    for (int ieta = 0; ieta < nEta; ++ieta)
    {
        for (int iphi = 0; iphi < nPhi; ++iphi)
        {

            c->cd(0);
            const int padEta = etaAscending ? ieta : (nEta - 1 - ieta);
            const int padPhi = phiAscending ? iphi : (nPhi - 1 - iphi);

            auto *pad = dynamic_cast<TPad *>(c->FindObject(TString::Format("pad_%d_%d", padEta, padPhi).Data()));
            if (!pad)
                continue;

            pad->cd();
            pad->SetTicks(1, 1);

            float xFactor = pad00->GetAbsWNDC() / pad->GetAbsWNDC();
            float yFactor = pad00->GetAbsHNDC() / pad->GetAbsHNDC();

            auto *hSig = GetCellHist(fin, "dR_Sig", prefix, mode, "sig", ieta, iphi, hGrid);
            auto *hBkg = GetCellHist(fin, "dR_Bkg", prefix, mode, "bkg", ieta, iphi, hGrid);
            if (!hSig && !hBkg)
                continue;

            // Clone axis template from whichever histogram is available.
            TH1D *baseHist = hSig ? hSig : hBkg;
            auto *frame = (TH1D *)baseHist->Clone(Form("hframe_ieta%d_iphi%d", ieta, iphi));
            frame->SetDirectory(nullptr);
            frame->Reset("ICES");
            frames.push_back(frame); // keep alive until after SaveAs

            // frame->SetMinimum(0.);
            frame->SetMaximum(shared_ymax);
            frame->SetTitle("");
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetYaxis()->SetNdivisions(505);
            if (use_custom_xrange)
            {
                frame->GetXaxis()->SetRangeUser(plot_xmin, plot_xmax);
            }
            frame->GetXaxis()->SetLabelFont(43);
            frame->GetXaxis()->SetLabelSize(34);
            frame->GetYaxis()->SetLabelFont(43);
            frame->GetYaxis()->SetLabelSize(34);
            frame->GetXaxis()->SetTickLength(yFactor * 0.06f / xFactor);
            frame->GetYaxis()->SetTickLength(xFactor * 0.04f / yFactor);
            frame->GetXaxis()->SetTitle("");
            frame->GetYaxis()->SetTitle("");
            frame->Draw("AXIS");

            createAndDrawCutLine(raw_count_xmin, 1);
            createAndDrawCutLine(raw_count_xmax, 1);

            if (hSig)
                styleHist(hSig, kFullCircle, kBlack, 0.0f);
            if (hBkg)
                styleHist(hBkg, kOpenCircle, kRed + 1, 0.0f);

            double nraw = h2_nraw->GetBinContent(ieta + 1, iphi + 1);
            double nraw_err = h2_nraw->GetBinError(ieta + 1, iphi + 1);
            if (!use_direct_raw_count_grid)
            {
                nraw = IntegralInRange(hSig, raw_count_xmin, raw_count_xmax) - IntegralInRange(hBkg, raw_count_xmin, raw_count_xmax);
                nraw_err = (nraw > 0.0) ? std::sqrt(nraw) : 0.0;
                h2_nraw->SetBinContent(ieta + 1, iphi + 1, nraw);
                h2_nraw->SetBinError(ieta + 1, iphi + 1, nraw_err);
            }

            const double etaLow = hGrid->GetXaxis()->GetBinLowEdge(ieta + 1);
            const double etaHigh = hGrid->GetXaxis()->GetBinUpEdge(ieta + 1);
            const double phiLow = hGrid->GetYaxis()->GetBinLowEdge(iphi + 1);
            const double phiHigh = hGrid->GetYaxis()->GetBinUpEdge(iphi + 1);

            TLatex lbl;
            lbl.SetNDC();
            lbl.SetTextFont(43);
            lbl.SetTextSize(28);
            lbl.DrawLatex(XtoPad(0.05), YtoPad(0.90), Form("#eta=[%s,%s]", to_string_with_precision(etaLow).c_str(), to_string_with_precision(etaHigh).c_str()));
            lbl.DrawLatex(XtoPad(0.05), YtoPad(0.81), Form("#phi=[%s,%s]", to_string_with_precision(phiLow).c_str(), to_string_with_precision(phiHigh).c_str()));
            lbl.DrawLatex(XtoPad(0.05), YtoPad(0.70), use_direct_raw_count_grid ? Form("N_{raw}^{stored}=%.1f", nraw) : Form("N_{raw}(%.3f<|#Delta#phi|<%.3f)=%.1f", raw_count_xmin, raw_count_xmax, nraw));

            // SaveSingleCellPlot(!is_simulation, hSig, hBkg, ieta, iphi, etaLow, etaHigh, phiLow, phiHigh, nraw, modeString, std::make_pair(raw_count_xmin, raw_count_xmax));

            pad->Modified();
            ++nDrawn;
        }
    }

    // Global axis titles and header drawn on the main canvas.
    c->cd();
    TLatex tx;
    tx.SetNDC();
    tx.SetTextFont(43);

    tx.SetTextSize(50);
    tx.SetTextAlign(kHAlignCenter + kVAlignTop);
    tx.DrawLatex(0.5, bM / 2., "|#Delta#phi| [radian]");

    tx.SetTextAngle(90);
    tx.SetTextAlign(kHAlignCenter + kVAlignBottom);
    tx.SetTextSize(50);
    tx.DrawLatex(lM / 2., 0.5, yAxisTitle.c_str());

    std::cout << "[plot_tklCombinatoric] Pads with histograms: " << nDrawn << " / " << nEta * nPhi << "\n";

    auto build_and_draw_ratio = [&](const TH2D *hNum,                   //
                                    const TH2D *hDen,                   //
                                    const std::string &ratio_hist_name, //
                                    const std::string &ztitle,          //
                                    const std::string &out_suffix,      //
                                    int round_zmax_mode = -1,           //
                                    const std::string &text_format = "") -> TH2D *
    {
        double ratio_zmax = -1;
        auto *h2_ratio = dynamic_cast<TH2D *>(h2_nraw->Clone(ratio_hist_name.c_str()));
        h2_ratio->SetDirectory(nullptr);
        h2_ratio->Reset("ICES");

        for (int ibx = 1; ibx <= h2_ratio->GetNbinsX(); ++ibx)
        {
            for (int iby = 1; iby <= h2_ratio->GetNbinsY(); ++iby)
            {
                const double num = hNum->GetBinContent(ibx, iby);
                const double den = hDen->GetBinContent(ibx, iby);
                const double ratio = (den > 0.0) ? (num / den) : 0.0;
                if (ratio > ratio_zmax)
                    ratio_zmax = ratio;

                // Counting errors with sqrt(N), propagated for ratio r = num/den.
                double ratio_err = 0.0;
                if (den > 0.0 && num > 0.0)
                {
                    const double sigma_num = std::sqrt(num);
                    const double sigma_den = std::sqrt(den);
                    const double rel_num = sigma_num / num;
                    const double rel_den = sigma_den / den;
                    ratio_err = ratio * std::sqrt(rel_num * rel_num + rel_den * rel_den);
                }
                h2_ratio->SetBinContent(ibx, iby, ratio);
                h2_ratio->SetBinError(ibx, iby, ratio_err);
            }
        }

        // The maximum should be the nearest number of 5, 10, 15, ... above the histogram maximum.
        const bool round_zmax = (round_zmax_mode < 0) ? use_rounded_zmax : (round_zmax_mode != 0);
        double ratio_zmax_plot = ratio_zmax;
        if (round_zmax)
        {
            ratio_zmax_plot = std::ceil(ratio_zmax / 5.0) * 5.0;
            if (ratio_zmax_plot < ratio_zmax)
                ratio_zmax_plot += 5.0;
        }
        const std::string ratio_paint_fmt = text_format.empty() ? ratio_text_format : text_format;

        draw2Dhistogram_adjZaxis(h2_ratio,        //
                                 false,           //
                                 "#eta",          //
                                 "#phi [radian]", //
                                 ztitle,          //
                                 std::pair<double, double>(0, ratio_zmax_plot), {sample_label},
                                 ratio_draw_option, //
                                 ratio_paint_fmt,   //
                                 true,              //
                                 false,             //
                                 plotdir + "/" + outstem + out_suffix);
        return h2_ratio;
    };

    auto draw_input_2d = [&](TH2D *h2,                                    //
                             const std::string &ztitle,                   //
                             const std::string &out_suffix,               //
                             const std::string &drawoption = "colz",      //
                             bool show_bin_value_with_error = false,      //
                             bool show_bin_value_only = false,            //
                             const std::string &xtitle = "#eta",          //
                             const std::string &ytitle = "#phi [radian]", //
                             bool logz = false)
    {
        if (!h2)
            return;

        // The maximum should be the nearest number of 5, 10, 15, ... above the histogram maximum.
        double zmax = h2->GetMaximum();
        double zmax_plot = zmax;
        if (use_rounded_zmax)
        {
            zmax_plot = std::ceil(zmax / 5.0) * 5.0;
            if (zmax_plot < zmax)
                zmax_plot += 5.0;
        }

        draw2Dhistogram_adjZaxis(h2,                                      //
                                 logz,                                    //
                                 xtitle,                                  //
                                 ytitle,                                  //
                                 ztitle,                                  //
                                 std::pair<double, double>(0, zmax_plot), //
                                 {sample_label},                          //
                                 drawoption,                              //
                                 ratio_text_format,                       //
                                 show_bin_value_with_error,               //
                                 show_bin_value_only,                     //
                                 plotdir + "/" + outstem + out_suffix);
    };

    auto make_projection = [&](TH2D *h2, bool project_x, const std::string &hist_name) -> TH1D *
    {
        if (!h2)
            return nullptr;

        TH1D *hproj = project_x ? h2->ProjectionX(hist_name.c_str(), 1, h2->GetNbinsY()) : h2->ProjectionY(hist_name.c_str(), 1, h2->GetNbinsX());
        if (!hproj)
            return nullptr;
        hproj->SetDirectory(nullptr);
        return hproj;
    };

    auto draw_projection = [&](TH2D *h2, bool project_x, const std::string &axis_title, const std::string &ytitle, const std::string &out_suffix)
    {
        if (!h2)
            return;

        TH1D *hproj = make_projection(h2, project_x, std::string(h2->GetName()) + (project_x ? "_projX" : "_projY"));
        if (!hproj)
            return;
        hproj->GetXaxis()->SetTitle(axis_title.c_str());
        hproj->GetYaxis()->SetTitle(ytitle.c_str());
        draw1Dhistogram(hproj, false, axis_title, ytitle, {sample_label}, "hist e1", plotdir + "/" + outstem + out_suffix);
        delete hproj;
    };

    auto build_projection_ratio = [&](std::pair<TH2D *, std::string> hNum2D, std::pair<TH2D *, std::string> hDen2D, bool project_x, const std::string &axis_title, const std::string &ratio_ytitle, const std::string &out_suffix, const std::string &ratio_hist_name) -> TH1D *
    {
        if (!hNum2D.first || !hDen2D.first)
            return nullptr;

        // Projection, do NOT include overflow and underflow bins for now
        TH1D *hNumProj = project_x ? hNum2D.first->ProjectionX((ratio_hist_name + "_num").c_str(), 1, hNum2D.first->GetNbinsY()) : hNum2D.first->ProjectionY((ratio_hist_name + "_num").c_str(), 1, hNum2D.first->GetNbinsX());
        TH1D *hDenProj = project_x ? hDen2D.first->ProjectionX((ratio_hist_name + "_den").c_str(), 1, hDen2D.first->GetNbinsY()) : hDen2D.first->ProjectionY((ratio_hist_name + "_den").c_str(), 1, hDen2D.first->GetNbinsX());
        if (!hNumProj || !hDenProj)
        {
            delete hNumProj;
            delete hDenProj;
            return nullptr;
        }

        hNumProj->SetDirectory(nullptr);
        hDenProj->SetDirectory(nullptr);
        drawComparisonWithRatio(std::make_pair(hDenProj, hDen2D.second), std::make_pair(hNumProj, hNum2D.second), false, axis_title, "Counts", ratio_ytitle, {sample_label}, plotdir + "/" + outstem + out_suffix + "_components");

        TH1D *hRatio = dynamic_cast<TH1D *>(hNumProj->Clone(ratio_hist_name.c_str()));
        if (hRatio)
        {
            hRatio->SetDirectory(nullptr);
            hRatio->Reset("ICES");
            hRatio->Divide(hNumProj, hDenProj, 1.0, 1.0);
            draw1Dhistogram(hRatio, false, axis_title, ratio_ytitle, {sample_label}, "hist e1", plotdir + "/" + outstem + out_suffix);
        }

        delete hNumProj;
        delete hDenProj;
        return hRatio;
    };

    auto draw_projection_stack = [&](TH2D *hInclusive, const std::vector<std::pair<TH2D *, std::string>> &components, bool project_x, const std::string &axis_title, const std::string &out_suffix)
    {
        if (!hInclusive)
            return;

        TH1D *hInclusiveProj = make_projection(hInclusive, project_x, std::string(hInclusive->GetName()) + (project_x ? "_projX_stack_inclusive" : "_projY_stack_inclusive"));
        if (!hInclusiveProj)
            return;

        const float alpha = 0.5f;
        const std::vector<std::string> fill_colors = {"#ffa90e", "#6CA651", "#832db6", "#bd1f01", "#3f90da"};

        TCanvas *cstack = new TCanvas((std::string("c") + out_suffix).c_str(), (std::string("c") + out_suffix).c_str(), 800, 700);
        cstack->cd();
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.08);

        THStack *hs = new THStack((std::string(hInclusive->GetName()) + (project_x ? "_projX_stack" : "_projY_stack")).c_str(), (std::string(";") + axis_title + ";Counts").c_str());
        std::vector<TH1D *> component_proj;
        std::vector<std::string> component_labels;
        component_proj.reserve(components.size());
        component_labels.reserve(components.size());

        for (size_t i = 0; i < components.size(); ++i)
        {
            TH2D *h2 = components[i].first;
            if (!h2)
                continue;

            TH1D *hproj = make_projection(h2, project_x, std::string(h2->GetName()) + (project_x ? "_projX_stack_component" : "_projY_stack_component"));
            if (!hproj)
                continue;

            const Color_t color = TColor::GetColor(fill_colors[i % fill_colors.size()].c_str());
            hproj->SetFillColorAlpha(color, alpha);
            hproj->SetLineColor(color);
            hproj->SetLineWidth(1);
            hs->Add(hproj);
            component_proj.push_back(hproj);
            component_labels.push_back(components[i].second);
        }

        hs->Draw("hist");
        hs->GetXaxis()->SetTitle(axis_title.c_str());
        hs->GetYaxis()->SetTitle("Counts");
        hs->GetYaxis()->SetTitleOffset(1.4);
        hs->SetMaximum(hInclusiveProj->GetMaximum() * 1.5);
        gPad->Modified();
        gPad->Update();

        hInclusiveProj->SetLineColor(kBlack);
        hInclusiveProj->SetLineWidth(2);
        hInclusiveProj->Draw("hist SAME");

        TLegend *leg = new TLegend(gPad->GetLeftMargin() + 0.02, 1 - gPad->GetTopMargin() - 0.24, gPad->GetLeftMargin() + 0.34, 1 - gPad->GetTopMargin() - 0.04);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(hInclusiveProj, "Inclusive", "l");
        for (size_t i = 0; i < component_proj.size(); ++i)
            leg->AddEntry(component_proj[i], component_labels[i].c_str(), "f");
        leg->Draw();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextSize(0.035);
        latex.SetTextAlign(kHAlignRight + kVAlignBottom);
        latex.DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, sample_label.c_str());

        cstack->SaveAs((plotdir + "/" + outstem + out_suffix + ".png").c_str());
        cstack->SaveAs((plotdir + "/" + outstem + out_suffix + ".pdf").c_str());

        delete leg;
        for (auto *h : component_proj)
            delete h;
        delete hs;
        delete hInclusiveProj;
        delete cstack;
    };

    auto draw_projection_stack_fractional = [&](TH2D *hInclusive, const std::vector<std::pair<TH2D *, std::string>> &components, bool project_x, const std::string &axis_title, const std::string &out_suffix)
    {
        if (!hInclusive)
            return;

        TH1D *hInclusiveProj = make_projection(hInclusive, project_x, std::string(hInclusive->GetName()) + (project_x ? "_projX_frac_inclusive" : "_projY_frac_inclusive"));
        if (!hInclusiveProj)
            return;

        const float alpha = 0.5f;
        const std::vector<std::string> fill_colors = {"#ffa90e", "#6CA651", "#832db6", "#bd1f01", "#3f90da"};

        TCanvas *cFrac = new TCanvas((std::string("c") + out_suffix).c_str(), (std::string("c") + out_suffix).c_str(), 800, 700);
        cFrac->cd();
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.23);

        THStack *hsFrac = new THStack((std::string(hInclusive->GetName()) + (project_x ? "_projX_stack_fractional" : "_projY_stack_fractional")).c_str(), (std::string(";") + axis_title + ";Fractional contribution").c_str());
        std::vector<TH1D *> component_frac;
        std::vector<std::string> component_labels;
        component_frac.reserve(components.size());
        component_labels.reserve(components.size());

        for (size_t i = 0; i < components.size(); ++i)
        {
            TH2D *h2 = components[i].first;
            if (!h2)
                continue;

            TH1D *hproj = make_projection(h2, project_x, std::string(h2->GetName()) + (project_x ? "_projX_frac_component" : "_projY_frac_component"));
            if (!hproj)
                continue;

            TH1D *hFrac = MakeFractionalHist(hproj, hInclusiveProj, std::string(hproj->GetName()) + "_fractional");
            delete hproj;
            if (!hFrac)
                continue;

            const Color_t color = TColor::GetColor(fill_colors[i % fill_colors.size()].c_str());
            hFrac->SetFillColorAlpha(color, alpha);
            hFrac->SetLineColor(color);
            hFrac->SetLineWidth(2);
            hsFrac->Add(hFrac);
            component_frac.push_back(hFrac);
            component_labels.push_back(components[i].second);
        }

        hsFrac->Draw("hist");
        hsFrac->GetXaxis()->SetTitle(axis_title.c_str());
        hsFrac->GetYaxis()->SetTitle("Fractional contribution");
        hsFrac->GetYaxis()->SetTitleOffset(1.5);
        hsFrac->SetMaximum(1.0);
        gPad->Modified();
        gPad->Update();

        TLine *line1 = new TLine(hInclusiveProj->GetXaxis()->GetXmin(), 1.0, hInclusiveProj->GetXaxis()->GetXmax(), 1.0);
        line1->SetLineColor(kBlack);
        line1->SetLineWidth(2);

        TLegend *legFrac = new TLegend(gPad->GetLeftMargin() + 0.02, 1 - gPad->GetTopMargin() + 0.01, gPad->GetLeftMargin() + 0.34, 0.99);
        legFrac->SetBorderSize(0);
        legFrac->SetFillStyle(0);
        legFrac->SetTextSize(0.035);
        legFrac->AddEntry(line1, "Inclusive", "l");
        for (size_t i = 0; i < component_frac.size(); ++i)
            legFrac->AddEntry(component_frac[i], component_labels[i].c_str(), "f");
        legFrac->Draw();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextSize(0.03);
        latex.SetTextAlign(kHAlignRight + kVAlignBottom);
        latex.DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, sample_label.c_str());

        cFrac->SaveAs((plotdir + "/" + outstem + out_suffix + ".png").c_str());
        cFrac->SaveAs((plotdir + "/" + outstem + out_suffix + ".pdf").c_str());

        delete legFrac;
        delete line1;
        for (auto *h : component_frac)
            delete h;
        delete hsFrac;
        delete hInclusiveProj;
        delete cFrac;
    };

    draw_input_2d(h2_nraw, "N_{Uncorrected}^{INTT Doublets}", "_h2_nraw");
    draw_input_2d(hSilSeed, "N_{Silicon seeds}", "_hSilSeed", "colz text", false, true);
    draw_input_2d(hSilSeed_nMVTX3nINTT2, "N_{Silicon seeds}^{nMVTX3nINTT2}", "_hSilSeed_nMVTX3nINTT2", "colz text", false, true);
    draw_input_2d(hSilSeed_nMVTX3nINTT1, "N_{Silicon seeds}^{nMVTX3nINTT1}", "_hSilSeed_nMVTX3nINTT1", "colz text", false, true);
    draw_input_2d(hSilSeed_nMVTX2nINTT2, "N_{Silicon seeds}^{nMVTX2nINTT2}", "_hSilSeed_nMVTX2nINTT2", "colz text", false, true);
    draw_input_2d(hSilSeed_nMVTX2nINTT1, "N_{Silicon seeds}^{nMVTX2nINTT1}", "_hSilSeed_nMVTX2nINTT1", "colz text", false, true);
    draw_input_2d(hSilSeed_nMVTXnINTTOther, "N_{Silicon seeds}^{nMVTXnINTTOther}", "_hSilSeed_nMVTXnINTTOther", "colz text", false, true);
    draw_input_2d(hSilSeed_good, "N_{Silicon seeds}^{good}", "_hSilSeed_good", "colz text", false, true);
    draw_input_2d(hINTTClusterVsAssociatedSilSeed, "Entries", "_hINTTClusterVsAssociatedSilSeed", "colz", false, false, "N_{INTT clusters} (selected crossing)", "N_{Silicon seeds}^{assoc. w. vertex} (selected crossing)", true);
    draw_input_2d(hINTTClusterVsAllSilSeed, "Entries", "_hINTTClusterVsAllSilSeed", "colz", false, false, "N_{INTT clusters} (selected crossing)", "N_{Silicon seeds}^{all in crossing} (selected crossing)", true);
    drawComparisonWithRatio(std::make_pair(hSeedCrossing_selectedCrossings, "Selected crossings"),                                                                        //
                            std::make_pair(hSeedCrossing_selectedCrossings_largeDiffNClusterSeeds, "Selected crossings with N_{INTT clusters}>4#timesN_{Silicon seeds}"), //
                            true,                                                                                                                                         //
                            "Crossing",                                                                                                                                   //
                            "Entries",                                                                                                                                    //
                            "Ratio",                                                                                                                                      //
                            {sample_label},                                                                                                                               //
                            plotdir + "/" + outstem + "_hSeedCrossing_selectedCrossings_comparison"                                                                       //
    );
    draw_input_2d(hSeed, "N_{Cluster}^{Layer 3+4}", "_hSeed"); // this is cluster in inner layer, not "silicon seeds"
    if (is_simulation)
        draw_input_2d(hPrimaryCharged, "N_{sPHENIX primary charged}", "_hPrimaryCharged");

    draw_projection(hSilSeed, false, "#phi [radian]", "Counts", "_hSilSeed_phiProjection");
    draw_projection(hSilSeed, true, "#eta", "Counts", "_hSilSeed_etaProjection");
    draw_projection(hSilSeed_nMVTX3nINTT2, false, "#phi [radian]", "Counts", "_hSilSeed_nMVTX3nINTT2_phiProjection");
    draw_projection(hSilSeed_nMVTX3nINTT2, true, "#eta", "Counts", "_hSilSeed_nMVTX3nINTT2_etaProjection");
    draw_projection(hSilSeed_nMVTX3nINTT1, false, "#phi [radian]", "Counts", "_hSilSeed_nMVTX3nINTT1_phiProjection");
    draw_projection(hSilSeed_nMVTX3nINTT1, true, "#eta", "Counts", "_hSilSeed_nMVTX3nINTT1_etaProjection");
    draw_projection(hSilSeed_nMVTX2nINTT2, false, "#phi [radian]", "Counts", "_hSilSeed_nMVTX2nINTT2_phiProjection");
    draw_projection(hSilSeed_nMVTX2nINTT2, true, "#eta", "Counts", "_hSilSeed_nMVTX2nINTT2_etaProjection");
    draw_projection(hSilSeed_nMVTX2nINTT1, false, "#phi [radian]", "Counts", "_hSilSeed_nMVTX2nINTT1_phiProjection");
    draw_projection(hSilSeed_nMVTX2nINTT1, true, "#eta", "Counts", "_hSilSeed_nMVTX2nINTT1_etaProjection");
    draw_projection(hSilSeed_nMVTXnINTTOther, false, "#phi [radian]", "Counts", "_hSilSeed_nMVTXnINTTOther_phiProjection");
    draw_projection(hSilSeed_nMVTXnINTTOther, true, "#eta", "Counts", "_hSilSeed_nMVTXnINTTOther_etaProjection");
    draw_projection_stack(hSilSeed, //
                          {
                              //
                              {hSilSeed_nMVTX3nINTT2, "MVTX+INTT: 3+2"},      //
                              {hSilSeed_nMVTX3nINTT1, "MVTX+INTT: 3+1"},      //
                              {hSilSeed_nMVTX2nINTT2, "MVTX+INTT: 2+2"},      //
                              {hSilSeed_nMVTX2nINTT1, "MVTX+INTT: 2+1"},      //
                              {hSilSeed_nMVTXnINTTOther, "MVTX+INTT: Other"}, //
                          },                                                  //
                          false,                                              //
                          "#phi [radian]",                                   //
                          "_hSilSeed_phiProjection_stack"                    //
    );
    draw_projection_stack(hSilSeed, //
                          {
                              //
                              {hSilSeed_nMVTX3nINTT2, "MVTX+INTT: 3+2"},      //
                              {hSilSeed_nMVTX3nINTT1, "MVTX+INTT: 3+1"},      //
                              {hSilSeed_nMVTX2nINTT2, "MVTX+INTT: 2+2"},      //
                              {hSilSeed_nMVTX2nINTT1, "MVTX+INTT: 2+1"},      //
                              {hSilSeed_nMVTXnINTTOther, "MVTX+INTT: Other"}, //
                          },                                                  //
                          true,                                               //
                          "#eta",                                            //
                          "_hSilSeed_etaProjection_stack"                    //
    );
    draw_projection_stack_fractional(hSilSeed, //
                                     {
                                         //
                                         {hSilSeed_nMVTX3nINTT2, "MVTX+INTT: 3+2"},      //
                                         {hSilSeed_nMVTX3nINTT1, "MVTX+INTT: 3+1"},      //
                                         {hSilSeed_nMVTX2nINTT2, "MVTX+INTT: 2+2"},      //
                                         {hSilSeed_nMVTX2nINTT1, "MVTX+INTT: 2+1"},      //
                                         {hSilSeed_nMVTXnINTTOther, "MVTX+INTT: Other"}, //
                                     },                                                  //
                                     false,                                              //
                                     "#phi [radian]",                                   //
                                     "_hSilSeed_phiProjection_stack_fractional"         //
    );
    draw_projection_stack_fractional(hSilSeed, //
                                     {
                                         //
                                         {hSilSeed_nMVTX3nINTT2, "MVTX+INTT: 3+2"},      //
                                         {hSilSeed_nMVTX3nINTT1, "MVTX+INTT: 3+1"},      //
                                         {hSilSeed_nMVTX2nINTT2, "MVTX+INTT: 2+2"},      //
                                         {hSilSeed_nMVTX2nINTT1, "MVTX+INTT: 2+1"},      //
                                         {hSilSeed_nMVTXnINTTOther, "MVTX+INTT: Other"}, //
                                     },                                                  //
                                     true,                                               //
                                     "#eta",                                            //
                                     "_hSilSeed_etaProjection_stack_fractional"         //
    );

    if (hDCA3D_Sig && hDCA3D_Bkg)
    {
        TH1D *hDCA3D_Sub = dynamic_cast<TH1D *>(hDCA3D_Sig->Clone((prefix + "_DCA3D_Sub").c_str()));
        if (hDCA3D_Sub && hSilSeedDCA3D)
        {
            hDCA3D_Sub->SetDirectory(nullptr);
            hDCA3D_Sub->Add(hDCA3D_Bkg, -1.0);

            double dca_ymin = std::min({hDCA3D_Sig->GetMinimum(0), hDCA3D_Bkg->GetMinimum(0), hDCA3D_Sub->GetMinimum(0), hSilSeedDCA3D->GetMinimum(0)});
            double dca_ymax = std::max({hDCA3D_Sig->GetMaximum(), hDCA3D_Bkg->GetMaximum(), hDCA3D_Sub->GetMaximum(), hSilSeedDCA3D->GetMaximum()});
            if (dca_ymax <= dca_ymin)
                dca_ymax = dca_ymin + 1.0;
            const double dca_margin = 0.15 * (dca_ymax - dca_ymin);

            draw_comparison_alt({hDCA3D_Sig, hDCA3D_Bkg, hDCA3D_Sub},                                                                //
                                {"#1f77b4", "#67B2D8", "#000000"},                                                                       //
                                {24, 25, 20},                                                                                                   //
                                {0.5, 0.5, 1},                                                                                                   //
                                {"INTT doublets - Unrotated", "INTT doublets - Rotated", "INTT doublets - Background subtracted"}, //
                                "DCA_{w.r.t vertex} [cm]",                                                                                          //
                                "Counts",                                                                                                           //
                                {0, std::max({hDCA3D_Sig->GetMaximum(), hDCA3D_Bkg->GetMaximum(), hDCA3D_Sub->GetMaximum()})*1.3},                                                                                                           //
                                false,                                                                                                               //
                                false,                                                                                                              //
                                510,                                                                                                                //
                                {false, 0.f},                                                                                                       //
                                plotdir + "/" + outstem + "_hDCA3D_comparison"                                                              //
            );

            draw_comparison_alt({hDCA3D_Sig, hDCA3D_Bkg, hDCA3D_Sub, hSilSeedDCA3D},                                                                //
                                {"#1f77b4", "#67B2D8", "#000000", "#bd1f01"},                                                                       //
                                {24, 25, 20, 21},                                                                                                   //
                                {0.5, 0.5, 1, 1},                                                                                                   //
                                {"INTT doublets - Unrotated", "INTT doublets - Rotated", "INTT doublets - Background subtracted", "Silicon seeds"}, //
                                "DCA_{w.r.t vertex} [cm]",                                                                                          //
                                "Counts",                                                                                                           //
                                {-1, -1},                                                                                                           //
                                true,                                                                                                               //
                                false,                                                                                                              //
                                510,                                                                                                                //
                                {false, 0.f},                                                                                                       //
                                plotdir + "/" + outstem + "_hDCA3D_SilSeed_comparison"                                                              //
            );
            // delete hDCA3D_Sub;
        }
    }

    if (hSilSeedDCA3D)
        draw1Dhistogram(hSilSeedDCA3D, true, "DCA_{Silicon seeds}^{w.r.t vertex} [cm]", "Counts", {sample_label}, "hist e1", plotdir + "/" + outstem + "_hSilSeedDCA3D");

    if (hSilSeedDCA3D_truthVtxZ_2D)
    {
        const double heatmap_zmax = std::max(0.05, hSilSeedDCA3D_truthVtxZ_2D->GetMaximum() * 1.1);
        hSilSeedDCA3D_truthVtxZ_2D->GetYaxis()->SetRangeUser(0, 3);
        draw2Dhistogram_adjZaxis(hSilSeedDCA3D_truthVtxZ_2D,
                                 true,
                                 "Truth vertex z [cm]",
                                 "DCA_{Silicon seeds}^{w.r.t vertex} [cm]",
                                 "Normalized entries",
                                 std::pair<double, double>(0.0, heatmap_zmax),
                                 {sample_label},
                                 "colz",
                                 ".3g",
                                 false,
                                 false,
                                 plotdir + "/" + outstem + "_hSilSeedDCA3D_truthVtxZ_2D");
    }

    if (!hSilSeedDCA3D_truthVtxZ_absMerged.empty())
    {
        std::vector<TH1 *> hists;
        std::vector<std::string> labels;
        hists.reserve(hSilSeedDCA3D_truthVtxZ_absMerged.size());
        labels.reserve(hSilSeedDCA3D_truthVtxZ_absMerged.size());
        for (const auto &entry : hSilSeedDCA3D_truthVtxZ_absMerged)
        {
            hists.push_back(entry.first);
            labels.push_back(entry.second);
        }

        draw_comparison(hists,
                        {"#4477aa", "#228833", "#ccbb44", "#ee8866", "#cc3311"},
                        labels,
                        "DCA_{Silicon seeds}^{w.r.t vertex} [cm]",
                        "Normalized counts",
                        {-1, -1},
                        true,
                        true,
                        510,
                        std::make_pair(false, 0.0f),
                        plotdir + "/" + outstem + "_hSilSeedDCA3D_truthVtxZAbs_comparison");
    }

    // Always produce nraw / silicon-seed ratio.
    TH2D *h2_ratio_nraw_over_silseed = build_and_draw_ratio(h2_nraw, hSilSeed,                                     //
                                                            prefix + "_NRawOverSilSeedEtaPhi",                     //
                                                            "N_{Uncorrected}^{INTT Doublets} / N_{Silicon seeds}", //
                                                            "_ratio_nraw_over_silseed");
    TH2D *h2_ratio_nraw_over_silseed_good = build_and_draw_ratio(h2_nraw, hSilSeed_good,                                      //
                                                                 prefix + "_NRawOverSilSeedGoodEtaPhi",                             //
                                                                 "N_{Uncorrected}^{INTT Doublets} / N_{Silicon seeds}^{good}",     //
                                                                 "_ratio_nraw_over_silseed_good");

    ratio_outputs.eta_ratio = build_projection_ratio(std::make_pair(h2_nraw, "N_{Uncorrected}^{INTT Doublets}"), //
                                                     std::make_pair(hSilSeed, "N_{Silicon seeds}"),              //
                                                     true,                                                       //
                                                     "#eta",                                                     //
                                                     "N_{Uncorrected}^{INTT Doublets} / N_{Silicon seeds}",      //
                                                     "_ratio1D_nraw_over_silseed_eta",                           //
                                                     prefix + "_NRawOverSilSeedEtaProjection"                    //
    );
    ratio_outputs.phi_ratio = build_projection_ratio(std::make_pair(h2_nraw, "N_{Uncorrected}^{INTT Doublets}"), //
                                                     std::make_pair(hSilSeed, "N_{Silicon seeds}"),              //
                                                     false,                                                      //
                                                     "#phi [radian]",                                            //
                                                     "N_{Uncorrected}^{INTT Doublets} / N_{Silicon seeds}",      //
                                                     "_ratio1D_nraw_over_silseed_phi",                           //
                                                     prefix + "_NRawOverSilSeedPhiProjection"                    //
    );
    ratio_outputs.eta_ratio_good = build_projection_ratio(std::make_pair(h2_nraw, "N_{Uncorrected}^{INTT Doublets}"), //
                                                          std::make_pair(hSilSeed_good, "N_{Silicon seeds}^{good}"),   //
                                                          true,                                                        //
                                                          "#eta",                                                      //
                                                          "N_{Uncorrected}^{INTT Doublets} / N_{Silicon seeds}^{good}", //
                                                          "_ratio1D_nraw_over_silseed_good_eta",                       //
                                                          prefix + "_NRawOverSilSeedGoodEtaProjection"                 //
    );
    ratio_outputs.phi_ratio_good = build_projection_ratio(std::make_pair(h2_nraw, "N_{Uncorrected}^{INTT Doublets}"), //
                                                          std::make_pair(hSilSeed_good, "N_{Silicon seeds}^{good}"),   //
                                                          false,                                                       //
                                                          "#phi [radian]",                                             //
                                                          "N_{Uncorrected}^{INTT Doublets} / N_{Silicon seeds}^{good}", //
                                                          "_ratio1D_nraw_over_silseed_good_phi",                       //
                                                          prefix + "_NRawOverSilSeedGoodPhiProjection"                 //
    );


    // In simulation, add extra ratios: primary charged / nraw and primary charged / silicon-seed.
    TH2D *h2_ratio_primary_over_nraw = nullptr;
    TH2D *h2_ratio_primary_over_silseed = nullptr;
    if (is_simulation)
    {
        h2_ratio_primary_over_nraw = build_and_draw_ratio(hPrimaryCharged, h2_nraw,                                //
                                                          prefix + "_PrimaryChargedOverNRawEtaPhi",                //
                                                          "N_{sPHENIX primary charged} / N_{raw}^{INTT Doublets}", //
                                                          "_ratio_primarycharged_over_nraw");
        h2_ratio_primary_over_silseed = build_and_draw_ratio(hPrimaryCharged, hSilSeed,                         //
                                                             prefix + "_PrimaryChargedOverSilSeedEtaPhi",       //
                                                             "N_{sPHENIX primary charged} / N_{Silicon seeds}", //
                                                             "_ratio_primarycharged_over_silseed");
    }

    TH2D *h2_ratio_silseed_nMVTX3nINTT2_over_inclusive = build_and_draw_ratio(hSilSeed_nMVTX3nINTT2, hSilSeed,                        //
                                                                              prefix + "_SilSeed_nMVTX3nINTT2_OverInclusive",         //
                                                                              "N_{SilSeed}^{nMVTX3nINTT2} / N_{SilSeed}^{Inclusive}", //
                                                                              "_ratio_silseed_nMVTX3nINTT2_over_inclusive",           //
                                                                              false,                                                  //
                                                                              ".2g");
    TH2D *h2_ratio_silseed_nMVTX3nINTT1_over_inclusive = build_and_draw_ratio(hSilSeed_nMVTX3nINTT1, hSilSeed,                        //
                                                                              prefix + "_SilSeed_nMVTX3nINTT1_OverInclusive",         //
                                                                              "N_{SilSeed}^{nMVTX3nINTT1} / N_{SilSeed}^{Inclusive}", //
                                                                              "_ratio_silseed_nMVTX3nINTT1_over_inclusive",           //
                                                                              false,                                                  //
                                                                              ".2g");
    TH2D *h2_ratio_silseed_nMVTX2nINTT2_over_inclusive = build_and_draw_ratio(hSilSeed_nMVTX2nINTT2, hSilSeed,                        //
                                                                              prefix + "_SilSeed_nMVTX2nINTT2_OverInclusive",         //
                                                                              "N_{SilSeed}^{nMVTX2nINTT2} / N_{SilSeed}^{Inclusive}", //
                                                                              "_ratio_silseed_nMVTX2nINTT2_over_inclusive",           //
                                                                              false,                                                  //
                                                                              ".2g");
    TH2D *h2_ratio_silseed_nMVTX2nINTT1_over_inclusive = build_and_draw_ratio(hSilSeed_nMVTX2nINTT1, hSilSeed,                        //
                                                                              prefix + "_SilSeed_nMVTX2nINTT1_OverInclusive",         //
                                                                              "N_{SilSeed}^{nMVTX2nINTT1} / N_{SilSeed}^{Inclusive}", //
                                                                              "_ratio_silseed_nMVTX2nINTT1_over_inclusive",           //
                                                                              false,                                                  //
                                                                              ".2g");
    TH2D *h2_ratio_silseed_nMVTXnINTTOther_over_inclusive = build_and_draw_ratio(hSilSeed_nMVTXnINTTOther, hSilSeed,                        //
                                                                                 prefix + "_SilSeed_nMVTXnINTTOther_OverInclusive",         //
                                                                                 "N_{SilSeed}^{nMVTXnINTTOther} / N_{SilSeed}^{Inclusive}", //
                                                                                 "_ratio_silseed_nMVTXnINTTOther_over_inclusive",           //
                                                                                 false,                                                     //
                                                                                 ".2g");

    c->SaveAs((plotdir + "/" + outstem + ".pdf").c_str());
    // c->SaveAs((plotdir + "/" + outstem + ".png").c_str()); // somehow this take a really long time, skip for now

    // Safe to clean up now that the canvas has been saved.
    for (auto *f : frames)
        delete f;
    for (auto *ln : cutLines)
        delete ln;
    delete h2_ratio_nraw_over_silseed;
    delete h2_ratio_nraw_over_silseed_good;
    delete h2_ratio_silseed_nMVTX3nINTT2_over_inclusive;
    delete h2_ratio_silseed_nMVTX3nINTT1_over_inclusive;
    delete h2_ratio_silseed_nMVTX2nINTT2_over_inclusive;
    delete h2_ratio_silseed_nMVTX2nINTT1_over_inclusive;
    delete h2_ratio_silseed_nMVTXnINTTOther_over_inclusive;
    if (h2_ratio_primary_over_nraw)
        delete h2_ratio_primary_over_nraw;
    if (h2_ratio_primary_over_silseed)
        delete h2_ratio_primary_over_silseed;
    delete hSilSeedDCA3D_truthVtxZ_2D;
    for (auto &entry : hSilSeedDCA3D_truthVtxZ_absMerged)
        delete entry.first;
    delete hSilSeed_good;
    delete h2_nraw;
    fin->Close();
    delete c;
    return ratio_outputs;
}

// void plot_tklCombinatoric(const std::string &inputfile = "test_OO_combinatoric.root",     //
//                           const std::string &prefix = "tkl_Combinatoric",                 //
//                           const std::string &mode = "absdPhi",                            //
//                           bool is_simulation = false,                                     //
//                           double zvtx_cut_min = -10.0,                                    //
//                           double zvtx_cut_max = 10.0,                                     //
//                           double raw_count_xmin = 0.0,                                    //
//                           double raw_count_xmax = kDirectRawCountCut,                     //
//                           double plot_xmin = 0.0,                                         //
//                           double plot_xmax = 0.2,                                         //
//                           const std::string &outstem = "plot_tklCombinatoric",            //
//                           const std::string &plotdir = "./figure/figure-tklCombinatoric", //
//                           bool zero_x_error = true,                                       //
//                           const std::string &ratio_draw_option = "colz text",             //
//                           const std::string &ratio_text_format = ".3f",                   //
//                           bool use_rounded_zmax = true                                    //
// )
void plot_tklCombinatoric(const std::string &inputfile = "/sphenix/tg/tg01/hf/hjheng/ppg-dNdEta-OOpp/TrackletAna-HistOutput/simulation/ACTS/cotThetaMax2p9/merged_hist.root", //
                          const std::string &prefix = "tkl_Combinatoric",                                                                                      //
                          const std::string &mode = "absdPhi",                                                                                                 //
                          bool is_simulation = true,                                                                                                           //
                          double zvtx_cut_min = -10.0,                                                                                                         //
                          double zvtx_cut_max = 10.0,                                                                                                          //
                          double raw_count_xmin = 0.0,                                                                                                         //
                          double raw_count_xmax = kDirectRawCountCut,                                                                                          //
                          double plot_xmin = 0.0,                                                                                                              //
                          double plot_xmax = 0.2,                                                                                                              //
                          const std::string &outstem = "plot_tklCombinatoric",                                                                                 //
                          const std::string &plotdir = "./figure/figure-tklCombinatoric-simulation-cotThetaMax2p9",                                                           //
                          bool zero_x_error = true,                                                                                                            //
                          const std::string &ratio_draw_option = "colz text",                                                                                  //
                          const std::string &ratio_text_format = ".3f",                                                                                        //
                          bool use_rounded_zmax = true                                                                                                         //
)
{
    TFile *fin = TFile::Open(inputfile.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "[plot_tklCombinatoric] Cannot open: " << inputfile << "\n";
        delete fin;
        return;
    }

    const std::vector<std::string> prefixes = CollectAvailablePrefixes(fin, prefix);
    const auto percentile_boundaries = LoadPercentileBoundariesByBin(is_simulation);
    fin->Close();
    delete fin;

    const std::vector<std::string> colors = {"#4477aa", "#228833", "#ccbb44", "#ee8866", "#cc3311", "#aa3377"};
    std::vector<ProjectionRatioOutputs> ratio_sets;
    ratio_sets.reserve(prefixes.size());

    for (const auto &current_prefix : prefixes)
    {
        const std::string suffix = PrefixSuffixForOutput(prefix, current_prefix);
        ProjectionRatioOutputs outputs = plot_tklCombinatoric_single(inputfile, current_prefix, mode, is_simulation, zvtx_cut_min, zvtx_cut_max, raw_count_xmin, raw_count_xmax, plot_xmin, plot_xmax, outstem + suffix, plotdir, zero_x_error, ratio_draw_option, ratio_text_format, use_rounded_zmax, PrefixDisplayLabel(prefix, current_prefix, percentile_boundaries));
        if (outputs.eta_ratio || outputs.phi_ratio || outputs.eta_ratio_good || outputs.phi_ratio_good)
            ratio_sets.push_back(outputs);
    }

    std::vector<TH1 *> eta_hists;
    std::vector<TH1 *> phi_hists;
    std::vector<TH1 *> eta_hists_good;
    std::vector<TH1 *> phi_hists_good;
    std::vector<TH1 *> eta_hists_inverse;
    std::vector<TH1 *> phi_hists_inverse;
    std::vector<TH1 *> eta_hists_good_inverse;
    std::vector<TH1 *> phi_hists_good_inverse;
    std::vector<std::string> labels;
    std::vector<std::string> labels_good;
    std::vector<std::string> comparison_colors;
    std::vector<std::string> comparison_colors_good;
    for (size_t i = 0; i < ratio_sets.size(); ++i)
    {
        if (!ratio_sets[i].eta_ratio || !ratio_sets[i].phi_ratio)
            continue;
        eta_hists.push_back(ratio_sets[i].eta_ratio);
        phi_hists.push_back(ratio_sets[i].phi_ratio);
        if (TH1D *hInv = MakeInverseHist(ratio_sets[i].eta_ratio, ratio_sets[i].prefix + "_NRawOverSilSeedEtaProjection_inverse"))
            eta_hists_inverse.push_back(hInv);
        if (TH1D *hInv = MakeInverseHist(ratio_sets[i].phi_ratio, ratio_sets[i].prefix + "_NRawOverSilSeedPhiProjection_inverse"))
            phi_hists_inverse.push_back(hInv);
        if (ratio_sets[i].eta_ratio_good && ratio_sets[i].phi_ratio_good)
        {
            eta_hists_good.push_back(ratio_sets[i].eta_ratio_good);
            phi_hists_good.push_back(ratio_sets[i].phi_ratio_good);
            if (TH1D *hInv = MakeInverseHist(ratio_sets[i].eta_ratio_good, ratio_sets[i].prefix + "_NRawOverSilSeedGoodEtaProjection_inverse"))
                eta_hists_good_inverse.push_back(hInv);
            if (TH1D *hInv = MakeInverseHist(ratio_sets[i].phi_ratio_good, ratio_sets[i].prefix + "_NRawOverSilSeedGoodPhiProjection_inverse"))
                phi_hists_good_inverse.push_back(hInv);
            labels_good.push_back(ratio_sets[i].label);
            comparison_colors_good.push_back((i == 0) ? "#000000" : colors[(i - 1) % colors.size()]);
        }
        labels.push_back(ratio_sets[i].label);
        comparison_colors.push_back((i == 0) ? "#000000" : colors[(i - 1) % colors.size()]);
    }

    if (!eta_hists.empty())
    {
        draw_comparison(eta_hists,                                                         //
                        comparison_colors,                                                 //
                        labels,                                                            //
                        "#eta",                                                            //
                        "N_{Uncorrected}^{INTT Doublets}/N_{Silicon seeds}",               //
                        {0, (!is_simulation) ? 2 : 3},                                     //
                        false,                                                             //
                        false,                                                             //
                        510,                                                               //
                        std::make_pair(true, 1.0),                                         //
                        plotdir + "/" + outstem + "_ratio1D_nraw_over_seed_eta_comparison" //
        );

        draw_comparison(eta_hists_inverse,                                                  //
                        comparison_colors,                                                  //
                        labels,                                                             //
                        "#eta",                                                             //
                        "N_{Silicon seeds}/N_{Uncorrected}^{INTT Doublets}",               //
                        {0, (!is_simulation) ? 5 : 2},                                                           //
                        false,                                                              //
                        false,                                                              //
                        510,                                                                //
                        std::make_pair(true, 1.0),                                          //
                        plotdir + "/" + outstem + "_ratio1D_seed_over_nraw_eta_comparison" //
        );
    }

    if (!phi_hists.empty())
    {
        draw_comparison(phi_hists,                                                         //
                        comparison_colors,                                                 //
                        labels,                                                            //
                        "#phi [radian]",                                                  //
                        "N_{Uncorrected}^{INTT Doublets}/N_{Silicon seeds}",              //
                        {0, (!is_simulation) ? 5 : 1.7},                                   //
                        false,                                                             //
                        false,                                                             //
                        510,                                                               //
                        std::make_pair(true, 1.0),                                         //
                        plotdir + "/" + outstem + "_ratio1D_nraw_over_seed_phi_comparison" //
        );

        draw_comparison(phi_hists_inverse,                                                  //
                        comparison_colors,                                                  //
                        labels,                                                             //
                        "#phi [radian]",                                                  //
                        "N_{Silicon seeds}/N_{Uncorrected}^{INTT Doublets}",              //
                        {0, (!is_simulation) ? 5 : 2},                                                           //
                        false,                                                              //
                        false,                                                              //
                        510,                                                                //
                        std::make_pair(true, 1.0),                                          //
                        plotdir + "/" + outstem + "_ratio1D_seed_over_nraw_phi_comparison" //
        );
    }

    if (!eta_hists_good.empty())
    {
        draw_comparison(eta_hists_good,                                                          //
                        comparison_colors_good,                                                   //
                        labels_good,                                                              //
                        "#eta",                                                                  //
                        "N_{Uncorrected}^{INTT Doublets}/N_{Silicon seeds}^{good}",             //
                        {0, (!is_simulation) ? 2 : 3},                                            //
                        false,                                                                    //
                        false,                                                                    //
                        510,                                                                      //
                        std::make_pair(true, 1.0),                                                //
                        plotdir + "/" + outstem + "_ratio1D_nraw_over_seed_good_eta_comparison" //
        );

        draw_comparison(eta_hists_good_inverse,                                                    //
                        comparison_colors_good,                                                    //
                        labels_good,                                                               //
                        "#eta",                                                                   //
                        "N_{Silicon seeds}^{good}/N_{Uncorrected}^{INTT Doublets}",              //
                        {0, (!is_simulation) ? 2 : 2},                                                                  //
                        false,                                                                     //
                        false,                                                                     //
                        510,                                                                       //
                        std::make_pair(true, 1.0),                                                 //
                        plotdir + "/" + outstem + "_ratio1D_seed_good_over_nraw_eta_comparison" //
        );
    }

    if (!phi_hists_good.empty())
    {
        draw_comparison(phi_hists_good,                                                          //
                        comparison_colors_good,                                                   //
                        labels_good,                                                              //
                        "#phi [radian]",                                                        //
                        "N_{Uncorrected}^{INTT Doublets}/N_{Silicon seeds}^{good}",             //
                        {0, (!is_simulation) ? 5 : 1.7},                                          //
                        false,                                                                    //
                        false,                                                                    //
                        510,                                                                      //
                        std::make_pair(true, 1.0),                                                //
                        plotdir + "/" + outstem + "_ratio1D_nraw_over_seed_good_phi_comparison" //
        );

        draw_comparison(phi_hists_good_inverse,                                                    //
                        comparison_colors_good,                                                    //
                        labels_good,                                                               //
                        "#phi [radian]",                                                         //
                        "N_{Silicon seeds}^{good}/N_{Uncorrected}^{INTT Doublets}",              //
                        {0, (!is_simulation) ? 5 : 2},                                                                  //
                        false,                                                                     //
                        false,                                                                     //
                        510,                                                                       //
                        std::make_pair(true, 1.0),                                                 //
                        plotdir + "/" + outstem + "_ratio1D_seed_good_over_nraw_phi_comparison" //
        );
    }

    for (auto &outputs : ratio_sets)
    {
        delete outputs.eta_ratio;
        delete outputs.phi_ratio;
        delete outputs.eta_ratio_good;
        delete outputs.phi_ratio_good;
    }

    for (auto *hist : eta_hists_inverse)
        delete hist;
    for (auto *hist : phi_hists_inverse)
        delete hist;
    for (auto *hist : eta_hists_good_inverse)
        delete hist;
    for (auto *hist : phi_hists_good_inverse)
        delete hist;
}
