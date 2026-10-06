#include "analysisHelper.h"

// ─────────────────────────────────────────────────────────────────────────────
//  makeVtxWeights.C
//
//  Forms the vtx_z reweighting factors needed to correct the MC vtx_z
//  distribution to match data.
//
//  Inputs:
//    mcFile   — merged MC file (output of hadd over all kVtx response files)
//               must contain hVtxMC
//    dataFile — merged data file (output of hadd over all fillData.C outputs)
//               must contain hVtxData
//
//  Procedure:
//    1. Load hVtxMC and hVtxData from the two merged files.
//    2. Normalize both to unit area within [-10, 10) cm.
//    3. Divide: hVtxWeight = hVtxData_norm / hVtxMC_norm.
//    4. Print the bin contents as a C++ array to paste into analysisHelper.h.
//    5. Write all histograms and a comparison canvas to outFile.
//
//  Usage:
//    root -l 'makeVtxWeights.C("vtxWeights.root",
//                               "mc_vtx_merged.root",
//                               "data_merged.root")'
//
//  After running, replace the vtxWeights[] array in analysisHelper.h with
//  the printed output, then rerun fillResponse.C in kData mode.
// ─────────────────────────────────────────────────────────────────────────────

void makeVtxWeights(const char* outDir  = "vtxWeights.root")
{

    std::string MCFile = std::format("{}/response-all-vtx.root",outDir);
    std::string dataFile = std::format("{}/data_measured-all.root",outDir);
    std::string outFile = std::format("{}/vtxWeights.root",outDir);

    // ── load histograms ───────────────────────────────────────────────────
    TFile* fMC = TFile::Open(MCFile.c_str(), "READ");
    if (!fMC || fMC->IsZombie()) {
        std::cerr << "Cannot open MC file: " << MCFile << "\n"; return; }

    TH1D* hMC = (TH1D*)fMC->Get("hVtxMC");
    if (!hMC) {
        std::cerr << "hVtxMC not found in " << MCFile << "\n"; return; }
    hMC = (TH1D*)hMC->Clone("hVtxMC_sum"); hMC->SetDirectory(0);
    fMC->Close(); delete fMC;

    TFile* fData = TFile::Open(dataFile.c_str(), "READ");
    if (!fData || fData->IsZombie()) {
        std::cerr << "Cannot open data file: " << dataFile << "\n"; return; }

    TH1D* hData = (TH1D*)fData->Get("hVtxData");
    if (!hData) {
        std::cerr << "hVtxData not found in " << dataFile << "\n"; return; }
    hData = (TH1D*)hData->Clone("hVtxData_sum"); hData->SetDirectory(0);
    fData->Close(); delete fData;

    if (hMC  ->Integral() <= 0) { std::cerr << "MC vtx histogram is empty\n";   return; }
    if (hData->Integral() <= 0) { std::cerr << "Data vtx histogram is empty\n"; return; }

    // ── normalize both to unit area within [-10, 10) cm ──────────────────
    // Use Integral() which sums all bins — since both histograms are already
    // restricted to this range, this is the full integral.
    TH1D* hMC_norm   = (TH1D*)hMC  ->Clone("hVtxMC_norm");
    TH1D* hData_norm = (TH1D*)hData->Clone("hVtxData_norm");
    hMC_norm  ->Scale(1.0 / hMC_norm  ->Integral());
    hData_norm->Scale(1.0 / hData_norm->Integral());

    // ── form ratio ────────────────────────────────────────────────────────
    TH1D* hWeight = (TH1D*)hData_norm->Clone("hVtxWeight");
    hWeight->SetTitle("vtx_z reweight: data/MC;v_{tz} (cm);data/MC");

    for (int i = 1; i <= nVtxBins; ++i) {
        double mc  = hMC_norm  ->GetBinContent(i);
        double dat = hData_norm->GetBinContent(i);
        double mc_err  = hMC_norm  ->GetBinError(i);
        double dat_err = hData_norm->GetBinError(i);
        if (mc > 0) {
            double w   = dat / mc;
            double err = w * std::sqrt(std::pow(dat_err/dat, 2) +
                                       std::pow(mc_err /mc,  2));
            hWeight->SetBinContent(i, w);
            hWeight->SetBinError  (i, err);
        } else {
            // MC has no entries in this bin — weight undefined, set to 1
            // (should not occur within [-10,10) cm for a large MC sample).
            std::cerr << std::format(
                "Warning: MC vtx bin {} ([{:.0f},{:.0f}) cm) is empty — weight set to 1\n",
                i, vtxZMin + (i-1)*vtxBinWidth, vtxZMin + i*vtxBinWidth);
            hWeight->SetBinContent(i, 1.0);
            hWeight->SetBinError  (i, 0.0);
        }
    }

    // ── print array for analysisHelper.h ─────────────────────────────────
    std::cout << "\n";
    std::cout << "// ── paste into analysisHelper.h vtxWeights[] ──────────\n";
    std::cout << "constexpr double vtxWeights[nVtxBins] = {\n";
    for (int i = 1; i <= nVtxBins; ++i) {
        double w   = hWeight->GetBinContent(i);
        double lo  = vtxZMin + (i-1) * vtxBinWidth;
        double hi  = lo + vtxBinWidth;
        std::cout << std::format("    {:.6f},  // bin {:2d}: [{:.0f}, {:.0f}) cm\n",
            w, i-1, lo, hi);
    }
    std::cout << "};\n\n";

    // ── draw comparison ───────────────────────────────────────────────────
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("cVtx", "vtx_z distributions", 1200, 500);
    c->Divide(3, 1);

    c->cd(1);
    hMC_norm  ->SetLineColor(kBlue+1);   hMC_norm  ->SetLineWidth(2);
    hData_norm->SetLineColor(kRed+1);    hData_norm->SetLineWidth(2);
    double ymax = std::max(hMC_norm->GetMaximum(), hData_norm->GetMaximum()) * 1.15;
    hMC_norm->GetYaxis()->SetRangeUser(0, ymax);
    hMC_norm->SetTitle("vtx_z: MC vs data (normalized);v_{tz} (cm);a.u.");
    hMC_norm->Draw("HIST");
    hData_norm->Draw("E SAME");
    TLegend* leg = new TLegend(0.15, 0.75, 0.55, 0.90);
    leg->AddEntry(hMC_norm,   "MC (xsec weighted)", "l");
    leg->AddEntry(hData_norm, "Data",               "lep");
    leg->Draw();

    c->cd(2);
    hWeight->SetLineColor(kBlack); hWeight->SetLineWidth(2);
    hWeight->SetMarkerStyle(20);
    hWeight->GetYaxis()->SetRangeUser(0, 2);
    hWeight->Draw("E");
    TLine* unity = new TLine(vtxZMin, 1, vtxZMax, 1);
    unity->SetLineStyle(2); unity->SetLineColor(kGray+1); unity->Draw();

    c->cd(3);
    // Show MC before and after reweighting as a closure check
    TH1D* hMC_rw = (TH1D*)hMC_norm->Clone("hVtxMC_reweighted");
    for (int i = 1; i <= nVtxBins; ++i)
        hMC_rw->SetBinContent(i,
            hMC_norm->GetBinContent(i) * hWeight->GetBinContent(i));
    hMC_rw->Scale(1.0 / hMC_rw->Integral());
    hMC_rw->SetLineColor(kGreen+2); hMC_rw->SetLineWidth(2);
    hMC_rw->SetTitle("Closure: reweighted MC vs data;v_{tz} (cm);a.u.");
    hData_norm->GetYaxis()->SetRangeUser(0, ymax);
    hData_norm->Draw("E");
    hMC_rw->Draw("HIST SAME");
    TLegend* leg2 = new TLegend(0.15, 0.75, 0.65, 0.90);
    leg2->AddEntry(hData_norm, "Data",              "lep");
    leg2->AddEntry(hMC_rw,    "MC (reweighted)",   "l");
    leg2->Draw();

    // ── write output ──────────────────────────────────────────────────────
    TFile* fOut = new TFile(outFile.c_str(), "RECREATE");
    hMC      ->Write("hVtxMC");
    hData    ->Write("hVtxData");
    hMC_norm ->Write("hVtxMC_norm");
    hData_norm->Write("hVtxData_norm");
    hWeight  ->Write("hVtxWeight");
    hMC_rw   ->Write("hVtxMC_reweighted");
    c->Write("cVtxComparison");
    fOut->Close();

    std::cout << "Written: " << outFile << "\n";
}
