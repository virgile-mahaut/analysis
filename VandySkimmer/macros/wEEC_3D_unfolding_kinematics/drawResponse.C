#include "analysisHelper.h"

// ─────────────────────────────────────────────────────────────────────────────
//  drawResponse.C
//
//  Visualises one ΔΦ-bin response matrix at four successive zoom levels.
//
//  Level 0 — full response matrix (nRecoFlat3D × nTrueFlat3D)
//  Level 1 — fixed truth lead pT bin (iLsel): all subl pT × pair weight bins
//  Level 2 — fixed truth (lead, subl) bin (iLsel, iSsel): all pair weight bins
//  Level 3 — fixed truth (lead, subl, pairweight) cell: column + row projections
//
//  Follows the same calling convention as drawClosure.C:
//    outDir   — working directory; plots written to outDir/Plots/
//    respFile — merged response ROOT file (output of merge macro)
//    mode     — Mode::kFull / kHalf / kData (selects ModeLabel for filenames)
//
//  Usage:
//    root -l 'drawResponse.C(".", "resp_merged.root", Mode::kFull,
//                             4,   // kDphi  — ΔΦ bin index
//                             2,   // iLsel  — lead pT bin  (0-based)
//                             3,   // iSsel  — subl pT bin  (0-based)
//                             5)'  // iPwsel — pair weight bin (0-based)
//
//  Outputs (in outDir/Plots/):
//    response_level0-{label}.pdf  — full matrix
//    response_level1-{label}.pdf  — lead pT slice
//    response_level2-{label}.pdf  — dijet pT slice
//    response_level3-{label}.pdf  — single-cell projections
// ─────────────────────────────────────────────────────────────────────────────

// ── helpers ──────────────────────────────────────────────────────────────────

// Build a sub-matrix TH2D from the full response TH2D by selecting a subset
// of reco (X) and truth (Y) bins.  recoSel and trueSel are 1-based bin indices
// into the full response histogram.
TH2D* SubMatrix(TH2D* hFull,
                const std::vector<int>& recoSel,
                const std::vector<int>& trueSel,
                const char* name, const char* title,
                const char* xlbl, const char* ylbl)
{
    int nX = (int)recoSel.size();
    int nY = (int)trueSel.size();
    TH2D* hSub = new TH2D(name, title, nX, 0, nX, nY, 0, nY);
    hSub->GetXaxis()->SetTitle(xlbl);
    hSub->GetYaxis()->SetTitle(ylbl);
    hSub->SetDirectory(0);
    for (int ix = 0; ix < nX; ++ix)
    for (int iy = 0; iy < nY; ++iy)
        hSub->SetBinContent(ix+1, iy+1,
            hFull->GetBinContent(recoSel[ix], trueSel[iy]));
    return hSub;
}

// Draw a TH2D on the current pad with a COLZ palette, log-z if any bin > 0.
// A dashed red diagonal is drawn whenever both axes have the same number of bins.
void DrawMatrix(TH2D* h, const char* ztitle = "cross-section weight")
{
    h->GetZaxis()->SetTitle(ztitle);
    gPad->SetLogz(h->GetMaximum() > 0);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.20);
    gPad->SetLeftMargin(0.20);
    h->Draw("COLZ");
    gPad->Update();  // establish axis transforms before any overlaid lines
    if (h->GetNbinsX() == h->GetNbinsY()) {
        double xlo = gPad->GetUxmin(), xhi = gPad->GetUxmax();
        double ylo = gPad->GetUymin(), yhi = gPad->GetUymax();
        TLine* diag = new TLine(xlo, ylo, xhi, yhi);
        diag->SetLineColor(kRed); diag->SetLineWidth(1);
        diag->SetLineStyle(2); diag->Draw();
    }
}

// ── main macro ───────────────────────────────────────────────────────────────
void drawResponse(const char* outDir   = ".",
                  const char* respFile = nullptr,
                  Mode        mode     = Mode::kFull,
                  int  kDphi   = 0,    // ΔΦ bin index
                  int  iLsel   = 2,    // lead pT bin index   (0-based)
                  int  iSsel   = 1,    // subl pT bin index   (0-based)
                  int  iPwsel  = 5)    // pair weight bin index (0-based)
{
    // ── sanity checks ────────────────────────────────────────────────────
    if (iLsel  < 0 || iLsel  >= nTrueLead) {
        std::cerr << "iLsel out of range [0," << nTrueLead-1 << "]\n"; return; }
    if (iSsel  < 0 || iSsel  >= nTrueSubl) {
        std::cerr << "iSsel out of range [0," << nTrueSubl-1 << "]\n"; return; }
    if (iPwsel < 0 || iPwsel >= nPairWeight) {
        std::cerr << "iPwsel out of range [0," << nPairWeight-1 << "]\n"; return; }
    if (TrueFlatIndex(iLsel, iSsel) < 0) {
        std::cerr << "iLsel/iSsel is a kinematically forbidden (iL,iS) pair\n"; return; }

    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    const std::string label   = ModeLabel(mode);
    const std::string plotDir = std::format("{}/Plots", outDir);

    // ── open file ────────────────────────────────────────────────────────
    if (!respFile) { std::cerr << "respFile is required\n"; return; }
    TFile* fIn = TFile::Open(respFile, "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "Cannot open " << respFile << "\n"; return; }

    RooUnfoldResponse* resp = (RooUnfoldResponse*)fIn->Get(
        std::format("response_wEEC3D_{}", kDphi).c_str());
    if (!resp) {
        std::cerr << "Cannot find response_wEEC3D_" << kDphi << "\n"; return; }

    TH2D* hFull = (TH2D*)resp->Hresponse();
    if (!hFull) {
        std::cerr << "Response has no internal TH2D\n"; return; }
    hFull = (TH2D*)hFull->Clone("hFull"); hFull->SetDirectory(0);

    // ── human-readable labels ────────────────────────────────────────────
    auto leadLbl = [](int iL) {
        return std::format("{:.0f}-{:.0f} GeV",
            trueLeadPtBins[iL], trueLeadPtBins[iL+1]); };
    auto sublLbl = [](int iS) {
        return std::format("{:.0f}-{:.0f} GeV",
            trueSublPtBins[iS], trueSublPtBins[iS+1]); };
    auto pwLbl = [](int iPw) {
        return std::format("{:.3g}-{:.3g}",
            pairWeightBins[iPw], pairWeightBins[iPw+1]); };

    const std::string dphiLbl = std::format(
        "#Delta#phi bin {} ({:.2f}-{:.2f})",
        kDphi, dPhiBins[kDphi], dPhiBins[kDphi+1]);

    // ══════════════════════════════════════════════════════════════════════
    //  Build bin-index selection vectors for each zoom level
    // ══════════════════════════════════════════════════════════════════════

    // ── Level 1: all (subl, pairweight) bins for lead == iLsel ───────────
    std::vector<int> recoSel1, trueSel1;
    for (int iF = 0; iF < nRecoFlat3D(); ++iF) {
        int iL, iS, iPw; RecoFlat3DToIJK(iF, iL, iS, iPw);
        if (iL == iLsel) recoSel1.push_back(iF + 1);  // 1-based
    }
    for (int iF = 0; iF < nTrueFlat3D(); ++iF) {
        int iL, iS, iPw; TrueFlat3DToIJK(iF, iL, iS, iPw);
        if (iL == iLsel) trueSel1.push_back(iF + 1);
    }

    // ── Level 2: all pair weight bins for (lead, subl) == (iLsel, iSsel) ─
    std::vector<int> recoSel2, trueSel2;
    for (int iPw = 0; iPw < nPairWeight; ++iPw) {
        int iF = RecoFlat3DIndex(iLsel, iSsel, iPw);
        if (iF >= 0) recoSel2.push_back(iF + 1);
    }
    for (int iPw = 0; iPw < nPairWeight; ++iPw) {
        int iF = TrueFlat3DIndex(iLsel, iSsel, iPw);
        if (iF >= 0) trueSel2.push_back(iF + 1);
    }

    // ── Level 3: single (lead, subl, pairweight) cell ────────────────────
    // Show the full reco column (all reco bins for this truth cell) and
    // full truth row (all truth bins for this reco cell) as 1D projections.
    int truthCell = TrueFlat3DIndex(iLsel, iSsel, iPwsel) + 1;  // 1-based
    int recoCell  = RecoFlat3DIndex(iLsel, iSsel, iPwsel) + 1;

    // Column projection: truth bin = truthCell, vary reco
    TH1D* hCol = new TH1D("hCol",
        std::format("Truth ({},{},{}) #rightarrow reco migration;reco flat3D bin;weight",
            iLsel, iSsel, iPwsel).c_str(),
        hFull->GetNbinsX(), 0, hFull->GetNbinsX());
    hCol->SetDirectory(0);
    for (int ib = 1; ib <= hFull->GetNbinsX(); ++ib)
        hCol->SetBinContent(ib, hFull->GetBinContent(ib, truthCell));

    // Row projection: reco bin = recoCell, vary truth
    TH1D* hRow = new TH1D("hRow",
        std::format("Reco ({},{},{}) #leftarrow truth sources;truth flat3D bin;weight",
            iLsel, iSsel, iPwsel).c_str(),
        hFull->GetNbinsY(), 0, hFull->GetNbinsY());
    hRow->SetDirectory(0);
    for (int ib = 1; ib <= hFull->GetNbinsY(); ++ib)
        hRow->SetBinContent(ib, hFull->GetBinContent(recoCell, ib));

    // ── derive output filenames ───────────────────────────────────────────
    auto outName = [&](int lvl) {
        return std::format("{}/response_level{}-{}.pdf", plotDir, lvl, label); };

    // ══════════════════════════════════════════════════════════════════════
    //  Level 0 — full response matrix
    // ══════════════════════════════════════════════════════════════════════
    {
        TCanvas* c0 = new TCanvas("cL0", "Level 0: full response matrix", 1000, 900);
        c0->SetRightMargin(0.15);
        c0->SetBottomMargin(0.15);
        c0->SetLeftMargin(0.15);

        TH2D* hL0 = (TH2D*)hFull->Clone("hL0");
        hL0->SetTitle(
            std::format("Full response matrix | {};reco flat3D;truth flat3D",
                dphiLbl).c_str());
        hL0->GetXaxis()->SetTitle("reco flat3D  [lead #times subl #times pair weight]");
        hL0->GetYaxis()->SetTitle("truth flat3D [lead #times subl #times pair weight]");
        DrawMatrix(hL0);

        // Draw lines at each lead pT bin boundary to aid orientation.
        // The boundary falls every nPairWeight * (number of valid subl bins
        // for that lead bin) reco/truth flat3D bins.
        auto drawBoundaries = [&](bool isTruth) {
            int nBins = isTruth ? nTrueFlat3D() : nRecoFlat3D();
            int nTotal = isTruth ? nTrueFlat3D() : nRecoFlat3D();
            int prev = -1, curL = -1;
            for (int iF = 0; iF < nBins; ++iF) {
                int iL, iS, iPw;
                if (isTruth) TrueFlat3DToIJK(iF, iL, iS, iPw);
                else         RecoFlat3DToIJK(iF, iL, iS, iPw);
                if (iL != curL) {
                    if (curL >= 0) {
                        // draw boundary line at iF
                        TLine* l = isTruth
                            ? new TLine(0, iF, nRecoFlat3D(), iF)
                            : new TLine(iF, 0, iF, nTrueFlat3D());
                        l->SetLineColor(kWhite); l->SetLineWidth(1);
                        l->SetLineStyle(2); l->Draw();
                    }
                    curL = iL;
                }
            }
        };
        drawBoundaries(false);  // reco (X) boundaries
        drawBoundaries(true);   // truth (Y) boundaries

        // Label each lead pT block center on both axes
        {
            // Build a map: iL -> first flat3D index for that lead bin
            auto labelLeadBins = [&](TAxis* ax, bool isTruth) {
                int nBins = isTruth ? nTrueFlat3D() : nRecoFlat3D();
                int prevL = -1, blockStart = 0;
                for (int iF = 0; iF <= nBins; ++iF) {
                    int iL = -1;
                    if (iF < nBins) {
                        int iS, iPw;
                        if (isTruth) TrueFlat3DToIJK(iF, iL, iS, iPw);
                        else         RecoFlat3DToIJK(iF, iL, iS, iPw);
                    }
                    if (iL != prevL) {
                        if (prevL >= 0) {
                            int mid = (blockStart + iF) / 2;
                            ax->SetBinLabel(mid + 1,
                                std::format("{:.0f}-{:.0f}",
                                    trueLeadPtBins[prevL],
                                    trueLeadPtBins[prevL+1]).c_str());
                        }
                        blockStart = iF;
                        prevL = iL;
                    }
                }
            };
            labelLeadBins(hL0->GetXaxis(), false);
            labelLeadBins(hL0->GetYaxis(), true);
            hL0->GetXaxis()->LabelsOption("h");
            hL0->GetYaxis()->LabelsOption("h");
            hL0->GetXaxis()->SetLabelSize(0.03);
            hL0->GetYaxis()->SetLabelSize(0.03);
        }

        /*
        // Dashed black lines at each lead pT block boundary
        {
            double xlo = gPad->GetUxmin(), xhi = gPad->GetUxmax();
            double ylo = gPad->GetUymin(), yhi = gPad->GetUymax();
            auto drawLeadBoundaries = [&](bool isTruth) {
                int nBins   = isTruth ? nTrueFlat3D() : nRecoFlat3D();
                int nOther  = isTruth ? nRecoFlat3D() : nTrueFlat3D();
                double span = isTruth ? (xhi - xlo) : (yhi - ylo);
                double base = isTruth ? xlo : ylo;
                double otherLo = isTruth ? ylo : xlo;
                double otherHi = isTruth ? yhi : xhi;
                int curL = -1;
                for (int iF = 0; iF < nBins; ++iF) {
                    int iL, iS, iPw;
                    if (isTruth) TrueFlat3DToIJK(iF, iL, iS, iPw);
                    else         RecoFlat3DToIJK(iF, iL, iS, iPw);
                    if (iL != curL && curL >= 0) {
                        double pos = base + span * iF / nBins;
                        TLine* l = isTruth
                            ? new TLine(otherLo, pos, otherHi, pos)
                            : new TLine(pos, otherLo, pos, otherHi);
                        l->SetLineColor(kBlack); l->SetLineWidth(1);
                        l->SetLineStyle(2); l->Draw();
                    }
                    curL = iL;
                }
            };
            drawLeadBoundaries(false);
            drawLeadBoundaries(true);
        }
        */

        gPad->Update();  // flush COLZ before printing
        c0->Print(outName(0).c_str());
        std::cout << "Written: " << outName(0) << "\n";
        delete c0;
    }

    // ══════════════════════════════════════════════════════════════════════
    //  Level 1 — lead pT slice
    // ══════════════════════════════════════════════════════════════════════
    {
        TCanvas* c1 = new TCanvas("cL1", "Level 1: lead pT slice", 900, 800);
        c1->SetRightMargin(0.15);
        c1->SetBottomMargin(0.22);
        c1->SetLeftMargin(0.22);

        TH2D* hL1 = SubMatrix(hFull, recoSel1, trueSel1, "hL1",
            std::format("Lead pT = {} | {}",
                leadLbl(iLsel), dphiLbl).c_str(),
            "", "");
        hL1->GetXaxis()->SetTitle("reco subl p_{T} bin");
        hL1->GetYaxis()->SetTitle("truth subl p_{T} bin");
        DrawMatrix(hL1);

        // Label each subl pT block center on both axes — same approach as
        // Level 0's lead pT block labeling.
        auto labelSublBins = [&](TAxis* ax, const std::vector<int>& sel, bool isTruth) {
            int n = (int)sel.size();
            int prevS = -1, blockStart = 0;
            for (int i = 0; i <= n; ++i) {
                int curS = -1;
                if (i < n) {
                    int iL, iS, iPw;
                    if (isTruth) TrueFlat3DToIJK(sel[i]-1, iL, iS, iPw);
                    else         RecoFlat3DToIJK(sel[i]-1, iL, iS, iPw);
                    curS = iS;
                }
                if (curS != prevS) {
                    if (prevS >= 0) {
                        int mid = (blockStart + i) / 2;
                        ax->SetBinLabel(mid + 1,
                            std::format("{:.0f}-{:.0f}",
                                trueSublPtBins[prevS],
                                trueSublPtBins[prevS+1]).c_str());
                    }
                    blockStart = i;
                    prevS = curS;
                }
            }
            ax->LabelsOption("h");
            ax->SetLabelSize(0.03);
        };
        labelSublBins(hL1->GetXaxis(), recoSel1, false);
        labelSublBins(hL1->GetYaxis(), trueSel1, true);

        /*
        // Dashed black lines at each subl pT block boundary
        {
            double xlo = gPad->GetUxmin(), xhi = gPad->GetUxmax();
            double ylo = gPad->GetUymin(), yhi = gPad->GetUymax();
            int nX = (int)recoSel1.size(), nY = (int)trueSel1.size();
            int prevS = -1;
            for (int ix = 0; ix < nX; ++ix) {
                int iL, iS, iPw; RecoFlat3DToIJK(recoSel1[ix]-1, iL, iS, iPw);
                if (iS != prevS && prevS >= 0) {
                    double pos = xlo + (xhi - xlo) * ix / nX;
                    TLine* l = new TLine(pos, ylo, pos, yhi);
                    l->SetLineColor(kBlack); l->SetLineWidth(1);
                    l->SetLineStyle(2); l->Draw();
                }
                prevS = iS;
            }
            prevS = -1;
            for (int iy = 0; iy < nY; ++iy) {
                int iL, iS, iPw; TrueFlat3DToIJK(trueSel1[iy]-1, iL, iS, iPw);
                if (iS != prevS && prevS >= 0) {
                    double pos = ylo + (yhi - ylo) * iy / nY;
                    TLine* l = new TLine(xlo, pos, xhi, pos);
                    l->SetLineColor(kBlack); l->SetLineWidth(1);
                    l->SetLineStyle(2); l->Draw();
                }
                prevS = iS;
            }
        }
        */

        c1->Print(outName(1).c_str());
        std::cout << "Written: " << outName(1) << "\n";
        delete c1;
    }

    // ══════════════════════════════════════════════════════════════════════
    //  Level 2 — dijet pT slice (lead + subl fixed, all pair weights)
    // ══════════════════════════════════════════════════════════════════════
    {
        TCanvas* c2 = new TCanvas("cL2", "Level 2: dijet pT slice", 800, 750);
        c2->SetRightMargin(0.15);
        c2->SetBottomMargin(0.28);
        c2->SetLeftMargin(0.28);

        TH2D* hL2 = SubMatrix(hFull, recoSel2, trueSel2, "hL2",
            std::format("Lead {} / Subl {} | {}",
                leadLbl(iLsel), sublLbl(iSsel), dphiLbl).c_str(),
            "", "");
        hL2->GetXaxis()->SetTitle("reco pair weight bin");
        hL2->GetYaxis()->SetTitle("truth pair weight bin");

        // Label each bin with its low and high pair weight edge.
        auto labelPWBins = [&](TAxis* ax, const std::vector<int>& sel, bool isTruth) {
            ax->SetLabelSize(0.025);
            for (int i = 0; i < (int)sel.size(); ++i) {
                int iL, iS, iPw;
                if (isTruth) TrueFlat3DToIJK(sel[i]-1, iL, iS, iPw);
                else         RecoFlat3DToIJK(sel[i]-1, iL, iS, iPw);
                ax->SetBinLabel(i + 1,
                    std::format("{:.2g}-{:.2g}",
                        pairWeightBins[iPw], pairWeightBins[iPw+1]).c_str());
            }
            ax->LabelsOption("v");
        };
        labelPWBins(hL2->GetXaxis(), recoSel2, false);
        labelPWBins(hL2->GetYaxis(), trueSel2, true);
        DrawMatrix(hL2);

        c2->Print(outName(2).c_str());
        std::cout << "Written: " << outName(2) << "\n";
        delete c2;
    }

    // ══════════════════════════════════════════════════════════════════════
    //  Level 3 — single cell: column and row projections side by side
    // ══════════════════════════════════════════════════════════════════════
    {
        TCanvas* c3 = new TCanvas("cL3", "Level 3: single cell", 1400, 600);
        c3->Divide(2, 1);

        // ── left: column (truth→reco migration) ──────────────────────────
        c3->cd(1);
        gPad->SetLogy(hCol->GetMaximum() > 0);
        gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.12);
        hCol->SetFillColor(kAzure+7); hCol->SetLineColor(kAzure+9);
        hCol->Draw("HIST");
        TLatex tx3; tx3.SetNDC(); tx3.SetTextSize(0.038);
        tx3.DrawLatex(0.17, 0.93,
            std::format("Truth ({},{},{}) #rightarrow reco",
                iLsel, iSsel, iPwsel).c_str());
        tx3.DrawLatex(0.17, 0.87,
            std::format("Lead {}  Subl {}  PW {}",
                leadLbl(iLsel), sublLbl(iSsel), pwLbl(iPwsel)).c_str());
        tx3.DrawLatex(0.17, 0.81, dphiLbl.c_str());

        // ── right: row (reco←truth sources) ──────────────────────────────
        c3->cd(2);
        gPad->SetLogy(hRow->GetMaximum() > 0);
        gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.12);
        hRow->SetFillColor(kOrange-3); hRow->SetLineColor(kOrange+2);
        hRow->Draw("HIST");
        TLatex tx4; tx4.SetNDC(); tx4.SetTextSize(0.038);
        tx4.DrawLatex(0.17, 0.93,
            std::format("Reco ({},{},{}) #leftarrow truth",
                iLsel, iSsel, iPwsel).c_str());
        tx4.DrawLatex(0.17, 0.87,
            std::format("Lead {}  Subl {}  PW {}",
                leadLbl(iLsel), sublLbl(iSsel), pwLbl(iPwsel)).c_str());
        tx4.DrawLatex(0.17, 0.81, dphiLbl.c_str());

        c3->Print(outName(3).c_str());
        std::cout << "Written: " << outName(3) << "\n";
        delete c3;
    }

    fIn->Close(); delete fIn;
}
