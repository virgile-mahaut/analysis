#include <algorithm>
#include <cmath>
#include <glob.h>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

int nfiles = 100; // number of simulation ntuples to read

bool HasWildcard(const std::string &path) { return path.find('*') != std::string::npos || path.find('?') != std::string::npos || path.find('[') != std::string::npos; }

std::vector<std::string> ResolveInputFiles(const std::string &inputPattern)
{
    std::vector<std::string> files;
    if (inputPattern.empty())
    {
        return files;
    }

    if (!HasWildcard(inputPattern))
    {
        files.push_back(inputPattern);
        return files;
    }

    glob_t globResult;
    const int status = glob(inputPattern.c_str(), 0, nullptr, &globResult);
    if (status != 0)
    {
        globfree(&globResult);
        return files;
    }

    files.reserve(globResult.gl_pathc);
    for (std::size_t i = 0; i < globResult.gl_pathc; ++i)
    {
        files.emplace_back(globResult.gl_pathv[i]);
    }
    globfree(&globResult);
    std::sort(files.begin(), files.end());

    return files;
}

std::vector<double> DerivePercentileBoundaries(const std::vector<int> &values, int nPercentileBins = 10)
{
    std::vector<double> boundaries;
    if (values.empty() || nPercentileBins < 2)
    {
        return boundaries;
    }

    std::vector<int> sorted_values = values;
    std::sort(sorted_values.begin(), sorted_values.end());

    const double last_index = static_cast<double>(sorted_values.size() - 1);
    for (int i = 1; i < nPercentileBins; ++i)
    {
        const double fraction = static_cast<double>(i) / nPercentileBins;
        const double position = fraction * last_index;
        const auto low_index = static_cast<std::size_t>(position);
        const auto high_index = std::min(low_index + 1, sorted_values.size() - 1);
        const double weight = position - static_cast<double>(low_index);
        const double value = (1.0 - weight) * sorted_values[low_index] + weight * sorted_values[high_index];
        boundaries.push_back(value);
    }

    return boundaries;
}

std::vector<std::pair<int, int>> BuildPercentileIntervals(const std::vector<double> &boundaries)
{
    std::vector<std::pair<int, int>> intervals;
    if (boundaries.empty())
    {
        return intervals;
    }

    std::vector<int> integer_boundaries;
    integer_boundaries.reserve(boundaries.size());
    for (double boundary : boundaries)
    {
        integer_boundaries.push_back(static_cast<int>(std::ceil(boundary)));
    }

    intervals.emplace_back(integer_boundaries.back(), std::numeric_limits<int>::max());
    for (int i = static_cast<int>(integer_boundaries.size()) - 1; i > 0; --i)
    {
        intervals.emplace_back(integer_boundaries[i - 1], integer_boundaries[i]);
    }
    intervals.emplace_back(0, integer_boundaries.front());

    return intervals;
}

std::vector<std::string> BuildPercentileLabels(int nPercentileBins)
{
    std::vector<std::string> labels;
    labels.reserve(nPercentileBins);
    for (int i = 0; i < nPercentileBins; ++i)
    {
        std::ostringstream out;
        out << i * 100 / nPercentileBins << "-" << (i + 1) * 100 / nPercentileBins << "%";
        labels.push_back(out.str());
    }
    return labels;
}

std::vector<std::string> BuildIntervalLabels(const std::vector<std::pair<int, int>> &intervals)
{
    std::vector<std::string> labels;
    labels.reserve(intervals.size());
    for (const auto &[low, high] : intervals)
    {
        std::ostringstream out;
        out << "[" << low << ", ";
        if (high == std::numeric_limits<int>::max())
        {
            out << "INT_MAX";
        }
        else
        {
            out << high;
        }
        out << ")";
        labels.push_back(out.str());
    }
    return labels;
}

int FindPercentileBinIndex(const std::vector<std::pair<int, int>> &intervals, int value)
{
    for (std::size_t i = 0; i < intervals.size(); ++i)
    {
        if (value >= intervals[i].first && value < intervals[i].second)
        {
            return static_cast<int>(i);
        }
    }
    return -1;
}

void Draw2DHist(TH2 *hist, bool logz, const std::string &exttext, const std::string &filename)
{
    if (!hist)
    {
        return;
    }

    TCanvas *c = new TCanvas(Form("c_%s", hist->GetName()), hist->GetName(), 800, 700);
    c->cd();
    c->SetLogz(logz);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.08);
    hist->GetZaxis()->SetTitleOffset(1.5);
    hist->Draw("colz");

    TLatex text;
    text.SetNDC();
    text.SetTextFont(43);
    text.SetTextSize(22);
    if (!exttext.empty())
    {
        text.DrawLatex(0.18, 0.86, exttext.c_str());
    }

    c->SaveAs(Form("%s.pdf", filename.c_str()));
    c->SaveAs(Form("%s.png", filename.c_str()));
    delete c;
}

bool FitGaussianAndDraw(TH1D *hist, const std::string &label, const std::string &filename, double &mean, double &meanerr, double &sigma, double &sigmaerr)
{
    mean = 0.0;
    meanerr = 0.0;
    sigma = 0.0;
    sigmaerr = 0.0;

    if (!hist || hist->GetEntries() < 10)
    {
        return false;
    }

    TCanvas *c = new TCanvas(Form("c_%s", hist->GetName()), hist->GetName(), 800, 700);
    c->SetLogy();
    c->cd();
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.08);
    gPad->SetTopMargin(0.05);
    // gPad->SetBottomMargin(0.12);

    hist->SetLineColor(1);
    hist->SetLineWidth(2);
    hist->GetYaxis()->SetRangeUser(std::max(hist->GetMinimum(0.0) * 0.5, 1e-2), hist->GetMaximum() * 80.0);
    hist->Draw("hist");

    TF1 *f1 = new TF1(Form("f_%s", hist->GetName()), "gaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    f1->SetParameter(1, hist->GetMean());
    f1->SetParLimits(1, -1.0, 1.0);
    f1->SetParameter(2, hist->GetRMS());
    f1->SetParLimits(2, 0.0, 100.0);
    f1->SetLineColor(TColor::GetColor("#F54748"));
    hist->Fit(f1, "BQ");
    f1->Draw("same");

    mean = f1->GetParameter(1);
    meanerr = f1->GetParError(1);
    sigma = f1->GetParameter(2);
    sigmaerr = f1->GetParError(2);

    TLatex latex;
    latex.SetNDC();
    // latex.SetTextFont(43);
    latex.SetTextSize(0.03);
    latex.SetTextAlign(kHAlignLeft + kVAlignCenter);
    latex.DrawLatex(gPad->GetLeftMargin()+0.04, (1-gPad->GetTopMargin()) - 0.04, "#it{#bf{sPHENIX}} Simulation");
    latex.DrawLatex(gPad->GetLeftMargin()+0.04, (1-gPad->GetTopMargin()) - 0.08, label.c_str());
    latex.DrawLatex(gPad->GetLeftMargin()+0.04, (1-gPad->GetTopMargin()) - 0.12, Form("#mu=%.3g#pm%.3g [cm]", mean, meanerr));
    latex.DrawLatex(gPad->GetLeftMargin()+0.04, (1-gPad->GetTopMargin()) - 0.16, Form("#sigma=%.3g#pm%.3g [cm]", sigma, sigmaerr));

    c->SaveAs(Form("%s.pdf", filename.c_str()));
    c->SaveAs(Form("%s.png", filename.c_str()));

    delete f1;
    delete c;
    return true;
}

void DrawSummaryHistogram(TH1D *hist,                                //
                          const std::vector<std::string> &binlabels, //
                          const std::string &ytitle,                 //
                          const std::pair<float, float> &yrange,     //
                          const std::string &filename                //
)
{
    if (!hist)
    {
        return;
    }

    TCanvas *c = new TCanvas(Form("c_%s", hist->GetName()), hist->GetName(), 900, 700);
    c->cd();
    gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.08);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.14);

    hist->SetLineWidth(2);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.2);
    hist->GetXaxis()->SetTitle("INTT cluster percentile");
    hist->GetYaxis()->SetTitle(ytitle.c_str());
    hist->GetYaxis()->SetTitleOffset(1.7);
    hist->GetYaxis()->SetRangeUser(yrange.first, yrange.second);

    for (int i = 0; i < hist->GetNbinsX() && i < static_cast<int>(binlabels.size()); ++i)
    {
        hist->GetXaxis()->SetBinLabel(i + 1, binlabels[i].c_str());
    }
    hist->Draw("PE1");

    c->SaveAs(Form("%s.pdf", filename.c_str()));
    c->SaveAs(Form("%s.png", filename.c_str()));
    delete c;
}

void SavePercentileIntervals(TFile *outfile, const char *tree_name, const std::vector<std::pair<int, int>> &intervals, int nPercentileBins)
{
    if (!outfile)
    {
        return;
    }

    TTree *tree = new TTree(tree_name, tree_name);
    int percentile_bin = -1;
    double percentile_low = 0.0;
    double percentile_high = 0.0;
    int cluster_low = 0;
    int cluster_high = 0;

    tree->Branch("percentile_bin", &percentile_bin, "percentile_bin/I");
    tree->Branch("percentile_low", &percentile_low, "percentile_low/D");
    tree->Branch("percentile_high", &percentile_high, "percentile_high/D");
    tree->Branch("cluster_low", &cluster_low, "cluster_low/I");
    tree->Branch("cluster_high", &cluster_high, "cluster_high/I");

    const double percentile_width = 100.0 / static_cast<double>(nPercentileBins);
    for (std::size_t i = 0; i < intervals.size(); ++i)
    {
        percentile_bin = static_cast<int>(i);
        percentile_low = static_cast<double>(i) * percentile_width;
        percentile_high = static_cast<double>(i + 1) * percentile_width;
        cluster_low = intervals[i].first;
        cluster_high = intervals[i].second;
        tree->Fill();
    }

    outfile->cd();
    tree->Write();
}

void VertexInSimulation(const std::string &inputPattern = "/sphenix/tg/tg01/hf/hjheng/ppg-dNdEta-OOpp/simulation-nopileup-ntuple/ACTS/ntuple_*.root", int nPercentileBins = 10)
{
    std::string plotdir = "./figure/figure-VertexInSimulation";
    gSystem->mkdir(plotdir.c_str(), true);
    gStyle->SetOptStat(0);

    TChain *t = new TChain("VTX");
    const auto inputFiles = ResolveInputFiles(inputPattern);
    if (inputFiles.empty())
    {
        std::cout << "[VertexInSimulation] No files matched: " << inputPattern << std::endl;
        delete t;
        return;
    }

    const std::size_t maxFiles = std::min(static_cast<std::size_t>(nfiles), inputFiles.size());
    int nAddedFiles = 0;
    for (std::size_t i = 0; i < maxFiles; ++i)
    {
        if (t->Add(inputFiles[i].c_str()) > 0)
        {
            ++nAddedFiles;
        }
    }

    if (nAddedFiles == 0)
    {
        std::cout << "[VertexInSimulation] Failed to add any input files from pattern: " << inputPattern << std::endl;
        delete t;
        return;
    }

    std::cout << "[VertexInSimulation] Added " << nAddedFiles << " files to TChain." << std::endl;

    int counter = 0;
    int nTruthVertex = 0;
    std::vector<float> *TruthVertexZ = nullptr;
    std::vector<float> *trackerVertexZ = nullptr;
    std::vector<short> *trackerVertexCrossing = nullptr;
    std::vector<int> *cluster_timeBucketID = nullptr;

    t->SetBranchAddress("counter", &counter);
    t->SetBranchAddress("nTruthVertex", &nTruthVertex);
    t->SetBranchAddress("TruthVertexZ", &TruthVertexZ);
    t->SetBranchAddress("trackerVertexZ", &trackerVertexZ);
    t->SetBranchAddress("trackerVertexCrossing", &trackerVertexCrossing);
    t->SetBranchAddress("cluster_timeBucketID", &cluster_timeBucketID);

    std::vector<int> nInttClustersSelectedCrossings;
    const Long64_t nEntries = t->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        t->GetEntry(i);
        if (!trackerVertexZ || !trackerVertexCrossing || !cluster_timeBucketID || !TruthVertexZ)
        {
            continue;
        }

        std::map<int, std::vector<int>> vertexCrossingIdxMap;
        for (int ivtx = 0; ivtx < static_cast<int>(trackerVertexCrossing->size()); ++ivtx)
        {
            vertexCrossingIdxMap[trackerVertexCrossing->at(ivtx)].push_back(ivtx);
        }

        for (const auto &[vertex_crossing, vertex_indices] : vertexCrossingIdxMap)
        {
            if (vertex_crossing == std::numeric_limits<short int>::max() || vertex_indices.size() != 1)
            {
                continue;
            }

            int nClustersInCrossing = 0;
            for (int time_bucket : *cluster_timeBucketID)
            {
                if (time_bucket == vertex_crossing)
                {
                    ++nClustersInCrossing;
                }
            }
            nInttClustersSelectedCrossings.push_back(nClustersInCrossing);
        }
    }

    if (nInttClustersSelectedCrossings.empty())
    {
        std::cout << "[VertexInSimulation] No eligible crossings found after requiring one silicon vertex per unique crossing." << std::endl;
        delete t;
        return;
    }

    const auto percentile_boundaries = DerivePercentileBoundaries(nInttClustersSelectedCrossings, nPercentileBins);
    const auto percentile_intervals = BuildPercentileIntervals(percentile_boundaries);
    const auto percentile_labels = BuildPercentileLabels(nPercentileBins);
    const auto interval_labels = BuildIntervalLabels(percentile_intervals);

    std::cout << "[VertexInSimulation] Selected-crossing percentile boundaries:";
    for (double boundary : percentile_boundaries)
    {
        std::cout << " " << boundary;
    }
    std::cout << std::endl;

    for (std::size_t i = 0; i < percentile_intervals.size(); ++i)
    {
        std::cout << "  bin " << i << " (" << percentile_labels[i] << "): " << interval_labels[i] << std::endl;
    }

    TH2D *hM_TruthVtxZ_RecoTruthDiff = new TH2D("hM_TruthVtxZ_RecoTruthDiff", ";vtx_{Z}^{Truth} [cm];vtx_{Z}^{Reco}-vtx_{Z}^{Truth} [cm];Entries", 200, -25, 25, 200, -5, 5);
    std::vector<TH2D *> h2_percentile;
    std::vector<TH1D *> h1_percentile;
    h2_percentile.reserve(nPercentileBins);
    h1_percentile.reserve(nPercentileBins);

    for (int i = 0; i < nPercentileBins; ++i)
    {
        h2_percentile.push_back(new TH2D(Form("hM_TruthVtxZ_RecoTruthDiff_Percentile%d", i), ";vtx_{Z}^{Truth} [cm];vtx_{Z}^{Reco}-vtx_{Z}^{Truth} [cm];Entries", 200, -25, 25, 200, -5, 5));

        const double diffMax = (i < 7) ? 2.0 : 5.0;
        h1_percentile.push_back(new TH1D(Form("hM_DiffVtxZ_Percentile%d", i), ";#Deltaz(vtx_{Reco}, vtx_{Truth}) [cm];Entries", 100, -diffMax, diffMax));
    }

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        t->GetEntry(i);
        if (!trackerVertexZ || !trackerVertexCrossing || !cluster_timeBucketID || !TruthVertexZ)
        {
            continue;
        }
        if (nTruthVertex != 1 || TruthVertexZ->size() != 1)
        {
            continue;
        }

        std::map<int, std::vector<int>> vertexCrossingIdxMap;
        for (int ivtx = 0; ivtx < static_cast<int>(trackerVertexCrossing->size()); ++ivtx)
        {
            vertexCrossingIdxMap[trackerVertexCrossing->at(ivtx)].push_back(ivtx);
        }

        for (const auto &[vertex_crossing, vertex_indices] : vertexCrossingIdxMap)
        {
            if (vertex_crossing == std::numeric_limits<short int>::max() || vertex_indices.size() != 1)
            {
                continue;
            }

            int nClustersInCrossing = 0;
            for (int time_bucket : *cluster_timeBucketID)
            {
                if (time_bucket == vertex_crossing)
                {
                    ++nClustersInCrossing;
                }
            }

            const int percentile_bin = FindPercentileBinIndex(percentile_intervals, nClustersInCrossing);
            if (percentile_bin < 0)
            {
                continue;
            }

            const int vtx_idx = vertex_indices.front();
            const double truthVtxZ = TruthVertexZ->at(0);
            const double recoTruthDiff = trackerVertexZ->at(vtx_idx) - truthVtxZ;
            hM_TruthVtxZ_RecoTruthDiff->Fill(truthVtxZ, recoTruthDiff);
            h2_percentile[percentile_bin]->Fill(truthVtxZ, recoTruthDiff);
            h1_percentile[percentile_bin]->Fill(recoTruthDiff);
        }
    }

    Draw2DHist(hM_TruthVtxZ_RecoTruthDiff, true, "All selected crossings", Form("%s/TruthVtxZ_RecoTruthDiff_Sim", plotdir.c_str()));
    for (int i = 0; i < nPercentileBins; ++i)
    {
        const std::string label = percentile_labels[i] + ", N_{INTT} " + interval_labels[i];
        Draw2DHist(h2_percentile[i], true, label, Form("%s/TruthVtxZ_RecoTruthDiff_Percentile%d_Sim", plotdir.c_str(), i));
    }

    TH1D *h_bias = new TH1D("h_bias_percentile", ";INTT cluster percentile;#mu_{#Deltaz}^{Gaussian} [cm]", nPercentileBins, 0, nPercentileBins);
    TH1D *h_resolution = new TH1D("h_resolution_percentile", ";INTT cluster percentile;#sigma_{#Deltaz}^{Gaussian} [cm]", nPercentileBins, 0, nPercentileBins);

    for (int i = 0; i < nPercentileBins; ++i)
    {
        double mean = 0.0;
        double meanerr = 0.0;
        double sigma = 0.0;
        double sigmaerr = 0.0;
        const std::string label = percentile_labels[i] + ", N_{INTT} " + interval_labels[i];
        if (FitGaussianAndDraw(h1_percentile[i], label, Form("%s/DiffVtxZ_Percentile%d_Sim", plotdir.c_str(), i), mean, meanerr, sigma, sigmaerr))
        {
            h_bias->SetBinContent(i + 1, mean);
            h_bias->SetBinError(i + 1, meanerr);
            h_resolution->SetBinContent(i + 1, sigma);
            h_resolution->SetBinError(i + 1, sigmaerr);
        }
    }

    DrawSummaryHistogram(h_bias, percentile_labels, "#mu_{#Deltaz}^{Gaussian} [cm]", std::make_pair(-0.0075, 0.0075), Form("%s/VtxZBias_Percentile", plotdir.c_str()));
    DrawSummaryHistogram(h_resolution, percentile_labels, "#sigma_{#Deltaz}^{Gaussian} [cm]", std::make_pair(0.0, 0.08), Form("%s/VtxZResolution_Percentile", plotdir.c_str()));

    TFile *outfile = TFile::Open(Form("%s/VertexInSimulation.root", plotdir.c_str()), "RECREATE");
    SavePercentileIntervals(outfile, "selected_crossing_percentile_intervals", percentile_intervals, nPercentileBins);
    hM_TruthVtxZ_RecoTruthDiff->Write();
    h_bias->Write();
    h_resolution->Write();
    for (TH2D *hist : h2_percentile)
    {
        hist->Write();
    }
    for (TH1D *hist : h1_percentile)
    {
        hist->Write();
    }
    outfile->Close();
    delete outfile;

    delete hM_TruthVtxZ_RecoTruthDiff;
    delete h_bias;
    delete h_resolution;
    for (TH2D *hist : h2_percentile)
    {
        delete hist;
    }
    for (TH1D *hist : h1_percentile)
    {
        delete hist;
    }
    delete t;
}
