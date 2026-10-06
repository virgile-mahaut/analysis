#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <glob.h>
#include <iostream>
#include <limits>
#include <sstream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

bool HasWildcard(const std::string &path)
{
    return path.find('*') != std::string::npos || path.find('?') != std::string::npos || path.find('[') != std::string::npos;
}

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

std::vector<double> DerivePercentileBoundaries(const std::vector<int> &values, int nPercentileBins = 5)
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

std::vector<std::string> BuildPercentileBoundaryLabels(int nPercentileBins)
{
    std::vector<std::string> labels;
    if (nPercentileBins < 2)
    {
        return labels;
    }

    const double percentile_width = 100.0 / static_cast<double>(nPercentileBins);
    for (int i = 1; i < nPercentileBins; ++i)
    {
        const double low = 100.0 - (i + 1) * percentile_width;
        const double high = 100.0 - i * percentile_width;

        auto format_percent = [](double value)
        {
            std::ostringstream out;
            if (std::fabs(value - std::round(value)) < 1e-6)
            {
                out << static_cast<int>(std::round(value));
            }
            else
            {
                out.setf(std::ios::fixed);
                out.precision(1);
                out << value;
            }
            return out.str();
        };

        labels.push_back(format_percent(low) + "-" + format_percent(high) + "%");
    }

    return labels;
}

std::vector<std::pair<int, int>> BuildPercentileIntervals(const std::vector<double> &boundaries)
{
    std::vector<std::pair<int, int>> intervals;
    const int nPercentileBins = static_cast<int>(boundaries.size()) + 1;
    if (nPercentileBins < 2)
    {
        return intervals;
    }

    std::vector<int> integer_boundaries;
    integer_boundaries.reserve(boundaries.size());
    for (double boundary : boundaries)
    {
        integer_boundaries.push_back(static_cast<int>(std::ceil(boundary)));
    }

    intervals.reserve(nPercentileBins);
    intervals.emplace_back(integer_boundaries.empty() ? 0 : integer_boundaries.back(), std::numeric_limits<int>::max());
    for (int i = static_cast<int>(integer_boundaries.size()) - 1; i > 0; --i)
    {
        intervals.emplace_back(integer_boundaries[i - 1], integer_boundaries[i]);
    }
    intervals.emplace_back(0, integer_boundaries.empty() ? std::numeric_limits<int>::max() : integer_boundaries.front());

    return intervals;
}

std::vector<std::string> BuildPercentileIntervalLabels(const std::vector<std::pair<int, int>> &intervals)
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

std::vector<double> BuildIntervalLabelPositions(const std::vector<double> &boundaries, double xMin, double xMax)
{
    std::vector<double> positions;
    if (xMax <= xMin)
    {
        return positions;
    }

    const double left_anchor = std::max(xMin, 0.0);
    if (boundaries.empty())
    {
        positions.push_back(0.5 * (left_anchor + xMax));
        return positions;
    }

    positions.reserve(boundaries.size() + 1);
    positions.push_back(0.5 * (boundaries.back() + xMax));
    for (int i = static_cast<int>(boundaries.size()) - 1; i > 0; --i)
    {
        positions.push_back(0.5 * (boundaries[i - 1] + boundaries[i]));
    }
    positions.push_back(0.5 * (left_anchor + boundaries.front()));

    return positions;
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

void Draw1DHistogramsWithPercentiles(const std::vector<TH1 *> &histograms,                                                             //
                                     const std::vector<std::vector<double>> &percentile_boundaries,                                    //
                                     const std::string plotinfo,                                                                       // single line info, as the legend header
                                     const std::vector<std::string> &hist_labels,                                                      //
                                     const std::string &x_axis_title,                                                                  //
                                     const std::string &outname,                                                                       //
                                     const std::vector<std::string> &colors = {"#3f90da", "#bd1f01", "#ffa90e", "#6CA651", "#832db6"}, //
                                     const std::vector<double> &label_boundaries = {},                                                 //
                                     const std::vector<std::string> &custom_boundary_labels = {},                                      //
                                     bool normalize = true,                                                                            //
                                     bool logy = true                                                                                  //
)
{
    if (histograms.empty())
    {
        std::cout << "[Draw1DHistogramsWithPercentiles] No histograms provided.\n";
        return;
    }
    if (percentile_boundaries.size() != histograms.size())
    {
        std::cout << "[Draw1DHistogramsWithPercentiles] Boundary set count does not match histogram count.\n";
        return;
    }

    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c_NInttClusterCrossing", "c_NInttClusterCrossing", 900, 700);
    c->cd();
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.08);
    gPad->SetTopMargin(0.08);
    if (logy)
    {
        c->SetLogy();
    }

    std::vector<TH1 *> draw_hists;
    draw_hists.reserve(histograms.size());

    double y_max = 0.0;
    double min_positive = 0.0;
    for (std::size_t i = 0; i < histograms.size(); ++i)
    {
        if (!histograms[i])
        {
            continue;
        }

        TH1 *hist = static_cast<TH1 *>(histograms[i]->Clone(Form("%s_draw_%zu", histograms[i]->GetName(), i)));
        hist->SetDirectory(nullptr);
        if (normalize && hist->Integral(0, hist->GetNbinsX() + 1) > 0.0)
        {
            hist->Scale(1.0 / hist->Integral(0, hist->GetNbinsX() + 1));
        }

        hist->SetLineColor(TColor::GetColor(colors[i % colors.size()].c_str()));
        hist->SetMarkerColor(TColor::GetColor(colors[i % colors.size()].c_str()));
        hist->SetLineWidth(2);
        hist->GetXaxis()->SetTitle(x_axis_title.c_str());
        hist->GetYaxis()->SetTitle(normalize ? "Normalized entries" : "Entries");

        y_max = std::max(y_max, hist->GetMaximum());
        const double hist_min_positive = hist->GetMinimum(0.0);
        if (hist_min_positive > 0.0)
        {
            min_positive = (min_positive > 0.0) ? std::min(min_positive, hist_min_positive) : hist_min_positive;
        }

        draw_hists.push_back(hist);
    }

    if (draw_hists.empty())
    {
        std::cout << "[Draw1DHistogramsWithPercentiles] All histograms were null.\n";
        delete c;
        return;
    }

    draw_hists.front()->GetYaxis()->SetRangeUser(logy ? std::max(min_positive * 0.5, 1e-6) : 0.0, y_max * (logy ? 15.0 : 1.3));
    draw_hists.front()->Draw("hist");
    for (std::size_t i = 1; i < draw_hists.size(); ++i)
    {
        draw_hists[i]->Draw("hist same");
    }

    TLegend *leg = new TLegend((1 - gPad->GetRightMargin()) - 0.3,                                                     //
                               1 - gPad->GetTopMargin() - 0.03 - 0.045 * (static_cast<double>(draw_hists.size()) + 1), //
                               (1 - gPad->GetRightMargin()) - 0.15,                                                    //
                               1 - gPad->GetTopMargin() - 0.03                                                         //
    );
    leg->SetHeader(plotinfo.c_str());
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);
    for (std::size_t i = 0; i < draw_hists.size(); ++i)
    {
        const std::string legend_label = (i < hist_labels.size() && !hist_labels[i].empty()) ? hist_labels[i] : draw_hists[i]->GetName();
        leg->AddEntry(draw_hists[i], legend_label.c_str(), "l");
    }
    leg->Draw();

    // Draw percentile boundary lines and labels, but only for the trigger crossing which is the first histogram in the list
    for (std::size_t i = 0; i < 1; ++i)
    {
        const auto &boundaries = percentile_boundaries[i];
        const double label_y = y_max * (logy ? (0.60 - 0.12 * static_cast<double>(i)) : (0.85 - 0.08 * static_cast<double>(i)));

        for (std::size_t j = 0; j < boundaries.size(); ++j)
        {
            TLine *line = new TLine(boundaries[j], draw_hists.front()->GetYaxis()->GetXmin(), boundaries[j], label_y);
            line->SetLineColor(TColor::GetColor(colors[i % colors.size()].c_str()));
            line->SetLineStyle(kDashed);
            line->SetLineWidth(2);
            line->Draw("same");
        }
    }

    const std::vector<double> &boundaries_for_labels = label_boundaries.empty() ? percentile_boundaries.front() : label_boundaries;
    const std::vector<std::string> labels_for_boundaries = custom_boundary_labels.empty() ? BuildPercentileBoundaryLabels(static_cast<int>(boundaries_for_labels.size()) + 1) : custom_boundary_labels;
    const std::vector<double> label_x_positions = (labels_for_boundaries.size() == boundaries_for_labels.size() + 1) ? BuildIntervalLabelPositions(boundaries_for_labels, draw_hists.front()->GetXaxis()->GetXmin(), draw_hists.front()->GetXaxis()->GetXmax()) : boundaries_for_labels;
    const double label_y = y_max * (logy ? 0.24 : 0.72);

    for (std::size_t i = 0; i < label_x_positions.size() && i < labels_for_boundaries.size(); ++i)
    {
        TLatex *label = new TLatex(label_x_positions[i], label_y, labels_for_boundaries[i].c_str());
        label->SetTextAngle(90);
        label->SetTextAlign(13);
        label->SetTextSize(0.025);
        label->SetTextColor(TColor::GetColor(colors[0].c_str()));
        label->Draw("same");
    }

    c->SaveAs(Form("%s.pdf", outname.c_str()));
    c->SaveAs(Form("%s.png", outname.c_str()));

    for (TH1 *hist : draw_hists)
    {
        delete hist;
    }
    delete leg;
    delete c;
}

// void NInttClusterCrossing(const std::string inputfile = "/sphenix/user/hjheng/sPHENIXRepo/TrackingAnalysis/Silicon_MBD_Vertexing/Silicon_MBD_Comparisons/VertexCompare_run_82405/files/outputVTX_Acts_Default.root", //
//                           const bool isMC = false,                                                                                                                                                                   //
//                           int nPercentileBins = 5                                                                                                                                                                    //
// )
void NInttClusterCrossing(const std::string inputfile = "/sphenix/tg/tg01/hf/hjheng/ppg-dNdEta-OOpp/simulation-nopileup-ntuple/ACTS/ntuple_*.root", //
                          const bool isMC = true,                                                                                                                                                                   //
                          int nPercentileBins = 5                                                                                                                                                                    //
)
{
    if (inputfile.empty())
    {
        std::cout << "[NInttClusterCrossing] Please provide an input ROOT file.\n";
        return;
    }

    std::string plotdir = "./figure/figure-NInttClusterCrossing";
    system(("mkdir -p " + plotdir).c_str());
    const std::string outputPrefix = isMC ? "MC_" : "Data_";

    std::vector<int> NInttCluster_TrgCrossing;
    std::vector<int> NInttCluster_NonTrgCrossing;

    TH1 *hM_NInttCluster_TrgCrossing = new TH1D("hM_NInttCluster_TrgCrossing", "N_{INTT} clusters in trigger crossing;N_{INTT} clusters;Counts", 200, -0.5, 999.5);
    TH1 *hM_NInttCluster_NonTrgCrossing = new TH1D("hM_NInttCluster_NonTrgCrossing", "N_{INTT} clusters in non-trigger crossing;N_{INTT} clusters;Counts", 200, -0.5, 999.5);

    TChain *t = new TChain("VTX");
    const auto inputFiles = ResolveInputFiles(inputfile);
    if (inputFiles.empty())
    {
        std::cout << "[NInttClusterCrossing] No files matched: " << inputfile << "\n";
        delete t;
        return;
    }

    const std::size_t maxFiles = isMC ? 100 : inputFiles.size();
    int nAddedFiles = 0;
    for (std::size_t i = 0; i < inputFiles.size() && i < maxFiles; ++i)
    {
        if (t->Add(inputFiles[i].c_str()) > 0)
        {
            ++nAddedFiles;
        }
    }

    if (isMC && inputFiles.size() > maxFiles)
    {
        std::cout << "[NInttClusterCrossing] MC mode: using first " << maxFiles << " files out of " << inputFiles.size() << " matched files.\n";
    }

    if (nAddedFiles == 0)
    {
        std::cout << "[NInttClusterCrossing] Failed to add input files to TChain.\n";
        delete t;
        return;
    }

    std::cout << "[NInttClusterCrossing] Added " << nAddedFiles << " files to TChain.\n";

    int counter = 0;
    std::vector<int> *firedTriggers = nullptr;
    std::vector<unsigned int> *cluster_layer = nullptr;
    std::vector<int> *cluster_timeBucketID = nullptr;
    t->SetBranchAddress("counter", &counter);
    t->SetBranchAddress("firedTriggers", &firedTriggers);
    t->SetBranchAddress("cluster_layer", &cluster_layer);
    t->SetBranchAddress("cluster_timeBucketID", &cluster_timeBucketID);

    const Long64_t nEntries = t->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        t->GetEntry(i);

        std::set<int> unique_time_buckets(cluster_timeBucketID->begin(), cluster_timeBucketID->end());
        for (int time_bucket : unique_time_buckets)
        {
            int n_clusters_in_bucket = 0;
            for (std::size_t j = 0; j < cluster_timeBucketID->size(); ++j)
            {
                if (cluster_timeBucketID->at(j) == time_bucket)
                {
                    ++n_clusters_in_bucket;
                }
            }

            if (time_bucket == 0)
            {
                hM_NInttCluster_TrgCrossing->Fill(n_clusters_in_bucket);
                NInttCluster_TrgCrossing.push_back(n_clusters_in_bucket);
            }
            else
            {
                hM_NInttCluster_NonTrgCrossing->Fill(n_clusters_in_bucket);
                NInttCluster_NonTrgCrossing.push_back(n_clusters_in_bucket);
            }
        }
    }

    const auto trg_percentile_boundaries = DerivePercentileBoundaries(NInttCluster_TrgCrossing, nPercentileBins);
    const auto nontrg_percentile_boundaries = DerivePercentileBoundaries(NInttCluster_NonTrgCrossing, nPercentileBins);
    const auto trg_percentile_intervals = BuildPercentileIntervals(trg_percentile_boundaries);
    const auto trg_percentile_interval_labels = BuildPercentileIntervalLabels(trg_percentile_intervals);

    std::cout << "[NInttClusterCrossing] Trigger crossing percentile boundaries:";
    for (double boundary : trg_percentile_boundaries)
    {
        std::cout << " " << boundary;
    }
    std::cout << "\n";

    std::cout << "[NInttClusterCrossing] Non-trigger crossing percentile boundaries:";
    for (double boundary : nontrg_percentile_boundaries)
    {
        std::cout << " " << boundary;
    }
    std::cout << "\n";

    std::cout << "[NInttClusterCrossing] Trigger crossing percentile intervals:\n";
    for (std::size_t i = 0; i < trg_percentile_interval_labels.size(); ++i)
    {
        std::cout << "  bin " << i << ": " << trg_percentile_interval_labels[i] << "\n";
    }

    Draw1DHistogramsWithPercentiles({hM_NInttCluster_TrgCrossing, hM_NInttCluster_NonTrgCrossing},                                                          //
                                    {trg_percentile_boundaries, nontrg_percentile_boundaries},                                                              //
                                    (isMC ? "Simulation" : "Data, O+O 82405"),                                                                              //
                                    {"Trigger crossing", "Non-trigger crossing"},                                                                           //
                                    "N_{INTT clusters}",                                                                                                     //
                                    Form("%s/%shM_NInttCluster_CrossingComparison_wPercentile", plotdir.c_str(), outputPrefix.c_str()),                    //
                                    {"#3f90da", "#bd1f01"},                                                                                                //
                                    trg_percentile_boundaries, BuildPercentileBoundaryLabels(nPercentileBins)                                               //
    );

    TFile *fout = TFile::Open(Form("%s/%sNInttClusterPercentileBoundaries.root", plotdir.c_str(), outputPrefix.c_str()), "RECREATE");
    SavePercentileIntervals(fout, "trigger_percentile_intervals", trg_percentile_intervals, nPercentileBins);
    fout->Close();
    delete fout;

    delete t;
    delete hM_NInttCluster_TrgCrossing;
    delete hM_NInttCluster_NonTrgCrossing;
}
