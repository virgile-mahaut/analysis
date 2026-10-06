#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH2F.h"
#include "TString.h"

struct SeedDisplayTrack
{
    TGraph *clus_xy = nullptr;
    TGraph *clus_zr = nullptr;
    TGraph *pca_xy = nullptr;
    TGraph *pca_zr = nullptr;
};

static bool GraphHasPoints(TGraph *graph) { return graph && graph->GetN() > 0; }

static TGraph *FirstNonEmptyGraph(const std::vector<TGraph *> &graphs)
{
    for (TGraph *graph : graphs)
    {
        if (GraphHasPoints(graph))
            return graph;
    }
    return nullptr;
}

static std::vector<std::string> CollectGraphNames(TDirectory *dir, const std::string &prefix)
{
    std::vector<std::string> names;
    if (!dir)
        return names;

    TIter next(dir->GetListOfKeys());
    while (TKey *key = static_cast<TKey *>(next()))
    {
        const std::string key_name = key->GetName();
        if (key_name.find(prefix) == 0)
            names.push_back(key_name);
    }

    std::sort(names.begin(), names.end());
    return names;
}

static std::string BuildCrossingLabel(const std::string &crossing_key, bool is_simulation)
{
    const size_t split_pos = crossing_key.rfind('_');
    if (split_pos == std::string::npos)
        return crossing_key;

    const std::string id_str = crossing_key.substr(0, split_pos);
    const std::string crossing_str = crossing_key.substr(split_pos + 1);

    try
    {
        const int crossing = std::stoi(crossing_str);
        if (is_simulation)
        {
            const int counter = std::stoi(id_str);
            return Form("Simulation, counter=%d, crossing=%d", counter, crossing);
        }

        const long long gl1bco = std::stoll(id_str);
        return Form("Data, GL1 BCO=%lld, crossing=%d", gl1bco, crossing);
    }
    catch (const std::exception &)
    {
        return crossing_key;
    }
}

static void DrawSeedClusterDisplay(TGraph *gr_clus_xy_all,                               //
                                   TGraph *gr_clus_zr_all,                               //
                                   TGraph *gr_primary_intt_xy,                           //
                                   TGraph *gr_primary_intt_zr,                           //
                                   TGraph *gr_nonprimary_intt_xy,                        //
                                   TGraph *gr_nonprimary_intt_zr,                        //
                                   TGraph *gr_pca_xy_all,                                //
                                   TGraph *gr_pca_zr_all,                                //
                                   TGraph *gr_reco_vertex_xy,                            //
                                   TGraph *gr_reco_vertex_zr,                            //
                                   TGraph *gr_truth_vertex_xy,                           //
                                   TGraph *gr_truth_vertex_zr,                           //
                                   const std::vector<SeedDisplayTrack> &tracks,          //
                                   const std::vector<TGraph *> &nonassociated_xy_graphs, //
                                   const std::vector<TGraph *> &nonassociated_zr_graphs, //
                                   const std::vector<TGraph *> &truth_xy_graphs,         //
                                   const std::vector<TGraph *> &truth_zr_graphs,         //
                                   const std::vector<TGraph *> &reco_xy_graphs,          //
                                   const std::vector<TGraph *> &reco_zr_graphs,          //
                                   const std::vector<std::string> &info,                 //
                                   const std::string &outname,                           //
                                   const bool is_simulation_flag                         //
)
{
    if (!is_simulation_flag && (!gr_clus_xy_all || !gr_clus_zr_all))
        return;
    if (is_simulation_flag && !GraphHasPoints(gr_primary_intt_xy) && !GraphHasPoints(gr_nonprimary_intt_xy) && !GraphHasPoints(gr_primary_intt_zr) && !GraphHasPoints(gr_nonprimary_intt_zr))
        return;

    TCanvas *c1 = new TCanvas("c_seeddisplay", "Seed-cluster display", 1150, 700);
    c1->Divide(2, 1);

    c1->cd(1);
    gPad->SetTopMargin(0.23);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    TH2F *frame_xy = new TH2F("frame_xy", ";x [cm];y [cm]", 10, -12.0, 12.0, 10, -12.0, 12.0);
    frame_xy->SetStats(0);
    frame_xy->GetXaxis()->SetTitleOffset(1.2);
    frame_xy->GetYaxis()->SetTitleOffset(1.4);
    frame_xy->Draw();

    if (is_simulation_flag)
    {
        for (size_t i = 0; i < truth_xy_graphs.size(); i++)
        {
            TGraph *gr_truth_xy = truth_xy_graphs[i];
            if (GraphHasPoints(gr_truth_xy))
            {
                gr_truth_xy->SetMarkerStyle(25);
                gr_truth_xy->SetMarkerSize(0.8);
                gr_truth_xy->SetMarkerColor(TColor::GetColor("#63C8FF"));
                gr_truth_xy->SetLineColor(TColor::GetColor("#63C8FF"));
                gr_truth_xy->SetLineStyle(3);
                gr_truth_xy->SetLineWidth(1);
                gr_truth_xy->Draw("PL SAME");
            }
        }
        for (size_t i = 0; i < reco_xy_graphs.size(); i++)
        {
            TGraph *gr_reco_xy = reco_xy_graphs[i];
            if (GraphHasPoints(gr_reco_xy))
            {
                gr_reco_xy->SetMarkerStyle(21);
                gr_reco_xy->SetMarkerSize(0.7);
                gr_reco_xy->SetMarkerColor(TColor::GetColor("#001BB7"));
                gr_reco_xy->SetLineColor(TColor::GetColor("#001BB7"));
                gr_reco_xy->SetLineStyle(3);
                gr_reco_xy->SetLineWidth(1);
                gr_reco_xy->Draw("PL SAME");
            }
        }
    }

    if (is_simulation_flag)
    {
        if (GraphHasPoints(gr_nonprimary_intt_xy))
        {
            gr_nonprimary_intt_xy->SetMarkerStyle(24);
            gr_nonprimary_intt_xy->SetMarkerSize(0.8);
            gr_nonprimary_intt_xy->SetMarkerColor(kBlack);
            gr_nonprimary_intt_xy->Draw("P SAME");
        }
        if (GraphHasPoints(gr_primary_intt_xy))
        {
            gr_primary_intt_xy->SetMarkerStyle(25);
            gr_primary_intt_xy->SetMarkerSize(1.4);
            gr_primary_intt_xy->SetMarkerColorAlpha(kBlack, 1.0);
            gr_primary_intt_xy->Draw("P SAME");
        }
    }
    else
    {
        gr_clus_xy_all->SetMarkerStyle(20);
        gr_clus_xy_all->SetMarkerSize(0.8);
        gr_clus_xy_all->SetMarkerColorAlpha(kBlack, 0.3);
        gr_clus_xy_all->Draw("P SAME");
    }

    for (const auto &track : tracks)
    {
        if (GraphHasPoints(track.clus_xy))
        {
            track.clus_xy->SetMarkerStyle(20);
            track.clus_xy->SetMarkerSize(0.7);
            track.clus_xy->SetMarkerColor(TColor::GetColor("#cc3311"));
            track.clus_xy->SetLineColor(TColor::GetColor("#cc3311"));
            track.clus_xy->SetLineStyle(3);
            track.clus_xy->SetLineWidth(1);
            track.clus_xy->Draw("PL SAME");
        }
    }
    if (GraphHasPoints(gr_pca_xy_all))
    {
        gr_pca_xy_all->SetMarkerStyle(47);
        gr_pca_xy_all->SetMarkerSize(0.6);
        gr_pca_xy_all->SetMarkerColor(TColor::GetColor("#ee8866"));
        gr_pca_xy_all->Draw("P SAME");
    }
    for (size_t i = 0; i < nonassociated_xy_graphs.size(); i++)
    {
        TGraph *gr_nonassociated_xy = nonassociated_xy_graphs[i];
        if (GraphHasPoints(gr_nonassociated_xy))
        {
            gr_nonassociated_xy->SetMarkerStyle(24);
            gr_nonassociated_xy->SetMarkerSize(0.7);
            gr_nonassociated_xy->SetMarkerColor(TColor::GetColor("#FF88BA"));
            gr_nonassociated_xy->SetLineColor(TColor::GetColor("#FF88BA"));
            gr_nonassociated_xy->SetLineStyle(3);
            gr_nonassociated_xy->SetLineWidth(1);
            gr_nonassociated_xy->Draw("PL SAME");
        }
    }

    if (GraphHasPoints(gr_reco_vertex_xy))
    {
        gr_reco_vertex_xy->SetMarkerStyle(34);
        gr_reco_vertex_xy->SetMarkerSize(1.2);
        gr_reco_vertex_xy->SetMarkerColorAlpha(TColor::GetColor("#346739"), 1.0);
        gr_reco_vertex_xy->Draw("P SAME");
    }

    if (is_simulation_flag && GraphHasPoints(gr_truth_vertex_xy))
    {
        gr_truth_vertex_xy->SetMarkerStyle(28);
        gr_truth_vertex_xy->SetMarkerSize(1.2);
        gr_truth_vertex_xy->SetMarkerColorAlpha(TColor::GetColor("#91D06C"), 1.0);
        gr_truth_vertex_xy->Draw("P SAME");
    }

    c1->cd(2);
    gPad->SetTopMargin(0.23);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    TH2F *frame_zr = new TH2F("frame_zr", ";z [cm];r [cm]", 10, -29.0, 29.0, 10, -13.0, 13.0);
    frame_zr->SetStats(0);
    frame_zr->GetXaxis()->SetTitleOffset(1.2);
    frame_zr->GetYaxis()->SetTitleOffset(1.4);
    frame_zr->Draw();

    if (is_simulation_flag)
    {
        for (size_t i = 0; i < truth_zr_graphs.size(); i++)
        {
            TGraph *gr_truth_zr = truth_zr_graphs[i];
            if (GraphHasPoints(gr_truth_zr))
            {
                gr_truth_zr->SetMarkerStyle(25);
                gr_truth_zr->SetMarkerSize(0.8);
                gr_truth_zr->SetMarkerColor(TColor::GetColor("#63C8FF"));
                gr_truth_zr->SetLineColor(TColor::GetColor("#63C8FF"));
                gr_truth_zr->SetLineStyle(3);
                gr_truth_zr->SetLineWidth(1);
                gr_truth_zr->Draw("PL SAME");
            }
        }
        for (size_t i = 0; i < reco_zr_graphs.size(); i++)
        {
            TGraph *gr_reco_zr = reco_zr_graphs[i];
            if (GraphHasPoints(gr_reco_zr))
            {
                gr_reco_zr->SetMarkerStyle(21);
                gr_reco_zr->SetMarkerSize(0.7);
                gr_reco_zr->SetMarkerColor(TColor::GetColor("#001BB7"));
                gr_reco_zr->SetLineColor(TColor::GetColor("#001BB7"));
                gr_reco_zr->SetLineStyle(3);
                gr_reco_zr->SetLineWidth(1);
                gr_reco_zr->Draw("PL SAME");
            }
        }
    }

    if (is_simulation_flag)
    {
        if (GraphHasPoints(gr_nonprimary_intt_zr))
        {
            gr_nonprimary_intt_zr->SetMarkerStyle(24);
            gr_nonprimary_intt_zr->SetMarkerSize(0.8);
            gr_nonprimary_intt_zr->SetMarkerColor(kBlack);
            gr_nonprimary_intt_zr->Draw("P SAME");
        }
        if (GraphHasPoints(gr_primary_intt_zr))
        {
            gr_primary_intt_zr->SetMarkerStyle(25);
            gr_primary_intt_zr->SetMarkerSize(1.4);
            gr_primary_intt_zr->SetMarkerColorAlpha(kBlack, 1.0);
            gr_primary_intt_zr->Draw("P SAME");
        }
    }
    else
    {
        gr_clus_zr_all->SetMarkerStyle(20);
        gr_clus_zr_all->SetMarkerSize(0.8);
        gr_clus_zr_all->SetMarkerColorAlpha(kBlack, 0.3);
        gr_clus_zr_all->Draw("P SAME");
    }

    for (const auto &track : tracks)
    {
        if (GraphHasPoints(track.clus_zr))
        {
            track.clus_zr->SetMarkerStyle(20);
            track.clus_zr->SetMarkerSize(0.7);
            track.clus_zr->SetMarkerColor(TColor::GetColor("#cc3311"));
            track.clus_zr->SetLineColor(TColor::GetColor("#cc3311"));
            track.clus_zr->SetLineStyle(3);
            track.clus_zr->SetLineWidth(1);
            track.clus_zr->Draw("PL SAME");
        }
    }
    if (GraphHasPoints(gr_pca_zr_all))
    {
        gr_pca_zr_all->SetMarkerStyle(47);
        gr_pca_zr_all->SetMarkerSize(0.6);
        gr_pca_zr_all->SetMarkerColor(TColor::GetColor("#ee8866"));
        gr_pca_zr_all->Draw("P SAME");
    }
    for (size_t i = 0; i < nonassociated_zr_graphs.size(); i++)
    {
        TGraph *gr_nonassociated_zr = nonassociated_zr_graphs[i];
        if (GraphHasPoints(gr_nonassociated_zr))
        {
            gr_nonassociated_zr->SetMarkerStyle(24);
            gr_nonassociated_zr->SetMarkerSize(0.7);
            gr_nonassociated_zr->SetMarkerColor(TColor::GetColor("#FF88BA"));
            gr_nonassociated_zr->SetLineColor(TColor::GetColor("#FF88BA"));
            gr_nonassociated_zr->SetLineStyle(3);
            gr_nonassociated_zr->SetLineWidth(1);
            gr_nonassociated_zr->Draw("PL SAME");
        }
    }

    if (GraphHasPoints(gr_reco_vertex_zr))
    {
        gr_reco_vertex_zr->SetMarkerStyle(34);
        gr_reco_vertex_zr->SetMarkerSize(1.2);
        gr_reco_vertex_zr->SetMarkerColorAlpha(TColor::GetColor("#346739"), 1.0);
        gr_reco_vertex_zr->Draw("P SAME");
    }

    if (is_simulation_flag && GraphHasPoints(gr_truth_vertex_zr))
    {
        gr_truth_vertex_zr->SetMarkerStyle(28);
        gr_truth_vertex_zr->SetMarkerSize(1.2);
        gr_truth_vertex_zr->SetMarkerColorAlpha(TColor::GetColor("#91D06C"), 1.0);
        gr_truth_vertex_zr->Draw("P SAME");
    }

    c1->cd(1);
    TLegend *legend_xy = new TLegend(gPad->GetLeftMargin(), 0.99, gPad->GetLeftMargin() + 0.29, 1 - gPad->GetTopMargin() + 0.01);
    for (const auto &line : info)
        legend_xy->AddEntry((TObject *)nullptr, line.c_str(), "");
    if (is_simulation_flag)
    {
        if (GraphHasPoints(gr_primary_intt_xy))
            legend_xy->AddEntry(gr_primary_intt_xy, "INTT clusters from sPHENIX primary particles", "p");
        if (GraphHasPoints(gr_nonprimary_intt_xy))
            legend_xy->AddEntry(gr_nonprimary_intt_xy, "INTT clusters not from sPHENIX primary particles", "p");
    }
    else
    {
        legend_xy->AddEntry(gr_clus_xy_all, "All INTT clusters of the crossing", "p");
    }
    if (!tracks.empty() && GraphHasPoints(tracks.front().clus_xy))
    {
        legend_xy->AddEntry(tracks.front().clus_xy, "Seed clusters of the crossing assoc. to vertex", "p");
    }
    if (TGraph *nonassociated_legend_graph = FirstNonEmptyGraph(nonassociated_xy_graphs))
    {
        legend_xy->AddEntry(nonassociated_legend_graph, "Seed clusters of the crossing not associated to vertex", "p");
    }
    // if (GraphHasPoints(gr_pca_xy_all))
    // {
    //     legend_xy->AddEntry(gr_pca_xy_all, "Seed PCA of the crossing", "p");
    // }
    legend_xy->SetTextSize(0.03);
    legend_xy->SetFillStyle(0);
    legend_xy->SetFillColor(0);
    legend_xy->SetBorderSize(0);
    legend_xy->Draw();

    c1->cd(2);
    TLegend *legend_zr = new TLegend(gPad->GetLeftMargin(), 0.99, gPad->GetLeftMargin() + 0.30, 1 - gPad->GetTopMargin() + 0.01);
    // for (const auto &line : info)
    legend_zr->AddEntry((TObject *)nullptr, "", "");
    if (GraphHasPoints(gr_pca_xy_all))
    {
        legend_zr->AddEntry(gr_pca_xy_all, "Seed PCA of the crossing", "p");
    }
    if (GraphHasPoints(gr_reco_vertex_xy))
    {
        legend_zr->AddEntry(gr_reco_vertex_xy, "Reconstructed silicon vertex", "p");
    }
    if (is_simulation_flag && GraphHasPoints(gr_truth_vertex_xy))
    {
        legend_zr->AddEntry(gr_truth_vertex_xy, "Truth vertex", "p");
    }
    if (is_simulation_flag)
    {
        if (TGraph *truth_legend_graph = FirstNonEmptyGraph(truth_xy_graphs))
        {
            legend_zr->AddEntry(truth_legend_graph, "Truth clusters of sPHENIX primary charged particles", "p");
        }
        if (TGraph *reco_legend_graph = FirstNonEmptyGraph(reco_xy_graphs))
        {
            legend_zr->AddEntry(reco_legend_graph, "Reco. clusters of sPHENIX primary charged particles", "p");
        }
    }
    legend_zr->SetTextSize(0.03);
    legend_zr->SetFillStyle(0);
    legend_zr->SetFillColor(0);
    legend_zr->SetBorderSize(0);
    legend_zr->Draw();

    c1->SaveAs((outname + ".png").c_str());
    c1->SaveAs((outname + ".pdf").c_str());

    delete legend_xy;
    delete legend_zr;
    delete frame_xy;
    delete frame_zr;
    delete c1;
}

// void plot_seeddisplay(const std::string &inputfile = "./test_OO_combinatoric.root",    //
//                       const bool is_simulation_flag = false,                           //
//                       const std::string &plotdir_base = "./figure/figure-seeddisplay", //
//                       const std::string &crossing_filter = ""                          //
// )
void plot_seeddisplay(const std::string &inputfile = "/sphenix/tg/tg01/hf/hjheng/ppg-dNdEta-OOpp/TrackletAna-HistOutput/simulation/ACTS/cotThetaMax4p0//hist_0.root", //
                      const bool is_simulation_flag = true,                                                                                                           //
                      const std::string &plotdir_base = "./figure/figure-seeddisplay-cotThetaMax4p0",                                                                 //
                      const std::string &crossing_filter = ""                                                                                                         //
)
{
    gStyle->SetOptStat(0);

    const std::string mode_dir = is_simulation_flag ? "simulation" : "data";
    const std::string plotdir = plotdir_base + "/" + mode_dir;
    gSystem->mkdir(plotdir.c_str(), true);

    TFile *fin = TFile::Open(inputfile.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "[plot_seeddisplay] Cannot open: " << inputfile << std::endl;
        return;
    }

    TDirectory *crossing_display_dir = fin->GetDirectory("CrossingDisplays");
    if (!crossing_display_dir)
    {
        std::cerr << "[plot_seeddisplay] Missing CrossingDisplays directory in " << inputfile << std::endl;
        fin->Close();
        return;
    }

    int n_drawn = 0;
    TIter next_crossing(crossing_display_dir->GetListOfKeys());
    while (TKey *key = static_cast<TKey *>(next_crossing()))
    {
        TObject *obj = key->ReadObj();
        TDirectory *crossing_dir = dynamic_cast<TDirectory *>(obj);
        if (!crossing_dir)
            continue;

        const std::string crossing_key = crossing_dir->GetName();
        if (!crossing_filter.empty() && crossing_key != crossing_filter)
            continue;

        TGraph *gr_all_xy = dynamic_cast<TGraph *>(crossing_dir->Get("gr_all_intt_clusters_xy"));
        TGraph *gr_all_zr = dynamic_cast<TGraph *>(crossing_dir->Get("gr_all_intt_clusters_zr"));
        TGraph *gr_primary_intt_xy = dynamic_cast<TGraph *>(crossing_dir->Get("gr_primary_intt_clusters_xy"));
        TGraph *gr_primary_intt_zr = dynamic_cast<TGraph *>(crossing_dir->Get("gr_primary_intt_clusters_zr"));
        TGraph *gr_nonprimary_intt_xy = dynamic_cast<TGraph *>(crossing_dir->Get("gr_nonprimary_intt_clusters_xy"));
        TGraph *gr_nonprimary_intt_zr = dynamic_cast<TGraph *>(crossing_dir->Get("gr_nonprimary_intt_clusters_zr"));
        TGraph *gr_pca_xy = dynamic_cast<TGraph *>(crossing_dir->Get("gr_associated_silseed_pca_xy"));
        TGraph *gr_pca_zr = dynamic_cast<TGraph *>(crossing_dir->Get("gr_associated_silseed_pca_zr"));
        TGraph *gr_reco_vertex_xy = dynamic_cast<TGraph *>(crossing_dir->Get("gr_reco_vertex_xy"));
        TGraph *gr_reco_vertex_zr = dynamic_cast<TGraph *>(crossing_dir->Get("gr_reco_vertex_zr"));
        TGraph *gr_truth_vertex_xy = dynamic_cast<TGraph *>(crossing_dir->Get("gr_truth_vertex_xy"));
        TGraph *gr_truth_vertex_zr = dynamic_cast<TGraph *>(crossing_dir->Get("gr_truth_vertex_zr"));

        if (!is_simulation_flag && (!gr_all_xy || !gr_all_zr))
        {
            std::cerr << "[plot_seeddisplay] Skipping " << crossing_key << " because required all-cluster graphs are missing." << std::endl;
            continue;
        }
        if (is_simulation_flag && !gr_primary_intt_xy && !gr_primary_intt_zr && !gr_nonprimary_intt_xy && !gr_nonprimary_intt_zr)
        {
            std::cerr << "[plot_seeddisplay] Skipping " << crossing_key << " because required primary/non-primary INTT cluster graphs are missing." << std::endl;
            continue;
        }

        std::vector<std::string> cluster_xy_names = CollectGraphNames(crossing_dir, "gr_associated_silseed_clusters_xy_");
        std::vector<std::string> cluster_zr_names = CollectGraphNames(crossing_dir, "gr_associated_silseed_clusters_zr_");
        std::vector<std::string> nonassociated_xy_names = CollectGraphNames(crossing_dir, "gr_nonassociated_silseed_clusters_xy_");
        std::vector<std::string> nonassociated_zr_names = CollectGraphNames(crossing_dir, "gr_nonassociated_silseed_clusters_zr_");
        std::vector<std::string> truth_xy_names;
        std::vector<std::string> truth_zr_names;
        std::vector<std::string> reco_xy_names;
        std::vector<std::string> reco_zr_names;
        if (is_simulation_flag)
        {
            truth_xy_names = CollectGraphNames(crossing_dir, "gr_sPHENIXPrimary_truthcluster_xy_");
            truth_zr_names = CollectGraphNames(crossing_dir, "gr_sPHENIXPrimary_truthcluster_zr_");
            reco_xy_names = CollectGraphNames(crossing_dir, "gr_sPHENIXPrimary_recocluster_xy_");
            reco_zr_names = CollectGraphNames(crossing_dir, "gr_sPHENIXPrimary_recocluster_zr_");
        }
        const size_t ntracks = std::min(cluster_xy_names.size(), cluster_zr_names.size());
        const size_t ntruth = std::min(truth_xy_names.size(), truth_zr_names.size());
        const size_t nreco = std::min(reco_xy_names.size(), reco_zr_names.size());

        std::vector<SeedDisplayTrack> tracks;
        tracks.reserve(ntracks);
        for (size_t i = 0; i < ntracks; ++i)
        {
            SeedDisplayTrack track;
            track.clus_xy = dynamic_cast<TGraph *>(crossing_dir->Get(cluster_xy_names[i].c_str()));
            track.clus_zr = dynamic_cast<TGraph *>(crossing_dir->Get(cluster_zr_names[i].c_str()));
            tracks.push_back(track);
        }

        std::vector<TGraph *> nonassociated_xy_graphs;
        std::vector<TGraph *> nonassociated_zr_graphs;
        const size_t nnonassociated = std::min(nonassociated_xy_names.size(), nonassociated_zr_names.size());
        nonassociated_xy_graphs.reserve(nnonassociated);
        nonassociated_zr_graphs.reserve(nnonassociated);
        for (size_t i = 0; i < nnonassociated; ++i)
        {
            nonassociated_xy_graphs.push_back(dynamic_cast<TGraph *>(crossing_dir->Get(nonassociated_xy_names[i].c_str())));
            nonassociated_zr_graphs.push_back(dynamic_cast<TGraph *>(crossing_dir->Get(nonassociated_zr_names[i].c_str())));
        }

        std::vector<TGraph *> truth_xy_graphs;
        std::vector<TGraph *> truth_zr_graphs;
        truth_xy_graphs.reserve(ntruth);
        truth_zr_graphs.reserve(ntruth);
        for (size_t i = 0; i < ntruth; ++i)
        {
            truth_xy_graphs.push_back(dynamic_cast<TGraph *>(crossing_dir->Get(truth_xy_names[i].c_str())));
            truth_zr_graphs.push_back(dynamic_cast<TGraph *>(crossing_dir->Get(truth_zr_names[i].c_str())));
        }

        std::vector<TGraph *> reco_xy_graphs;
        std::vector<TGraph *> reco_zr_graphs;
        reco_xy_graphs.reserve(nreco);
        reco_zr_graphs.reserve(nreco);
        for (size_t i = 0; i < nreco; ++i)
        {
            reco_xy_graphs.push_back(dynamic_cast<TGraph *>(crossing_dir->Get(reco_xy_names[i].c_str())));
            reco_zr_graphs.push_back(dynamic_cast<TGraph *>(crossing_dir->Get(reco_zr_names[i].c_str())));
        }

        const int n_intt_clusters = is_simulation_flag ? ((gr_primary_intt_xy ? gr_primary_intt_xy->GetN() : 0) + (gr_nonprimary_intt_xy ? gr_nonprimary_intt_xy->GetN() : 0)) : (gr_all_xy ? gr_all_xy->GetN() : 0);

        std::vector<std::string> info;
        info.push_back(BuildCrossingLabel(crossing_key, is_simulation_flag));
        info.push_back(Form("N_{INTT clusters}=%d, N_{associated seeds}=%d", n_intt_clusters, gr_pca_xy ? gr_pca_xy->GetN() : static_cast<int>(ntracks)));

        const std::string outname = plotdir + "/display_" + crossing_key;
        DrawSeedClusterDisplay(gr_all_xy, gr_all_zr, gr_primary_intt_xy, gr_primary_intt_zr, gr_nonprimary_intt_xy, gr_nonprimary_intt_zr, gr_pca_xy, gr_pca_zr, gr_reco_vertex_xy, gr_reco_vertex_zr, gr_truth_vertex_xy, gr_truth_vertex_zr, tracks, nonassociated_xy_graphs, nonassociated_zr_graphs, truth_xy_graphs, truth_zr_graphs, reco_xy_graphs, reco_zr_graphs, info, outname, is_simulation_flag);
        ++n_drawn;
    }

    std::cout << "[plot_seeddisplay] Wrote " << n_drawn << " display plots to " << plotdir << std::endl;
    fin->Close();

    // go to plotdir and merge individual PDFs into one
    system(Form("cd %s && gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=merge-display.pdf -dBATCH display_*.pdf", plotdir.c_str()));
}