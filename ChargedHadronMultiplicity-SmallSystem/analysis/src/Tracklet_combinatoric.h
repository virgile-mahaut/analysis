#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TSystem.h>
#include <TTree.h>

#include "GenHadron.h"
#include "Tracklet.h"
#include "Hit.h"
#include "LightHelixFit.h"
#include "Utilities.h"

static inline int logicalInttLayerIndex(int physical_layer)
{
    if (physical_layer == 3 || physical_layer == 4)
        return 0; // seed layer
    if (physical_layer == 5 || physical_layer == 6)
        return 1; // search layer
    return -1;
}

// static inline double wrapPhi_mpi_pi(double phi)
// {
//     const double pi = TMath::Pi();
//     const double twopi = 2.0 * pi;
//     while (phi <= -pi)
//         phi += twopi;
//     while (phi > pi)
//         phi -= twopi;
//     return phi;
// }

// static inline double wrapPhi_0_2pi(double phi)
// {
//     const double twopi = 2.0 * TMath::Pi();
//     while (phi < 0)
//         phi += twopi;
//     while (phi >= twopi)
//         phi -= twopi;
//     return phi;
// }

// static inline double rotatePhiByPi(double phi, bool phiAxisIs_0_2pi)
// {
//     // PHOBOS rotated-layer technique: phi -> phi + pi (with wrap)
//     const double phi_rot = phi + TMath::Pi();
//     return phiAxisIs_0_2pi ? wrapPhi_0_2pi(phi_rot) : wrapPhi_mpi_pi(phi_rot);
// }

// Binning helper: support uniform or explicit edges (both axes must match mode)
struct AxisBinning
{
    bool use_edges = false;
    std::vector<double> edges; // size = nbins+1 when use_edges=true
    double min = 0.0;
    double max = 0.0;
    int nbins = 0;

    AxisBinning() {}

    AxisBinning(double mn, double mx, int n)
        : use_edges(false)
        , min(mn)
        , max(mx)
        , nbins(n)
    {
    }

    AxisBinning(const std::vector<double> &e)
        : use_edges(true)
        , edges(e)
    {
        nbins = (edges.size() >= 2) ? int(edges.size()) - 1 : 0;
        min = (edges.size() ? edges.front() : 0.0);
        max = (edges.size() ? edges.back() : 0.0);
    }
};

struct TrackletCombiDRConfig
{
    enum class PairCutMode
    {
        kDeltaR,
        kDeltaPhi,
        kAbsDeltaPhi
    };

    // Pair cut mode and threshold.
    // kDeltaPhi mode still uses a symmetric threshold (-pair_cut, +pair_cut).
    PairCutMode pair_cut_mode = PairCutMode::kDeltaPhi;
    double pair_cut = 0.016; // default PHOBOS-like scale; tune as needed

    // Output grid binning (\eta, \phi) for uncorrected counts after subtraction
    AxisBinning etaAxis = AxisBinning(-1.0, 1.0, 20);
    AxisBinning phiAxis = AxisBinning(-3.2, 3.2, 32); // radians

    // Per-cell pair-variable histogram settings
    int dR_nbins = 200;
    double dR_min = 0.0;
    double dR_max = 0.5;

    // Performance: reduce the all-pairs scan using an η window (only meaningful for kDeltaR mode)
    bool enable_eta_preselect = true;

    // Rotate which logical layer for background estimation:
    // True = rotate search layer (logical 1) by pi; False = rotate seed layer (logical 0).
    bool rotate_search_layer = true;

    // Background scaling (keep configurable even if default is 1.0)
    double bkg_scale = 1.0;
    // Independent |delta-phi| counting threshold for a directly mergeable background-subtracted
    // (eta, phi) yield map. This is filled regardless of the main pair_cut_mode selection.
    double abs_dphi_count_cut = 0.15;

    // Optional z vertex cut (applied to silicon vertex z); enabled by default with +-10 cm window for silicon
    bool enable_zvtx_cut = true;
    double zvtx_cut_min = -10.0; // cm
    double zvtx_cut_max = 10.0;  // cm
    // option for simulation to use simulation truth vertex cut
    bool use_truth_vertex_for_zvtx_cut = false;

    // Optional crossing-level cut on the number of silicon seeds associated to the selected vertex.
    // The vertex finder already implies >=2 seeds, so this only acts when set to >=3.
    int min_vertex_associated_silseeds = 2;

    // Truth-vertex-z bin edges used for simulation-only silicon-seed DCA histograms.
    // Open-ended bins are implied below the first edge and above the last edge.
    std::vector<double> truth_vtxz_bin_edges = {-30.0, -20.0, -10.0, -5.0, 0.0, 5.0, 10.0, 20.0, 30.0};
};

// Keep the same formatter behavior as plotutil.h.
template <typename T> static inline std::string to_string_with_precision(const T a_value, const int n = 2)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

static inline std::string ToLowerCopy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

static inline TrackletCombiDRConfig::PairCutMode ParsePairCutMode(const std::string &mode_in)
{
    const std::string mode = ToLowerCopy(mode_in);
    if (mode == "dr" || mode == "deltar")
        return TrackletCombiDRConfig::PairCutMode::kDeltaR;
    if (mode == "dphi" || mode == "deltaphi")
        return TrackletCombiDRConfig::PairCutMode::kDeltaPhi;
    if (mode == "absdphi" || mode == "abs_dphi" || mode == "absolutedeltaphi")
        return TrackletCombiDRConfig::PairCutMode::kAbsDeltaPhi;
    throw std::invalid_argument("Unknown pair cut mode: " + mode_in + ". Use one of: dR, dPhi, absdPhi.");
}

static inline const char *PairCutModeLabel(TrackletCombiDRConfig::PairCutMode mode)
{
    switch (mode)
    {
    case TrackletCombiDRConfig::PairCutMode::kDeltaR:
        return "dR";
    case TrackletCombiDRConfig::PairCutMode::kDeltaPhi:
        return "dPhi";
    case TrackletCombiDRConfig::PairCutMode::kAbsDeltaPhi:
        return "absdPhi";
    }
    return "dR";
}

struct TrackletCombiDROutput
{
    struct CrossingDisplayGraphs
    {
        std::string crossing_key;
        TGraph *grAllINTTClustersXY = nullptr;
        TGraph *grAllINTTClustersZR = nullptr;
        TGraph *grPrimaryINTTClustersXY = nullptr;
        TGraph *grPrimaryINTTClustersZR = nullptr;
        TGraph *grNonPrimaryINTTClustersXY = nullptr;
        TGraph *grNonPrimaryINTTClustersZR = nullptr;
        TGraph *grAssociatedSilSeedPCAXY = nullptr;
        TGraph *grAssociatedSilSeedPCAZR = nullptr;
        TGraph *grRecoVertexXY = nullptr;
        TGraph *grRecoVertexZR = nullptr;
        TGraph *grTruthVertexXY = nullptr;
        TGraph *grTruthVertexZR = nullptr;
        std::vector<TGraph *> grAssociatedSilSeedClustersXY;
        std::vector<TGraph *> grAssociatedSilSeedClustersZR;
        std::vector<TGraph *> grNonAssociatedSilSeedClustersXY;
        std::vector<TGraph *> grNonAssociatedSilSeedClustersZR;
        std::vector<TGraph *> grPrimaryTruthClustersXY;
        std::vector<TGraph *> grPrimaryTruthClustersZR;
        std::vector<TGraph *> grPrimaryRecoClustersXY;
        std::vector<TGraph *> grPrimaryRecoClustersZR;
    };

    TH2D *hSig = nullptr;                // unrotated all-pairs passing deltaR
    TH2D *hBkg = nullptr;                // rotated all-pairs passing deltaR
    TH2D *hSub = nullptr;                // hSig - scale*hBkg
    TH2D *hSub_AbsDphiCount = nullptr;   // directly filled background-subtracted yield with |delta-phi| < abs_dphi_count_cut
    TH1D *hNSelectedCrossings = nullptr; // number of crossings passing (only-1-vertex && crossing != short-int-max)
    long long nSelectedCrossings = 0;    // authoritative counter; written into hNSelectedCrossings at output time
    TH1D *hSeedPhi = nullptr;            // diagnostic phi(seed hit) w.r.t. the silicon vertex
    TH1D *hDCA3D_Sig = nullptr;               // integrated 3D DCA of unrotated signal doublets to selected silicon vertex
    TH1D *hDCA3D_Bkg = nullptr;               // integrated 3D DCA of pi-rotated background doublets to selected silicon vertex
    TH1D *hDCA3DHelix_Sig = nullptr;          // integrated 3D DCA of unrotated signal doublets from lightweight helix fit
    TH1D *hDCA3DHelix_Bkg = nullptr;          // integrated 3D DCA of pi-rotated background doublets from lightweight helix fit
    TH1D *hSilSeedDCA3D = nullptr;            // integrated PCA-to-associated-vertex distance for silicon seeds
    std::vector<TH1D *> hSilSeedDCA3D_truthVtxZ; // simulation-only PCA-to-associated-vertex distance for silicon seeds by truth vertex z range
    TH1D *hSilSeed_Eta = nullptr;        // selected-crossing silicon-seed eta w.r.t. vertex
    TH1D *hSilSeed_Eta_nMVTX3nINTT2 = nullptr;
    TH1D *hSilSeed_Eta_nMVTX3nINTT1 = nullptr;
    TH1D *hSilSeed_Eta_nMVTX2nINTT2 = nullptr;
    TH1D *hSilSeed_Eta_nMVTX2nINTT1 = nullptr;
    TH1D *hSilSeed_Eta_nMVTXnINTTOther = nullptr;
    TH1D *hAllSilSeed_Eta = nullptr;
    TH1D *hAllSilSeed_Eta_nMVTX3nINTT2 = nullptr;
    TH1D *hAllSilSeed_Eta_nMVTX3nINTT1 = nullptr;
    TH1D *hAllSilSeed_Eta_nMVTX2nINTT2 = nullptr;
    TH1D *hAllSilSeed_Eta_nMVTX2nINTT1 = nullptr;
    TH1D *hAllSilSeed_Eta_nMVTXnINTTOther = nullptr;
    TH2D *hSeed_Eta_Phi = nullptr;       // diagnostic 2D distribution of seed hit eta vs phi (unrotated)
    TH2D *hSilSeed_Eta_Phi = nullptr;    // silicon-seed eta/phi from ntuple (matched by trackerVertexTrackIDs -> silseed_id)
    TH2D *hSilSeed_Eta_Phi_nMVTX3nINTT2 = nullptr;
    TH2D *hSilSeed_Eta_Phi_nMVTX3nINTT1 = nullptr;
    TH2D *hSilSeed_Eta_Phi_nMVTX2nINTT2 = nullptr;
    TH2D *hSilSeed_Eta_Phi_nMVTX2nINTT1 = nullptr;
    TH2D *hSilSeed_Eta_Phi_nMVTXnINTTOther = nullptr;
    TH2D *hINTTClusterVsAssociatedSilSeed = nullptr; // selected-crossing diagnostic: N(INTT clusters in crossing) vs N(associated silicon seeds)
    TH2D *hINTTClusterVsAllSilSeed = nullptr;        // selected-crossing diagnostic: N(INTT clusters in crossing) vs N(all silicon seeds in crossing)
    TH2D *hPrimaryCharged_Eta_Phi = nullptr;         // simulation-only: sPHENIX primary charged hadron eta/phi
    // seed crossing histogram
    TH1D *hSeedCrossing_selectedCrossings = nullptr;                        // the crossing number for crossing with only-1-vertex && crossing != short-int-max
    TH1D *hSeedCrossing_selectedCrossings_largeDiffNClusterSeeds = nullptr; // the crossing number for crossing with only-1-vertex && crossing != short-int-max where the INTT clusters are much more than silicon seeds associated to vertex

    // Per-(eta,phi) bin deltaR distributions (lazy-created; nullptr if never filled)
    std::vector<TH1D *> hDR_Sig;
    std::vector<TH1D *> hDR_Bkg;
    std::vector<TH1D *> hDR_Sub;

    int nEta = 0;
    int nPhi = 0;

    // Used to decide how to wrap rotated phi for binning consistency
    bool phiAxisIs_0_2pi = false;

    // Crossing-level display objects for seed/cluster scatter visualization.
    std::vector<CrossingDisplayGraphs> crossingDisplays;
};

struct MultiplicityPercentileBin
{
    int percentile_bin = -1;
    int cluster_low = 0;
    int cluster_high = std::numeric_limits<int>::max();
};

static inline size_t flatIndex(int ieta, int iphi, int nPhi) { return size_t(ieta) * size_t(nPhi) + size_t(iphi); }

static inline std::pair<double, double> GetAxisBinRange(const AxisBinning &axis, int ibin)
{
    if (axis.use_edges)
    {
        return {axis.edges.at(ibin), axis.edges.at(ibin + 1)};
    }
    const double width = (axis.max - axis.min) / axis.nbins;
    const double low = axis.min + ibin * width;
    return {low, low + width};
}

static inline std::string EtaPhiRangeNameSuffix(const TrackletCombiDRConfig &cfg, int ieta, int iphi)
{
    const auto etaRange = GetAxisBinRange(cfg.etaAxis, ieta);
    const auto phiRange = GetAxisBinRange(cfg.phiAxis, iphi);
    std::ostringstream nm;
    nm << "eta" << to_string_with_precision(etaRange.first) << "-" << to_string_with_precision(etaRange.second) << "_phi" << to_string_with_precision(phiRange.first) << "-" << to_string_with_precision(phiRange.second);
    return nm.str();
}

static inline int FindTruthVtxZBin(const std::vector<double> &edges, double z)
{
    return static_cast<int>(std::upper_bound(edges.begin(), edges.end(), z) - edges.begin());
}

static inline std::pair<double, double> GetTruthVtxZBinRange(const std::vector<double> &edges, int ibin)
{
    const double low = (ibin == 0) ? -std::numeric_limits<double>::infinity() : edges.at(ibin - 1);
    const double high = (ibin == static_cast<int>(edges.size())) ? std::numeric_limits<double>::infinity() : edges.at(ibin);
    return {low, high};
}

static inline std::string TruthVtxZRangeNameSuffix(const std::vector<double> &edges, int ibin)
{
    const auto range = GetTruthVtxZBinRange(edges, ibin);

    auto format_bound = [](double value) -> std::string
    {
        if (!std::isfinite(value))
            return (value < 0.0) ? "mInf" : "pInf";

        std::string s = to_string_with_precision(value);
        std::replace(s.begin(), s.end(), '-', 'm');
        std::replace(s.begin(), s.end(), '.', 'p');
        return s;
    };

    std::ostringstream nm;
    nm << "truthVtxZ_bin" << std::setw(2) << std::setfill('0') << ibin << "_" << format_bound(range.first) << "_" << format_bound(range.second);
    return nm.str();
}

static inline std::string TruthVtxZRangeTitle(const std::vector<double> &edges, int ibin)
{
    const auto range = GetTruthVtxZBinRange(edges, ibin);

    auto format_bound = [](double value) -> std::string
    {
        if (!std::isfinite(value))
            return (value < 0.0) ? "-inf" : "inf";
        return to_string_with_precision(value, 1);
    };

    const bool is_first_bin = (ibin == 0);
    const bool is_last_bin = (ibin == static_cast<int>(edges.size()));

    std::ostringstream title;
    title << (is_first_bin ? "(" : "[")
          << format_bound(range.first) << ", " << format_bound(range.second)
          << (is_last_bin ? "]" : ")");
    return title.str();
}

static inline TH1D *getOrMakeDRHist(std::vector<TH1D *> &vec, int ieta, int iphi, int nPhi, const TrackletCombiDRConfig &cfg, const std::string &name_prefix, const std::string &tag)
{
    const size_t idx = flatIndex(ieta, iphi, nPhi);
    if (vec[idx])
        return vec[idx];

    std::ostringstream nm;
    nm << name_prefix << "_" << PairCutModeLabel(cfg.pair_cut_mode) << "_" << tag << "_" << EtaPhiRangeNameSuffix(cfg, ieta, iphi);

    vec[idx] = new TH1D(nm.str().c_str(), nm.str().c_str(), cfg.dR_nbins, cfg.dR_min, cfg.dR_max);
    vec[idx]->SetDirectory(nullptr);
    vec[idx]->Sumw2();
    return vec[idx];
}

static inline TH1D *MakeTH1D(const std::string &name, int nbins, double xmin, double xmax)
{
    TH1D *hist = new TH1D(name.c_str(), name.c_str(), nbins, xmin, xmax);
    hist->SetDirectory(nullptr);
    hist->Sumw2();
    return hist;
}

static inline TH1D *MakeTH1D(const std::string &name, const AxisBinning &axis)
{
    TH1D *hist = axis.use_edges ? new TH1D(name.c_str(), name.c_str(), axis.nbins, axis.edges.data()) : new TH1D(name.c_str(), name.c_str(), axis.nbins, axis.min, axis.max);
    hist->SetDirectory(nullptr);
    hist->Sumw2();
    return hist;
}

static inline TH2D *MakeEtaPhiTH2D(const std::string &name, const AxisBinning &eta_axis, const AxisBinning &phi_axis)
{
    TH2D *hist = eta_axis.use_edges ? new TH2D(name.c_str(), name.c_str(), eta_axis.nbins, eta_axis.edges.data(), phi_axis.nbins, phi_axis.edges.data()) : new TH2D(name.c_str(), name.c_str(), eta_axis.nbins, eta_axis.min, eta_axis.max, phi_axis.nbins, phi_axis.min, phi_axis.max);
    hist->SetDirectory(nullptr);
    hist->Sumw2();
    return hist;
}

static inline void SetAxisTitles(TH1 *hist, const char *xtitle, const char *ytitle)
{
    hist->GetXaxis()->SetTitle(xtitle);
    hist->GetYaxis()->SetTitle(ytitle);
}

static inline void SetAxisTitles(TH2 *hist, const char *xtitle, const char *ytitle, const char *ztitle)
{
    hist->GetXaxis()->SetTitle(xtitle);
    hist->GetYaxis()->SetTitle(ytitle);
    if (ztitle)
        hist->GetZaxis()->SetTitle(ztitle);
}

template <class TObjectType> static inline void WriteIfValid(TObjectType *obj)
{
    if (!obj)
        return;
    obj->Write(obj->GetName(), TObject::kOverwrite);
}

static inline std::string TrackletCombiDROutputPrefix(const TrackletCombiDROutput &out)
{
    if (!out.hSig)
        return "tklCombiDR";

    const std::string sig_name = out.hSig->GetName();
    const std::string suffix = "_Sig";
    if (sig_name.size() > suffix.size() && sig_name.compare(sig_name.size() - suffix.size(), suffix.size(), suffix) == 0)
        return sig_name.substr(0, sig_name.size() - suffix.size());

    return sig_name;
}

static inline std::string MultiplicityPercentileBinPrefix(const std::string &base_prefix, int percentile_bin) { return base_prefix + "_multPercentileBin" + std::to_string(percentile_bin); }

static inline bool IsInMultiplicityPercentileBin(const MultiplicityPercentileBin &bin, int n_intt_clusters)
{
    if (n_intt_clusters < bin.cluster_low)
        return false;
    if (bin.cluster_high == std::numeric_limits<int>::max())
        return true;
    return n_intt_clusters < bin.cluster_high;
}

static inline int FindMultiplicityPercentileBinIndex(const std::vector<MultiplicityPercentileBin> &bins, int n_intt_clusters)
{
    for (size_t i = 0; i < bins.size(); ++i)
    {
        if (IsInMultiplicityPercentileBin(bins[i], n_intt_clusters))
            return static_cast<int>(i);
    }
    return -1;
}

static inline std::vector<MultiplicityPercentileBin> LoadMultiplicityPercentileBins(const std::string &input_path, const std::string &tree_name = "trigger_percentile_intervals")
{
    if (input_path.empty())
        throw std::runtime_error("Multiplicity percentile boundary file path is empty.");

    TFile fin(input_path.c_str(), "READ");
    if (fin.IsZombie())
        throw std::runtime_error("Unable to open NIntt cluster percentile boundary ROOT file: " + input_path);

    TTree *tree = dynamic_cast<TTree *>(fin.Get(tree_name.c_str()));
    if (!tree)
        throw std::runtime_error("Unable to find percentile interval tree '" + tree_name + "' in " + input_path);

    int percentile_bin = -1;
    int cluster_low = 0;
    int cluster_high = 0;
    tree->SetBranchAddress("percentile_bin", &percentile_bin);
    tree->SetBranchAddress("cluster_low", &cluster_low);
    tree->SetBranchAddress("cluster_high", &cluster_high);

    std::vector<MultiplicityPercentileBin> bins;
    const Long64_t nentries = tree->GetEntries();
    bins.reserve(static_cast<size_t>(nentries));
    for (Long64_t i = 0; i < nentries; ++i)
    {
        tree->GetEntry(i);
        bins.push_back({percentile_bin, cluster_low, cluster_high});
    }

    std::sort(bins.begin(), bins.end(), [](const MultiplicityPercentileBin &lhs, const MultiplicityPercentileBin &rhs) { return lhs.percentile_bin < rhs.percentile_bin; });

    if (bins.empty())
        throw std::runtime_error("Percentile interval tree is empty in " + input_path);

    return bins;
}

static inline void InitTrackletCombiDROutput(TrackletCombiDROutput &out,                    //
                                             const TrackletCombiDRConfig &cfg,              //
                                             const std::string &name_prefix = "tklCombiDR", //
                                             bool make_primary_charged_hist = false,        //
                                             bool make_seed_eta_1d_hists = true             //
)
{
    const double twopi = 2.0 * TMath::Pi();
    out.phiAxisIs_0_2pi = (cfg.phiAxis.min >= 0.0 && cfg.phiAxis.max <= twopi + 1e-6);

    out.nEta = cfg.etaAxis.nbins;
    out.nPhi = cfg.phiAxis.nbins;

    const std::string nSig = name_prefix + "_Sig";
    const std::string nBkg = name_prefix + "_Bkg";
    const std::string nSub = name_prefix + "_Sub";
    const std::string nSubAbsDphiCount = name_prefix + "_Sub_AbsDphiCount";
    const std::string nNSelectedCrossings = name_prefix + "_NSelectedCrossings";
    const std::string nSeedPhi = name_prefix + "_SeedPhi";
    const std::string nDCA3DSig = name_prefix + "_DCA3D_Sig";
    const std::string nDCA3DBkg = name_prefix + "_DCA3D_Bkg";
    const std::string nDCA3DHelixSig = name_prefix + "_DCA3DHelix_Sig";
    const std::string nDCA3DHelixBkg = name_prefix + "_DCA3DHelix_Bkg";
    const std::string nSilSeedDCA3D = name_prefix + "_SilSeedDCA3D";
    const std::string nSilSeedEta = name_prefix + "_SilSeedEta";
    const std::string nSilSeedEta_nMVTX3nINTT2 = name_prefix + "_SilSeedEta_nMVTX3nINTT2";
    const std::string nSilSeedEta_nMVTX3nINTT1 = name_prefix + "_SilSeedEta_nMVTX3nINTT1";
    const std::string nSilSeedEta_nMVTX2nINTT2 = name_prefix + "_SilSeedEta_nMVTX2nINTT2";
    const std::string nSilSeedEta_nMVTX2nINTT1 = name_prefix + "_SilSeedEta_nMVTX2nINTT1";
    const std::string nSilSeedEta_nMVTXnINTTOther = name_prefix + "_SilSeedEta_nMVTXnINTTOther";
    const std::string nAllSilSeedEta = name_prefix + "_AllSilSeedEta";
    const std::string nAllSilSeedEta_nMVTX3nINTT2 = name_prefix + "_AllSilSeedEta_nMVTX3nINTT2";
    const std::string nAllSilSeedEta_nMVTX3nINTT1 = name_prefix + "_AllSilSeedEta_nMVTX3nINTT1";
    const std::string nAllSilSeedEta_nMVTX2nINTT2 = name_prefix + "_AllSilSeedEta_nMVTX2nINTT2";
    const std::string nAllSilSeedEta_nMVTX2nINTT1 = name_prefix + "_AllSilSeedEta_nMVTX2nINTT1";
    const std::string nAllSilSeedEta_nMVTXnINTTOther = name_prefix + "_AllSilSeedEta_nMVTXnINTTOther";
    const std::string nSilSeedEtaPhi = name_prefix + "_SilSeedEtaPhi";
    const std::string nSilSeedEtaPhi_nMVTX3nINTT2 = name_prefix + "_SilSeedEtaPhi_nMVTX3nINTT2";
    const std::string nSilSeedEtaPhi_nMVTX3nINTT1 = name_prefix + "_SilSeedEtaPhi_nMVTX3nINTT1";
    const std::string nSilSeedEtaPhi_nMVTX2nINTT2 = name_prefix + "_SilSeedEtaPhi_nMVTX2nINTT2";
    const std::string nSilSeedEtaPhi_nMVTX2nINTT1 = name_prefix + "_SilSeedEtaPhi_nMVTX2nINTT1";
    const std::string nSilSeedEtaPhi_nMVTXnINTTOther = name_prefix + "_SilSeedEtaPhi_nMVTXnINTTOther";
    const std::string nINTTClusterVsAssociatedSilSeed = name_prefix + "_INTTClusterVsAssociatedSilSeed";
    const std::string nINTTClusterVsAllSilSeed = name_prefix + "_INTTClusterVsAllSilSeed";
    const std::string nPrimaryChargedEtaPhi = name_prefix + "_PrimaryChargedEtaPhi";
    const std::string nSeedCrossingSelectedCrossings = name_prefix + "_SeedCrossing_selectedCrossings";
    const std::string nSeedCrossingSelectedCrossingsLargeDiffNClusterSeeds = name_prefix + "_SeedCrossing_selectedCrossings_largeDiffNClusterSeeds";

    if (cfg.etaAxis.use_edges != cfg.phiAxis.use_edges)
    {
        // Keeping this strict avoids subtle axis-mode bugs.
        throw std::runtime_error("InitTrackletCombiDROutput: eta/phi axis must both be uniform or both be edges.");
    }

    out.hSig = MakeEtaPhiTH2D(nSig, cfg.etaAxis, cfg.phiAxis);
    out.hBkg = MakeEtaPhiTH2D(nBkg, cfg.etaAxis, cfg.phiAxis);
    out.hSub = MakeEtaPhiTH2D(nSub, cfg.etaAxis, cfg.phiAxis);
    out.hSub_AbsDphiCount = MakeEtaPhiTH2D(nSubAbsDphiCount, cfg.etaAxis, cfg.phiAxis);

    out.hNSelectedCrossings = MakeTH1D(nNSelectedCrossings, 1, 0.5, 1.5);
    out.hNSelectedCrossings->GetXaxis()->SetBinLabel(1, "selected_crossings");
    out.hNSelectedCrossings->GetYaxis()->SetTitle("Counts");
    out.nSelectedCrossings = 0;

    out.hSeedPhi = MakeTH1D(nSeedPhi, cfg.phiAxis);
    out.hDCA3D_Sig = MakeTH1D(nDCA3DSig, 200, 0.0, 10.0);
    out.hDCA3D_Bkg = MakeTH1D(nDCA3DBkg, 200, 0.0, 10.0);
    out.hDCA3DHelix_Sig = MakeTH1D(nDCA3DHelixSig, 200, 0.0, 10.0);
    out.hDCA3DHelix_Bkg = MakeTH1D(nDCA3DHelixBkg, 200, 0.0, 10.0);
    out.hSilSeedDCA3D = MakeTH1D(nSilSeedDCA3D, 200, 0.0, 10.0);
    out.hSilSeedDCA3D_truthVtxZ.clear();
    if (make_primary_charged_hist)
    {
        out.hSilSeedDCA3D_truthVtxZ.reserve(cfg.truth_vtxz_bin_edges.size() + 1);
        for (int ibin = 0; ibin <= static_cast<int>(cfg.truth_vtxz_bin_edges.size()); ++ibin)
        {
            const std::string hist_name = nSilSeedDCA3D + "_" + TruthVtxZRangeNameSuffix(cfg.truth_vtxz_bin_edges, ibin);
            TH1D *hist = MakeTH1D(hist_name, 200, 0.0, 10.0);
            hist->SetTitle((hist_name + " " + TruthVtxZRangeTitle(cfg.truth_vtxz_bin_edges, ibin)).c_str());
            SetAxisTitles(hist, "|PCA - associated silicon vertex| [cm]", "Counts");
            out.hSilSeedDCA3D_truthVtxZ.push_back(hist);
        }
    }
    if (make_seed_eta_1d_hists)
    {
        out.hSilSeed_Eta = MakeTH1D(nSilSeedEta, 70, -3.5, 3.5);
        out.hSilSeed_Eta_nMVTX3nINTT2 = MakeTH1D(nSilSeedEta_nMVTX3nINTT2, 70, -3.5, 3.5);
        out.hSilSeed_Eta_nMVTX3nINTT1 = MakeTH1D(nSilSeedEta_nMVTX3nINTT1, 70, -3.5, 3.5);
        out.hSilSeed_Eta_nMVTX2nINTT2 = MakeTH1D(nSilSeedEta_nMVTX2nINTT2, 70, -3.5, 3.5);
        out.hSilSeed_Eta_nMVTX2nINTT1 = MakeTH1D(nSilSeedEta_nMVTX2nINTT1, 70, -3.5, 3.5);
        out.hSilSeed_Eta_nMVTXnINTTOther = MakeTH1D(nSilSeedEta_nMVTXnINTTOther, 70, -3.5, 3.5);
        out.hAllSilSeed_Eta = MakeTH1D(nAllSilSeedEta, 70, -3.5, 3.5);
        out.hAllSilSeed_Eta_nMVTX3nINTT2 = MakeTH1D(nAllSilSeedEta_nMVTX3nINTT2, 70, -3.5, 3.5);
        out.hAllSilSeed_Eta_nMVTX3nINTT1 = MakeTH1D(nAllSilSeedEta_nMVTX3nINTT1, 70, -3.5, 3.5);
        out.hAllSilSeed_Eta_nMVTX2nINTT2 = MakeTH1D(nAllSilSeedEta_nMVTX2nINTT2, 70, -3.5, 3.5);
        out.hAllSilSeed_Eta_nMVTX2nINTT1 = MakeTH1D(nAllSilSeedEta_nMVTX2nINTT1, 70, -3.5, 3.5);
        out.hAllSilSeed_Eta_nMVTXnINTTOther = MakeTH1D(nAllSilSeedEta_nMVTXnINTTOther, 70, -3.5, 3.5);
    }
    out.hSeed_Eta_Phi = MakeEtaPhiTH2D(name_prefix + "_SeedEtaPhi", cfg.etaAxis, cfg.phiAxis);
    out.hSilSeed_Eta_Phi = MakeEtaPhiTH2D(nSilSeedEtaPhi, cfg.etaAxis, cfg.phiAxis);
    out.hSilSeed_Eta_Phi_nMVTX3nINTT2 = MakeEtaPhiTH2D(nSilSeedEtaPhi_nMVTX3nINTT2, cfg.etaAxis, cfg.phiAxis);
    out.hSilSeed_Eta_Phi_nMVTX3nINTT1 = MakeEtaPhiTH2D(nSilSeedEtaPhi_nMVTX3nINTT1, cfg.etaAxis, cfg.phiAxis);
    out.hSilSeed_Eta_Phi_nMVTX2nINTT2 = MakeEtaPhiTH2D(nSilSeedEtaPhi_nMVTX2nINTT2, cfg.etaAxis, cfg.phiAxis);
    out.hSilSeed_Eta_Phi_nMVTX2nINTT1 = MakeEtaPhiTH2D(nSilSeedEtaPhi_nMVTX2nINTT1, cfg.etaAxis, cfg.phiAxis);
    out.hSilSeed_Eta_Phi_nMVTXnINTTOther = MakeEtaPhiTH2D(nSilSeedEtaPhi_nMVTXnINTTOther, cfg.etaAxis, cfg.phiAxis);
    out.hINTTClusterVsAssociatedSilSeed = new TH2D(nINTTClusterVsAssociatedSilSeed.c_str(), nINTTClusterVsAssociatedSilSeed.c_str(), 150, -0.5, 299.5, 150, -0.5, 299.5);
    out.hINTTClusterVsAssociatedSilSeed->SetDirectory(nullptr);
    out.hINTTClusterVsAssociatedSilSeed->Sumw2();
    out.hINTTClusterVsAllSilSeed = new TH2D(nINTTClusterVsAllSilSeed.c_str(), nINTTClusterVsAllSilSeed.c_str(), 150, -0.5, 299.5, 150, -0.5, 299.5);
    out.hINTTClusterVsAllSilSeed->SetDirectory(nullptr);
    out.hINTTClusterVsAllSilSeed->Sumw2();
    out.hSeedCrossing_selectedCrossings = MakeTH1D(nSeedCrossingSelectedCrossings, 680, -130.5, 549.5);
    out.hSeedCrossing_selectedCrossings_largeDiffNClusterSeeds = MakeTH1D(nSeedCrossingSelectedCrossingsLargeDiffNClusterSeeds, 680, -130.5, 549.5);

    SetAxisTitles(out.hSeedPhi, "#phi_{seed} (rad)", "Counts");
    SetAxisTitles(out.hDCA3D_Sig, "3D DCA to selected silicon vertex [cm]", "Counts");
    SetAxisTitles(out.hDCA3D_Bkg, "3D DCA to selected silicon vertex [cm]", "Counts");
    SetAxisTitles(out.hDCA3DHelix_Sig, "3D DCA from lightweight helix fit to selected silicon vertex [cm]", "Counts");
    SetAxisTitles(out.hDCA3DHelix_Bkg, "3D DCA from lightweight helix fit to selected silicon vertex [cm]", "Counts");
    SetAxisTitles(out.hSilSeedDCA3D, "|PCA - associated silicon vertex| [cm]", "Counts");
    if (make_seed_eta_1d_hists)
    {
        SetAxisTitles(out.hSilSeed_Eta, "#eta_{silseed,vtx}", "Counts");
        SetAxisTitles(out.hSilSeed_Eta_nMVTX3nINTT2, "#eta_{silseed,vtx}", "Counts");
        SetAxisTitles(out.hSilSeed_Eta_nMVTX3nINTT1, "#eta_{silseed,vtx}", "Counts");
        SetAxisTitles(out.hSilSeed_Eta_nMVTX2nINTT2, "#eta_{silseed,vtx}", "Counts");
        SetAxisTitles(out.hSilSeed_Eta_nMVTX2nINTT1, "#eta_{silseed,vtx}", "Counts");
        SetAxisTitles(out.hSilSeed_Eta_nMVTXnINTTOther, "#eta_{silseed,vtx}", "Counts");
        SetAxisTitles(out.hAllSilSeed_Eta, "#eta_{silseed,vtx}", "Counts");
        SetAxisTitles(out.hAllSilSeed_Eta_nMVTX3nINTT2, "#eta_{silseed,vtx}", "Counts");
        SetAxisTitles(out.hAllSilSeed_Eta_nMVTX3nINTT1, "#eta_{silseed,vtx}", "Counts");
        SetAxisTitles(out.hAllSilSeed_Eta_nMVTX2nINTT2, "#eta_{silseed,vtx}", "Counts");
        SetAxisTitles(out.hAllSilSeed_Eta_nMVTX2nINTT1, "#eta_{silseed,vtx}", "Counts");
        SetAxisTitles(out.hAllSilSeed_Eta_nMVTXnINTTOther, "#eta_{silseed,vtx}", "Counts");
    }
    SetAxisTitles(out.hSeed_Eta_Phi, "#eta_{seed}", "#phi_{seed} (rad)", "Counts");
    SetAxisTitles(out.hSilSeed_Eta_Phi, "#eta_{silseed,vtx}", "#phi_{silseed,vtx} (rad)", "Counts");
    SetAxisTitles(out.hSilSeed_Eta_Phi_nMVTX3nINTT2, "#eta_{silseed,vtx}", "#phi_{silseed,vtx} (rad)", "Counts");
    SetAxisTitles(out.hSilSeed_Eta_Phi_nMVTX3nINTT1, "#eta_{silseed,vtx}", "#phi_{silseed,vtx} (rad)", "Counts");
    SetAxisTitles(out.hSilSeed_Eta_Phi_nMVTX2nINTT2, "#eta_{silseed,vtx}", "#phi_{silseed,vtx} (rad)", "Counts");
    SetAxisTitles(out.hSilSeed_Eta_Phi_nMVTX2nINTT1, "#eta_{silseed,vtx}", "#phi_{silseed,vtx} (rad)", "Counts");
    SetAxisTitles(out.hSilSeed_Eta_Phi_nMVTXnINTTOther, "#eta_{silseed,vtx}", "#phi_{silseed,vtx} (rad)", "Counts");
    SetAxisTitles(out.hINTTClusterVsAssociatedSilSeed, "N_{INTT clusters}^{crossing}", "N_{silicon seeds}^{associated to vertex}", "Selected crossings");
    SetAxisTitles(out.hINTTClusterVsAllSilSeed, "N_{INTT clusters}^{crossing}", "N_{silicon seeds}^{all in crossing}", "Selected crossings");
    SetAxisTitles(out.hSeedCrossing_selectedCrossings, "Crossing", "Selected crossings");
    SetAxisTitles(out.hSeedCrossing_selectedCrossings_largeDiffNClusterSeeds, "Crossing", "Selected crossings");

    if (make_primary_charged_hist)
    {
        out.hPrimaryCharged_Eta_Phi = MakeEtaPhiTH2D(nPrimaryChargedEtaPhi, cfg.etaAxis, cfg.phiAxis);
        SetAxisTitles(out.hPrimaryCharged_Eta_Phi, "#eta_{primary,chg}", "#phi_{primary,chg} (rad)", "Counts");
    }

    SetAxisTitles(out.hSig, "#eta_{doublet}", "#phi_{doublet} (rad)", nullptr);
    SetAxisTitles(out.hBkg, "#eta_{doublet}", "#phi_{doublet} (rad)", nullptr);
    SetAxisTitles(out.hSub, "#eta_{doublet}", "#phi_{doublet} (rad)", nullptr);
    SetAxisTitles(out.hSub_AbsDphiCount, "#eta_{doublet}", "#phi_{doublet} (rad)", nullptr);

    const size_t n = size_t(out.nEta) * size_t(out.nPhi);
    out.hDR_Sig.assign(n, nullptr);
    out.hDR_Bkg.assign(n, nullptr);
    out.hDR_Sub.assign(n, nullptr);
}

static inline bool PairZRangeCoversSiliconVertex(Hit *seedHit,                            //
                                                 Hit *searchHit,                          //
                                                 const std::array<double, 3> &silicon_vtx //
)
{
    if (!seedHit || !searchHit)
        return false;

    const double vtx_x = silicon_vtx[0];
    const double vtx_y = silicon_vtx[1];
    const double vtx_z = silicon_vtx[2];

    // Project the seed-search segment to rho=0 around the silicon vertex, using cluster z-edges to build
    // an allowed zvtx segment (same construction as in INTTVtxZ.cxx).
    const double rho_seed = std::hypot(seedHit->posX() - vtx_x, seedHit->posY() - vtx_y);
    const double rho_search = std::hypot(searchHit->posX() - vtx_x, searchHit->posY() - vtx_y);
    const double drho = rho_search - rho_seed;
    if (std::abs(drho) < 1e-9)
        return false;

    const double edge1 = seedHit->Edge().first - (searchHit->Edge().second - seedHit->Edge().first) / drho * rho_seed;
    const double edge2 = seedHit->Edge().second - (searchHit->Edge().first - seedHit->Edge().second) / drho * rho_seed;

    const double zmin = std::min(edge1, edge2);
    const double zmax = std::max(edge1, edge2);
    return (vtx_z >= zmin && vtx_z <= zmax);
}

// Core pair-filling: all pairs between seed and search are considered; deltaR cut applied.
// INTT doublet coordinate for output grid uses pair-average (eta, phi).
static inline void FillCombinatoricPairs_DR(const std::vector<Hit *> &seedHits,       //
                                            const std::vector<Hit *> &searchHits,     //
                                            const TrackletCombiDRConfig &cfg,         //
                                            const std::array<double, 3> &silicon_vtx, //
                                            TrackletCombiDROutput &out,               //
                                            TH2D *hGrid,                              //
                                            TH2D *hAbsDphiCountGrid,                  //
                                            double abs_dphi_count_weight,             //
                                            TH1D *hDCA3D,                             //
                                            TH1D *hDCA3DHelix,                        //
                                            std::vector<TH1D *> &hDR,                 //
                                            const std::string &name_prefix,           //
                                            const std::string &tag                    //
)
{
    if (!hGrid)
        return;
    if (seedHits.empty() || searchHits.empty())
        return;

    const bool do_eta_preselect = (cfg.enable_eta_preselect && cfg.pair_cut_mode == TrackletCombiDRConfig::PairCutMode::kDeltaR);

    // For eta preselection we sort search hits by η once per call.
    std::vector<Hit *> search_sorted = searchHits;
    if (do_eta_preselect)
    {
        std::sort(search_sorted.begin(), search_sorted.end(), [](Hit *a, Hit *b) { return a->Eta() < b->Eta(); });
    }

    const double etaWin = cfg.pair_cut; // because |delta\eta| <= deltaR for pairs passing the cut

    auto *axEta = hGrid->GetXaxis();
    auto *axPhi = hGrid->GetYaxis();

    // a counter for debuging purpose
    int pairCount = 0;

    for (auto *h1 : seedHits)
    {
        const double eta_seed = h1->Eta();

        double phi_seed = h1->Phi();
        phi_seed = out.phiAxisIs_0_2pi ? wrapPhi_0_2pi(phi_seed) : wrapPhi_mpi_pi(phi_seed);

        auto lo = search_sorted.begin();
        auto hi = search_sorted.end();

        if (do_eta_preselect)
        {
            lo = std::lower_bound(search_sorted.begin(), search_sorted.end(), eta_seed - etaWin, [](Hit *h, double val) { return h->Eta() < val; });
            hi = std::upper_bound(search_sorted.begin(), search_sorted.end(), eta_seed + etaWin, [](double val, Hit *h) { return val < h->Eta(); });
        }

        for (auto it = lo; it != hi; ++it)
        {
            Hit *h2 = *it;

            const double eta_search = h2->Eta();
            double phi_search = h2->Phi();

            const double deta = eta_search - eta_seed;
            const double dphi = deltaPhi(phi_search, phi_seed); // must wrap; provided by your Utilities.h

            const double dR = std::sqrt(dphi * dphi + deta * deta);

            const double eta_doublet = 0.5 * (eta_seed + eta_search);
            // Use wrapped delta-phi to compute an angular mean robustly across branch boundaries.
            const double phi_doublet_raw = phi_seed + 0.5 * dphi;
            const double phi_doublet = out.phiAxisIs_0_2pi ? wrapPhi_0_2pi(phi_doublet_raw) : wrapPhi_mpi_pi(phi_doublet_raw);

            const int bEta = axEta->FindBin(eta_doublet);
            const int bPhi = axPhi->FindBin(phi_doublet);
            if (bEta < 1 || bEta > out.nEta)
                continue;
            if (bPhi < 1 || bPhi > out.nPhi)
                continue;

            const int ieta = bEta - 1;
            const int iphi = bPhi - 1;

            double pair_value = dR;
            bool pass_pair_cut = (dR < cfg.pair_cut);
            if (cfg.pair_cut_mode == TrackletCombiDRConfig::PairCutMode::kDeltaPhi)
            {
                pair_value = dphi;
                pass_pair_cut = (std::abs(dphi) < cfg.pair_cut);
            }
            else if (cfg.pair_cut_mode == TrackletCombiDRConfig::PairCutMode::kAbsDeltaPhi)
            {
                pair_value = std::abs(dphi);
                pass_pair_cut = (pair_value < cfg.pair_cut);
            }

            // Keep only pairs whose projected zvtx segment covers this crossing's silicon vertex z.
            if (!PairZRangeCoversSiliconVertex(h1, h2, silicon_vtx))
                continue;

            if (hAbsDphiCountGrid && std::abs(dphi) < cfg.abs_dphi_count_cut)
                hAbsDphiCountGrid->Fill(eta_doublet, phi_doublet, abs_dphi_count_weight);

            if (!pass_pair_cut)
                continue;

            hGrid->Fill(eta_doublet, phi_doublet, 1.0);

            if (hDCA3D)
            {
                Point3D p1 = {h1->posX(), h1->posY(), h1->posZ()};
                Point3D p2 = {h2->posX(), h2->posY(), h2->posZ()};
                const Point3D p3 = {silicon_vtx[0], silicon_vtx[1], silicon_vtx[2]};
                const float dca3d = computeDCA(p1, p2, p3);
                if (std::isfinite(dca3d))
                    hDCA3D->Fill(dca3d, 1.0);
            }

            if (hDCA3DHelix)
            {
                const LightHelixFit::Point3D helix_vertex = {silicon_vtx[0], silicon_vtx[1], silicon_vtx[2]};
                const auto helix_fit = LightHelixFit::FitFromVertexAndDoublet(helix_vertex, h1, h2);
                if (helix_fit.success && std::isfinite(helix_fit.dca3d))
                    hDCA3DHelix->Fill(helix_fit.dca3d, 1.0);
            }

            TH1D *hdr = getOrMakeDRHist(hDR, ieta, iphi, out.nPhi, cfg, name_prefix, tag);
            hdr->Fill(pair_value, 1.0);
            pairCount++;
        }
    }

    {
        std::cout << "Number of pairs passing " << PairCutModeLabel(cfg.pair_cut_mode) << " cut in this event: " << pairCount << std::endl;
    }
}

// Fill signal + rotated background for one event.
// Expects tkldata.layers[0] and [1] to already be grouped as:
//   layers[0] = logical seed = {3,4}
//   layers[1] = logical search = {5,6}
static inline void Tracklets_CombinatoricDR_FillEvent(TrackletData &tkldata,                        //
                                                      const TrackletCombiDRConfig &cfg,             //
                                                      const std::array<double, 3> &silicon_vtx,     //
                                                      TrackletCombiDROutput &out,                   //
                                                      const std::string &name_prefix = "tklCombiDR" //
)
{
    const std::vector<Hit *> &L0 = tkldata.layers[0]; // seed ({3,4})
    const std::vector<Hit *> &L1 = tkldata.layers[1]; // search ({5,6})

    // Seed-phi diagnostic is a pure seed-hit observable (no signal/background split).
    // if (out.hSeedPhi)
    {
        for (auto *hseed : L0)
        {
            double phi_seed = hseed->Phi();
            phi_seed = out.phiAxisIs_0_2pi ? wrapPhi_0_2pi(phi_seed) : wrapPhi_mpi_pi(phi_seed);
            out.hSeedPhi->Fill(phi_seed, 1.0);
            out.hSeed_Eta_Phi->Fill(hseed->Eta(), phi_seed, 1.0);
        }
    }

    // Signal (unrotated)
    FillCombinatoricPairs_DR(L0, L1, cfg, silicon_vtx, out, out.hSig, out.hSub_AbsDphiCount, 1.0, out.hDCA3D_Sig, out.hDCA3DHelix_Sig, out.hDR_Sig, name_prefix, "sig");

    // Background (rotate by pi)
    if (cfg.rotate_search_layer)
    {
        // Rotate the search hits in place for the background estimate, then restore them so any
        // subsequent fill on the same tkldata sees the original unrotated geometry.
        for (auto *h : L1)
        {
            h->RotatePhiByPi();
        }
        FillCombinatoricPairs_DR(L0, L1, cfg, silicon_vtx, out, out.hBkg, out.hSub_AbsDphiCount, -cfg.bkg_scale, out.hDCA3D_Bkg, out.hDCA3DHelix_Bkg, out.hDR_Bkg, name_prefix, "bkg");
        for (auto *h : L1)
        {
            h->RotatePhiByPi();
        }
    }
    else
    {
        // Rotate the seed hits in place for the background estimate, then restore them so any
        // subsequent fill on the same tkldata sees the original unrotated geometry.
        for (auto *h : L0)
        {
            h->RotatePhiByPi();
        }
        FillCombinatoricPairs_DR(L0, L1, cfg, silicon_vtx, out, out.hBkg, out.hSub_AbsDphiCount, -cfg.bkg_scale, out.hDCA3D_Bkg, out.hDCA3DHelix_Bkg, out.hDR_Bkg, name_prefix, "bkg");
        for (auto *h : L0)
        {
            h->RotatePhiByPi();
        }
    }

    // for debuging purpose, print out the entries of histograms after each event
    std::cout << "Event combinatoric pairs passing " << PairCutModeLabel(cfg.pair_cut_mode) << " cut: Sig = " << out.hSig->Integral() << ", Bkg = " << out.hBkg->Integral() << std::endl;
}

// After event loop: build subtracted grids & per-bin deltaR hists.
// Uncorrected = Sig - scale * Bkg
static inline void Tracklets_CombinatoricDR_Finalize(const TrackletCombiDRConfig &cfg, TrackletCombiDROutput &out, const std::string &name_prefix = "tklCombiDR")
{
    if (!out.hSig || !out.hBkg || !out.hSub)
        return;

    out.hSub->Reset("ICES");
    out.hSub->Add(out.hSig, 1.0);
    out.hSub->Add(out.hBkg, -cfg.bkg_scale);

    for (int ieta = 0; ieta < out.nEta; ++ieta)
    {
        for (int iphi = 0; iphi < out.nPhi; ++iphi)
        {
            const size_t idx = flatIndex(ieta, iphi, out.nPhi);
            TH1D *hs = out.hDR_Sig[idx];
            TH1D *hb = out.hDR_Bkg[idx];
            if (!hs && !hb)
                continue;

            if (!out.hDR_Sub[idx])
            {
                std::ostringstream nm;
                nm << name_prefix << "_" << PairCutModeLabel(cfg.pair_cut_mode) << "_sub_" << EtaPhiRangeNameSuffix(cfg, ieta, iphi);
                out.hDR_Sub[idx] = new TH1D(nm.str().c_str(), nm.str().c_str(), cfg.dR_nbins, cfg.dR_min, cfg.dR_max);
                out.hDR_Sub[idx]->SetDirectory(nullptr);
                out.hDR_Sub[idx]->Sumw2();
            }

            out.hDR_Sub[idx]->Reset("ICES");
            if (hs)
                out.hDR_Sub[idx]->Add(hs, 1.0);
            if (hb)
                out.hDR_Sub[idx]->Add(hb, -cfg.bkg_scale);
        }
    }
}

// Write all combinatoric DR outputs to a ROOT file.
static inline void WriteTracklets_CombinatoricDR_OutputToDirectory(const TrackletCombiDROutput &out, TDirectory *base_dir, bool preserve_legacy_dirnames = false)
{
    if (!base_dir)
        return;

    base_dir->cd();
    WriteIfValid(out.hSig);
    WriteIfValid(out.hBkg);
    WriteIfValid(out.hSub);

    if (out.hNSelectedCrossings)
    {
        out.hNSelectedCrossings->SetBinContent(1, static_cast<double>(out.nSelectedCrossings));
        out.hNSelectedCrossings->SetBinError(1, 0.0);
        out.hNSelectedCrossings->Write(out.hNSelectedCrossings->GetName(), TObject::kOverwrite);
    }
    else
    {
        std::string nNSelectedCrossings = "tklCombiDR_NSelectedCrossings";
        if (out.hSig)
        {
            const std::string sigName = out.hSig->GetName();
            const std::string sigSuffix = "_Sig";
            if (sigName.size() > sigSuffix.size() && sigName.compare(sigName.size() - sigSuffix.size(), sigSuffix.size(), sigSuffix) == 0)
                nNSelectedCrossings = sigName.substr(0, sigName.size() - sigSuffix.size()) + "_NSelectedCrossings";
        }
        TH1D hTmp(nNSelectedCrossings.c_str(), nNSelectedCrossings.c_str(), 1, 0.5, 1.5);
        hTmp.GetXaxis()->SetBinLabel(1, "selected_crossings");
        hTmp.GetYaxis()->SetTitle("Counts");
        hTmp.SetBinContent(1, static_cast<double>(out.nSelectedCrossings));
        hTmp.SetBinError(1, 0.0);
        hTmp.Write(hTmp.GetName(), TObject::kOverwrite);
    }

    WriteIfValid(out.hSub_AbsDphiCount);
    WriteIfValid(out.hSeedPhi);
    WriteIfValid(out.hDCA3D_Sig);
    WriteIfValid(out.hDCA3D_Bkg);
    WriteIfValid(out.hDCA3DHelix_Sig);
    WriteIfValid(out.hDCA3DHelix_Bkg);
    WriteIfValid(out.hSilSeedDCA3D);
    for (auto *h : out.hSilSeedDCA3D_truthVtxZ)
        WriteIfValid(h);
    WriteIfValid(out.hSilSeed_Eta);
    WriteIfValid(out.hSilSeed_Eta_nMVTX3nINTT2);
    WriteIfValid(out.hSilSeed_Eta_nMVTX3nINTT1);
    WriteIfValid(out.hSilSeed_Eta_nMVTX2nINTT2);
    WriteIfValid(out.hSilSeed_Eta_nMVTX2nINTT1);
    WriteIfValid(out.hSilSeed_Eta_nMVTXnINTTOther);
    WriteIfValid(out.hAllSilSeed_Eta);
    WriteIfValid(out.hAllSilSeed_Eta_nMVTX3nINTT2);
    WriteIfValid(out.hAllSilSeed_Eta_nMVTX3nINTT1);
    WriteIfValid(out.hAllSilSeed_Eta_nMVTX2nINTT2);
    WriteIfValid(out.hAllSilSeed_Eta_nMVTX2nINTT1);
    WriteIfValid(out.hAllSilSeed_Eta_nMVTXnINTTOther);
    WriteIfValid(out.hSeed_Eta_Phi);
    WriteIfValid(out.hSilSeed_Eta_Phi);
    WriteIfValid(out.hSilSeed_Eta_Phi_nMVTX3nINTT2);
    WriteIfValid(out.hSilSeed_Eta_Phi_nMVTX3nINTT1);
    WriteIfValid(out.hSilSeed_Eta_Phi_nMVTX2nINTT2);
    WriteIfValid(out.hSilSeed_Eta_Phi_nMVTX2nINTT1);
    WriteIfValid(out.hSilSeed_Eta_Phi_nMVTXnINTTOther);
    WriteIfValid(out.hINTTClusterVsAssociatedSilSeed);
    WriteIfValid(out.hINTTClusterVsAllSilSeed);
    WriteIfValid(out.hSeedCrossing_selectedCrossings);
    WriteIfValid(out.hSeedCrossing_selectedCrossings_largeDiffNClusterSeeds);
    WriteIfValid(out.hPrimaryCharged_Eta_Phi);

    const std::string output_prefix = TrackletCombiDROutputPrefix(out);
    const std::string dir_prefix = preserve_legacy_dirnames ? "" : (output_prefix + "_");

    auto write_h1_collection = [base_dir](const std::vector<TH1D *> &hists, const std::string &dirname)
    {
        base_dir->cd();
        TDirectory *dir = base_dir->GetDirectory(dirname.c_str());
        if (!dir)
            dir = base_dir->mkdir(dirname.c_str());
        if (!dir)
            return;
        dir->cd();
        for (auto *h : hists)
        {
            if (!h)
                continue;
            h->Write(h->GetName(), TObject::kOverwrite);
        }
        base_dir->cd();
    };

    write_h1_collection(out.hDR_Sig, dir_prefix + "dR_Sig");
    write_h1_collection(out.hDR_Bkg, dir_prefix + "dR_Bkg");
    write_h1_collection(out.hDR_Sub, dir_prefix + "dR_Sub");

    if (!out.crossingDisplays.empty())
    {
        const std::string display_dir_name = preserve_legacy_dirnames ? "CrossingDisplays" : (output_prefix + "_CrossingDisplays");
        TDirectory *crossingDisplayDir = base_dir->GetDirectory(display_dir_name.c_str());
        if (!crossingDisplayDir)
            crossingDisplayDir = base_dir->mkdir(display_dir_name.c_str());
        if (crossingDisplayDir)
        {
            for (const auto &display : out.crossingDisplays)
            {
                if (display.crossing_key.empty())
                    continue;

                TDirectory *crossingDir = crossingDisplayDir->GetDirectory(display.crossing_key.c_str());
                if (!crossingDir)
                    crossingDir = crossingDisplayDir->mkdir(display.crossing_key.c_str());
                if (!crossingDir)
                    continue;

                crossingDir->cd();

                if (display.grAllINTTClustersXY)
                    display.grAllINTTClustersXY->Write(display.grAllINTTClustersXY->GetName(), TObject::kOverwrite);
                if (display.grAllINTTClustersZR)
                    display.grAllINTTClustersZR->Write(display.grAllINTTClustersZR->GetName(), TObject::kOverwrite);
                if (display.grPrimaryINTTClustersXY)
                    display.grPrimaryINTTClustersXY->Write(display.grPrimaryINTTClustersXY->GetName(), TObject::kOverwrite);
                if (display.grPrimaryINTTClustersZR)
                    display.grPrimaryINTTClustersZR->Write(display.grPrimaryINTTClustersZR->GetName(), TObject::kOverwrite);
                if (display.grNonPrimaryINTTClustersXY)
                    display.grNonPrimaryINTTClustersXY->Write(display.grNonPrimaryINTTClustersXY->GetName(), TObject::kOverwrite);
                if (display.grNonPrimaryINTTClustersZR)
                    display.grNonPrimaryINTTClustersZR->Write(display.grNonPrimaryINTTClustersZR->GetName(), TObject::kOverwrite);
                if (display.grAssociatedSilSeedPCAXY)
                    display.grAssociatedSilSeedPCAXY->Write(display.grAssociatedSilSeedPCAXY->GetName(), TObject::kOverwrite);
                if (display.grAssociatedSilSeedPCAZR)
                    display.grAssociatedSilSeedPCAZR->Write(display.grAssociatedSilSeedPCAZR->GetName(), TObject::kOverwrite);
                if (display.grRecoVertexXY)
                    display.grRecoVertexXY->Write(display.grRecoVertexXY->GetName(), TObject::kOverwrite);
                if (display.grRecoVertexZR)
                    display.grRecoVertexZR->Write(display.grRecoVertexZR->GetName(), TObject::kOverwrite);
                if (display.grTruthVertexXY)
                    display.grTruthVertexXY->Write(display.grTruthVertexXY->GetName(), TObject::kOverwrite);
                if (display.grTruthVertexZR)
                    display.grTruthVertexZR->Write(display.grTruthVertexZR->GetName(), TObject::kOverwrite);

                for (auto *g : display.grAssociatedSilSeedClustersXY)
                {
                    if (g)
                        g->Write(g->GetName(), TObject::kOverwrite);
                }
                for (auto *g : display.grAssociatedSilSeedClustersZR)
                {
                    if (g)
                        g->Write(g->GetName(), TObject::kOverwrite);
                }
                for (auto *g : display.grNonAssociatedSilSeedClustersXY)
                {
                    if (g)
                        g->Write(g->GetName(), TObject::kOverwrite);
                }
                for (auto *g : display.grNonAssociatedSilSeedClustersZR)
                {
                    if (g)
                        g->Write(g->GetName(), TObject::kOverwrite);
                }
                for (auto *g : display.grPrimaryTruthClustersXY)
                {
                    if (g)
                        g->Write(g->GetName(), TObject::kOverwrite);
                }
                for (auto *g : display.grPrimaryTruthClustersZR)
                {
                    if (g)
                        g->Write(g->GetName(), TObject::kOverwrite);
                }
                for (auto *g : display.grPrimaryRecoClustersXY)
                {
                    if (g)
                        g->Write(g->GetName(), TObject::kOverwrite);
                }
                for (auto *g : display.grPrimaryRecoClustersZR)
                {
                    if (g)
                        g->Write(g->GetName(), TObject::kOverwrite);
                }

                base_dir->cd();
                crossingDisplayDir->cd();
            }
            base_dir->cd();
        }
    }
}

static inline bool WriteTracklets_CombinatoricDR_Output(const TrackletCombiDROutput &out, const std::string &outfilename)
{
    TFile fout(outfilename.c_str(), "RECREATE");
    if (fout.IsZombie())
        return false;

    WriteTracklets_CombinatoricDR_OutputToDirectory(out, &fout, true);
    fout.Write();
    fout.Close();
    return true;
}

static inline bool WriteTracklets_CombinatoricDR_Output(const TrackletCombiDROutput &inclusive_out, const std::vector<TrackletCombiDROutput> &multiplicity_outputs, const std::string &outfilename)
{
    TFile fout(outfilename.c_str(), "RECREATE");
    if (fout.IsZombie())
        return false;

    WriteTracklets_CombinatoricDR_OutputToDirectory(inclusive_out, &fout, true);
    for (const auto &out : multiplicity_outputs)
    {
        WriteTracklets_CombinatoricDR_OutputToDirectory(out, &fout, false);
    }

    fout.Write();
    fout.Close();
    return true;
}
