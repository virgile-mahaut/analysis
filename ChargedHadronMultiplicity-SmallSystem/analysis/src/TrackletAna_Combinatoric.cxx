//

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TMath.h>
#include <TObjString.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TTreeIndex.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <array>
#include <limits>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "GenHadron.h"
#include "LightHelixFit.h"
#include "Tracklet_combinatoric.h"
#include "Vertex.h"
#include "pdgidfunc.h"

bool verbose_debug = true;
bool savedetail = true;
bool dotruthmatching = false;
bool use_sPHENIX_primary_definition = true; // false: use original AuAu dNdEta primary definition

// ---------------------------------------------------------------------------
// Branch registration helpers — replaces repetitive SetBranchAddress calls.
// ---------------------------------------------------------------------------

// T* covers both scalar branches (&scalar) and vector branches (&vec_ptr where vec_ptr is T*).
// SetBranchAddress takes void* internally so one indirection level is sufficient for both.
template<typename T>
void RegisterBranch(TTree* t, const char* name, T* ptr)
{
    t->SetBranchAddress(name, ptr);
}

// Registers only if the branch exists. Returns true on success.
// Use for truly optional single branches (e.g. silseed_crossing).
// For all-or-nothing groups, check existence explicitly then call RegisterBranch.
template<typename T>
bool TryRegisterBranch(TTree* t, const char* name, T* ptr)
{
    if (!t->GetBranch(name)) return false;
    t->SetBranchAddress(name, ptr);
    return true;
}

int main(int argc, char *argv[])
{
    if (argc < 8 || argc > 10)
    {
        cout << "Usage: ./TrackletAna_Combinatoric [isdata] [infile] [outfile] [NevtToRun] [pairCut] [clusadccutset] [clusphisizecutset] [pairCutMode(optional: dR|dPhi|absdPhi)] [minCrossingSilSeeds(optional)]" << endl;
        exit(1);
    }

    std::cout << "Primary definition = " << (use_sPHENIX_primary_definition ? "sPHENIXPrimary" : "PrimaryPHG4Ptcl (AuAu dNdEta)") << std::endl;

    for (int i = 0; i < argc; i++)
        cout << "argv[" << i << "] = " << argv[i] << endl;

    bool IsData          = (TString(argv[1]).Atoi() == 1);
    TString infilename   = TString(argv[2]);
    TString outfilename  = TString(argv[3]);
    int NevtToRun_       = TString(argv[4]).Atoi();
    float pairCut        = TString(argv[5]).Atof();
    int clusadccutset    = TString(argv[6]).Atoi();
    int clusphisizecutset = TString(argv[7]).Atoi();

    TString idxstr = IsData ? "gl1bco" : "counter";

    TrackletData tkldata = {};

    // Setup for combinatoric tracklet method
    TrackletCombiDRConfig cfg;
    cfg.pair_cut = pairCut;

    auto parse_int_arg = [](const char *arg, const char *name) -> int
    {
        char *endptr = nullptr;
        const long value = std::strtol(arg, &endptr, 10);
        if (!endptr || *endptr != '\0')
            throw std::invalid_argument(std::string("Invalid integer for ") + name + ": " + arg);
        return static_cast<int>(value);
    };

    if (argc >= 9)
    {
        const std::string arg8 = argv[8];
        try
        {
            cfg.pair_cut_mode = ParsePairCutMode(arg8);
        }
        catch (const std::exception &e)
        {
            if (argc == 9)
            {
                try
                {
                    cfg.min_vertex_associated_silseeds = parse_int_arg(argv[8], "minCrossingSilSeeds");
                }
                catch (const std::exception &seed_e)
                {
                    std::cerr << e.what() << std::endl;
                    std::cerr << seed_e.what() << std::endl;
                    return 1;
                }
            }
            else
            {
                std::cerr << e.what() << std::endl;
                return 1;
            }
        }
    }
    if (argc == 10)
    {
        try
        {
            cfg.min_vertex_associated_silseeds = parse_int_arg(argv[9], "minCrossingSilSeeds");
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << std::endl;
            return 1;
        }
    }

    // Output (eta,phi) grid binning: uniform or edges
    cfg.etaAxis = AxisBinning(-1.1, 1.1, 11);
    cfg.phiAxis = AxisBinning(-3.2, 3.2, 16);
    // Background rotation settings
    cfg.enable_eta_preselect = false;
    cfg.rotate_search_layer  = true; // rotate logical search layer ({5,6})
    cfg.bkg_scale            = 1.0;
    cfg.abs_dphi_count_cut   = 0.15; // default PHOBOS-like scale; tune as needed
    // Pair-variable hist settings per (eta,phi) cell
    cfg.dR_nbins = 120;
    if (cfg.pair_cut_mode == TrackletCombiDRConfig::PairCutMode::kDeltaPhi)
    {
        cfg.dR_min = -cfg.pair_cut;
        cfg.dR_max =  cfg.pair_cut;
    }
    else
    {
        cfg.dR_min = 0.0;
        cfg.dR_max = cfg.pair_cut;
    }
    cfg.enable_zvtx_cut               = true;
    cfg.zvtx_cut_min                  = -10.0;
    cfg.zvtx_cut_max                  =  10.0;
    cfg.use_truth_vertex_for_zvtx_cut = true; // set to true for simulation study with truth vertex only for now
    if (cfg.use_truth_vertex_for_zvtx_cut)
        std::cout << "Using truth vertex for z vertex cut with range [" << cfg.zvtx_cut_min << ", " << cfg.zvtx_cut_max << "] cm" << std::endl;

    if (IsData && cfg.use_truth_vertex_for_zvtx_cut)
    {
        std::cout << "[WARNING] use_truth_vertex_for_zvtx_cut currently set to true but running with data. Setting to false" << std::endl;
        cfg.use_truth_vertex_for_zvtx_cut = false;
    }

    const std::string output_prefix                     = "tkl_Combinatoric";
    const std::string multiplicity_boundary_prefix      = IsData ? "Data_" : "MC_";
    const std::string multiplicity_percentile_boundary_file = "figure/figure-NInttClusterCrossing/" + multiplicity_boundary_prefix + "NInttClusterPercentileBoundaries.root";
    const std::vector<MultiplicityPercentileBin> multiplicity_percentile_bins = LoadMultiplicityPercentileBins(multiplicity_percentile_boundary_file);

    std::cout << "Loaded " << multiplicity_percentile_bins.size() << " multiplicity percentile bins from " << multiplicity_percentile_boundary_file << std::endl;
    for (const auto &bin : multiplicity_percentile_bins)
    {
        std::cout << "  multPercentileBin" << bin.percentile_bin << ": [" << bin.cluster_low << ", ";
        if (bin.cluster_high == std::numeric_limits<int>::max())
            std::cout << "INT_MAX";
        else
            std::cout << bin.cluster_high;
        std::cout << ") INTT clusters" << std::endl;
    }

    TrackletCombiDROutput out;
    InitTrackletCombiDROutput(out, cfg, output_prefix, !IsData, true);

    std::vector<TrackletCombiDROutput> multiplicity_percentile_outputs(multiplicity_percentile_bins.size());
    for (size_t i = 0; i < multiplicity_percentile_bins.size(); ++i)
        InitTrackletCombiDROutput(multiplicity_percentile_outputs[i], cfg, MultiplicityPercentileBinPrefix(output_prefix, multiplicity_percentile_bins[i].percentile_bin), !IsData, false);

    std::cout << "Pair cut mode = " << PairCutModeLabel(cfg.pair_cut_mode) << ", threshold = " << cfg.pair_cut << std::endl;
    std::cout << "Direct |delta-phi| count cut for merged (eta, phi) subtracted histogram = " << cfg.abs_dphi_count_cut << std::endl;
    std::cout << "Minimum silicon seeds associated to selected crossing vertex = " << cfg.min_vertex_associated_silseeds;
    if (cfg.min_vertex_associated_silseeds >= 3)
        std::cout << " (active)";
    else
        std::cout << " (inactive; vertex finder already requires >=2)";
    std::cout << std::endl;

    // cluster adc cut
    int clusadccut       = ConstADCCut(clusadccutset);
    float clusphisizecut = ConstClusPhiCut(clusphisizecutset);
    std::cout << "Cluster ADC cut (greater than): " << clusadccut << "; Cluster Phi-size cut (smaller than): " << clusphisizecut << std::endl;

    // -----------------------------------------------------------------------
    // Display graph helpers — used exclusively in the crossing display section.
    // -----------------------------------------------------------------------

    // Allocates a TGraph with name and title set in one call.
    auto MakeGraph = [](const char* name, const char* title) -> TGraph*
    {
        auto* g = new TGraph();
        g->SetName(name);
        g->SetTitle(title);
        return g;
    };

    // Adds a single point to an (x-y) and a (z, r_signed) display TGraph pair.
    // r_signed = ±hypot(x,y), sign follows y.
    auto AddPoint = [](TGraph* grXY, TGraph* grZR, double x, double y, double z)
    {
        const double r_signed = (y < 0.0) ? -std::hypot(x, y) : std::hypot(x, y);
        grXY->SetPoint(grXY->GetN(), x, y);
        grZR->SetPoint(grZR->GetN(), z, r_signed);
    };

    // Iterates parallel coordinate vectors and fills two display TGraph pairs.
    // Skips NaN entries (uniform treatment across truth and reco clusters).
    // Returns the number of points actually added.
    auto FillXYZR = [](TGraph* grXY, TGraph* grZR,
                       const std::vector<float>& xs,
                       const std::vector<float>& ys,
                       const std::vector<float>& zs) -> int
    {
        const size_t n = std::min({xs.size(), ys.size(), zs.size()});
        int added = 0;
        for (size_t ic = 0; ic < n; ++ic)
        {
            const double x = xs[ic], y = ys[ic], z = zs[ic];
            if (std::isnan(x) || std::isnan(y) || std::isnan(z)) continue;
            const double r_signed = (y < 0.0) ? -std::hypot(x, y) : std::hypot(x, y);
            grXY->SetPoint(grXY->GetN(), x, y);
            grZR->SetPoint(grZR->GetN(), z, r_signed);
            ++added;
        }
        return added;
    };

    // -----------------------------------------------------------------------
    // Open input file and tree
    // -----------------------------------------------------------------------
    TFile *f = new TFile(infilename, "READ");
    TTree *t = (TTree *)f->Get("VTX");
    t->BuildIndex(idxstr); // Reference: https://root-forum.cern.ch/t/sort-ttree-entries/13138
    TTreeIndex *index = (TTreeIndex *)t->GetTreeIndex();
    const Long64_t nIndexed = index ? index->GetN() : 0;
    if (nIndexed <= 0)
    {
        std::cerr << "No entries available in tree index." << std::endl;
        f->Close(); delete f; return 1;
    }
    if (NevtToRun_ > nIndexed)
    {
        std::cout << "Requested NevtToRun (" << NevtToRun_ << ") exceeds indexed entries (" << nIndexed << "), clamping." << std::endl;
        NevtToRun_ = static_cast<int>(nIndexed);
    }

    // -----------------------------------------------------------------------
    // Branch variable declarations
    // -----------------------------------------------------------------------
    int event, nTruthVertex;
    std::vector<int> *firedTriggers = 0;
    uint64_t gl1bco, bcotr;
    bool is_min_bias;
    bool InttBco_IsToBeRemoved = false; // For this study we don't remove any BCO
    float mbd_south_charge_sum, mbd_north_charge_sum, mbd_charge_sum, mbd_charge_asymm, mbd_z_vtx;

    // Vertex information
    int nSvtxVertices;
    std::vector<int>       *trackerVertexId       = 0;
    std::vector<float>     *trackerVertexX        = 0, *trackerVertexY        = 0, *trackerVertexZ        = 0;
    std::vector<int>       *trackerVertexNTracks   = 0;
    std::vector<std::vector<int>> *trackerVertexTrackIDs = 0;
    std::vector<short int> *trackerVertexCrossing  = 0;

    // Truth vertex (simulation only)
    std::vector<float> *TruthVertexX = 0, *TruthVertexY = 0, *TruthVertexZ = 0;

    // Silicon seeds
    std::vector<float> *silseed_x = 0, *silseed_y = 0, *silseed_z = 0;
    std::vector<float> *silseed_eta_vtx = 0, *silseed_phi_vtx = 0;
    std::vector<int>   *silseed_nMvtx = 0, *silseed_nIntt = 0, *silseed_crossing = 0;
    std::vector<std::vector<float>> *silseed_cluster_globalX = 0, *silseed_cluster_globalY = 0, *silseed_cluster_globalZ = 0;

    // Clusters
    std::vector<int>          *ClusLayer = 0;
    std::vector<float>        *ClusX = 0, *ClusY = 0, *ClusZ = 0, *ClusPhiSize = 0, *ClusZSize = 0;
    std::vector<unsigned int> *ClusAdc = 0;
    std::vector<int>          *cluster_timeBucketID = 0, *ClusLadderZId = 0;

    // Cluster-to-G4Particle truth matching
    std::vector<int>   *cluster_matchedG4P_trackID = 0, *cluster_matchedG4P_PID = 0;
    std::vector<float> *cluster_matchedG4P_E = 0, *cluster_matchedG4P_pT = 0, *cluster_matchedG4P_eta = 0, *cluster_matchedG4P_phi = 0;

    // sPHENIX primary definition
    std::vector<int>    *sPHENIXPrimary_PID = 0, *sPHENIXPrimary_trackID = 0;
    std::vector<float>  *sPHENIXPrimary_pT = 0, *sPHENIXPrimary_eta = 0, *sPHENIXPrimary_phi = 0, *sPHENIXPrimary_E = 0;
    std::vector<double> *sPHENIXPrimary_charge = 0;
    std::vector<bool>   *sPHENIXPrimary_isChargedHadron = 0;

    // Original AuAu dNdEta primary definition
    std::vector<int>    *PrimaryPHG4Ptcl_PID = 0, *PrimaryPHG4Ptcl_trackID = 0;
    std::vector<float>  *PrimaryPHG4Ptcl_pT = 0, *PrimaryPHG4Ptcl_eta = 0, *PrimaryPHG4Ptcl_phi = 0, *PrimaryPHG4Ptcl_E = 0;
    std::vector<double> *PrimaryPHG4Ptcl_charge = 0;
    std::vector<bool>   *PrimaryPHG4Ptcl_isChargedHadron = 0;

    // sPHENIX primary truth/reco cluster matching
    std::vector<std::vector<float>> *sPHENIXPrimary_truthcluster_X = 0, *sPHENIXPrimary_truthcluster_Y = 0, *sPHENIXPrimary_truthcluster_Z = 0;
    std::vector<std::vector<float>> *sPHENIXPrimary_recocluster_globalX = 0, *sPHENIXPrimary_recocluster_globalY = 0, *sPHENIXPrimary_recocluster_globalZ = 0;

    std::vector<int>   *G4P_trackID = 0;
    std::vector<float> *G4P_Pt = 0, *G4P_Eta = 0, *G4P_Phi = 0, *G4P_E = 0;

    // -----------------------------------------------------------------------
    // Branch registration
    // -----------------------------------------------------------------------
    RegisterBranch(t, "counter",               &event);
    RegisterBranch(t, "gl1bco",                &gl1bco);
    RegisterBranch(t, "bcotr",                 &bcotr);
    RegisterBranch(t, "nSvtxVertices",         &nSvtxVertices);
    RegisterBranch(t, "trackerVertexId",       &trackerVertexId);
    RegisterBranch(t, "trackerVertexX",        &trackerVertexX);
    RegisterBranch(t, "trackerVertexY",        &trackerVertexY);
    RegisterBranch(t, "trackerVertexZ",        &trackerVertexZ);
    RegisterBranch(t, "trackerVertexNTracks",  &trackerVertexNTracks);
    RegisterBranch(t, "trackerVertexTrackIDs", &trackerVertexTrackIDs);
    RegisterBranch(t, "trackerVertexCrossing", &trackerVertexCrossing);
    RegisterBranch(t, "MBD_charge_sum",        &mbd_charge_sum);
    RegisterBranch(t, "cluster_layer",         &ClusLayer);
    RegisterBranch(t, "cluster_globalX",       &ClusX);
    RegisterBranch(t, "cluster_globalY",       &ClusY);
    RegisterBranch(t, "cluster_globalZ",       &ClusZ);
    RegisterBranch(t, "cluster_phiSize",       &ClusPhiSize);
    RegisterBranch(t, "cluster_zSize",         &ClusZSize);
    RegisterBranch(t, "cluster_adc",           &ClusAdc);
    RegisterBranch(t, "cluster_timeBucketID",  &cluster_timeBucketID);
    RegisterBranch(t, "cluster_ladderZId",     &ClusLadderZId);
    RegisterBranch(t, "silseed_eta_vtx",       &silseed_eta_vtx);
    RegisterBranch(t, "silseed_phi_vtx",       &silseed_phi_vtx);
    RegisterBranch(t, "silseed_nMvtx",         &silseed_nMvtx);
    RegisterBranch(t, "silseed_nIntt",         &silseed_nIntt);

    if (!IsData)
    {
        RegisterBranch(t, "nTruthVertex", &nTruthVertex);
        RegisterBranch(t, "TruthVertexX", &TruthVertexX);
        RegisterBranch(t, "TruthVertexY", &TruthVertexY);
        RegisterBranch(t, "TruthVertexZ", &TruthVertexZ);

        if (use_sPHENIX_primary_definition)
        {
            RegisterBranch(t, "sPHENIXPrimary_PID",             &sPHENIXPrimary_PID);
            RegisterBranch(t, "sPHENIXPrimary_isChargedHadron", &sPHENIXPrimary_isChargedHadron);
            RegisterBranch(t, "sPHENIXPrimary_trackID",         &sPHENIXPrimary_trackID);
            RegisterBranch(t, "sPHENIXPrimary_pT",              &sPHENIXPrimary_pT);
            RegisterBranch(t, "sPHENIXPrimary_eta",             &sPHENIXPrimary_eta);
            RegisterBranch(t, "sPHENIXPrimary_phi",             &sPHENIXPrimary_phi);
            RegisterBranch(t, "sPHENIXPrimary_E",               &sPHENIXPrimary_E);
        }
        else
        {
            RegisterBranch(t, "PrimaryPHG4Ptcl_PID",             &PrimaryPHG4Ptcl_PID);
            RegisterBranch(t, "PrimaryPHG4Ptcl_isChargedHadron", &PrimaryPHG4Ptcl_isChargedHadron);
            RegisterBranch(t, "PrimaryPHG4Ptcl_trackID",         &PrimaryPHG4Ptcl_trackID);
            RegisterBranch(t, "PrimaryPHG4Ptcl_pT",              &PrimaryPHG4Ptcl_pT);
            RegisterBranch(t, "PrimaryPHG4Ptcl_eta",             &PrimaryPHG4Ptcl_eta);
            RegisterBranch(t, "PrimaryPHG4Ptcl_phi",             &PrimaryPHG4Ptcl_phi);
            RegisterBranch(t, "PrimaryPHG4Ptcl_E",               &PrimaryPHG4Ptcl_E);
            RegisterBranch(t, "PrimaryPHG4Ptcl_charge",          &PrimaryPHG4Ptcl_charge);
        }

        // if (dotruthmatching)
        {
            RegisterBranch(t, "cluster_matchedG4P_trackID", &cluster_matchedG4P_trackID);
            RegisterBranch(t, "cluster_matchedG4P_PID",     &cluster_matchedG4P_PID);
            RegisterBranch(t, "cluster_matchedG4P_E",       &cluster_matchedG4P_E);
            RegisterBranch(t, "cluster_matchedG4P_pT",      &cluster_matchedG4P_pT);
            RegisterBranch(t, "cluster_matchedG4P_eta",     &cluster_matchedG4P_eta);
            RegisterBranch(t, "cluster_matchedG4P_phi",     &cluster_matchedG4P_phi);
        }

        RegisterBranch(t, "sPHENIXPrimary_charge",              &sPHENIXPrimary_charge);
        RegisterBranch(t, "sPHENIXPrimary_truthcluster_X",      &sPHENIXPrimary_truthcluster_X);
        RegisterBranch(t, "sPHENIXPrimary_truthcluster_Y",      &sPHENIXPrimary_truthcluster_Y);
        RegisterBranch(t, "sPHENIXPrimary_truthcluster_Z",      &sPHENIXPrimary_truthcluster_Z);
        RegisterBranch(t, "sPHENIXPrimary_recocluster_globalX", &sPHENIXPrimary_recocluster_globalX);
        RegisterBranch(t, "sPHENIXPrimary_recocluster_globalY", &sPHENIXPrimary_recocluster_globalY);
        RegisterBranch(t, "sPHENIXPrimary_recocluster_globalZ", &sPHENIXPrimary_recocluster_globalZ);

        InttBco_IsToBeRemoved = false;
    }
    if (IsData)
        RegisterBranch(t, "firedTriggers", &firedTriggers);

    // Optional branches — register only if present; presence flag drives downstream logic.
    const bool has_silseed_crossing_branch = TryRegisterBranch(t, "silseed_crossing", &silseed_crossing);
    if (!has_silseed_crossing_branch)
        std::cout << "[INFO] silseed_crossing branch is not available; all-seed crossing histograms will remain empty." << std::endl;

    // XYZ branch groups are all-or-nothing: check all three exist before registering any.
    const bool has_silseed_xyz_branches =
        t->GetBranch("silseed_x") && t->GetBranch("silseed_y") && t->GetBranch("silseed_z");
    if (has_silseed_xyz_branches)
    {
        RegisterBranch(t, "silseed_x", &silseed_x);
        RegisterBranch(t, "silseed_y", &silseed_y);
        RegisterBranch(t, "silseed_z", &silseed_z);
    }
    else
    {
        std::cout << "[INFO] silseed_x/y/z branches are not available; seed-PCA display graphs will be skipped." << std::endl;
    }

    const bool has_silseed_cluster_xyz_branches =
        t->GetBranch("silseed_cluster_globalX") && t->GetBranch("silseed_cluster_globalY") && t->GetBranch("silseed_cluster_globalZ");
    if (has_silseed_cluster_xyz_branches)
    {
        RegisterBranch(t, "silseed_cluster_globalX", &silseed_cluster_globalX);
        RegisterBranch(t, "silseed_cluster_globalY", &silseed_cluster_globalY);
        RegisterBranch(t, "silseed_cluster_globalZ", &silseed_cluster_globalZ);
    }
    else
    {
        std::cout << "[INFO] silseed_cluster_globalX/Y/Z branches are not available; per-seed cluster display graphs will be skipped." << std::endl;
    }

    // -----------------------------------------------------------------------
    // Event loop
    // -----------------------------------------------------------------------
    std::unordered_map<int, std::vector<int>> vertexCrossing_idx_map; // key: vertex crossing, value: vertex indices
    std::unordered_map<int, std::vector<size_t>> silseedCrossing_idx_map; // key: silicon-seed crossing, value: seed indices in this event

    for (int i = 0; i < NevtToRun_; i++)
    {
        Long64_t local = t->LoadTree(index->GetIndex()[i]);
        t->GetEntry(local);

        std::vector<float> *selectedPrimaryEta             = use_sPHENIX_primary_definition ? sPHENIXPrimary_eta             : PrimaryPHG4Ptcl_eta;
        std::vector<float> *selectedPrimaryPhi             = use_sPHENIX_primary_definition ? sPHENIXPrimary_phi             : PrimaryPHG4Ptcl_phi;
        std::vector<bool>  *selectedPrimaryIsChargedHadron = use_sPHENIX_primary_definition ? sPHENIXPrimary_isChargedHadron : PrimaryPHG4Ptcl_isChargedHadron;
        bool filledPrimaryChargedAfterZvtxCut = false;
        std::vector<char> filledPrimaryChargedAfterZvtxCutByMultBin(multiplicity_percentile_outputs.size(), 0);

        silseedCrossing_idx_map.clear();
        if (has_silseed_crossing_branch && silseed_crossing)
        {
            for (size_t iseed = 0; iseed < silseed_crossing->size(); ++iseed)
                silseedCrossing_idx_map[silseed_crossing->at(iseed)].push_back(iseed);
        }

        auto fill_selected_crossing_histograms = [&](TrackletCombiDROutput &target_out, auto &filled_primary_charged, int crossing_value, int vtx_idx_value, const std::vector<size_t> &used_seed_indices_value, int n_associated_silseeds_value, int n_intt_clusters_in_crossing_value, int n_all_silseeds_in_crossing_value, double truth_vtxz_value, bool has_truth_vtxz_value)
        {
            if (target_out.hSeedCrossing_selectedCrossings)
                target_out.hSeedCrossing_selectedCrossings->Fill(crossing_value, 1.0);

            if (target_out.hSeedCrossing_selectedCrossings_largeDiffNClusterSeeds && n_intt_clusters_in_crossing_value > 4 * n_associated_silseeds_value)
                target_out.hSeedCrossing_selectedCrossings_largeDiffNClusterSeeds->Fill(crossing_value, 1.0);

            if (!filled_primary_charged && !IsData && target_out.hPrimaryCharged_Eta_Phi && selectedPrimaryEta && selectedPrimaryPhi && selectedPrimaryIsChargedHadron)
            {
                const size_t nprim = std::min({selectedPrimaryEta->size(), selectedPrimaryPhi->size(), selectedPrimaryIsChargedHadron->size()});
                for (size_t iprim = 0; iprim < nprim; ++iprim)
                {
                    if (!selectedPrimaryIsChargedHadron->at(iprim)) continue;
                    const double eta_primary = selectedPrimaryEta->at(iprim);
                    double phi_primary = selectedPrimaryPhi->at(iprim);
                    phi_primary = target_out.phiAxisIs_0_2pi ? wrapPhi_0_2pi(phi_primary) : wrapPhi_mpi_pi(phi_primary);
                    target_out.hPrimaryCharged_Eta_Phi->Fill(eta_primary, phi_primary, 1.0);
                }
                filled_primary_charged = true;
            }

            if (target_out.hSilSeed_Eta_Phi && trackerVertexTrackIDs && vtx_idx_value >= 0 && vtx_idx_value < static_cast<int>(trackerVertexTrackIDs->size()) && silseed_eta_vtx && silseed_phi_vtx)
            {
                for (const size_t iseed : used_seed_indices_value)
                {
                    if (iseed >= silseed_eta_vtx->size() || iseed >= silseed_phi_vtx->size()) continue;

                    const double eta_seed_vtx = silseed_eta_vtx->at(iseed);
                    const double phi_seed_vtx_raw = silseed_phi_vtx->at(iseed);
                    if (!std::isfinite(eta_seed_vtx) || !std::isfinite(phi_seed_vtx_raw) || std::abs(phi_seed_vtx_raw) > 10.0 * TMath::Pi())
                    {
                        continue;
                    }
                    double phi_seed_vtx = target_out.phiAxisIs_0_2pi ? wrapPhi_0_2pi(phi_seed_vtx_raw) : wrapPhi_mpi_pi(phi_seed_vtx_raw);
                    target_out.hSilSeed_Eta_Phi->Fill(eta_seed_vtx, phi_seed_vtx, 1.0);
                    if (target_out.hSilSeed_Eta)
                        target_out.hSilSeed_Eta->Fill(eta_seed_vtx, 1.0);

                    if (target_out.hSilSeedDCA3D && has_silseed_xyz_branches && silseed_x && silseed_y && silseed_z && iseed < silseed_x->size() && iseed < silseed_y->size() && iseed < silseed_z->size())
                    {
                        const double dx = silseed_x->at(iseed) - trackerVertexX->at(vtx_idx_value);
                        const double dy = silseed_y->at(iseed) - trackerVertexY->at(vtx_idx_value);
                        const double dz = silseed_z->at(iseed) - trackerVertexZ->at(vtx_idx_value);
                        const double silseed_dca3d = std::sqrt(dx * dx + dy * dy + dz * dz);
                        target_out.hSilSeedDCA3D->Fill(silseed_dca3d, 1.0);

                        if (has_truth_vtxz_value && !target_out.hSilSeedDCA3D_truthVtxZ.empty())
                        {
                            const int truth_vtxz_bin = FindTruthVtxZBin(cfg.truth_vtxz_bin_edges, truth_vtxz_value);
                            if (truth_vtxz_bin >= 0 && truth_vtxz_bin < static_cast<int>(target_out.hSilSeedDCA3D_truthVtxZ.size()) && target_out.hSilSeedDCA3D_truthVtxZ[truth_vtxz_bin])
                                target_out.hSilSeedDCA3D_truthVtxZ[truth_vtxz_bin]->Fill(silseed_dca3d, 1.0);
                        }
                    }

                    if (silseed_nMvtx && silseed_nIntt && iseed < silseed_nMvtx->size() && iseed < silseed_nIntt->size())
                    {
                        const int nMvtx = silseed_nMvtx->at(iseed);
                        const int nIntt = silseed_nIntt->at(iseed);

                        if      (nMvtx == 3 && nIntt == 2)
                        {
                            if (target_out.hSilSeed_Eta_Phi_nMVTX3nINTT2) target_out.hSilSeed_Eta_Phi_nMVTX3nINTT2->Fill(eta_seed_vtx, phi_seed_vtx, 1.0);
                            if (target_out.hSilSeed_Eta_nMVTX3nINTT2)     target_out.hSilSeed_Eta_nMVTX3nINTT2->Fill(eta_seed_vtx, 1.0);
                        }
                        else if (nMvtx == 3 && nIntt == 1)
                        {
                            if (target_out.hSilSeed_Eta_Phi_nMVTX3nINTT1) target_out.hSilSeed_Eta_Phi_nMVTX3nINTT1->Fill(eta_seed_vtx, phi_seed_vtx, 1.0);
                            if (target_out.hSilSeed_Eta_nMVTX3nINTT1)     target_out.hSilSeed_Eta_nMVTX3nINTT1->Fill(eta_seed_vtx, 1.0);
                        }
                        else if (nMvtx == 2 && nIntt == 2)
                        {
                            if (target_out.hSilSeed_Eta_Phi_nMVTX2nINTT2) target_out.hSilSeed_Eta_Phi_nMVTX2nINTT2->Fill(eta_seed_vtx, phi_seed_vtx, 1.0);
                            if (target_out.hSilSeed_Eta_nMVTX2nINTT2)     target_out.hSilSeed_Eta_nMVTX2nINTT2->Fill(eta_seed_vtx, 1.0);
                        }
                        else if (nMvtx == 2 && nIntt == 1)
                        {
                            if (target_out.hSilSeed_Eta_Phi_nMVTX2nINTT1) target_out.hSilSeed_Eta_Phi_nMVTX2nINTT1->Fill(eta_seed_vtx, phi_seed_vtx, 1.0);
                            if (target_out.hSilSeed_Eta_nMVTX2nINTT1)     target_out.hSilSeed_Eta_nMVTX2nINTT1->Fill(eta_seed_vtx, 1.0);
                        }
                        else
                        {
                            if (target_out.hSilSeed_Eta_Phi_nMVTXnINTTOther) target_out.hSilSeed_Eta_Phi_nMVTXnINTTOther->Fill(eta_seed_vtx, phi_seed_vtx, 1.0);
                            if (target_out.hSilSeed_Eta_nMVTXnINTTOther)     target_out.hSilSeed_Eta_nMVTXnINTTOther->Fill(eta_seed_vtx, 1.0);
                        }
                    }
                }
            }

            if (silseed_eta_vtx && silseed_nMvtx && silseed_nIntt)
            {
                const auto crossing_it = silseedCrossing_idx_map.find(crossing_value);
                if (crossing_it != silseedCrossing_idx_map.end())
                {
                    for (const size_t iseed : crossing_it->second)
                    {
                        if (iseed >= silseed_eta_vtx->size() || iseed >= silseed_nMvtx->size() || iseed >= silseed_nIntt->size()) continue;

                        const double eta_seed_vtx = silseed_eta_vtx->at(iseed);
                        const int nMvtx = silseed_nMvtx->at(iseed);
                        const int nIntt = silseed_nIntt->at(iseed);

                        if (target_out.hAllSilSeed_Eta)
                            target_out.hAllSilSeed_Eta->Fill(eta_seed_vtx, 1.0);

                        if      (nMvtx == 3 && nIntt == 2)
                        {
                            if (target_out.hAllSilSeed_Eta_nMVTX3nINTT2) target_out.hAllSilSeed_Eta_nMVTX3nINTT2->Fill(eta_seed_vtx, 1.0);
                        }
                        else if (nMvtx == 3 && nIntt == 1)
                        {
                            if (target_out.hAllSilSeed_Eta_nMVTX3nINTT1) target_out.hAllSilSeed_Eta_nMVTX3nINTT1->Fill(eta_seed_vtx, 1.0);
                        }
                        else if (nMvtx == 2 && nIntt == 2)
                        {
                            if (target_out.hAllSilSeed_Eta_nMVTX2nINTT2) target_out.hAllSilSeed_Eta_nMVTX2nINTT2->Fill(eta_seed_vtx, 1.0);
                        }
                        else if (nMvtx == 2 && nIntt == 1)
                        {
                            if (target_out.hAllSilSeed_Eta_nMVTX2nINTT1) target_out.hAllSilSeed_Eta_nMVTX2nINTT1->Fill(eta_seed_vtx, 1.0);
                        }
                        else
                        {
                            if (target_out.hAllSilSeed_Eta_nMVTXnINTTOther) target_out.hAllSilSeed_Eta_nMVTXnINTTOther->Fill(eta_seed_vtx, 1.0);
                        }
                    }
                }
            }

            if (target_out.hINTTClusterVsAssociatedSilSeed)
                target_out.hINTTClusterVsAssociatedSilSeed->Fill(n_intt_clusters_in_crossing_value, n_associated_silseeds_value, 1.0);

            if (target_out.hINTTClusterVsAllSilSeed && n_all_silseeds_in_crossing_value >= 0)
                target_out.hINTTClusterVsAllSilSeed->Fill(n_intt_clusters_in_crossing_value, n_all_silseeds_in_crossing_value, 1.0);
        };

        cout << "i = " << i << " event = " << event << " GL1 BCO = " << gl1bco << " NevtToRun_ = " << NevtToRun_ << " local = " << local << endl;

        // Build unique vertex crossing map for this trigger frame
        for (int j = 0; j < (int)trackerVertexCrossing->size(); j++)
            vertexCrossing_idx_map[trackerVertexCrossing->at(j)].push_back(j);
        cout << "Number of unique vertex crossings in this trigger frame = " << vertexCrossing_idx_map.size() << endl;

        // Only consider crossings with exactly 1 silicon vertex to avoid in-time pileup contamination.
        // Early-continue on the multi-vertex / sentinel case to reduce nesting.
        for (auto &pair : vertexCrossing_idx_map)
        {
            auto vertex_crossing = pair.first;
            if (pair.second.size() != 1 || vertex_crossing == std::numeric_limits<short int>::max())
            {
                std::cout << "Skipping crossing " << pair.first << " with multiple vertices" << std::endl;
                continue;
            }

            out.nSelectedCrossings++;

            int vtx_idx = pair.second[0];
            const bool reco_vertex_vectors_ok =
                trackerVertexNTracks && trackerVertexX && trackerVertexY && trackerVertexZ &&
                vtx_idx >= 0 &&
                vtx_idx < static_cast<int>(trackerVertexNTracks->size()) &&
                vtx_idx < static_cast<int>(trackerVertexX->size()) &&
                vtx_idx < static_cast<int>(trackerVertexY->size()) &&
                vtx_idx < static_cast<int>(trackerVertexZ->size());
            if (!reco_vertex_vectors_ok)
            {
                std::cerr << "[ERROR] Reconstructed vertex index out of range for crossing " << vertex_crossing
                          << ": vtx_idx=" << vtx_idx
                          << ", trackerVertexNTracks size=" << (trackerVertexNTracks ? trackerVertexNTracks->size() : 0)
                          << ", trackerVertexX size=" << (trackerVertexX ? trackerVertexX->size() : 0)
                          << ", trackerVertexY size=" << (trackerVertexY ? trackerVertexY->size() : 0)
                          << ", trackerVertexZ size=" << (trackerVertexZ ? trackerVertexZ->size() : 0)
                          << std::endl;
                continue;
            }
            std::cout << "Vertex crossing " << vertex_crossing << " has 1 vertex with index " << vtx_idx << " and Ntracks = " << trackerVertexNTracks->at(vtx_idx) << std::endl;

            // Truth vertex validity check
            const bool truth_vertex_vectors_ok =
                TruthVertexX && TruthVertexY && TruthVertexZ &&
                vtx_idx >= 0 &&
                vtx_idx < static_cast<int>(TruthVertexX->size()) &&
                vtx_idx < static_cast<int>(TruthVertexY->size()) &&
                vtx_idx < static_cast<int>(TruthVertexZ->size());
            if (cfg.use_truth_vertex_for_zvtx_cut && !truth_vertex_vectors_ok)
            {
                std::cerr << "[ERROR] Truth vertex index out of range for crossing " << vertex_crossing
                          << ": vtx_idx=" << vtx_idx
                          << ", TruthVertexX size=" << (TruthVertexX ? TruthVertexX->size() : 0)
                          << ", TruthVertexY size=" << (TruthVertexY ? TruthVertexY->size() : 0)
                          << ", TruthVertexZ size=" << (TruthVertexZ ? TruthVertexZ->size() : 0)
                          << ". Falling back to reconstructed vertex for z-vertex selection." << std::endl;
            }

            const bool use_truth_vertex_here = cfg.use_truth_vertex_for_zvtx_cut && truth_vertex_vectors_ok;
            float silvtx_x = use_truth_vertex_here ? TruthVertexX->at(vtx_idx) : trackerVertexX->at(vtx_idx);
            float silvtx_y = use_truth_vertex_here ? TruthVertexY->at(vtx_idx) : trackerVertexY->at(vtx_idx);
            float silvtx_z = use_truth_vertex_here ? TruthVertexZ->at(vtx_idx) : trackerVertexZ->at(vtx_idx);
            const std::array<double, 3> silicon_vtx = {silvtx_x, silvtx_y, silvtx_z};

            // Collect seed indices associated to this vertex via trackerVertexTrackIDs
            std::vector<size_t> used_seed_indices;
            if (trackerVertexTrackIDs && vtx_idx >= 0 && vtx_idx < static_cast<int>(trackerVertexTrackIDs->size()) && silseed_eta_vtx && silseed_phi_vtx)
            {
                const auto &vtx_track_ids = trackerVertexTrackIDs->at(vtx_idx);
                std::cout << "There are " << vtx_track_ids.size() << " Track IDs associated to this vertex." << std::endl;
                for (const int track_id : vtx_track_ids)
                    used_seed_indices.push_back(track_id);
            }
            std::cout << "Vertex crossing " << vertex_crossing << " has " << used_seed_indices.size() << " associated silicon seeds." << std::endl;

            // Z-vertex cut
            if (cfg.enable_zvtx_cut && (silvtx_z < cfg.zvtx_cut_min || silvtx_z > cfg.zvtx_cut_max))
            {
                std::cout << "Skipping vertex crossing " << vertex_crossing << " due to z vertex cut: " << silvtx_z << " not in [" << cfg.zvtx_cut_min << ", " << cfg.zvtx_cut_max << "]" << std::endl;
                continue;
            }

            // Collect INTT cluster indices belonging to this crossing
            std::vector<size_t> crossing_cluster_indices;
            crossing_cluster_indices.reserve(ClusLayer->size());
            for (size_t ihit = 0; ihit < ClusLayer->size(); ++ihit)
            {
                if (cluster_timeBucketID->at(ihit) != vertex_crossing) continue;
                if (logicalInttLayerIndex(ClusLayer->at(ihit)) < 0)   continue;
                crossing_cluster_indices.push_back(ihit);
            }

            const int n_associated_silseeds       = static_cast<int>(used_seed_indices.size());
            const int n_intt_clusters_in_crossing = static_cast<int>(crossing_cluster_indices.size());
            const int n_all_silseeds_in_crossing  = has_silseed_crossing_branch
                ? static_cast<int>(silseedCrossing_idx_map[vertex_crossing].size()) : -1;
            const bool   has_truth_vtxz = (!IsData && TruthVertexZ && vtx_idx >= 0 && vtx_idx < static_cast<int>(TruthVertexZ->size()));
            const double truth_vtxz     = has_truth_vtxz ? TruthVertexZ->at(vtx_idx) : 0.0;
            const int multiplicity_percentile_bin_index = FindMultiplicityPercentileBinIndex(multiplicity_percentile_bins, n_intt_clusters_in_crossing);

            if (multiplicity_percentile_bin_index >= 0)
                multiplicity_percentile_outputs[multiplicity_percentile_bin_index].nSelectedCrossings++;
            else
                std::cout << "Skipping multiplicity-binned fill for crossing " << vertex_crossing << " because N_{INTT clusters} = " << n_intt_clusters_in_crossing << " did not match any configured percentile interval." << std::endl;

            // Bind the 8 common call-site args so each fill site only specifies its target output + flag.
            auto fill_crossing = [&](TrackletCombiDROutput &target_out, auto &filled_flag)
            {
                fill_selected_crossing_histograms(target_out, filled_flag,
                    vertex_crossing, vtx_idx, used_seed_indices, n_associated_silseeds,
                    n_intt_clusters_in_crossing, n_all_silseeds_in_crossing, truth_vtxz, has_truth_vtxz);
            };
            fill_crossing(out, filledPrimaryChargedAfterZvtxCut);
            if (multiplicity_percentile_bin_index >= 0)
            {
                fill_crossing(multiplicity_percentile_outputs[multiplicity_percentile_bin_index],
                              filledPrimaryChargedAfterZvtxCutByMultBin[multiplicity_percentile_bin_index]);
            }

            // ---------------------------------------------------------------
            // Crossing display 
            // const bool make_crossing_display = (n_associated_silseeds >= 4 && n_intt_clusters_in_crossing > 4 * n_associated_silseeds);
            // const bool make_crossing_display = true; // enable for full debug display
            const bool make_crossing_display = false;
            // ---------------------------------------------------------------

            // Fill tkldata with clusters from this crossing — display is handled separately below
            for (const size_t ihit : crossing_cluster_indices)
            {
                const int logical_layer = logicalInttLayerIndex(ClusLayer->at(ihit));
                Hit *hit = new Hit(ClusX->at(ihit), ClusY->at(ihit), ClusZ->at(ihit), silvtx_x, silvtx_y, silvtx_z, logical_layer, ClusPhiSize->at(ihit), ClusAdc->at(ihit));
                hit->SetEdgeFromZId(ClusLadderZId->at(ihit));
                tkldata.layers[logical_layer].push_back(hit);
            }

            std::cout << "Number of clusters in this crossing (inner, outer) = (" << tkldata.layers[0].size() << ", " << tkldata.layers[1].size() << ")" << std::endl;

            // BuildCrossingDisplay owns the full display graph lifecycle in one contiguous block.
            // It re-iterates crossing_cluster_indices for per-cluster graphs rather than
            // interleaving display fills into the tkldata loop above.
            // display, nonassociated_seed_indices, and sphenix_primary_track_ids are local to it.
            auto BuildCrossingDisplay = [&]()
            {
                TrackletCombiDROutput::CrossingDisplayGraphs display;

                display.crossing_key = IsData
                    ? (std::to_string(gl1bco) + "_" + std::to_string(vertex_crossing))
                    : (std::to_string(event)  + "_" + std::to_string(vertex_crossing));

                // Identify seeds in this crossing that are NOT associated to the selected vertex
                std::vector<size_t> nonassociated_seed_indices;
                if (has_silseed_crossing_branch && silseed_crossing)
                {
                    std::vector<size_t> associated_sorted(used_seed_indices.begin(), used_seed_indices.end());
                    std::sort(associated_sorted.begin(), associated_sorted.end());
                    for (size_t iseed = 0; iseed < silseed_crossing->size(); ++iseed)
                    {
                        if (silseed_crossing->at(iseed) != static_cast<int>(vertex_crossing)) continue;
                        if (std::binary_search(associated_sorted.begin(), associated_sorted.end(), iseed)) continue;
                        nonassociated_seed_indices.push_back(iseed);
                    }
                }

                // Fixed cluster graphs: all / primary / non-primary
                display.grAllINTTClustersXY = MakeGraph("gr_all_intt_clusters_xy", "All INTT clusters in crossing (x-y)");
                display.grAllINTTClustersZR = MakeGraph("gr_all_intt_clusters_zr", "All INTT clusters in crossing (z-r)");

                std::unordered_set<int> sphenix_primary_track_ids;
                if (!IsData)
                {
                    display.grPrimaryINTTClustersXY    = MakeGraph("gr_primary_intt_clusters_xy",    "INTT clusters from sPHENIX primary particles in crossing (x-y)");
                    display.grPrimaryINTTClustersZR    = MakeGraph("gr_primary_intt_clusters_zr",    "INTT clusters from sPHENIX primary particles in crossing (z-r)");
                    display.grNonPrimaryINTTClustersXY = MakeGraph("gr_nonprimary_intt_clusters_xy", "INTT clusters not from sPHENIX primary particles in crossing (x-y)");
                    display.grNonPrimaryINTTClustersZR = MakeGraph("gr_nonprimary_intt_clusters_zr", "INTT clusters not from sPHENIX primary particles in crossing (z-r)");

                    if (sPHENIXPrimary_trackID)
                        sphenix_primary_track_ids.insert(sPHENIXPrimary_trackID->begin(), sPHENIXPrimary_trackID->end());
                }

                // Vertex position graphs
                display.grRecoVertexXY = MakeGraph("gr_reco_vertex_xy", "Selected reconstructed vertex (x-y)");
                display.grRecoVertexZR = MakeGraph("gr_reco_vertex_zr", "Selected reconstructed vertex (z-r)");
                AddPoint(display.grRecoVertexXY, display.grRecoVertexZR, silvtx_x, silvtx_y, silvtx_z);

                if (!IsData && TruthVertexX && TruthVertexY && TruthVertexZ &&
                    vtx_idx < static_cast<int>(TruthVertexX->size()) &&
                    vtx_idx < static_cast<int>(TruthVertexY->size()) &&
                    vtx_idx < static_cast<int>(TruthVertexZ->size()))
                {
                    display.grTruthVertexXY = MakeGraph("gr_truth_vertex_xy", "Truth vertex (x-y)");
                    display.grTruthVertexZR = MakeGraph("gr_truth_vertex_zr", "Truth vertex (z-r)");
                    AddPoint(display.grTruthVertexXY, display.grTruthVertexZR,
                             TruthVertexX->at(vtx_idx), TruthVertexY->at(vtx_idx), TruthVertexZ->at(vtx_idx));
                }

                // Per-cluster display graphs — re-iterates crossing_cluster_indices
                for (const size_t ihit : crossing_cluster_indices)
                {
                    const double x = ClusX->at(ihit), y = ClusY->at(ihit), z = ClusZ->at(ihit);
                    AddPoint(display.grAllINTTClustersXY, display.grAllINTTClustersZR, x, y, z);

                    if (!IsData && cluster_matchedG4P_trackID && ihit < cluster_matchedG4P_trackID->size())
                    {
                        const bool is_primary = (cluster_matchedG4P_trackID->at(ihit) != 0) &&
                                                sphenix_primary_track_ids.count(cluster_matchedG4P_trackID->at(ihit));
                        if (is_primary)
                        {
                            if (display.grPrimaryINTTClustersXY)    AddPoint(display.grPrimaryINTTClustersXY,    display.grPrimaryINTTClustersZR,    x, y, z);
                        }
                        else
                        {
                            if (display.grNonPrimaryINTTClustersXY) AddPoint(display.grNonPrimaryINTTClustersXY, display.grNonPrimaryINTTClustersZR, x, y, z);
                        }
                    }
                }

                // Seed PCA positions
                if (has_silseed_xyz_branches && silseed_x && silseed_y && silseed_z)
                {
                    display.grAssociatedSilSeedPCAXY = MakeGraph("gr_associated_silseed_pca_xy", "Associated silicon seed PCA points (x-y)");
                    display.grAssociatedSilSeedPCAZR = MakeGraph("gr_associated_silseed_pca_zr", "Associated silicon seed PCA points (z-r)");
                    for (const size_t iseed : used_seed_indices)
                    {
                        if (iseed >= silseed_x->size() || iseed >= silseed_y->size() || iseed >= silseed_z->size()) continue;
                        AddPoint(display.grAssociatedSilSeedPCAXY, display.grAssociatedSilSeedPCAZR,
                                 silseed_x->at(iseed), silseed_y->at(iseed), silseed_z->at(iseed));
                    }
                }

                // Per-seed cluster graphs — associated seeds
                if (has_silseed_cluster_xyz_branches && silseed_cluster_globalX && silseed_cluster_globalY && silseed_cluster_globalZ)
                {
                    std::vector<size_t> seed_indices_sorted(used_seed_indices.begin(), used_seed_indices.end());
                    std::sort(seed_indices_sorted.begin(), seed_indices_sorted.end());

                    for (const size_t iseed : seed_indices_sorted)
                    {
                        if (iseed >= silseed_cluster_globalX->size() || iseed >= silseed_cluster_globalY->size() || iseed >= silseed_cluster_globalZ->size()) continue;

                        TGraph *gr_xy = MakeGraph(Form("gr_associated_silseed_clusters_xy_idx%zu", iseed), Form("Associated silicon seed idx=%zu clusters (x-y)", iseed));
                        TGraph *gr_zr = MakeGraph(Form("gr_associated_silseed_clusters_zr_idx%zu", iseed), Form("Associated silicon seed idx=%zu clusters (z-r)", iseed));
                        FillXYZR(gr_xy, gr_zr, silseed_cluster_globalX->at(iseed), silseed_cluster_globalY->at(iseed), silseed_cluster_globalZ->at(iseed));
                        display.grAssociatedSilSeedClustersXY.push_back(gr_xy);
                        display.grAssociatedSilSeedClustersZR.push_back(gr_zr);
                    }
                }

                // Per-seed cluster graphs — non-associated seeds
                if (has_silseed_cluster_xyz_branches && silseed_cluster_globalX && silseed_cluster_globalY && silseed_cluster_globalZ)
                {
                    for (const size_t iseed : nonassociated_seed_indices)
                    {
                        if (iseed >= silseed_cluster_globalX->size() || iseed >= silseed_cluster_globalY->size() || iseed >= silseed_cluster_globalZ->size()) continue;

                        TGraph *gr_xy = MakeGraph(Form("gr_nonassociated_silseed_clusters_xy_idx%zu", iseed), Form("Non-associated silicon seed idx=%zu clusters (x-y)", iseed));
                        TGraph *gr_zr = MakeGraph(Form("gr_nonassociated_silseed_clusters_zr_idx%zu", iseed), Form("Non-associated silicon seed idx=%zu clusters (z-r)", iseed));
                        FillXYZR(gr_xy, gr_zr, silseed_cluster_globalX->at(iseed), silseed_cluster_globalY->at(iseed), silseed_cluster_globalZ->at(iseed));
                        display.grNonAssociatedSilSeedClustersXY.push_back(gr_xy);
                        display.grNonAssociatedSilSeedClustersZR.push_back(gr_zr);
                    }
                }

                // sPHENIX primary truth clusters
                if (!IsData && sPHENIXPrimary_truthcluster_X && sPHENIXPrimary_truthcluster_Y && sPHENIXPrimary_truthcluster_Z)
                {
                    const size_t nprimary_truth = std::min({sPHENIXPrimary_truthcluster_X->size(), sPHENIXPrimary_truthcluster_Y->size(), sPHENIXPrimary_truthcluster_Z->size()});
                    for (size_t iprim = 0; iprim < nprimary_truth; ++iprim)
                    {
                        if (!sPHENIXPrimary_charge || iprim >= sPHENIXPrimary_charge->size() || sPHENIXPrimary_charge->at(iprim) == 0.0) continue;
                        const auto &truth_cl_x = sPHENIXPrimary_truthcluster_X->at(iprim);
                        const auto &truth_cl_y = sPHENIXPrimary_truthcluster_Y->at(iprim);
                        const auto &truth_cl_z = sPHENIXPrimary_truthcluster_Z->at(iprim);

                        // Guard before graph allocation to avoid leaking TGraphs for empty primaries
                        if (std::min({truth_cl_x.size(), truth_cl_y.size(), truth_cl_z.size()}) == 0) continue;

                        TGraph *gr_xy = MakeGraph(Form("gr_sPHENIXPrimary_truthcluster_xy_idx%zu", iprim), Form("sPHENIX primary idx=%zu truth clusters (x-y)", iprim));
                        TGraph *gr_zr = MakeGraph(Form("gr_sPHENIXPrimary_truthcluster_zr_idx%zu", iprim), Form("sPHENIX primary idx=%zu truth clusters (z-r)", iprim));
                        FillXYZR(gr_xy, gr_zr, truth_cl_x, truth_cl_y, truth_cl_z);
                        display.grPrimaryTruthClustersXY.push_back(gr_xy);
                        display.grPrimaryTruthClustersZR.push_back(gr_zr);
                    }
                }

                // sPHENIX primary reconstructed clusters
                if (!IsData && sPHENIXPrimary_recocluster_globalX && sPHENIXPrimary_recocluster_globalY && sPHENIXPrimary_recocluster_globalZ)
                {
                    const size_t nprimary_reco = std::min({sPHENIXPrimary_recocluster_globalX->size(), sPHENIXPrimary_recocluster_globalY->size(), sPHENIXPrimary_recocluster_globalZ->size()});
                    for (size_t iprim = 0; iprim < nprimary_reco; ++iprim)
                    {
                        if (!sPHENIXPrimary_charge || iprim >= sPHENIXPrimary_charge->size() || sPHENIXPrimary_charge->at(iprim) == 0.0) continue;
                        const auto &reco_cl_x = sPHENIXPrimary_recocluster_globalX->at(iprim);
                        const auto &reco_cl_y = sPHENIXPrimary_recocluster_globalY->at(iprim);
                        const auto &reco_cl_z = sPHENIXPrimary_recocluster_globalZ->at(iprim);

                        TGraph *gr_xy = MakeGraph(Form("gr_sPHENIXPrimary_recocluster_xy_idx%zu", iprim), Form("sPHENIX primary idx=%zu reconstructed clusters (x-y)", iprim));
                        TGraph *gr_zr = MakeGraph(Form("gr_sPHENIXPrimary_recocluster_zr_idx%zu", iprim), Form("sPHENIX primary idx=%zu reconstructed clusters (z-r)", iprim));

                        const int points_added = FillXYZR(gr_xy, gr_zr, reco_cl_x, reco_cl_y, reco_cl_z);
                        if (points_added == 0) { delete gr_xy; delete gr_zr; continue; }

                        display.grPrimaryRecoClustersXY.push_back(gr_xy);
                        display.grPrimaryRecoClustersZR.push_back(gr_zr);
                    }
                }

                out.crossingDisplays.push_back(display);
                std::cout << "Prepared crossing display graphs for crossing key " << display.crossing_key << " (associated seeds = " << n_associated_silseeds << ", INTT clusters = " << n_intt_clusters_in_crossing << ")" << std::endl;
            };

            if (make_crossing_display)
                BuildCrossingDisplay();

            // Minimum associated silicon seed cut
            if (cfg.min_vertex_associated_silseeds >= 3 && n_associated_silseeds < cfg.min_vertex_associated_silseeds)
            {
                std::cout << "Skipping vertex crossing " << vertex_crossing << " because it has only " << used_seed_indices.size() << " associated silicon seeds, below the requested minimum of " << cfg.min_vertex_associated_silseeds << std::endl;
                continue;
            }

            // Fill running tracklet histograms for this crossing
            Tracklets_CombinatoricDR_FillEvent(tkldata, cfg, silicon_vtx, out, output_prefix);
            if (multiplicity_percentile_bin_index >= 0)
            {
                Tracklets_CombinatoricDR_FillEvent(tkldata, cfg, silicon_vtx,
                    multiplicity_percentile_outputs[multiplicity_percentile_bin_index],
                    MultiplicityPercentileBinPrefix(output_prefix, multiplicity_percentile_bins[multiplicity_percentile_bin_index].percentile_bin));
            }
            ResetVec(tkldata);

        } // end crossing loop
        vertexCrossing_idx_map.clear();
        ResetVec(tkldata);

    } // end event loop

    // -----------------------------------------------------------------------
    // Finalize and write output
    // -----------------------------------------------------------------------
    Tracklets_CombinatoricDR_Finalize(cfg, out, output_prefix);
    for (size_t i = 0; i < multiplicity_percentile_outputs.size(); ++i)
        Tracklets_CombinatoricDR_Finalize(cfg, multiplicity_percentile_outputs[i], MultiplicityPercentileBinPrefix(output_prefix, multiplicity_percentile_bins[i].percentile_bin));

    const bool write_ok = WriteTracklets_CombinatoricDR_Output(out, multiplicity_percentile_outputs, outfilename.Data());
    if (!write_ok)
    {
        std::cerr << "Failed to create output ROOT file: " << outfilename << std::endl;
        f->Close(); delete f; return 1;
    }

    std::cout << "Wrote combinatoric histograms to " << outfilename << std::endl;
    ResetVec(tkldata);
    f->Close();
    delete f;

    return 0;
}