#ifndef __MC_CALO_PYTHIA_PIO_H__
#define __MC_CALO_PYTHIA_PIO_H__

#include <fun4all/SubsysReco.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertexMap.h>
#include <calobase/RawClusterContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>

#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include <Math/Vector4D.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>

// Jet info needed for MC sample stitching.
#include <jetbase/Jetv1.h>
#include <jetbase/Jetv2.h>
#include <jetbase/JetContainer.h>

class MC_calo_pythia_pi0 : public SubsysReco
{
 public:
  //! Constructor
  MC_calo_pythia_pi0(const std::string &name = "MC_calo_pythia_pi0",
                     const std::string &outputfilename = "analysis_mc.root",
                     const int seednb = 0);

  //! Destructor
  virtual ~MC_calo_pythia_pi0();

  // Full Initialization
  int Init(PHCompositeNode *);

  int InitRun(PHCompositeNode *);
  
  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  //! Determine whether the simulation is MB or Photon:
  void set_pythia_photon(bool val)
  {
    if (val)
    {
      m_mb_pythia = false;
      m_photon_pythia = true;
    } else {
      m_mb_pythia = true;
      m_photon_pythia = false;
    }
  }

  void smear_photon(float& eta, float& phi, float& ecore, int smear_idx);

  void smear_ultra(float& ecore);

  void set_smearing(bool val, bool val_ecore = true, bool val_eta = true, bool val_phi = true)
  {
    do_smearing = val;
    do_smear_ecore = val_ecore;
    do_smear_eta = val_eta;
    do_smear_phi = val_phi;
  }

  void set_ultra_smearing(bool val)
  {
    do_ultra_smearing = val; 
  }

  //! Check diphoton cut
  bool diphoton_cut(ROOT::Math::PtEtaPhiMVector p1,
                    ROOT::Math::PtEtaPhiMVector p2,
                    ROOT::Math::PtEtaPhiMVector ppair);

  //! Gives the angle value between -pi and + pi
  void WrapAngle(float& phi);

  //! Absolute angle difference with wrapping
  float WrapAngleDifference(const float& phi1, const float& phi2);

  // Struct for pi0/eta decay photons
  struct TruthParticleInfo
  {
    int pid = 0;
    int track_id = -1;
    int parent_id = -1;
    ROOT::Math::PxPyPzEVector p4;
  };

  enum TruthMatchType
  {
    None = 0,
    DecayPhoton,
    NonDecayPhoton,
    NonPhoton
  };

  enum MassType
  {
    Pi0 = 0,
    Eta,
    Combinatorial,
    Unmatched,
    NMassTypes
  };

  // reco cluster + truth matching index
  struct customCluster {
    int match_truth_track_id = -1;
    int match_truth_pid = 0; // Can be non-photon
    TruthMatchType match_truth_type = TruthMatchType::None;
    int match_parent_track_id = -1;
    int match_parent_pid = 0;
    double match_truth_dR = -1;
    double match_truth_Eratio = -1;
    ROOT::Math::PtEtaPhiMVector p4;
  };

  void set_vertex_reweighting(bool val, float sigma)
  {
    do_vertex_reweighting = val;
    m_sigma_vertex = sigma;
  }

  void set_photon_trigger_efficiency_threshold(float efficiency_target)
  {
    efficiency_index = FindClosestValFromVector(efficiency_target, efficiency_thresholds, nThresholds);
    energy_threshold_3 = array_energy_threshold_3[efficiency_index];
    energy_threshold_4 = array_energy_threshold_4[efficiency_index];
    std::cout << "efficiency_index = " << efficiency_index << std::endl;
    std::cout << "efficiency threshold is " << efficiency_thresholds[efficiency_index] << std::endl;
    std::cout << "min energies = (" << energy_threshold_3 << ", " << energy_threshold_4 << ")" << std::endl;
  }

  // helper method for trigger efficiency threshold
  int FindClosestValFromVector(float val, const float* binValues, int nBins)
  {
    const float *it = std::upper_bound(binValues, binValues + nBins, val);
    int icenter = static_cast<int>(it - binValues) - 1;
    int iminus = std::max(0, icenter - 1);
    int imaxi = std::min(nBins - 1, icenter + 1);
    float valminus = *(binValues + iminus);
    float valcenter = *(binValues + icenter);
    float valmaxi = *(binValues + imaxi);
    float diffminus = std::abs(valminus - val);
    float diffcenter = std::abs(valcenter - val);
    float diffmaxi = std::abs(valmaxi - val);
    std::vector<float> diffs = {diffminus, diffcenter, diffmaxi};
    auto pos = std::min_element(diffs.begin(), diffs.end());
    int ipos = int(pos - diffs.begin());
    if (ipos == 0) return iminus;
    else if (ipos == 1) return icenter;
    else if (ipos == 2) return imaxi;
    return 0;
  }

  // Proxy Trigger Matching
  bool trigger_efficiency_matching(const ROOT::Math::PtEtaPhiMVector& photon1,
                                   const ROOT::Math::PtEtaPhiMVector& photon2,
                                   const ROOT::Math::PtEtaPhiMVector& diphoton);

  float get_max_energy(std::vector<customCluster>& reco_clusters);

  void match_truth_index(std::vector<customCluster>& reco_clusters);

  void book_histograms();

  void fill_truth_histograms();
  
  void fill_histograms(const customCluster& photon1,
                       const customCluster& photon2,
                       const ROOT::Math::PtEtaPhiMVector& diphoton);

  int FindBinBinary(float val, const float* binEdges, int nBins)
  {
    const float *it = std::upper_bound(binEdges, binEdges + nBins, val);
    int ibin = static_cast<int>(it - binEdges) -1;
    return ibin;
  }

  int FindBinDirect(float val, float valMin, float valMax, int nBins)
  {
    if (val < valMin) { // Underflow
      return -1;
    }
    else if (val > valMax) { // Overflow
      return nBins;
    }
    int ibin = static_cast<int>((val - valMin) * nBins / (valMax - valMin));
    if (ibin == nBins) ibin = nBins - 1;
    else if (ibin < 0) ibin = 0;
    return ibin;
  }

  int global_index(int iKin, int iPhi)
  {
    return iKin * nPhiBins + iPhi;
  }

  // To match the photon-triggered analysis, mimic the proxy trigger matching
  void set_apply_photon_reco(bool val)
  {
    require_efficiency_matching = val;
  }

  //! Set diphoton pT cut in pi0 mass range
  void set_ptcut(float pTMin = 1.0, float pTMax = 1000.0)
  {
    pTCutMin = pTMin;
    pTCutMax = pTMax;
  }

  void set_chi2_cut(const float chi2)
  {
    clus_chisq_cut = chi2;
  }

  void set_E_cut(const float E)
  {
    clus_E_cut = E;
  }
  
  enum MCSample
  {
    MB = 0,
    Jet5,
    Jet8,
    Jet12,
    Jet20,
    NMCSamples
  };

  void set_stitch_index(int ind) { _mc_index = ind; }
  
  bool stitch_mc_sample_truth_pt();

  bool stitch_mc_sample_diphoton(float diphoton_pt);

  void set_smear_index(int si) { _smear_index = si; }

  void set_scale_variation(float scale_diff=1.5)
  {
    m_scale_diff = scale_diff;
  }

 private:

  // event counter
  int _eventcounter;

  // Choice of pythia set (default is Minimum Bias)
  bool m_mb_pythia = true;
  bool m_photon_pythia = false;

  // Stitching values in pT for different MC samples
  static constexpr int nMCSamples = MCSample::NMCSamples;
  float max_truth_jet_pt = 0;
  const float _truth_pt_min[nMCSamples] = {0, 7, 9, 14, 21};
  const float _truth_pt_max[nMCSamples] = {7, 9, 14, 21, 32};
  const float _reco_pt_max[nMCSamples] = {14, 24, 1000, 1000, 1000};
  const float _diphoton_pt_min[nMCSamples] = {0, 3.5, 4.5, 7, 10.5};
  const float _diphoton_pt_max[nMCSamples] = {3.5, 4.5, 7, 10.5, 16};
  int _mc_index = -1;
  // In principle the upper limit for Jet 12 GeV sample is 21 GeV,
  // but we don't consider any MC sample above it in this study. 

  // Input nodes
  PHG4TruthInfoContainer *truthinfo = nullptr;
  JetContainer *_truth_jets = nullptr;
  JetContainer *_reco_jets = nullptr;
  GlobalVertexMap *vertexmap = nullptr;
  RawClusterContainer *clusterContainer = nullptr;

  // Truth info
  const int m_pi0_pid = 111; // pi0 in MC numbering scheme
  const float Vs = 200.0; // GeV
  // std::vector<TruthParticleInfo> decaying_mesons;
  // std::vector<TruthParticleInfo> decay_photons;
  std::unordered_map<int, TruthParticleInfo> photon_by_track_id;
  std::unordered_map<int, TruthParticleInfo> meson_by_track_id;
  std::unordered_map<int, std::vector<TruthParticleInfo>> photons_by_parent;
  std::unordered_map<int, TruthParticleInfo> nondecay_photon_by_track_id;
  std::unordered_map<int, TruthParticleInfo> nonphoton_by_track_id;

  // Vertices
  float vertex_z_true = 1000;
  float vertex_z_reco = 1000;

  // Optional vertex reweighting
  float m_event_weight = 1;
  bool do_vertex_reweighting = false;
  float m_sigma_vertex = 65; // cm
  const float m_sigma_reference = 65; // cm

  // Analysis cut
  float alphaCut = 0.7;
  float pTCutMin = 1.0;
  float pTCutMax = 1000.0;
  float clus_E_cut = 1.0; // GeV
  float clus_chisq_cut = 1000.0; // Essentially no chi2 cut

  // Photon info container (per event)
  std::vector<customCluster> good_photons;
  int num_photons = 0;

  // trigger efficiency matching (Cluster energy above 70% efficiency threshold of trigger)
  const std::string trigger_efficiency_path = "/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/turnon_efficiency/FitFunctions_Photon_plus_MBD_NS_geq_1.root";
  TFile *trigger_efficiency_file = nullptr;
  TF1 *trigger_turnon_curve = nullptr;
  
  bool require_efficiency_matching = false;
  bool efficiency_match = false;
  int efficiency_index = 6; // 70% threshold
  float energy_threshold_3 = 3.0;
  float energy_threshold_4 = 4.0;
  // 95 -> 4.2/5.3
  // 70 -> 3.5/4.3
  static constexpr int nThresholds = 10;
  const float efficiency_thresholds[nThresholds] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 95}; // efficiency thresholds
  const float array_energy_threshold_3[nThresholds] = {2.4, 2.7, 2.9, 3.1, 3.2, 3.4, 3.5, 3.7, 4.0, 4.3}; // energy thresholds for the photon 3 GeV trigger
  const float array_energy_threshold_4[nThresholds] = {2.9, 3.3, 3.6, 3.7, 3.9, 4.1, 4.3, 4.5, 4.9, 5.3}; // energy thresholds for the photon 4 GeV trigger
  static constexpr float delta_eta_threshold = 0.192; // pseudo-rapidity distance between 8 towers
  static constexpr float delta_phi_threshold = 0.192; // azimuthal distance between 8 towers

  // Output histogram file
  std::string outfilename;
  TFile *outfile = nullptr;
  
  // Seed Number (for smearing)
  int segnumber;
  int seednumber;

  // Smearing info for EMCal Resolution
  bool do_smearing = false;
  bool do_smear_ecore = false;
  bool do_smear_eta = false;
  bool do_smear_phi = false;
  TRandom3 *rnd = nullptr;

  int _smear_index = 2;
  bool do_ultra_smearing = false;
  TF1 *f_energy_ultra = nullptr;
  
  //! Quadrature diff graphs from function_compare_wide.root; x = E (GeV), Eval = width for smearing)
  static constexpr int nDifferences = 3;
  std::array<TF1*, nDifferences> f_energy_smear{};
  std::array<TF1*, nDifferences> f_position_smear{};

  // QA (Check that kinematics are actually smeared)
  TH1F *h_smear_eta = nullptr;
  TH1F *h_smear_phi = nullptr;
  TH1F *h_smear_E = nullptr;
  TH2F *h_smear_E_deta = nullptr;
  TH2F *h_smear_E_dphi = nullptr;
  TH2F *h_smear_E_dE = nullptr;

  // Energy scale variation
  float m_scale_diff = 0;

  // pT bins, same as those used in PHENIX 2021 Asymmetries (+ low pT)
  // static constexpr int nPtBins = 9;
  // const float pTBins[nPtBins + 1] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 20};
  static constexpr int nPtBins = 8;
  const float pTBins[nPtBins + 1] = {1, 2, 3, 4, 5, 6, 8, 10, 20};

  // For energy resolution studies
  static constexpr int nEBins = 7;
  const float EBins[nEBins + 1] = {1, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 5.0};

  // New binning -> more equally distributed
  static constexpr int nEtaBins = 8;
  const float etaBins[nEtaBins + 1] = {-2.00, -1.05, -0.86, -0.61, 0.0, 0.61, 0.86, 1.05, 2.0};

  // New binning -> more equally distributed
  static constexpr int nXfBins = 8;
  const float xfBins[nXfBins + 1] = {-0.200, -0.048, -0.035, -0.022, 0.0, 0.022, 0.035, 0.048, 0.200};
  
  static constexpr int nZvtxBins = 7;
  const float zvtxBins[nZvtxBins + 1] = {0, 10, 30, 50, 70, 100, 150, 200};

  static constexpr int nPhiBins = 12;
  const float phiBins[nPhiBins + 1] = {
    -M_PI, -5.0 * M_PI / 6.0, -2.0 * M_PI / 3.0, -M_PI / 2.0, -M_PI / 3.0, -M_PI / 6.0,
    0, M_PI / 6.0, M_PI / 3.0, M_PI / 2.0, 2.0 * M_PI / 3.0, 5.0 * M_PI / 6.0, M_PI
  };

  //2D -> 1D mapping for two-dimensional response matrices

  // Histograms

  // Vertices
  TH1F *h_true_zvtx = nullptr;
  TH1F *h_reco_zvtx = nullptr;
  TH1F *h_reco_true_zvtx = nullptr;

  // Invariant mass decomposition
  static constexpr int nMassTypes = MassType::NMassTypes;
  const std::string type_suffix[nMassTypes] = {"pi0", "eta", "comb", "unmatched"};
  TH1F *h_pair_mass_total = nullptr;
  TH1F *h_pair_mass[nMassTypes] = {nullptr};
  TH1F *h_pair_mass_zvtx[nMassTypes][nZvtxBins] = {nullptr};
  TH1F *h_pair_mass_eta[nMassTypes][nEtaBins] = {nullptr};
  TH1F *h_pair_mass_E[nMassTypes][nEBins] = {nullptr};
  TH1F *h_pair_mass_pt[nMassTypes][nPtBins] = {nullptr};
  TH1F *h_pair_mass_xf[nMassTypes][nXfBins] = {nullptr};

  // 1D true distributions
  TH1F *h_reco_E[nMassTypes] = {nullptr};
  TH1F *h_reco_pt[nMassTypes] = {nullptr};
  TH1F *h_reco_eta[nMassTypes] = {nullptr};
  TH1F *h_reco_xf[nMassTypes] = {nullptr};
  TH1F *h_reco_phi[nMassTypes] = {nullptr};

  // 1D reco distributions
  TH1F *h_true_E[nMassTypes - 1] = {nullptr};
  TH1F *h_true_pt[nMassTypes - 1] = {nullptr};
  TH1F *h_true_eta[nMassTypes - 1] = {nullptr};
  TH1F *h_true_xf[nMassTypes - 1] = {nullptr};
  TH1F *h_true_phi[nMassTypes - 1] = {nullptr};

  // Distributions with truth-matched cluster -> nMassTypes - 1

  // Differences Reco/True
  TH1F *h_reco_true_E[nMassTypes - 1] = {nullptr};
  TH1F *h_reco_true_pt[nMassTypes - 1] = {nullptr};
  TH1F *h_reco_true_eta[nMassTypes - 1] = {nullptr};
  TH1F *h_reco_true_xf[nMassTypes - 1] = {nullptr};
  TH1F *h_reco_true_phi[nMassTypes - 1] = {nullptr};

  // Matrix projection
  TH1F *h_true_coarse_pt[nMassTypes-1] = {nullptr};
  TH1F *h_true_coarse_eta[nMassTypes-1] = {nullptr};
  TH1F *h_true_coarse_xf[nMassTypes-1] = {nullptr};
  TH1F *h_true_coarse_phi[nMassTypes-1] = {nullptr};
  TH1F *h_true_coarse_pt_phi[nMassTypes-1] = {nullptr};
  TH1F *h_true_coarse_eta_phi[nMassTypes-1] = {nullptr};
  TH1F *h_true_coarse_xf_phi[nMassTypes-1] = {nullptr};
  TH1F *h_reco_coarse_pt[nMassTypes-1] = {nullptr};
  TH1F *h_reco_coarse_eta[nMassTypes-1] = {nullptr};
  TH1F *h_reco_coarse_xf[nMassTypes-1] = {nullptr};
  TH1F *h_reco_coarse_phi[nMassTypes-1] = {nullptr};
  TH1F *h_reco_coarse_pt_phi[nMassTypes-1] = {nullptr};
  TH1F *h_reco_coarse_eta_phi[nMassTypes-1] = {nullptr};
  TH1F *h_reco_coarse_xf_phi[nMassTypes-1] = {nullptr};

  // Response matrices for unfolding
  TH2F *h_response_pt[nMassTypes-1] = {nullptr};
  TH2F *h_response_eta[nMassTypes-1] = {nullptr};
  TH2F *h_response_xf[nMassTypes-1] = {nullptr};
  TH2F *h_response_phi[nMassTypes-1] = {nullptr};
  TH2F *h_response_pt_phi[nMassTypes-1] = {nullptr};
  TH2F *h_response_eta_phi[nMassTypes-1] = {nullptr};
  TH2F *h_response_xf_phi[nMassTypes-1] = {nullptr};
  TH2F *h_response_fine_pt[nMassTypes-1] = {nullptr};
  TH2F *h_response_fine_eta[nMassTypes-1] = {nullptr};
  TH2F *h_response_fine_xf[nMassTypes-1] = {nullptr};
  TH2F *h_response_fine_phi[nMassTypes-1] = {nullptr};

  // Matrices for meson reconstruction efficiency
  TH1F* h_true_pt_n[nMassTypes-1] = {nullptr};
  TH1F* h_true_eta_n[nMassTypes-1] = {nullptr};
  TH1F* h_true_xf_n[nMassTypes-1] = {nullptr};
  TH1F* h_true_pt_phi_n[nMassTypes-1] = {nullptr};
  TH1F* h_true_eta_phi_n[nMassTypes-1] = {nullptr};
  TH1F* h_true_xf_phi_n[nMassTypes-1] = {nullptr};
  TH1F* h_reco_pt_n_matched[nMassTypes-1] = {nullptr};
  TH1F* h_reco_eta_n_matched[nMassTypes-1] = {nullptr};
  TH1F* h_reco_xf_n_matched[nMassTypes-1] = {nullptr};
  TH1F* h_reco_pt_phi_n_matched[nMassTypes-1] = {nullptr};
  TH1F* h_reco_eta_phi_n_matched[nMassTypes-1] = {nullptr};
  TH1F* h_reco_xf_phi_n_matched[nMassTypes-1] = {nullptr};
};

#endif
