#ifndef __AN_NEUTRAL_MESON_MICRO_DST_H__
#define __AN_NEUTRAL_MESON_MICRO_DST_H__

#include <fun4all/SubsysReco.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <queue>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <globalvertex/GlobalVertexMap.h>
#include "ClusterSmallInfoContainer.h"

#include <calobase/RawTowerGeomContainer.h>

#include <ffaobjects/SyncDefs.h>
#include <ffaobjects/SyncObject.h>

// Trigger emulator
#include <triggeremulator/TriggerTile.h>

// MBD info
#include <pileupfilter/MBDSmallInfo.h>

// Forward declarations
class PHCompositeNode;
class TFile;
class TTree;
class TRandom3;
class TH1;
class TH1F;
class TH1I;
class TH2F;
class TGraph;
class TF1;
class LorentzVector;

void monitorMemoryUsage(const std::string& label="");

class AnNeutralMeson_micro_dst : public SubsysReco
{
 public:
  //! constructor
  AnNeutralMeson_micro_dst(const std::string &name = "AnNeutralMeson_micro_dst",
                           const int runnb = 48746,
                           const std::string &outputfilename = "analysis_per_run/analysis_48746.root",
                           const std::string &outputfiletreename = "analysis_per_run/diphoton_minimal_48746.root",
                           const int seednb = 0);

  //! destructor
  virtual ~AnNeutralMeson_micro_dst();

  //! full initialization
  int Init(PHCompositeNode *);

  int InitRun(PHCompositeNode *);
  
  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  int FindBinBinary(float val, const float* binEdges, int nBins)
  {
    const float *it = std::upper_bound(binEdges, binEdges + nBins, val);
    int ibin = static_cast<int>(it - binEdges) -1;
    return ibin;
  }

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

  // Determine whether MBD Trigger is fired for this event
  bool mbd_trigger_bit();

  // Determine whether Photon Trigger is fired for this event
  bool photon_trigger_bit();

  // Read input tile info from the trigger emulator
  void trigger_emulator_input();

  // Check if photon pair satsisfies the trigger efficiency matching using the proxy method
  bool trigger_efficiency_matching(const ROOT::Math::PtEtaPhiMVector& photon1,
                                   const ROOT::Math::PtEtaPhiMVector& photon2,
                                   const ROOT::Math::PtEtaPhiMVector& diphoton);

  void smear_photon(float& eta, float& phi, float& ecore);

  void event_mixing_mbd(PHCompositeNode *topNode);

  void event_mixing_photon();

  void set_sigma_number(const float val) {
    sigma_number = val;
    const int N = nParticles * (nRegions + 1) * 2;
    switch(sigma_number) {
    case 1:
      for (int j = 0; j < N; j++) {
        band_limits[j] = band_limits_15[j];
      }
      break;
    case 2:
      for (int j = 0; j < N; j++) {
        band_limits[j] = band_limits_2[j];
      }
      break;
    case 3:
      for (int j = 0; j < N; j++) {
        band_limits[j] = band_limits_3[j];
      }
      break;
    default:
      std::cout << "invalid sigma_number, take 3 sigma as default" << std::endl;
      break;
    }
  }

  void set_smearing(bool val, bool val_ecore = true, bool val_eta = true, bool val_phi = true)
  {
    do_smearing = val;
    do_smear_ecore = val_ecore;
    do_smear_eta = val_eta;
    do_smear_phi = val_phi;
  }

  void set_scale_variation(int scale_variation = -1, float scale_diff=2.6)
  {
    m_scale_diff = scale_diff;
    if (scale_variation == 1) {
      do_small_scale = true;
      do_high_scale = false;
      std::cout << "Reduce energy scale by " << m_scale_diff <<  "%" << std::endl;
    }
    else if (scale_variation == 2) {
      do_small_scale = false;
      do_high_scale = true;
      std::cout << "Increase energy scale by " << m_scale_diff << "%" << std::endl;
    }
    else {
      do_small_scale = false;
      do_high_scale = false;
      std::cout << "Apply no energy scale variation" << std::endl;
    }
  }

  void set_event_mixing(bool use) { use_event_mixing = use; }

  void set_mbd_trigger_bit_requirement(bool require, int bit_index = -1)
  {
    require_mbd_trigger_bit = require;

    if (bit_index == 1) {
      require_mbd_any_vtx = require;
    }
    else if (bit_index == 2) {
      require_mbd_vtx_10 = require;
    }
  }

  void set_photon_trigger_bit_requirement(bool require, int bit_index = -1)
  {
    require_photon_trigger_bit = require;

    if (bit_index == 1) {
      require_photon_3_any_vtx = require;
    }
    else if (bit_index == 2) {
      require_photon_3_vtx_10 = require;
    }
    else if (bit_index == 3) {
      require_photon_4_any_vtx = require;
    }
    else if (bit_index == 4) {
      require_photon_4_vtx_10 = require;
    }
  }

  void set_photon_trigger_emulator_matching_requirement(bool require) { require_emulator_matching = require; }

  void set_photon_trigger_efficiency_matching_requirement(bool require) { require_efficiency_matching = require; }

  void set_photon_trigger_efficiency_threshold(float efficiency_target) {
    efficiency_index = FindClosestValFromVector(efficiency_target, efficiency_thresholds, nThresholds);
    energy_threshold_3 = array_energy_threshold_3[efficiency_index];
    energy_threshold_4 = array_energy_threshold_4[efficiency_index];
    std::cout << "efficiency_index = " << efficiency_index << std::endl;
    std::cout << "efficiency threshold is " << efficiency_thresholds[efficiency_index] << std::endl;
    std::cout << "min energies = (" << energy_threshold_3 << ", " << energy_threshold_4 << ")" << std::endl;
  }

  void set_vertex_max(float vtx) { vertex_max = vtx; }

  //! Apply a cluster mask in the (eta, phi) phase space:
  void set_custom_mask(bool mask, float z_min_val=-1, float z_max_val=-1,
                       float phi_min_val=-1, float phi_max_val=-1)
  {
    do_custom_mask = mask;
    _z_min = z_min_val;
    _z_max = z_max_val;
    _phi_min = phi_min_val;
    _phi_max = phi_max_val;
  }

  //! Set cluster level chi2 cut. Only first cut is applied (the rest is for QA)
  void set_chi2cut(const std::vector<float>& chi2) { chi2_cuts = chi2; n_chi2_cuts = chi2.size(); } 

  //! Set cluster level ecore cut. Only first cut is applied (the rest is for QA)
  void set_ecorecut(const std::vector<float>& ecore) { ecore_cuts = ecore; n_ecore_cuts = ecore.size(); } 

  //! Set diphoton alpha cut
  void set_alphacut(float alpha) { alphaCut = alpha; }

  //! Set diphoton pT cut in pi0 mass range
  void set_ptcut(float pTMin = 1.0, float pTMax = 1000.0)
  {
    pTCutMin = pTMin;
    pTCutMax = pTMax;
  }

  //! Set diphoton pT threshold between MBD and photon trigger selection
  void set_pt_threshold(float pT = 0)
  {
    pTCutThreshold = pT;
  }

  //! Choose to record QA histograms or not
  void set_store_qa(bool val) { store_qa = val; }

  //! Produce QA histograms at different chi2 / energy cuts
  void cluster_cuts();

  //! Check diphoton cut
  bool diphoton_cut(ROOT::Math::PtEtaPhiMVector p1,
                    ROOT::Math::PtEtaPhiMVector p2,
                    ROOT::Math::PtEtaPhiMVector ppair);
  
  //! Check trigger matching
  bool trigger_matching(const ROOT::Math::PtEtaPhiMVector&, const ROOT::Math::PtEtaPhiMVector&, const ROOT::Math::PtEtaPhiMVector&);

  //! MBD RMS time cut for in-time pileup
  void set_mbd_timing_cut(bool decision, float rms_min, float rms_max)
  {
    require_mbd_timing = decision;
    mbd_rms_cut_min = rms_min;
    mbd_rms_cut_max = rms_max;
  }

  bool mbd_timing_cut();

  bool startswith(const std::string&, const std::string&);

  //! Gives the angle value between -pi and + pi
  void WrapAngle(float& phi);

  //! Absolute angle difference with wrapping
  float WrapAngleDifference(const float& phi1, const float& phi2);

  ROOT::Math::XYZVector EnergyWeightedAverageP3(const std::vector<ROOT::Math::PtEtaPhiMVector>& v4s);

  void set_store_tree(bool does_store_tree) { store_tree = does_store_tree; }

  struct Cluster {
    ROOT::Math::PtEtaPhiMVector p4;
    bool isTrigger = false;
  };


  // Minimal information for pooling event in event mixing
  struct EventInPool {
    float zvtx = 1000;
    float eta_lead = -10;
    float phi_lead = -10;
    std::vector<Cluster> clusters;
  };
 
 private:
  
  // run number
  int runnumber;

  // event counter
  int _eventcounter;

  // TTree file and object
  std::string infilename;
  TFile *infile;
  std::string treename;
  TTree *microDST;

  // Global information in clusters' node
  unsigned long live_trigger;
  unsigned long scaled_trigger;
  float vertex_z;
  int cluster_number;

  // Smearing info for EMCal Resolution
  bool do_smearing = false;
  bool do_smear_ecore = false;
  bool do_smear_eta = false;
  bool do_smear_phi = false;
  TRandom3 *rnd = nullptr;
  //! Quadrature diff graphs from function_compare_wide.root; x = E (GeV), Eval = width for smearing)
  static constexpr int nDifferences = 3;
  std::array<TF1*, nDifferences> f_energy_smear{};
  std::array<TF1*, nDifferences> f_position_smear{};

  // Energy scale variation
  float m_scale_diff = 2.6; // %
  bool do_high_scale = false;
  bool do_small_scale = false;

  // QA (Check that kinematics are actually smeared)
  TH1F *h_smear_eta = nullptr;
  TH1F *h_smear_phi = nullptr;
  TH1F *h_smear_E = nullptr;
  TH2F *h_smear_E_deta = nullptr;
  TH2F *h_smear_E_dphi = nullptr;
  TH2F *h_smear_E_dE = nullptr;

  // Output histogram file
  std::string outfilename;
  std::string outtreename;
  TFile *outfile = nullptr;
  TFile *outfile_tree = nullptr;
  bool store_tree = false;

  // Seed Number (for smearing)
  int seednumber;

  // Simplified output tree (for nano analysis)
  int diphoton_bunchnumber;
  float diphoton_vertex_z;
  float diphoton_mass;
  float diphoton_eta;
  float diphoton_pt;
  float diphoton_xf;
  float diphoton_phi;
  TTree *output_tree;

  // List of configurations
  static constexpr int nBeams = 2; // Yellow or blue beam
  const std::string beams[nBeams] = {"yellow", "blue"};
  static constexpr int nParticles = 2; // pi0 or eta
  const std::string particle[nParticles] = {"pi0", "eta"};
  static constexpr int nRegions = 2; // peak band or side_band invariant mass region
  const std::string regions[nRegions] = {"peak", "side"};
  static constexpr int nSpins = 2; // up or down spin
  const std::string spins[nSpins] = {"up", "down"};

  // // pT bins, same as those used in PHENIX 2021 Asymmetries
  // static constexpr int nPtBins = 9;
  // const float pTBins[nPtBins + 1] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 20};

  // pT bins, re-binned to limit migration effects
  static constexpr int nPtBins = 8;
  const float pTBins[nPtBins + 1] = {1, 2, 3, 4, 5, 6, 8, 10, 20};

  // // pT bins, re-binned to limit migration effects (more stringent)
  // static constexpr int nPtBins = 7;
  // const float pTBins[nPtBins + 1] = {1, 2, 3, 4, 5, 7, 10, 20};

  // New binning -> more equally distributed
  static constexpr int nEtaBins = 8;
  const float etaBins[nEtaBins + 1] = {-2.00, -1.05, -0.86, -0.61, 0.0, 0.61, 0.86, 1.05, 2.0};

  // New binning -> more equally distributed
  static constexpr int nXfBins = 8;
  const float xfBins[nXfBins + 1] = {-0.200, -0.048, -0.035, -0.022, 0.0, 0.022, 0.035, 0.048, 0.200};
  
  static constexpr int nZvtxBins = 7;
  const float zvtxBins[nZvtxBins + 1] = {0, 10, 30, 50, 70, 100, 150, 200};

  static constexpr int nEtaRegions = 2;
  // |eta| < thresh or |eta| > thresh -> low or high
  // thresh defines the PHENIX acceptance
  const std::string etaRegions[nEtaRegions] = {"low", "high"};
  static constexpr float Vs = 200; // 200 GeV
  static constexpr int nDirections = 2;
  // (eta * beamDirection) > 0 or (eta * beamDirection < 0) -> forward or backward
  const std::string directions[nDirections] = {"forward", "backward"};
  
  // Constants
  const double etaThreshold = 0.35; // PHENIX acceptance threshold
  const float beamDirection[nBeams] = {-1, 1}; // Yes, in that order for yellow and blue

  // List spin info
  static constexpr int nBunches = 120;
  int crossingshift;
  int beamspinpat[nBeams][nBunches];

  // Define the regions (in invariant mass) for pi0/eta peak/side (at 3 sigma)
  int sigma_number = 3;
  float band_limits_3[nParticles * (nRegions + 1) * 2] =
    {0.030, 0.070, // pi0 left side invariant mass range (in GeV/c^2)
     0.080, 0.199, // pi0 peak
     0.209, 0.249, // pi0 right side
     0.257, 0.371, // eta left side
     0.399,0.739, // eta peak
     0.767, 0.880}; // eta right side
  // 2 sigma window
  float band_limits_2[nParticles * (nRegions + 1) * 2] =
    {0.030, 0.070, // pi0 left side invariant mass range (in GeV/c^2)
     0.100, 0.180, // pi0 peak
     0.209, 0.249, // pi0 right side
     0.257, 0.371, // eta left side
     0.456,0.683, // eta peak
     0.767, 0.880}; // eta right side
  // 1.5 sigma window
  float band_limits_15[nParticles * (nRegions + 1) * 2] =
    {0.030, 0.070, // pi0 left side invariant mass range (in GeV/c^2)
     0.110, 0.170, // pi0 peak
     0.209, 0.249, // pi0 right side
     0.257, 0.371, // eta left side
     0.484,0.654, // eta peak
     0.767, 0.880}; // eta right side
  // Default window (defined to 3 sigma)
  float band_limits[nParticles * (nRegions + 1) * 2] =
    {0.030, 0.070, // pi0 left side invariant mass range (in GeV/c^2)
     0.080, 0.199, // pi0 peak
     0.209, 0.249, // pi0 right side
     0.257, 0.371, // eta left side
     0.399,0.739, // eta peak
     0.767, 0.880}; // eta right side
  

  // List of cuts

  // Custom mask
  bool do_custom_mask = false;
  float _z_min = -1;
  float _z_max = -1;
  float _phi_min = -1;
  float _phi_max = -1;

  // Event level cuts (aside from trigger)
  float vertex_max = 150; // cm

  // Activated trigger bit
  bool require_mbd_trigger_bit = false;
  bool require_mbd_any_vtx = false;
  bool require_mbd_vtx_10 = false;
  bool trigger_mbd_any_vtx = false; // Event with MinBias trigger (no vtx requirement)
  bool trigger_mbd_vtx_10 = false; // Event with MinBias trigger (|vertex| < 10 cm requirement)
  bool trigger_mbd = false; // Event with MinBias detector trigger (scaled)
  bool mbd_trigger_bit_event = false;

  bool require_photon_trigger_bit = false;
  bool require_photon_3_any_vtx = false;
  bool require_photon_3_vtx_10 = false;
  bool require_photon_4_any_vtx = false;
  bool require_photon_4_vtx_10 = false;
  bool trigger_mbd_photon_3_any_vtx = false;
  bool trigger_mbd_photon_4_any_vtx = false;
  bool trigger_mbd_photon_3_vtx_10 = false;
  bool trigger_mbd_photon_4_vtx_10 = false;
  bool trigger_mbd_photon_3 = false; // Event with photon 3 GeV trigger (scaled)
  bool trigger_mbd_photon_4 = false; // Event with photon 4 GeV trigger (scaled)
  bool photon_trigger_bit_event = false;
  
  // Cluster level cuts
  int n_chi2_cuts = 0;
  std::vector<float> chi2_cuts;
  std::vector<std::string> chi2_cuts_labels;
  int n_ecore_cuts = 0;
  std::vector<float> ecore_cuts;
  std::vector<std::string> ecore_cuts_labels;
  
  // diphoton level cuts
  float alphaCut = 1.0;
  float pTCutMin = 1.0;
  float pTCutMax = 1000.0;
  float pTCutThreshold = 0;

  // trigger selection (emulator)
  bool require_emulator_matching = false;
  bool emulator_match = false;
  int emulator_selection = 0;
  std::vector<int> fired_indices;
  int tileNumber = 0;
  int adc_threshold_3 = 13;
  int adc_threshold_4 = 17;
  std::map<int, std::vector<float>> emulator_energies;
  std::map<int, std::vector<int>> emulator_clusters;
  std::map<int, int> emulator_multiplicities;
  
  // trigger efficiency matching (Cluster energy above 70% efficiency threshold of trigger)
  bool require_efficiency_matching = false;
  bool efficiency_match = false;
  int efficiency_index = 6; // 70% threshold
  float energy_threshold_3 = 3.5;
  float energy_threshold_4 = 4.3;
  // 95 -> 4.2/5.3
  // 70 -> 3.5/4.3
  static constexpr int nThresholds = 10;
  const float efficiency_thresholds[nThresholds] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 95}; // efficiency thresholds
  const float array_energy_threshold_3[nThresholds] = {2.4, 2.7, 2.9, 3.1, 3.2, 3.4, 3.5, 3.7, 4.0, 4.3}; // energy thresholds for the photon 3 GeV trigger
  const float array_energy_threshold_4[nThresholds] = {2.9, 3.3, 3.6, 3.7, 3.9, 4.1, 4.3, 4.5, 4.9, 5.3}; // energy thresholds for the photon 4 GeV trigger
  static constexpr float delta_eta_threshold = 0.192; // pseudo-rapidity distance between 8 towers
  static constexpr float delta_phi_threshold = 0.192; // azimuthal distance between 8 towers
  
  // Histograms

  bool store_qa = false; // Store QA histograms in output ROOT file
  
  // Trigger histograms
  TH1I *h_emulator_selection;
  TH1F *h_trigger_live;
  TH1F *h_trigger_scaled;
  const std::map<unsigned int, std::string> trigger_names{
    {0, "Clock"}, {3, "ZDC Coincidence"}, {10, "MBD N&S >= 1"}, {12, "MBD N&S >=1, vtx < 10"}, {17, "Jet 8 GeV + MBD NS >= 1"}, {18, "Jet 10 GeV + MBD NS >= 1"}, {19, "Jet 12 GeV + MBD NS >= 1"}, {21, "Jet 8 GeV"}, {22, "Jet 10 GeV"}, {23, "Jet 12 GeV"}, {25, "Photon 3 GeV + MBD NS >= 1"}, {26, "Photon 4 GeV + MBD NS >= 1"}, {27, "Photon 5 GeV + MBD NS >= 1"}, {29, "Photon 3 Gev"}, {30, "Photon 4 Gev"}, {31, "Photon 5 Gev"}
  };

  // Cluster info node
  TriggerTile *emcaltiles;
  GlobalVertexMap *vertexmap;
  ClusterSmallInfoContainer *_smallclusters;

  // Photon info container (per event)
  std::vector<Cluster> good_photons;
  int num_photons = 0;

  // False positive histograms

  TH1I *h_efficiency_x_matching = nullptr;
  TH1I *h_efficiency_x_matching_pt[nPtBins] = {nullptr};

  // Misc QA histograms
  TH1I *h_cluster_multiplicity;
  TH1I *h_count_diphoton_nomatch;
  TH1F *h_event_zvtx;
  TH1F *h_triggered_event_zvtx;
  TH1F *h_triggered_event_photons;
  TH1F *h_analysis_event_zvtx;
  TH1F *h_photon_eta;
  TH1F *h_photon_phi;
  TH1F *h_photon_pt;
  TH1F *h_photon_zvtx;
  TH2F *h_photon_eta_phi;
  TH2F *h_photon_eta_phi_bin[nPtBins] = {nullptr};
  TH2F *h_photon_z_phi_bin[nPtBins] = {nullptr};
  TH2F *h_photon_eta_pt;
  TH2F *h_photon_eta_zvtx;
  TH2F *h_photon_phi_pt;
  TH2F *h_photon_phi_zvtx;
  TH2F *h_photon_pt_zvtx;
  TH1F *h_selected_photon_eta;
  TH1F *h_selected_photon_phi;
  TH1F *h_selected_photon_pt;
  TH1F *h_selected_photon_zvtx;
  TH2F *h_selected_photon_eta_phi;
  TH2F *h_selected_photon_eta_phi_bin[nPtBins] = {nullptr};
  TH2F *h_selected_photon_z_phi_bin[nPtBins] = {nullptr};
  TH2F *h_selected_photon_eta_pt;
  TH2F *h_selected_photon_eta_zvtx;
  TH2F *h_selected_photon_phi_pt;
  TH2F *h_selected_photon_phi_zvtx;
  TH2F *h_selected_photon_pt_zvtx;
  TH1F *h_pair_eta;
  TH1F *h_pair_phi;
  TH1F *h_pair_pt;
  TH1F *h_pair_zvtx;
  TH2F *h_pair_eta_phi;
  TH2F *h_pair_eta_pt;
  TH2F *h_pair_eta_zvtx;
  TH2F *h_pair_phi_pt;
  TH2F *h_pair_phi_zvtx;
  TH2F *h_pair_pt_zvtx;
  TH1F *h_pair_E1;
  TH1F *h_pair_E2;
  TH1F *h_pair_theta_12;
  TH1F *h_pair_DeltaR;
  TH1F *h_pair_alpha;
  TH2F *h_pair_alpha_pt;
  TH1F *h_pair_pi0_eta_pt_1;
  TH1F *h_pair_pi0_eta_pt_2;
  TH1F *h_pair_pi0_eta_pt_3;
  TH1F *h_pair_pi0_eta_pt_4;
  TH1F *h_pair_pi0_eta_pt_5;
  TH1F *h_pair_pi0_eta_pt_6;
  TH1F *h_pair_pi0_eta_pt_7;
  TH1F *h_pair_pi0_eta_pt_8;
  TH1F *h_pair_pi0_eta_pt_9;
  TH1F *h_pair_eta_eta_pt_1;
  TH1F *h_pair_eta_eta_pt_2;
  TH1F *h_pair_eta_eta_pt_3;
  TH1F *h_pair_eta_eta_pt_4;
  TH1F *h_pair_eta_eta_pt_5;
  TH1F *h_pair_eta_eta_pt_6;
  TH1F *h_pair_eta_eta_pt_7;
  TH1F *h_pair_eta_eta_pt_8;
  TH1F *h_pair_eta_eta_pt_9;
  TH1F *h_pair_pi0_xf_pt_1;
  TH1F *h_pair_pi0_xf_pt_2;
  TH1F *h_pair_pi0_xf_pt_3;
  TH1F *h_pair_pi0_xf_pt_4;
  TH1F *h_pair_pi0_xf_pt_5;
  TH1F *h_pair_pi0_xf_pt_6;
  TH1F *h_pair_pi0_xf_pt_7;
  TH1F *h_pair_pi0_xf_pt_8;
  TH1F *h_pair_pi0_xf_pt_9;
  TH1F *h_pair_eta_xf_pt_1;
  TH1F *h_pair_eta_xf_pt_2;
  TH1F *h_pair_eta_xf_pt_3;
  TH1F *h_pair_eta_xf_pt_4;
  TH1F *h_pair_eta_xf_pt_5;
  TH1F *h_pair_eta_xf_pt_6;
  TH1F *h_pair_eta_xf_pt_7;
  TH1F *h_pair_eta_xf_pt_8;
  TH1F *h_pair_eta_xf_pt_9;
  TH1F *h_meson_pi0_pt;
  TH1F *h_meson_pi0_E;
  TH1F *h_meson_eta_pt;
  TH1F *h_meson_eta_E;

  TH1F *h_photon_fired;
  TH1F *h_photon_fired_scale;

  // Cluster level cuts
  TH2F *h_cluster_level_cuts_total;
  TH2F *h_cluster_level_cuts_particle[nParticles][nRegions];
  TH2F *h_cluster_level_cuts_total_pt[nPtBins];
  TH2F *h_cluster_level_cuts_particle_pt[nParticles][nRegions][nPtBins];

  // Histograms for the average bin values
  TH1F* h_average_pt[nParticles];
  TH1F* h_average_eta[nParticles];
  TH1F* h_average_xf[nParticles];
  TH1F* h_norm_pt[nParticles];
  TH1F* h_norm_eta[nParticles];
  TH1F* h_norm_xf[nParticles];
  
  // Invariant mass histograms;
  TH1F *h_pair_mass;
  TH1F *h_pair_mass_pt[nPtBins];
  TH1F *h_pair_mass_zvtx[nZvtxBins];
  TH1F *h_pair_mass_eta[nEtaBins];
  TH1F *h_pair_mass_xf[nXfBins];
  TH1F *h_pair_mass_mixing;
  TH1F *h_pair_mass_mixing_pt[nPtBins];

  // Beam- spin- and kinematic-dependent yields -> pT dependent
  TH1F *h_yield_1[nBeams][nParticles][nRegions][nPtBins][nEtaRegions][nSpins]; // low vs high |eta|
  TH1F *h_yield_2[nBeams][nParticles][nRegions][nPtBins][nDirections][nSpins]; // forward vs backward
  TH1F *h_yield_3[nBeams][nParticles][nRegions][nPtBins][nEtaRegions][nDirections][nSpins]; // altogether
  
  // Beam- spin- and kinematic-dependent yields -> eta dependent
  TH1F *h_yield_eta[nBeams][nParticles][nRegions][nEtaBins][nSpins];

  // Beam- spin- and kinematic-dependent yields -> xf dependent
  TH1F *h_yield_xf[nBeams][nParticles][nRegions][nXfBins][nSpins];

  TH1F *h_multiplicity_efficiency;
  TH1F *h_multiplicity_emulator;
  TH1F *h_matching_consistency;

  // Geometry (for tile matching)
  RawTowerGeomContainer *towergeom;
  static constexpr double radius = 103; // cm
  
  int count_1 = 0;
  int count_2 = 0;
  int count_3 = 0;
  int count_4 = 0;
  int count_5 = 0;
  int count_yellow_pi0_peak_pt_0 = 0;

  // Scaledown list:
  int scaledown[64] = {0};
  double livetime[64] = {0};
  int prescale = 0;
  
  // For the SQL access
  std::string db_name = "daq";
  std::string user_name = "phnxrc";
  std::string table_name = "gl1_scaledown";


  // For event mixing
  bool use_event_mixing = false;
  static constexpr float dR_min = 0.034; // Minimum distance between two clusters in reconstruction (unused)
  static constexpr int nPoolZvtxBins = 10;
  const float poolZvtxBins[nPoolZvtxBins + 1] = {-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50};
  static constexpr int nMultiplicityBins = 9;
  const float multiplicityBins[nMultiplicityBins + 1] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  static constexpr int nPoolEtaBins = 12;
  const float poolEtaBins[nPoolEtaBins + 1] = {-1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2};
  static constexpr int Nmix = 10;
  std::deque<EventInPool> pool[nPoolZvtxBins][nMultiplicityBins][nPoolEtaBins];

  // Event mixing diagnosis
  TH1F *h_current_E = nullptr;
  TH1F *h_current_eta = nullptr;
  TH1F *h_current_phi = nullptr;
  TH1F *h_pool_E = nullptr;
  TH1F *h_pool_eta = nullptr;
  TH1F *h_pool_phi = nullptr;

  TH1F *h_same_delta_eta = nullptr;
  TH1F *h_same_delta_phi = nullptr;
  TH1F *h_same_delta_R = nullptr;
  TH1F *h_same_alpha = nullptr;
  TH1F *h_mixed_delta_eta = nullptr;
  TH1F *h_mixed_delta_phi = nullptr;
  TH1F *h_mixed_delta_R = nullptr;
  TH1F *h_mixed_alpha = nullptr;

  SyncObject *syncobject{nullptr};

  // MBD small info for timing
  MBDSmallInfo *m_mbd_small_info = nullptr;
  bool require_mbd_timing = false;
  float mbd_rms_cut_min = 0;
  float mbd_rms_cut_max = 30; // ns
};

#endif
