#ifndef __AN_NEUTRAL_MESON_NANO_H__
#define __AN_NEUTRAL_MESON_NANO_H__

#include <fun4all/SubsysReco.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>

// Forward declarations
class PHCompositeNode;

class AnNeutralMeson_nano : public SubsysReco
{
 public:
  
  //! constructor
  AnNeutralMeson_nano(const std::string &name = "AnNeutralMeson_nano", const std::string &inputlistname = "inputruns.txt", const std::string& inputfiletemplate = "analysis_diphoton_minimal_", const std::string &outputfilename = "analysis_0.root");

  //! destructor
  virtual ~AnNeutralMeson_nano();

  //! full initialization
  int Init(PHCompositeNode *);

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

  int FindBinDirect(float val, float valMin, float valMax, int nBins)
  {
    if (val < valMin || val > valMax) { 
      return -1;
    }
    int ibin = static_cast<int>((val - valMin) * nBins / (valMax - valMin));
    if (ibin == nBins) ibin = nBins - 1;
    else if (ibin < 0) ibin = 0;
    return ibin;
  }

  void set_sigma_number(const int val) {
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
    }
    std::cout << "band_limits = ";
    for (int j = 0; j < N; j++) {
      std::cout << band_limits[j] << " ";
    }
    std::cout << std::endl;
  }

  void shuffle_spin_pattern(const int irun);

  void set_store_bunch_yields(bool val,
                              const std::string& filename_template)
  {
    store_bunch_yields = val;
    outbunchtemplate = filename_template;
  }
  
  //! Set diphoton pT cut in pi0 mass range
  void set_ptcut(float pTMin = 1.0, float pTMax = 1000.0)
  {
    pTCutMin = pTMin;
    pTCutMax = pTMax;
  }

  //! Book histograms
  void BookHistos(const std::string &outputfilename = "");

  void SaveBunchYields(const std::string &outputfilename);

  // Give angle in the [-pi, pi] range
  float WrapAngle(const float);

  // Custom cuts
  void set_phenix_cut(bool val) { require_phenix_cut = val; } // |eta| < 0.35
  void set_high_xf_cut(bool val) { require_high_xf_cut = val; } // xF > 0.035
  void set_low_xf_cut(bool val) { require_low_xf_cut = val; } // xF < 0.035
  void set_low_vtx_cut(bool val, float thresh=30)
  {
    require_low_vtx_cut = val;
    vertex_min = thresh;
  } // |zvtx| < 30 cm
  void set_high_vtx_cut(bool val, float thresh=30)
  {
    require_high_vtx_cut = val;
    vertex_max = thresh;
  } // |zvtx| > 30 cm

  // Set trigger
  void set_trigger_mbd(bool val) { trigger_mbd = val; }
  void set_trigger_photon(bool val) { trigger_photon = val; }
 
 private:

  // Trigger selection -> changes the pT thresholds
  bool trigger_mbd = false;
  bool trigger_photon = false;
  
  // run number
  std::vector<int> runList;
  int nRuns = 0;

  // Seed number -> shuffle the spin pattern if non-zero
  int seednb = 0;

  // TTree file and object
  std::string inlistname;
  std::string infiletemplate;
  std::string treename;

  // tree branches
  float diphoton_vertex_z;
  int diphoton_bunchnumber;
  float diphoton_mass;
  float diphoton_eta;
  float diphoton_phi;
  float diphoton_pt;
  float diphoton_xf;

  // Special cuts
  bool require_phenix_cut = false;
  bool require_high_xf_cut = false;
  bool require_low_xf_cut = false;
  bool require_low_vtx_cut = false;
  bool require_high_vtx_cut = false;

  // Only when vertex cut is applied
  float vertex_min = 30; // cm
  float vertex_max = 30; // cm

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

  // New binning -> equally distributed
  static constexpr int nEtaBins = 8;
  const float etaBins[nEtaBins + 1] = {-2.00, -1.05, -0.86, -0.61, 0.0, 0.61, 0.86, 1.05, 2.0};

  // New binning -> equally distributed
  static constexpr int nXfBins = 8;
  const float xfBins[nXfBins + 1] = {-0.200, -0.048, -0.035, -0.022, 0.0, 0.022, 0.035, 0.048, 0.200};

  static constexpr int nZvtxBins = 11;
  const float zvtxBins[nZvtxBins + 1] = {0, 10, 20, 30, 40, 50, 60, 80, 100, 120, 150, 200};
  
  static constexpr int nDirections = 2;
  const std::string directions[nDirections] = {"forward", "backward"};
  static constexpr int nPhiBins = 12;
  const float phiMin = - M_PI;
  const float phiMax = + M_PI;
  
  // Constants
  const float phi_shift[nBeams] = {M_PI / 2, -M_PI / 2};
  const float beamDirection[nBeams] = {-1, 1}; // yellow / blue
  static constexpr int nEtaRegions = 2;
  static constexpr double etaThreshold = 0.35;

  // List spin info
  static constexpr int nBunches = 120;
  int crossingshift;
  int beamspinpat[nBeams][nBunches];

  // pT Cut
  float pTCutMin = 1.0;
  float pTCutMax = 1000.0;

  // Output histogram file
  std::string outfiletemplate;
  TFile *outfile = nullptr;

  // Invariant mass histograms
  TH1F *h_pair_mass;
  TH1F *h_pair_mass_pt[nPtBins];
  TH1F *h_pair_mass_eta[nEtaBins];
  TH1F *h_pair_mass_xf[nXfBins];
  TH1F *h_pair_mass_zvtx[nZvtxBins];

  // Histograms for the average bin values
  TH1F* h_average_pt[nParticles];
  TH1F* h_average_eta[nParticles];
  TH1F* h_average_xf[nParticles];
  TH1F* h_norm_pt[nParticles];
  TH1F* h_norm_eta[nParticles];
  TH1F* h_norm_xf[nParticles];

  // Kinematic correlations
  TH1F *h_pair_meson_zvtx[2] = {nullptr};
  TH1F *h_pair_meson_pt_eta[2][9] = {nullptr};
  TH1F *h_pair_meson_pt_xf[2][9] = {nullptr};
  TH1F *h_pair_meson_eta_pt[2][8] = {nullptr};
  TH1F *h_pair_meson_eta_xf[2][8] = {nullptr};
  TH1F *h_pair_meson_xf_pt[2][8] = {nullptr};
  TH1F *h_pair_meson_xf_eta[2][8] = {nullptr};
  TH2F *h_pair_meson_2D_xf_pt[2] = {nullptr};
  TH2F *h_pair_meson_2D_xf_eta[2] = {nullptr};
  TH2F *h_pair_meson_2D_xf_zvtx[2] = {nullptr};
  TH2F *h_pair_meson_2D_pt_eta[2] = {nullptr};
  TH2F *h_pair_meson_2D_pt_zvtx[2] = {nullptr};
  TH2F *h_pair_meson_2D_eta_zvtx[2] = {nullptr};

  // Output bunch yields, for bunch shuffling
  bool store_bunch_yields = false;
  std::string outbunchtemplate = "";
  uint32_t array_yield_pt[nBeams][nParticles][nRegions][nPtBins][nDirections][nSpins][nPhiBins][nBunches] = {0};
  uint32_t array_yield_eta[nBeams][nParticles][nRegions][nEtaBins][nSpins][nPhiBins][nBunches] = {0};
  uint32_t array_yield_xf[nBeams][nParticles][nRegions][nXfBins][nSpins][nPhiBins][nBunches] = {0};
  
  // Beam- spin- and kinematic-dependent yields -> pT dependent
  TH1F *h_yield_pt[nBeams][nParticles][nRegions][nPtBins][nDirections][nSpins]; // forward vs backward
  
  // Beam- spin- and kinematic-dependent yields -> eta dependent
  TH1F *h_yield_eta[nBeams][nParticles][nRegions][nEtaBins][nSpins];

  // Beam- spin- and kinematic-dependent yields -> xf dependent
  TH1F *h_yield_xf[nBeams][nParticles][nRegions][nXfBins][nSpins];

  // Define the regions (in invariant mass) for pi0/eta peak/side
  int sigma_number = 3;
  float band_limits_3[nParticles * (nRegions + 1) * 2] =
    {0.030, 0.070, // pi0 left side invariant mass range (in GeV/c^2)
     0.080, 0.199, // pi0 peak
     0.209, 0.249, // pi0 right side
     0.257, 0.371, // eta left side
     0.399,0.739, // eta peak
     0.767, 0.880}; // eta right side
  float band_limits_2[nParticles * (nRegions + 1) * 2] =
    {0.030, 0.070, // pi0 left side invariant mass range (in GeV/c^2)
     0.100, 0.180, // pi0 peak
     0.209, 0.249, // pi0 right side
     0.257, 0.371, // eta left side
     0.456,0.683, // eta peak
     0.767, 0.880}; // eta right side
  float band_limits_15[nParticles * (nRegions + 1) * 2] =
    {0.030, 0.070, // pi0 left side invariant mass range (in GeV/c^2)
     0.110, 0.170, // pi0 peak
     0.209, 0.249, // pi0 right side
     0.257, 0.371, // eta left side
     0.484,0.654, // eta peak
     0.767, 0.880}; // eta right side
  float band_limits[nParticles * (nRegions + 1) * 2] =
    {0.030, 0.070, // pi0 left side invariant mass range (in GeV/c^2)
     0.080, 0.199, // pi0 peak
     0.209, 0.249, // pi0 right side
     0.257, 0.371, // eta left side
     0.399,0.739, // eta peak
     0.767, 0.880}; // eta right side
};

#endif
