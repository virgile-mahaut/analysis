#ifndef UPCMESON_H
#define UPCMESON_H

#include <Rtypes.h>
#include <fun4all/SubsysReco.h>

/// Class declarations for use in the analysis module
class PHCompositeNode;
class TFile;
class TTree;
class TH1;
class TH2;
class PHCompositeNode;
class SvtxTrackMap;
class GlobalVertex;
class SvtxTrackEval;
class PHG4TruthInfoContainer;
class PHHepMCGenEventMap;
class PHHepMCGenEvent;
class CaloTriggerInfo;
class SvtxEvalStack;
class EventHeader;
class Gl1Packet;
//class TrackVertexCrossingAssoc;

/// Definition of this analysis module class
class UPCMeson : public SubsysReco
{
 public:
  /// Constructor
  UPCMeson(const std::string &name = "UPCMeson",
              const std::string &outputfname = "UPCMeson.root");

  // Destructor
  virtual ~UPCMeson();

  /// SubsysReco initialize processing method
  int Init(PHCompositeNode *) override;

  /// SubsysReco event processing method
  int process_event(PHCompositeNode *) override;

  /// SubsysReco end processing method
  int End(PHCompositeNode *) override;

  void SetGuessMass(const double m) { _mguess = m; }

  /// Set things to analyze
  void analyzeTracks(bool analyzeTracks) { m_analyzeTracks = analyzeTracks; }
  void analyzeTruth(bool analyzeTruth) { m_analyzeTruth = analyzeTruth; }

  int Reset(PHCompositeNode * /*topNode*/) override;

 private:
  /// String to contain the outfile name containing the trees
  std::string m_outfilename;

  /// A boolean for running over tracks
  bool m_analyzeTracks;

  /// A boolean for collecting hepmc information
  bool m_analyzeTruth;

  /// TFile to hold the following TTrees and histograms
  TFile *m_outfile{nullptr};
  TTree *m_tracktree{nullptr};
  TTree *m_hepmctree{nullptr};
  TTree *m_truthtree{nullptr};
  TTree *m_globaltree{nullptr};
  TTree *m_pairtree{nullptr};
 
  TH1 *h_phi[2]{nullptr,nullptr};       // [0]=opp. sign, [1]=like sign
  TH2 *h2_eta_phi[2]{nullptr,nullptr};
  TH1 *h_mass[2]{nullptr,nullptr};
  TH1 *h_pt[2]{nullptr,nullptr};
  TH1 *h_y[2]{nullptr,nullptr};
  TH1 *h_eta[2]{nullptr,nullptr};

  TH1 *h_trig{nullptr};
  TH1 *h_ntracks{nullptr};
  TH2 *h2_ntrksvsb{nullptr};
  TH1 *h_cross{nullptr};      //crossings with low num tracks
  TH1 *h_cross_evt{nullptr};  //crossings sampled
  TH1 *h_bunch{nullptr};      //bunches from cross 0

  TH1 *h_b_mb{nullptr};
  TH1 *h_npart_mb{nullptr};
  TH1 *h_ncoll_mb{nullptr};
  TH1 *h_b{nullptr};
  TH1 *h_npart{nullptr};
  TH1 *h_ncoll{nullptr};

  const double E_MASS = 0.000510998950;  // electron mass [Gev]
  const double MU_MASS = 0.000510998950;  // muon mass [Gev]
  const double PI_MASS = 0.13957039;  // pion mass [Gev]
  double _mguess{E_MASS};

  static inline constexpr int CROSS_OFFSET = 150; // used to map cross to a pos. array index

  /// Methods for grabbing the data
  int GetNodes(PHCompositeNode *topNode);
  int getTracks(PHCompositeNode *topNode);
  void getHEPMCTruth();
  void getPHG4Truth();

  /// Data
  EventHeader *_evthdr{nullptr};
  Gl1Packet* _gl1raw{nullptr};
  SvtxTrackMap *_trackmap{nullptr};
  PHG4TruthInfoContainer *_truthinfo{nullptr};
  PHHepMCGenEventMap *_genevent_map{nullptr};
  //TrackVertexCrossingAssoc* _track_vertex_crossing_map{nullptr};
  SvtxEvalStack *m_svtxEvalStack{nullptr};

  void initializeVariables();
  void initializeTrees();

  /**
   * Variables for the trees
   */

  // Global Info
  Int_t    m_run{ 0 };
  Int_t    m_evt{ 0 };
  Short_t  m_cross{ -1 };  // Streaming cross number
  Short_t  m_bunch{ 0 };  // RHIC bunch number
  uint64_t m_strig{ 0 };  // Scaled Trigger
  Int_t    m_npart_targ{ 0 };
  Int_t    m_npart_proj{ 0 };
  Int_t    m_npart{ 0 };
  Int_t    m_ncoll{ 0 };
  Int_t    m_ncoll_hard{ 0 };
  Float_t  m_bimpact{ -1. };
  Int_t    m_ntrks{ 0 };         // total reconstructed ntracks
  Int_t    m_ntrks_cross{ 0 };   // total ntracks in a crossing
  Int_t    m_ntrk_sphenix{ 0 };  // monte carlo tracks in sphenix acceptance (>0.4 pt, |eta|<1.1)
  Int_t    m_ntrk_mc{ 0 };       // total monte carlo charged tracks

  /// HEPMC Tree variables
  int m_partid1;
  int m_partid2;
  double m_x1;
  double m_x2;
  int m_mpi;
  int m_process_id;
  double m_truthenergy;
  double m_trutheta;
  double m_truthphi;
  double m_truthpx;
  double m_truthpy;
  double m_truthpz;
  double m_truthpt;
  double m_truthp;
  int m_numparticlesinevent;
  int m_truthpid;
  int m_truthcharge;


  /// Track variables
  double m_tr_px;
  double m_tr_py;
  double m_tr_pz;
  double m_tr_p;
  double m_tr_pt;
  double m_tr_phi;
  double m_tr_eta;
  int m_charge;
  double m_chisq;
  int m_ndf;
  double m_dca;
  double m_tr_x;
  double m_tr_y;
  double m_tr_z;
  int m_truth_is_primary;
  double m_truthtrackpx;
  double m_truthtrackpy;
  double m_truthtrackpz;
  double m_truthtrackp;
  double m_truthtracke;
  double m_truthtrackpt;
  double m_truthtrackphi;
  double m_truthtracketa;
  int m_truthtrackpid;

  /// Pair variables
  Float_t m_pm{ 0. };    // pair mass
  Float_t m_ppt{ 0. };
  Float_t m_pphi{ 0. };
  Float_t m_py{ 0. };
  Float_t m_peta{ 0. };
  Float_t m_pdphi{ 0. };
  Float_t m_ppt1{ 0. };
  Float_t m_ppz1{ 0. };
  Float_t m_pphi1{ 0. };
  Float_t m_peta1{ 0. };
  Float_t m_ppt2{ 0. };
  Float_t m_ppz2{ 0. };
  Float_t m_pphi2{ 0. };
  Float_t m_peta2{ 0. };
  Short_t  m_pq1{ 0 };
  Short_t  m_pq2{ 0 };
};

#endif  // UPCMESON_H
