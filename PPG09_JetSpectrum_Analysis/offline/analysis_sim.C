#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <tuple>
#include <algorithm>
#include "TRandom.h"
#include "TRandom3.h"
#include "unfold_Def.h"
#include "BootstrapGenerator/BootstrapGenerator.h"
#include "BootstrapGenerator/TH1DBootstrap.h"
#include "BootstrapGenerator/TH2DBootstrap.h"

R__LOAD_LIBRARY(libBootstrapGenerator.so)

void get_leading_subleading_jet(int& leadingjet_index, int& subleadingjet_index, std::vector<float>* jet_et);
bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi);
void get_recojet_filter(std::vector<bool>& jet_filter, std::vector<float>* jet_e, std::vector<float>* jet_pt, std::vector<float>* jet_eta, float zvertex, float jet_radius);
void get_truthjet_filter(std::vector<bool>& jet_filter, std::vector<float>* jet_e, std::vector<float>* jet_pt, std::vector<float>* jet_eta, float zvertex, float jet_radius);
void get_smear_value(std::vector<float>& smear_value, std::vector<bool>& jet_filter);
void get_calibjet(std::vector<float>& calibjet_pt, std::vector<float>& calibjet_eta, std::vector<float>& calibjet_phi, std::vector<bool>& jet_filter, std::vector<float>* jet_pt_calib, std::vector<float>* jet_eta, std::vector<float>* jet_phi, double jet_radius, float jes_para, float jer_para, std::vector<float>& smear_value, bool jet_background);
void get_truthjet(std::vector<float>& goodtruthjet_pt, std::vector<float>& goodtruthjet_eta, std::vector<float>& goodtruthjet_phi, std::vector<bool>& truthjet_filter, std::vector<float>* jet_pt, std::vector<float>* jet_eta, std::vector<float>* jet_phi);
void match_meas_truth(std::vector<float>& meas_eta, std::vector<float>& meas_phi, std::vector<float>& meas_matched, std::vector<float>& truth_eta, std::vector<float>& truth_phi, std::vector<float>& truth_matched, float jet_radius, bool has_zvertex);
void fill_response_matrix(TH1D*& h_truth, TH1D*& h_meas, TH2D*& h_resp, TH1D*& h_fake, TH1D*& h_miss, double scale, TF1* f_reweight, std::vector<float>& meas_pt, std::vector<float>& meas_matched, std::vector<float>& truth_pt, std::vector<float>& truth_matched);
void fill_response_matrix_altz(TH1D*& h_truth, TH1D*& h_meas, TH2D*& h_resp, TH1D*& h_fake, TH1D*& h_miss, double scale, TF1* f_reweight, std::vector<float>& meas_pt, std::vector<float>& meas_matched, std::vector<float>& truth_pt, std::vector<float>& truth_matched, bool is_mbd_zvertex60, bool is_truth_zvertex60, bool has_zvertex);
void fill_response_matrix_stat(TH1DBootstrap*& h_truth, TH1DBootstrap*& h_meas, TH2DBootstrap*& h_resp, TH1DBootstrap*& h_fake, TH1DBootstrap*& h_miss, double scale, TF1* f_reweight, std::vector<float>& meas_pt, std::vector<float>& meas_matched, std::vector<float>& truth_pt, std::vector<float>& truth_matched);
void fill_response_matrix_altz_stat(TH1DBootstrap*& h_truth, TH1DBootstrap*& h_meas, TH2DBootstrap*& h_resp, TH1DBootstrap*& h_fake, TH1DBootstrap*& h_miss, double scale, TF1* f_reweight, std::vector<float>& meas_pt, std::vector<float>& meas_matched, std::vector<float>& truth_pt, std::vector<float>& truth_matched, bool is_mbd_zvertex60, bool is_truth_zvertex60, bool has_zvertex);
void build_hist(std::string tag, TH1D*& h_truth, TH1D*& h_measure, TH2D*& h_respmatrix, TH1D*& h_fake, TH1D*& h_miss);
void build_hist_stat(std::string tag, TH1DBootstrap*& h_truth, TH1DBootstrap*& h_measure, TH2DBootstrap*& h_respmatrix, TH1DBootstrap*& h_fake, TH1DBootstrap*& h_miss, int nrep, BootstrapGenerator* bsgen);

TRandom3 randGen(1234);

////////////////////////////////////////// Main Function //////////////////////////////////////////
void analysis_sim(std::string runtype, int nseg, int iseg, double jet_radius, bool do_reweighting, int mrad, int interaction) {
  ////////// General Set up //////////
  int jet_radius_index = (int)(10 * jet_radius);
  double truthjet_pt_min = 0, truthjet_pt_max = 1000, calibjet_pt_max = 1000;
  if (runtype == "MB") {
    if (jet_radius == 0.2) {truthjet_pt_min = 0; truthjet_pt_max = 5; calibjet_pt_max = 17;}
    else if (jet_radius == 0.3) {truthjet_pt_min = 0; truthjet_pt_max = 6; calibjet_pt_max = 18;}
    else if (jet_radius == 0.4) {truthjet_pt_min = 0; truthjet_pt_max = 7; calibjet_pt_max = 19;}
    else if (jet_radius == 0.5) {truthjet_pt_min = 0; truthjet_pt_max = 10; calibjet_pt_max = 22;}
    else if (jet_radius == 0.6) {truthjet_pt_min = 0; truthjet_pt_max = 12; calibjet_pt_max = 24;}
    else if (jet_radius == 0.7) {truthjet_pt_min = 0; truthjet_pt_max = 14; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.8) {truthjet_pt_min = 0; truthjet_pt_max = 15; calibjet_pt_max = 26;}
  } else if (runtype == "Jet5GeV") {
    if (jet_radius == 0.2) {truthjet_pt_min = 5; truthjet_pt_max = 12; calibjet_pt_max = 29;}
    else if (jet_radius == 0.3) {truthjet_pt_min = 6; truthjet_pt_max = 13; calibjet_pt_max = 30;}
    else if (jet_radius == 0.4) {truthjet_pt_min = 7; truthjet_pt_max = 14; calibjet_pt_max = 30;}
    else if (jet_radius == 0.5) {truthjet_pt_min = 10; truthjet_pt_max = 19; calibjet_pt_max = 31;}
    else if (jet_radius == 0.6) {truthjet_pt_min = 11; truthjet_pt_max = 22; calibjet_pt_max = 31;}
    else if (jet_radius == 0.7) {truthjet_pt_min = 14; truthjet_pt_max = 24; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.8) {truthjet_pt_min = 15; truthjet_pt_max = 25; calibjet_pt_max = 33;}
  } else if (runtype == "Jet12GeV") {
    if (jet_radius == 0.2) {truthjet_pt_min = 12; truthjet_pt_max = 20; calibjet_pt_max = 47;}
    else if (jet_radius == 0.3) {truthjet_pt_min = 13; truthjet_pt_max = 21; calibjet_pt_max = 50;}
    else if (jet_radius == 0.4) {truthjet_pt_min = 14; truthjet_pt_max = 21; calibjet_pt_max = 49;}
    else if (jet_radius == 0.5) {truthjet_pt_min = 19; truthjet_pt_max = 27; calibjet_pt_max = 50;}
    else if (jet_radius == 0.6) {truthjet_pt_min = 22; truthjet_pt_max = 29; calibjet_pt_max = 52;}
    else if (jet_radius == 0.7) {truthjet_pt_min = 24; truthjet_pt_max = 32; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.8) {truthjet_pt_min = 25; truthjet_pt_max = 34; calibjet_pt_max = 52;}
  } else if (runtype == "Jet20GeV") {
    if (jet_radius == 0.2) {truthjet_pt_min = 20; truthjet_pt_max = 30; calibjet_pt_max = 64;}
    else if (jet_radius == 0.3) {truthjet_pt_min = 21; truthjet_pt_max = 31; calibjet_pt_max = 66;}
    else if (jet_radius == 0.4) {truthjet_pt_min = 21; truthjet_pt_max = 32; calibjet_pt_max = 69;}
    else if (jet_radius == 0.5) {truthjet_pt_min = 27; truthjet_pt_max = 38; calibjet_pt_max = 70;}
    else if (jet_radius == 0.6) {truthjet_pt_min = 29; truthjet_pt_max = 41; calibjet_pt_max = 72;}
    else if (jet_radius == 0.7) {truthjet_pt_min = 32; truthjet_pt_max = 45; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.8) {truthjet_pt_min = 34; truthjet_pt_max = 47; calibjet_pt_max = 70;}
  } else if (runtype == "Jet30GeV") {
    if (jet_radius == 0.2) {truthjet_pt_min = 30; truthjet_pt_max = 40; calibjet_pt_max = 92;}
    else if (jet_radius == 0.3) {truthjet_pt_min = 31; truthjet_pt_max = 41; calibjet_pt_max = 92;}
    else if (jet_radius == 0.4) {truthjet_pt_min = 32; truthjet_pt_max = 42; calibjet_pt_max = 90;}
    else if (jet_radius == 0.5) {truthjet_pt_min = 38; truthjet_pt_max = 49; calibjet_pt_max = 91;}
    else if (jet_radius == 0.6) {truthjet_pt_min = 41; truthjet_pt_max = 53; calibjet_pt_max = 92;}
    else if (jet_radius == 0.7) {truthjet_pt_min = 45; truthjet_pt_max = 57; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.8) {truthjet_pt_min = 47; truthjet_pt_max = 58; calibjet_pt_max = 86;}
  } else if (runtype == "Jet40GeV") {
    if (jet_radius == 0.2) {truthjet_pt_min = 40; truthjet_pt_max = 50; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.3) {truthjet_pt_min = 41; truthjet_pt_max = 51; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.4) {truthjet_pt_min = 42; truthjet_pt_max = 52; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.5) {truthjet_pt_min = 49; truthjet_pt_max = 58; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.6) {truthjet_pt_min = 53; truthjet_pt_max = 63; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.7) {truthjet_pt_min = 57; truthjet_pt_max = 66; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.8) {truthjet_pt_min = 58; truthjet_pt_max = 68; calibjet_pt_max = 3000;}
  } else if (runtype == "Jet50GeV") {
    if (jet_radius == 0.2) {truthjet_pt_min = 50; truthjet_pt_max = 60; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.3) {truthjet_pt_min = 51; truthjet_pt_max = 61; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.4) {truthjet_pt_min = 52; truthjet_pt_max = 62; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.5) {truthjet_pt_min = 58; truthjet_pt_max = 68; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.6) {truthjet_pt_min = 63; truthjet_pt_max = 72; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.7) {truthjet_pt_min = 66; truthjet_pt_max = 75; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.8) {truthjet_pt_min = 68; truthjet_pt_max = 76; calibjet_pt_max = 3000;}
  } else if (runtype == "Jet60GeV") {
    if (jet_radius == 0.2) {truthjet_pt_min = 60; truthjet_pt_max = 3000; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.3) {truthjet_pt_min = 61; truthjet_pt_max = 3000; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.4) {truthjet_pt_min = 62; truthjet_pt_max = 3000; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.5) {truthjet_pt_min = 68; truthjet_pt_max = 3000; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.6) {truthjet_pt_min = 72; truthjet_pt_max = 3000; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.7) {truthjet_pt_min = 75; truthjet_pt_max = 3000; calibjet_pt_max = 3000;}
    else if (jet_radius == 0.8) {truthjet_pt_min = 76; truthjet_pt_max = 3000; calibjet_pt_max = 3000;}
  }
  
  ////////// Files //////////
  std::string output_dir = "output_sim/output_sim";
  if (mrad == 0) output_dir += "_0mrad";
  else if (mrad == 1) output_dir += "_1p5mrad";
  else {
    std::cout << "Error: invalid mrad value" << std::endl;
    return;
  }
  if (interaction == 1) output_dir += "_single";
  else if (interaction == 2) output_dir += "_double";
  else {
    std::cout << "Error: invalid interaction value" << std::endl;
    return;
  }
  TFile *f_out = new TFile(Form("%s/output_r0%d_%s_%d_%d.root", output_dir.c_str(), jet_radius_index, runtype.c_str(), iseg, iseg+nseg), "RECREATE");
  TChain chain("ttree");
  for (int i = iseg; i < iseg + nseg; ++i) {
    if (interaction == 1) chain.Add(Form("/sphenix/user/hanpuj/Dataset/Run28/output_%s_%d.root", runtype.c_str(), i));
    else if (interaction == 2) chain.Add(Form("/sphenix/user/hanpuj/Dataset/Run28_double/output_%s_%d.root", runtype.c_str(), i));
  }
  chain.SetBranchStatus("*", 0);

  float zvertex; chain.SetBranchStatus("z_vertex", 1); chain.SetBranchAddress("z_vertex", &zvertex);
  float truth_zvertex; chain.SetBranchStatus("truth_z_vertex", 1); chain.SetBranchAddress("truth_z_vertex", &truth_zvertex);
  std::vector<float>* unsubjet_e = nullptr; chain.SetBranchStatus(Form("unsubjet0%d_e", jet_radius_index), 1); chain.SetBranchAddress(Form("unsubjet0%d_e", jet_radius_index), &unsubjet_e);
  std::vector<float>* unsubjet_pt = nullptr; chain.SetBranchStatus(Form("unsubjet0%d_pt", jet_radius_index), 1); chain.SetBranchAddress(Form("unsubjet0%d_pt", jet_radius_index), &unsubjet_pt);
  std::vector<float>* unsubjet_eta = nullptr; chain.SetBranchStatus(Form("unsubjet0%d_eta", jet_radius_index), 1); chain.SetBranchAddress(Form("unsubjet0%d_eta", jet_radius_index), &unsubjet_eta);
  std::vector<float>* unsubjet_phi = nullptr; chain.SetBranchStatus(Form("unsubjet0%d_phi", jet_radius_index), 1); chain.SetBranchAddress(Form("unsubjet0%d_phi", jet_radius_index), &unsubjet_phi);
  std::vector<float>* unsubjet_pt_calib = nullptr; chain.SetBranchStatus(Form("unsubjet0%d_pt_calib", jet_radius_index), 1); chain.SetBranchAddress(Form("unsubjet0%d_pt_calib", jet_radius_index), &unsubjet_pt_calib);

  std::vector<float>* truthjet_e = nullptr; chain.SetBranchStatus(Form("truthjet0%d_e", jet_radius_index), 1); chain.SetBranchAddress(Form("truthjet0%d_e", jet_radius_index), &truthjet_e);
  std::vector<float>* truthjet_pt = nullptr; chain.SetBranchStatus(Form("truthjet0%d_pt", jet_radius_index), 1); chain.SetBranchAddress(Form("truthjet0%d_pt", jet_radius_index), &truthjet_pt);
  std::vector<float>* truthjet_eta = nullptr; chain.SetBranchStatus(Form("truthjet0%d_eta", jet_radius_index), 1); chain.SetBranchAddress(Form("truthjet0%d_eta", jet_radius_index), &truthjet_eta);
  std::vector<float>* truthjet_phi = nullptr; chain.SetBranchStatus(Form("truthjet0%d_phi", jet_radius_index), 1); chain.SetBranchAddress(Form("truthjet0%d_phi", jet_radius_index), &truthjet_phi);

  // Read R=0.4 jet branches for event selection.
  std::vector<float>* unsubjet04_e_reader = nullptr, *unsubjet04_pt_reader = nullptr, *unsubjet04_eta_reader = nullptr, *unsubjet04_phi_reader = nullptr;
  if (jet_radius_index != 4) {
    chain.SetBranchStatus("unsubjet04_e", 1); chain.SetBranchAddress("unsubjet04_e", &unsubjet04_e_reader);
    chain.SetBranchStatus("unsubjet04_pt", 1); chain.SetBranchAddress("unsubjet04_pt", &unsubjet04_pt_reader);
    chain.SetBranchStatus("unsubjet04_eta", 1); chain.SetBranchAddress("unsubjet04_eta", &unsubjet04_eta_reader);
    chain.SetBranchStatus("unsubjet04_phi", 1); chain.SetBranchAddress("unsubjet04_phi", &unsubjet04_phi_reader);
  }
  std::vector<float>* &unsubjet04_e = (jet_radius_index != 4) ? unsubjet04_e_reader : unsubjet_e;
  std::vector<float>* &unsubjet04_pt = (jet_radius_index != 4) ? unsubjet04_pt_reader : unsubjet_pt;
  std::vector<float>* &unsubjet04_eta = (jet_radius_index != 4) ? unsubjet04_eta_reader : unsubjet_eta;
  std::vector<float>* &unsubjet04_phi = (jet_radius_index != 4) ? unsubjet04_phi_reader : unsubjet_phi;

  /////////////// Z-vertex Reweighting ///////////////
  std::string filename_zvertexreweight = "input_zvertexreweight";
  if (mrad == 0) filename_zvertexreweight += "_0mrad.root";
  else if (mrad == 1) filename_zvertexreweight += "_1p5mrad.root";
  else {
    std::cout << "Error: invalid mrad value" << std::endl;
    return;
  }
  TFile *f_zvertexreweight = new TFile(filename_zvertexreweight.c_str(), "READ");
  if (!f_zvertexreweight || f_zvertexreweight->IsZombie()) {
    std::cout << "Error: cannot open " << filename_zvertexreweight << std::endl;
    return;
  }
  TH1D *h_zvertex_reweight = (TH1D*)f_zvertexreweight->Get("h_zvertex_reweight");
  if (!h_zvertex_reweight) {
    std::cout << "Error: cannot open h_zvertex_reweight" << std::endl;
    return;
  }

  /////////////// Unfolding Reweighting ///////////////
  TFile *f_reweight;
  TF1 *f_reweightfunc_all, *f_reweightfunc_all_jesup, *f_reweightfunc_all_jesdown, *f_reweightfunc_all_jerup, *f_reweightfunc_all_jerdown, *f_reweightfunc_all_jetup, *f_reweightfunc_all_jetdown;
  TF1 *f_reweightfunc_zvertex60, *f_reweightfunc_zvertex60_jesup, *f_reweightfunc_zvertex60_jesdown, *f_reweightfunc_zvertex60_jerup, *f_reweightfunc_zvertex60_jerdown, *f_reweightfunc_zvertex60_jetup, *f_reweightfunc_zvertex60_jetdown;
  if (do_reweighting) {
    f_reweight = new TFile(Form("output_reweightfunction_r0%d.root", jet_radius_index), "READ");
    f_reweightfunc_all = (TF1*)f_reweight->Get(Form("reweightfunc0%d_all", jet_radius_index));
    f_reweightfunc_all_jesup = (TF1*)f_reweight->Get(Form("reweightfunc0%d_all_jesup", jet_radius_index));
    f_reweightfunc_all_jesdown = (TF1*)f_reweight->Get(Form("reweightfunc0%d_all_jesdown", jet_radius_index));
    f_reweightfunc_all_jerup = (TF1*)f_reweight->Get(Form("reweightfunc0%d_all_jerup", jet_radius_index));
    f_reweightfunc_all_jerdown = (TF1*)f_reweight->Get(Form("reweightfunc0%d_all_jerdown", jet_radius_index));
    f_reweightfunc_all_jetup = (TF1*)f_reweight->Get(Form("reweightfunc0%d_all_jetup", jet_radius_index));
    f_reweightfunc_all_jetdown = (TF1*)f_reweight->Get(Form("reweightfunc0%d_all_jetdown", jet_radius_index));
    f_reweightfunc_zvertex60 = (TF1*)f_reweight->Get(Form("reweightfunc0%d_zvertex60", jet_radius_index));
    f_reweightfunc_zvertex60_jesup = (TF1*)f_reweight->Get(Form("reweightfunc0%d_zvertex60_jesup", jet_radius_index));
    f_reweightfunc_zvertex60_jesdown = (TF1*)f_reweight->Get(Form("reweightfunc0%d_zvertex60_jesdown", jet_radius_index));
    f_reweightfunc_zvertex60_jerup = (TF1*)f_reweight->Get(Form("reweightfunc0%d_zvertex60_jerup", jet_radius_index));
    f_reweightfunc_zvertex60_jerdown = (TF1*)f_reweight->Get(Form("reweightfunc0%d_zvertex60_jerdown", jet_radius_index));
    f_reweightfunc_zvertex60_jetup = (TF1*)f_reweight->Get(Form("reweightfunc0%d_zvertex60_jetup", jet_radius_index));
    f_reweightfunc_zvertex60_jetdown = (TF1*)f_reweight->Get(Form("reweightfunc0%d_zvertex60_jetdown", jet_radius_index));
  } else {
    f_reweightfunc_all = new TF1("f_reweightfunc_all", "1", 0, 1000);
    f_reweightfunc_all_jesup = new TF1("f_reweightfunc_all_jesup", "1", 0, 1000);
    f_reweightfunc_all_jesdown = new TF1("f_reweightfunc_all_jesdown", "1", 0, 1000);
    f_reweightfunc_all_jerup = new TF1("f_reweightfunc_all_jerup", "1", 0, 1000);
    f_reweightfunc_all_jerdown = new TF1("f_reweightfunc_all_jerdown", "1", 0, 1000);
    f_reweightfunc_all_jetup = new TF1("f_reweightfunc_all_jetup", "1", 0, 1000);
    f_reweightfunc_all_jetdown = new TF1("f_reweightfunc_all_jetdown", "1", 0, 1000);
    f_reweightfunc_zvertex60 = new TF1("f_reweightfunc_zvertex60", "1", 0, 1000);
    f_reweightfunc_zvertex60_jesup = new TF1("f_reweightfunc_zvertex60_jesup", "1", 0, 1000);
    f_reweightfunc_zvertex60_jesdown = new TF1("f_reweightfunc_zvertex60_jesdown", "1", 0, 1000);
    f_reweightfunc_zvertex60_jerup = new TF1("f_reweightfunc_zvertex60_jerup", "1", 0, 1000);
    f_reweightfunc_zvertex60_jerdown = new TF1("f_reweightfunc_zvertex60_jerdown", "1", 0, 1000);
    f_reweightfunc_zvertex60_jetup = new TF1("f_reweightfunc_zvertex60_jetup", "1", 0, 1000);
    f_reweightfunc_zvertex60_jetdown = new TF1("f_reweightfunc_zvertex60_jetdown", "1", 0, 1000);
  }
  // Unweight functions (set to 1)
  TF1 *f_reweightfunc_unreweighted = new TF1("f_reweightfunc_unreweighted", "1", 0, 1000);

  ////////// Histograms //////////
  BootstrapGenerator* bsgen = new BootstrapGenerator("bsGen","bsGen",nrep,iseg);

  TH1D *h_event_all = new TH1D("h_event_all", ";Event Number", 1, 0, 1);
  TH1D *h_event_beforecut = new TH1D("h_event_beforecut", ";Event Number", 1, 0, 1);
  TH1D *h_event_passed = new TH1D("h_event_passed", ";Event Number", 1, 0, 1);
  TH1D *h_recojet_pt_record_nocut_all = new TH1D("h_recojet_pt_record_nocut_all", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_all = new TH1D("h_recojet_pt_record_all", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_nocut_zvertex60 = new TH1D("h_recojet_pt_record_nocut_zvertex60", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_zvertex60 = new TH1D("h_recojet_pt_record_zvertex60", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_truthjet_pt_record_all = new TH1D("h_truthjet_pt_record_all", ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);

  // Nominal histograms
  TH2D *h_respmatrix_all; TH1D *h_truth_all, *h_measure_all, *h_fake_all, *h_miss_all;
  build_hist("all", h_truth_all, h_measure_all, h_respmatrix_all, h_fake_all, h_miss_all);
  TH2D *h_respmatrix_zvertex60; TH1D *h_truth_zvertex60, *h_measure_zvertex60, *h_fake_zvertex60, *h_miss_zvertex60;
  build_hist("zvertex60", h_truth_zvertex60, h_measure_zvertex60, h_respmatrix_zvertex60, h_fake_zvertex60, h_miss_zvertex60);

  TH2DBootstrap *h_respmatrix_all_stat; TH1DBootstrap *h_truth_all_stat, *h_measure_all_stat, *h_fake_all_stat, *h_miss_all_stat;
  build_hist_stat("all_stat", h_truth_all_stat, h_measure_all_stat, h_respmatrix_all_stat, h_fake_all_stat, h_miss_all_stat, nrep, bsgen);
  TH2DBootstrap *h_respmatrix_zvertex60_stat; TH1DBootstrap *h_truth_zvertex60_stat, *h_measure_zvertex60_stat, *h_fake_zvertex60_stat, *h_miss_zvertex60_stat;
  build_hist_stat("zvertex60_stat", h_truth_zvertex60_stat, h_measure_zvertex60_stat, h_respmatrix_zvertex60_stat, h_fake_zvertex60_stat, h_miss_zvertex60_stat, nrep, bsgen);

  // Un-reweighted histograms (for uncertainty)
  TH2D *h_respmatrix_all_unreweighted; TH1D *h_truth_all_unreweighted, *h_measure_all_unreweighted, *h_fake_all_unreweighted, *h_miss_all_unreweighted;
  build_hist("all_unreweighted", h_truth_all_unreweighted, h_measure_all_unreweighted, h_respmatrix_all_unreweighted, h_fake_all_unreweighted, h_miss_all_unreweighted);
  TH2D *h_respmatrix_zvertex60_unreweighted; TH1D *h_truth_zvertex60_unreweighted, *h_measure_zvertex60_unreweighted, *h_fake_zvertex60_unreweighted, *h_miss_zvertex60_unreweighted;
  build_hist("zvertex60_unreweighted", h_truth_zvertex60_unreweighted, h_measure_zvertex60_unreweighted, h_respmatrix_zvertex60_unreweighted, h_fake_zvertex60_unreweighted, h_miss_zvertex60_unreweighted);

  // JES up/down histograms
  TH2D *h_respmatrix_all_jesup; TH1D *h_truth_all_jesup, *h_measure_all_jesup, *h_fake_all_jesup, *h_miss_all_jesup;
  build_hist("all_jesup", h_truth_all_jesup, h_measure_all_jesup, h_respmatrix_all_jesup, h_fake_all_jesup, h_miss_all_jesup);
  TH2D *h_respmatrix_all_jesdown; TH1D *h_truth_all_jesdown, *h_measure_all_jesdown, *h_fake_all_jesdown, *h_miss_all_jesdown;
  build_hist("all_jesdown", h_truth_all_jesdown, h_measure_all_jesdown, h_respmatrix_all_jesdown, h_fake_all_jesdown, h_miss_all_jesdown);
  TH2D *h_respmatrix_zvertex60_jesup; TH1D *h_truth_zvertex60_jesup, *h_measure_zvertex60_jesup, *h_fake_zvertex60_jesup, *h_miss_zvertex60_jesup;
  build_hist("zvertex60_jesup", h_truth_zvertex60_jesup, h_measure_zvertex60_jesup, h_respmatrix_zvertex60_jesup, h_fake_zvertex60_jesup, h_miss_zvertex60_jesup);
  TH2D *h_respmatrix_zvertex60_jesdown; TH1D *h_truth_zvertex60_jesdown, *h_measure_zvertex60_jesdown, *h_fake_zvertex60_jesdown, *h_miss_zvertex60_jesdown;
  build_hist("zvertex60_jesdown", h_truth_zvertex60_jesdown, h_measure_zvertex60_jesdown, h_respmatrix_zvertex60_jesdown, h_fake_zvertex60_jesdown, h_miss_zvertex60_jesdown);

  // JER up/down histograms
  TH2D *h_respmatrix_all_jerup; TH1D *h_truth_all_jerup, *h_measure_all_jerup, *h_fake_all_jerup, *h_miss_all_jerup;
  build_hist("all_jerup", h_truth_all_jerup, h_measure_all_jerup, h_respmatrix_all_jerup, h_fake_all_jerup, h_miss_all_jerup);
  TH2D *h_respmatrix_all_jerdown; TH1D *h_truth_all_jerdown, *h_measure_all_jerdown, *h_fake_all_jerdown, *h_miss_all_jerdown;
  build_hist("all_jerdown", h_truth_all_jerdown, h_measure_all_jerdown, h_respmatrix_all_jerdown, h_fake_all_jerdown, h_miss_all_jerdown);
  TH2D *h_respmatrix_zvertex60_jerup; TH1D *h_truth_zvertex60_jerup, *h_measure_zvertex60_jerup, *h_fake_zvertex60_jerup, *h_miss_zvertex60_jerup;
  build_hist("zvertex60_jerup", h_truth_zvertex60_jerup, h_measure_zvertex60_jerup, h_respmatrix_zvertex60_jerup, h_fake_zvertex60_jerup, h_miss_zvertex60_jerup);
  TH2D *h_respmatrix_zvertex60_jerdown; TH1D *h_truth_zvertex60_jerdown, *h_measure_zvertex60_jerdown, *h_fake_zvertex60_jerdown, *h_miss_zvertex60_jerdown;
  build_hist("zvertex60_jerdown", h_truth_zvertex60_jerdown, h_measure_zvertex60_jerdown, h_respmatrix_zvertex60_jerdown, h_fake_zvertex60_jerdown, h_miss_zvertex60_jerdown);

  // Jet trigger up/down histograms
  TH2D *h_respmatrix_all_jetup; TH1D *h_truth_all_jetup, *h_measure_all_jetup, *h_fake_all_jetup, *h_miss_all_jetup;
  build_hist("all_jetup", h_truth_all_jetup, h_measure_all_jetup, h_respmatrix_all_jetup, h_fake_all_jetup, h_miss_all_jetup);
  TH2D *h_respmatrix_all_jetdown; TH1D *h_truth_all_jetdown, *h_measure_all_jetdown, *h_fake_all_jetdown, *h_miss_all_jetdown;
  build_hist("all_jetdown", h_truth_all_jetdown, h_measure_all_jetdown, h_respmatrix_all_jetdown, h_fake_all_jetdown, h_miss_all_jetdown);
  TH2D *h_respmatrix_zvertex60_jetup; TH1D *h_truth_zvertex60_jetup, *h_measure_zvertex60_jetup, *h_fake_zvertex60_jetup, *h_miss_zvertex60_jetup;
  build_hist("zvertex60_jetup", h_truth_zvertex60_jetup, h_measure_zvertex60_jetup, h_respmatrix_zvertex60_jetup, h_fake_zvertex60_jetup, h_miss_zvertex60_jetup);
  TH2D *h_respmatrix_zvertex60_jetdown; TH1D *h_truth_zvertex60_jetdown, *h_measure_zvertex60_jetdown, *h_fake_zvertex60_jetdown, *h_miss_zvertex60_jetdown;
  build_hist("zvertex60_jetdown", h_truth_zvertex60_jetdown, h_measure_zvertex60_jetdown, h_respmatrix_zvertex60_jetdown, h_fake_zvertex60_jetdown, h_miss_zvertex60_jetdown);

  // Closure test histograms
  TH1D *h_fullclosure_measure_zvertex60 = new TH1D("h_fullclosure_measure_zvertex60", ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  TH1D *h_fullclosure_truth_zvertex60 = new TH1D("h_fullclosure_truth_zvertex60", ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  TH2D *h_fullclosure_respmatrix_zvertex60 = new TH2D("h_fullclosure_respmatrix_zvertex60", ";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", calibnpt, calibptbins, truthnpt, truthptbins);
  TH1D *h_fullclosure_fake_zvertex60 = new TH1D("h_fullclosure_fake_zvertex60", ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  TH1D *h_fullclosure_miss_zvertex60 = new TH1D("h_fullclosure_miss_zvertex60", ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);

  TH1D *h_halfclosure_inputmeasure_zvertex60 = new TH1D("h_halfclosure_inputmeasure_zvertex60", ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  TH1D *h_halfclosure_measure_zvertex60 = new TH1D("h_halfclosure_measure_zvertex60", ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  TH1D *h_halfclosure_truth_zvertex60 = new TH1D("h_halfclosure_truth_zvertex60", ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  TH2D *h_halfclosure_respmatrix_zvertex60 = new TH2D("h_halfclosure_respmatrix_zvertex60", ";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", calibnpt, calibptbins, truthnpt, truthptbins);
  TH1D *h_halfclosure_fake_zvertex60 = new TH1D("h_halfclosure_fake_zvertex60", ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  TH1D *h_halfclosure_miss_zvertex60 = new TH1D("h_halfclosure_miss_zvertex60", ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);

  ////////// Event Loop //////////
  std::cout << "Data analysis started." << std::endl;
  int n_events = chain.GetEntries();
  //n_events = 10;
  std::cout << "Total number of events: " << n_events << std::endl;
  // Event variables setup.
  int n_events_all = 0;
  bool unsubjet_background_nojet, unsubjet_background_dijet, unsubjet_background_njet;
  bool is_mbd_zvertex60, is_truth_zvertex60;
  std::vector<bool> jet_filter, truthjet_filter;
  std::vector<float> smear_value;
  std::vector<float> goodtruthjet_pt, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched;
  std::vector<float> calibjet_pt, calibjet_eta, calibjet_phi, calibjet_matched;
  std::vector<float> calibjet_pt_jesup, calibjet_eta_jesup, calibjet_phi_jesup, calibjet_matched_jesup;
  std::vector<float> calibjet_pt_jesdown, calibjet_eta_jesdown, calibjet_phi_jesdown, calibjet_matched_jesdown;
  std::vector<float> calibjet_pt_jerup, calibjet_eta_jerup, calibjet_phi_jerup, calibjet_matched_jerup;
  std::vector<float> calibjet_pt_jerdown, calibjet_eta_jerdown, calibjet_phi_jerdown, calibjet_matched_jerdown;
  for (int ie = 0; ie < n_events; ++ie) { // event loop start
    // Load event.
    if (ie % 1000 == 0) {
      std::cout << "Processing event " << ie << "..." << std::endl;
    }
    chain.GetEntry(ie);
    bsgen->Generate(iseg, ie);

    n_events_all++;
    // Get Z-vertex check.
    is_mbd_zvertex60 = (zvertex > -60 && zvertex < 60);
    is_truth_zvertex60 = (truth_zvertex > -60 && truth_zvertex < 60);
    double scale_zvertexreweight = h_zvertex_reweight->GetBinContent(h_zvertex_reweight->GetXaxis()->FindBin(truth_zvertex));

    bool has_zvertex = true;
    if (zvertex < -990) {
      has_zvertex = false;
      zvertex = 0;
    }

    // Do data set efficiency check.
    int leadingtruthjet_index = -1;
    float leadingtruthjet_pt = -9999;
    for (int ij = 0; ij < truthjet_pt->size(); ++ij) {
      if (truthjet_eta->at(ij) > 1.5 - jet_radius || truthjet_eta->at(ij) < -1.5 + jet_radius) continue;
      float jetpt = truthjet_pt->at(ij);
      if (jetpt > leadingtruthjet_pt) {
        leadingtruthjet_pt = jetpt;
        leadingtruthjet_index = ij;
      }
    }
    if (leadingtruthjet_index < 0) continue;
    if (leadingtruthjet_pt < truthjet_pt_min || leadingtruthjet_pt >= truthjet_pt_max) continue;

    // Do reco jet trimming.
    unsubjet_background_nojet = true;
    int leadingcalibjet_index = -1;
    float leadingcalibjet_pt = -9999;
    for (int ij = 0; ij < unsubjet_pt_calib->size(); ++ij) {
      float jetpt = unsubjet_pt_calib->at(ij);
      if (jetpt > leadingcalibjet_pt) {
        leadingcalibjet_pt = jetpt;
        leadingcalibjet_index = ij;
      }
    }
    if (leadingcalibjet_index >= 0) {
      unsubjet_background_nojet = false;
      if (leadingcalibjet_pt > calibjet_pt_max) continue;
    }

    for (int ij = 0; ij < truthjet_pt->size(); ++ij) {
      if (truthjet_eta->at(ij) > 1.1 - jet_radius || truthjet_eta->at(ij) < -1.1 + jet_radius) continue;
      h_truthjet_pt_record_all->Fill(truthjet_pt->at(ij));
    }

    // Get reco leading jet and subleading jet. Do basic jet cuts.
    unsubjet_background_dijet = true;
    unsubjet_background_njet = true;
    if (!unsubjet_background_nojet) {
      int leadingunsubjet04_index = -1;
      int subleadingunsubjet04_index = -1;
      get_leading_subleading_jet(leadingunsubjet04_index, subleadingunsubjet04_index, unsubjet04_e);
      if (leadingunsubjet04_index < 0) continue;
      // Dijet cut.
      bool match_dijet = false;
      if (subleadingunsubjet04_index >= 0 && unsubjet04_e->at(subleadingunsubjet04_index) / (float) unsubjet04_e->at(leadingunsubjet04_index) > 0.3) {
        match_dijet = match_leading_subleading_jet(unsubjet04_phi->at(leadingunsubjet04_index), unsubjet04_phi->at(subleadingunsubjet04_index));
      }
      if (match_dijet) unsubjet_background_dijet = false;
      // Njet cut.
      int n_jets = 0;
      for (int ij = 0; ij < unsubjet04_e->size(); ++ij) {
        if (unsubjet04_e->at(ij) > 5.0) {
          n_jets++;
        }
      }
      if (n_jets <= 9) unsubjet_background_njet = false;
    }

    // Do jet filter
    get_recojet_filter(jet_filter, unsubjet_e, unsubjet_pt, unsubjet_eta, zvertex, jet_radius);
    get_truthjet_filter(truthjet_filter, truthjet_e, truthjet_pt, truthjet_eta, zvertex, jet_radius);
    for (int ij = 0; ij < unsubjet_pt->size(); ++ij) {
      if (jet_filter[ij]) continue;

      h_recojet_pt_record_nocut_all->Fill(unsubjet_pt->at(ij));
      if (!unsubjet_background_dijet && !unsubjet_background_njet) h_recojet_pt_record_all->Fill(unsubjet_pt->at(ij));

      if (is_truth_zvertex60) {
        h_recojet_pt_record_nocut_zvertex60->Fill(unsubjet_pt->at(ij));
        if (!unsubjet_background_dijet && !unsubjet_background_njet) h_recojet_pt_record_zvertex60->Fill(unsubjet_pt->at(ij));
      }
    }

    // Get truth jet and calib jet.
    get_truthjet(goodtruthjet_pt, goodtruthjet_eta, goodtruthjet_phi, truthjet_filter, truthjet_pt, truthjet_eta, truthjet_phi);
    get_smear_value(smear_value, jet_filter);
    get_calibjet(calibjet_pt, calibjet_eta, calibjet_phi, jet_filter, unsubjet_pt_calib, unsubjet_eta, unsubjet_phi, jet_radius, 1, 0.1, smear_value, (unsubjet_background_dijet||unsubjet_background_njet));
    get_calibjet(calibjet_pt_jesup, calibjet_eta_jesup, calibjet_phi_jesup, jet_filter, unsubjet_pt_calib, unsubjet_eta, unsubjet_phi, jet_radius, 1.025, 0.1, smear_value, (unsubjet_background_dijet||unsubjet_background_njet));
    get_calibjet(calibjet_pt_jesdown, calibjet_eta_jesdown, calibjet_phi_jesdown, jet_filter, unsubjet_pt_calib, unsubjet_eta, unsubjet_phi, jet_radius, 0.975, 0.1, smear_value, (unsubjet_background_dijet||unsubjet_background_njet));
    get_calibjet(calibjet_pt_jerup, calibjet_eta_jerup, calibjet_phi_jerup, jet_filter, unsubjet_pt_calib, unsubjet_eta, unsubjet_phi, jet_radius, 1, 0.12, smear_value, (unsubjet_background_dijet||unsubjet_background_njet));
    get_calibjet(calibjet_pt_jerdown, calibjet_eta_jerdown, calibjet_phi_jerdown, jet_filter, unsubjet_pt_calib, unsubjet_eta, unsubjet_phi, jet_radius, 1, 0.08, smear_value, (unsubjet_background_dijet||unsubjet_background_njet));
    // Match truth jet and calib jet.
    match_meas_truth(calibjet_eta, calibjet_phi, calibjet_matched, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched, jet_radius, has_zvertex);
    fill_response_matrix(h_truth_all, h_measure_all, h_respmatrix_all, h_fake_all, h_miss_all,
                         scale_zvertexreweight, f_reweightfunc_all,
                         calibjet_pt, calibjet_matched,
                         goodtruthjet_pt, goodtruthjet_matched);
    fill_response_matrix_stat(h_truth_all_stat, h_measure_all_stat, h_respmatrix_all_stat, h_fake_all_stat, h_miss_all_stat,
                              scale_zvertexreweight, f_reweightfunc_all,
                              calibjet_pt, calibjet_matched,
                              goodtruthjet_pt, goodtruthjet_matched);
    fill_response_matrix(h_truth_all_unreweighted, h_measure_all_unreweighted, h_respmatrix_all_unreweighted, h_fake_all_unreweighted, h_miss_all_unreweighted,
                         scale_zvertexreweight, f_reweightfunc_unreweighted,
                         calibjet_pt, calibjet_matched,
                         goodtruthjet_pt, goodtruthjet_matched);
    fill_response_matrix(h_truth_all_jetup, h_measure_all_jetup, h_respmatrix_all_jetup, h_fake_all_jetup, h_miss_all_jetup,
                         scale_zvertexreweight, f_reweightfunc_all_jetup,
                         calibjet_pt, calibjet_matched,
                         goodtruthjet_pt, goodtruthjet_matched);
    fill_response_matrix(h_truth_all_jetdown, h_measure_all_jetdown, h_respmatrix_all_jetdown, h_fake_all_jetdown, h_miss_all_jetdown,
                         scale_zvertexreweight, f_reweightfunc_all_jetdown,
                         calibjet_pt, calibjet_matched,
                         goodtruthjet_pt, goodtruthjet_matched);
    if (is_truth_zvertex60 || is_mbd_zvertex60) {
      fill_response_matrix_altz(h_truth_zvertex60, h_measure_zvertex60, h_respmatrix_zvertex60, h_fake_zvertex60, h_miss_zvertex60,
                                scale_zvertexreweight, f_reweightfunc_zvertex60,
                                calibjet_pt, calibjet_matched,
                                goodtruthjet_pt, goodtruthjet_matched, is_mbd_zvertex60, is_truth_zvertex60, has_zvertex);
      fill_response_matrix_altz_stat(h_truth_zvertex60_stat, h_measure_zvertex60_stat, h_respmatrix_zvertex60_stat, h_fake_zvertex60_stat, h_miss_zvertex60_stat,
                                     scale_zvertexreweight, f_reweightfunc_zvertex60,
                                     calibjet_pt, calibjet_matched,
                                     goodtruthjet_pt, goodtruthjet_matched, is_mbd_zvertex60, is_truth_zvertex60, has_zvertex);
      fill_response_matrix_altz(h_truth_zvertex60_unreweighted, h_measure_zvertex60_unreweighted, h_respmatrix_zvertex60_unreweighted, h_fake_zvertex60_unreweighted, h_miss_zvertex60_unreweighted,
                                scale_zvertexreweight, f_reweightfunc_unreweighted,
                                calibjet_pt, calibjet_matched,
                                goodtruthjet_pt, goodtruthjet_matched, is_mbd_zvertex60, is_truth_zvertex60, has_zvertex);
      fill_response_matrix_altz(h_truth_zvertex60_jetup, h_measure_zvertex60_jetup, h_respmatrix_zvertex60_jetup, h_fake_zvertex60_jetup, h_miss_zvertex60_jetup,
                                scale_zvertexreweight, f_reweightfunc_zvertex60_jetup,
                                calibjet_pt, calibjet_matched,
                                goodtruthjet_pt, goodtruthjet_matched, is_mbd_zvertex60, is_truth_zvertex60, has_zvertex);
      fill_response_matrix_altz(h_truth_zvertex60_jetdown, h_measure_zvertex60_jetdown, h_respmatrix_zvertex60_jetdown, h_fake_zvertex60_jetdown, h_miss_zvertex60_jetdown,
                                scale_zvertexreweight, f_reweightfunc_zvertex60_jetdown,
                                calibjet_pt, calibjet_matched,
                                goodtruthjet_pt, goodtruthjet_matched, is_mbd_zvertex60, is_truth_zvertex60, has_zvertex);
    }
    //if (is_zvertex60) {
    //  for (int im = 0; im < calibjet_pt.size(); ++im) {
    //    h_fullclosure_measure_zvertex60->Fill(calibjet_pt[im]);
    //    if (ie % 2 == 0) h_halfclosure_measure_zvertex60->Fill(calibjet_pt[im]);
    //    else h_halfclosure_inputmeasure_zvertex60->Fill(calibjet_pt[im]);
    //    if (calibjet_matched[im] < 0) {
    //      h_fullclosure_fake_zvertex60->Fill(calibjet_pt[im]);
    //      if (ie % 2 == 0) h_halfclosure_fake_zvertex60->Fill(calibjet_pt[im]);
    //    } else {
    //      h_fullclosure_respmatrix_zvertex60->Fill(calibjet_pt[im], goodtruthjet_pt[calibjet_matched[im]]);
    //      if (ie % 2 == 0) h_halfclosure_respmatrix_zvertex60->Fill(calibjet_pt[im], goodtruthjet_pt[calibjet_matched[im]]);
    //    }
    //  }
    //  for (int it = 0; it < goodtruthjet_pt.size(); ++it) {
    //    h_fullclosure_truth_zvertex60->Fill(goodtruthjet_pt[it]);
    //    if (ie % 2 == 0) h_halfclosure_truth_zvertex60->Fill(goodtruthjet_pt[it]);
    //    if (goodtruthjet_matched[it] < 0) {
    //      h_fullclosure_miss_zvertex60->Fill(goodtruthjet_pt[it]);
    //      if (ie % 2 == 0) h_halfclosure_miss_zvertex60->Fill(goodtruthjet_pt[it]);
    //    }
    //  }
    //}

    match_meas_truth(calibjet_eta_jesup, calibjet_phi_jesup, calibjet_matched_jesup, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched, jet_radius, has_zvertex);
    fill_response_matrix(h_truth_all_jesup, h_measure_all_jesup, h_respmatrix_all_jesup, h_fake_all_jesup, h_miss_all_jesup,
                         scale_zvertexreweight, f_reweightfunc_all_jesup,
                         calibjet_pt_jesup, calibjet_matched_jesup,
                         goodtruthjet_pt, goodtruthjet_matched);
    if (is_truth_zvertex60 || is_mbd_zvertex60) fill_response_matrix_altz(h_truth_zvertex60_jesup, h_measure_zvertex60_jesup, h_respmatrix_zvertex60_jesup, h_fake_zvertex60_jesup, h_miss_zvertex60_jesup,
                                                                 scale_zvertexreweight, f_reweightfunc_zvertex60_jesup,
                                                                 calibjet_pt_jesup, calibjet_matched_jesup,
                                                                 goodtruthjet_pt, goodtruthjet_matched, is_mbd_zvertex60, is_truth_zvertex60, has_zvertex);
                                           
    match_meas_truth(calibjet_eta_jesdown, calibjet_phi_jesdown, calibjet_matched_jesdown, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched, jet_radius, has_zvertex);
    fill_response_matrix(h_truth_all_jesdown, h_measure_all_jesdown, h_respmatrix_all_jesdown, h_fake_all_jesdown, h_miss_all_jesdown,
                         scale_zvertexreweight, f_reweightfunc_all_jesdown,
                         calibjet_pt_jesdown, calibjet_matched_jesdown,
                         goodtruthjet_pt, goodtruthjet_matched);
    if (is_truth_zvertex60 || is_mbd_zvertex60) fill_response_matrix_altz(h_truth_zvertex60_jesdown, h_measure_zvertex60_jesdown, h_respmatrix_zvertex60_jesdown, h_fake_zvertex60_jesdown, h_miss_zvertex60_jesdown,
                                                                 scale_zvertexreweight, f_reweightfunc_zvertex60_jesdown,
                                                                 calibjet_pt_jesdown, calibjet_matched_jesdown,
                                                                 goodtruthjet_pt, goodtruthjet_matched, is_mbd_zvertex60, is_truth_zvertex60, has_zvertex);

    match_meas_truth(calibjet_eta_jerup, calibjet_phi_jerup, calibjet_matched_jerup, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched, jet_radius, has_zvertex);
    fill_response_matrix(h_truth_all_jerup, h_measure_all_jerup, h_respmatrix_all_jerup, h_fake_all_jerup, h_miss_all_jerup,
                         scale_zvertexreweight, f_reweightfunc_all_jerup,
                         calibjet_pt_jerup, calibjet_matched_jerup,
                         goodtruthjet_pt, goodtruthjet_matched);
    if (is_truth_zvertex60 || is_mbd_zvertex60) fill_response_matrix_altz(h_truth_zvertex60_jerup, h_measure_zvertex60_jerup, h_respmatrix_zvertex60_jerup, h_fake_zvertex60_jerup, h_miss_zvertex60_jerup,
                                                                 scale_zvertexreweight, f_reweightfunc_zvertex60_jerup,
                                                                 calibjet_pt_jerup, calibjet_matched_jerup,
                                                                 goodtruthjet_pt, goodtruthjet_matched, is_mbd_zvertex60, is_truth_zvertex60, has_zvertex);

    match_meas_truth(calibjet_eta_jerdown, calibjet_phi_jerdown, calibjet_matched_jerdown, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched, jet_radius, has_zvertex);
    fill_response_matrix(h_truth_all_jerdown, h_measure_all_jerdown, h_respmatrix_all_jerdown, h_fake_all_jerdown, h_miss_all_jerdown,
                         scale_zvertexreweight, f_reweightfunc_all_jerdown,
                         calibjet_pt_jerdown, calibjet_matched_jerdown,
                         goodtruthjet_pt, goodtruthjet_matched);
    if (is_truth_zvertex60 || is_mbd_zvertex60) fill_response_matrix_altz(h_truth_zvertex60_jerdown, h_measure_zvertex60_jerdown, h_respmatrix_zvertex60_jerdown, h_fake_zvertex60_jerdown, h_miss_zvertex60_jerdown,
                                                                 scale_zvertexreweight, f_reweightfunc_zvertex60_jerdown,
                                                                 calibjet_pt_jerdown, calibjet_matched_jerdown,
                                                                 goodtruthjet_pt, goodtruthjet_matched, is_mbd_zvertex60, is_truth_zvertex60, has_zvertex);
  } // event loop end
  // Fill event histograms.
  h_event_all->Fill(0.5, n_events_all);

  // Write histograms.
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();

  h_event_all->Write();
  h_event_beforecut->Write();
  h_event_passed->Write();
  h_recojet_pt_record_nocut_all->Write();
  h_recojet_pt_record_all->Write();
  h_recojet_pt_record_nocut_zvertex60->Write();
  h_recojet_pt_record_zvertex60->Write();
  h_truthjet_pt_record_all->Write();

  // Nominal histograms
  h_truth_all->Write(); h_measure_all->Write(); h_respmatrix_all->Write(); h_fake_all->Write(); h_miss_all->Write();
  h_truth_zvertex60->Write(); h_measure_zvertex60->Write(); h_respmatrix_zvertex60->Write(); h_fake_zvertex60->Write(); h_miss_zvertex60->Write();
  h_truth_all_stat->Write(); h_measure_all_stat->Write(); h_respmatrix_all_stat->Write(); h_fake_all_stat->Write(); h_miss_all_stat->Write();
  h_truth_zvertex60_stat->Write(); h_measure_zvertex60_stat->Write(); h_respmatrix_zvertex60_stat->Write(); h_fake_zvertex60_stat->Write(); h_miss_zvertex60_stat->Write();
  // Unreweighted histograms
  h_truth_all_unreweighted->Write(); h_measure_all_unreweighted->Write(); h_respmatrix_all_unreweighted->Write(); h_fake_all_unreweighted->Write(); h_miss_all_unreweighted->Write();
  h_truth_zvertex60_unreweighted->Write(); h_measure_zvertex60_unreweighted->Write(); h_respmatrix_zvertex60_unreweighted->Write(); h_fake_zvertex60_unreweighted->Write(); h_miss_zvertex60_unreweighted->Write();
  // JES up/down histograms
  h_truth_all_jesup->Write(); h_measure_all_jesup->Write(); h_respmatrix_all_jesup->Write(); h_fake_all_jesup->Write(); h_miss_all_jesup->Write();
  h_truth_all_jesdown->Write(); h_measure_all_jesdown->Write(); h_respmatrix_all_jesdown->Write(); h_fake_all_jesdown->Write(); h_miss_all_jesdown->Write();
  h_truth_zvertex60_jesup->Write(); h_measure_zvertex60_jesup->Write(); h_respmatrix_zvertex60_jesup->Write(); h_fake_zvertex60_jesup->Write(); h_miss_zvertex60_jesup->Write();
  h_truth_zvertex60_jesdown->Write(); h_measure_zvertex60_jesdown->Write(); h_respmatrix_zvertex60_jesdown->Write(); h_fake_zvertex60_jesdown->Write(); h_miss_zvertex60_jesdown->Write();
  // JER up/down histograms
  h_truth_all_jerup->Write(); h_measure_all_jerup->Write(); h_respmatrix_all_jerup->Write(); h_fake_all_jerup->Write(); h_miss_all_jerup->Write();
  h_truth_all_jerdown->Write(); h_measure_all_jerdown->Write(); h_respmatrix_all_jerdown->Write(); h_fake_all_jerdown->Write(); h_miss_all_jerdown->Write();
  h_truth_zvertex60_jerup->Write(); h_measure_zvertex60_jerup->Write(); h_respmatrix_zvertex60_jerup->Write(); h_fake_zvertex60_jerup->Write(); h_miss_zvertex60_jerup->Write();
  h_truth_zvertex60_jerdown->Write(); h_measure_zvertex60_jerdown->Write(); h_respmatrix_zvertex60_jerdown->Write(); h_fake_zvertex60_jerdown->Write(); h_miss_zvertex60_jerdown->Write();
  // Jet trigger up/down histograms
  h_truth_all_jetup->Write(); h_measure_all_jetup->Write(); h_respmatrix_all_jetup->Write(); h_fake_all_jetup->Write(); h_miss_all_jetup->Write();
  h_truth_all_jetdown->Write(); h_measure_all_jetdown->Write(); h_respmatrix_all_jetdown->Write(); h_fake_all_jetdown->Write(); h_miss_all_jetdown->Write();
  h_truth_zvertex60_jetup->Write(); h_measure_zvertex60_jetup->Write(); h_respmatrix_zvertex60_jetup->Write(); h_fake_zvertex60_jetup->Write(); h_miss_zvertex60_jetup->Write();
  h_truth_zvertex60_jetdown->Write(); h_measure_zvertex60_jetdown->Write(); h_respmatrix_zvertex60_jetdown->Write(); h_fake_zvertex60_jetdown->Write(); h_miss_zvertex60_jetdown->Write();

  // Full closure test histograms
  h_fullclosure_truth_zvertex60->Write();
  h_fullclosure_measure_zvertex60->Write();
  h_fullclosure_respmatrix_zvertex60->Write();
  h_fullclosure_fake_zvertex60->Write();
  h_fullclosure_miss_zvertex60->Write();

  h_halfclosure_truth_zvertex60->Write();
  h_halfclosure_inputmeasure_zvertex60->Write();
  h_halfclosure_measure_zvertex60->Write();
  h_halfclosure_respmatrix_zvertex60->Write();
  h_halfclosure_fake_zvertex60->Write();
  h_halfclosure_miss_zvertex60->Write();

  f_out->Close();
  std::cout << "All done!" << std::endl;
}

////////////////////////////////////////// Helper Functions //////////////////////////////////////////
float get_deta(float eta1, float eta2) {
  return eta1 - eta2;
}

float get_dphi(float phi1, float phi2) {
  float dphi1 = phi1 - phi2;
  float dphi2 = phi1 - phi2 + 2*TMath::Pi();
  float dphi3 = phi1 - phi2 - 2*TMath::Pi();
  if (fabs(dphi1) > fabs(dphi2)) {
    dphi1 = dphi2;
  }
  if (fabs(dphi1) > fabs(dphi3)) {
    dphi1 = dphi3;
  }
  return dphi1;
}

float get_dR(float eta1, float phi1, float eta2, float phi2) {
  float deta = get_deta(eta1, eta2);
  float dphi = get_dphi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}

float get_dR_ellipsoid(float eta1, float phi1, float eta2, float phi2) {
  float deta = get_deta(eta1, eta2) / 3.;
  float dphi = get_dphi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}

float get_emcal_mineta_zcorrected(float zvertex) {
  float minz_EM = -130.23;
  float radius_EM = 93.5;
  float z = minz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_emcal_maxeta_zcorrected(float zvertex) {
  float maxz_EM = 130.23;
  float radius_EM = 93.5;
  float z = maxz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_ihcal_mineta_zcorrected(float zvertex) {
  float minz_IH = -170.299;
  float radius_IH = 127.503;
  float z = minz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ihcal_maxeta_zcorrected(float zvertex) {
  float maxz_IH = 170.299;
  float radius_IH = 127.503;
  float z = maxz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ohcal_mineta_zcorrected(float zvertex) {
  float minz_OH = -301.683;
  float radius_OH = 225.87;
  float z = minz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

float get_ohcal_maxeta_zcorrected(float zvertex) {
  float maxz_OH = 301.683;
  float radius_OH = 225.87;
  float z = maxz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

bool check_bad_jet_eta(float jet_eta, float zvertex, float jet_radius) {
  float emcal_mineta = get_emcal_mineta_zcorrected(zvertex);
  float emcal_maxeta = get_emcal_maxeta_zcorrected(zvertex);
  float ihcal_mineta = get_ihcal_mineta_zcorrected(zvertex);
  float ihcal_maxeta = get_ihcal_maxeta_zcorrected(zvertex);
  float ohcal_mineta = get_ohcal_mineta_zcorrected(zvertex);
  float ohcal_maxeta = get_ohcal_maxeta_zcorrected(zvertex);
  float minlimit = emcal_mineta;
  if (ihcal_mineta > minlimit) minlimit = ihcal_mineta;
  if (ohcal_mineta > minlimit) minlimit = ohcal_mineta;
  float maxlimit = emcal_maxeta;
  if (ihcal_maxeta < maxlimit) maxlimit = ihcal_maxeta;
  if (ohcal_maxeta < maxlimit) maxlimit = ohcal_maxeta;
  minlimit += jet_radius;
  maxlimit -= jet_radius;
  return jet_eta < minlimit || jet_eta > maxlimit;
}

////////////////////////////////////////// Functions //////////////////////////////////////////
void get_leading_subleading_jet(int& leadingjet_index, int& subleadingjet_index, std::vector<float>* jet_et) {
  leadingjet_index = -1;
  subleadingjet_index = -1;
  float leadingjet_et = -9999;
  float subleadingjet_et = -9999;
  for (int ij = 0; ij < jet_et->size(); ++ij) {
    float jetet = jet_et->at(ij);
    if (jetet > leadingjet_et) {
      subleadingjet_et = leadingjet_et;
      subleadingjet_index = leadingjet_index;
      leadingjet_et = jetet;
      leadingjet_index = ij;
    } else if (jetet > subleadingjet_et) {
      subleadingjet_et = jetet;
      subleadingjet_index = ij;
    }
  }
}

bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi) {
  float dijet_min_phi = 3*TMath::Pi()/4.;
  float dphi = fabs(get_dphi(leadingjet_phi, subleadingjet_phi));
  return dphi > dijet_min_phi;
}

void get_recojet_filter(std::vector<bool>& jet_filter, std::vector<float>* jet_e, std::vector<float>* jet_pt, std::vector<float>* jet_eta, float zvertex, float jet_radius) {
  jet_filter.clear();
  int njet = jet_e->size();
  for (int ij = 0; ij < njet; ++ij) {
    jet_filter.push_back(jet_e->at(ij) < 0 || check_bad_jet_eta(jet_eta->at(ij), zvertex, jet_radius) || jet_eta->at(ij) > 1.1 - jet_radius || jet_eta->at(ij) < -1.1 + jet_radius);
  }
}

void get_truthjet_filter(std::vector<bool>& jet_filter, std::vector<float>* jet_e, std::vector<float>* jet_pt, std::vector<float>* jet_eta, float zvertex, float jet_radius) {
  jet_filter.clear();
  int njet = jet_e->size();
  for (int ij = 0; ij < njet; ++ij) {
    jet_filter.push_back(jet_e->at(ij) < 0 || jet_eta->at(ij) > 1.1 - jet_radius || jet_eta->at(ij) < -1.1 + jet_radius);
  }
}

void get_smear_value(std::vector<float>& smear_value, std::vector<bool>& jet_filter) {
  smear_value.clear();
  for (int ij = 0; ij < jet_filter.size(); ++ij) {
    smear_value.push_back(randGen.Gaus(0.0, 1.0));
  }
}

void get_calibjet(std::vector<float>& calibjet_pt, std::vector<float>& calibjet_eta, std::vector<float>& calibjet_phi, std::vector<bool>& jet_filter, std::vector<float>* jet_pt_calib, std::vector<float>* jet_eta, std::vector<float>* jet_phi, double jet_radius, float jes_para, float jer_para, std::vector<float>& smear_value, bool jet_background) {
  calibjet_pt.clear();
  calibjet_eta.clear();
  calibjet_phi.clear();
  for (int ij = 0; ij < jet_filter.size(); ++ij) {
    if (jet_filter.at(ij) || jet_background) continue;
    double calib_pt = jet_pt_calib->at(ij) * (1 + smear_value.at(ij) * jer_para) * jes_para;
    if (calib_pt < calibptbins[0] || calib_pt > calibptbins[calibnpt]) continue;
    calibjet_pt.push_back(calib_pt);
    calibjet_eta.push_back(jet_eta->at(ij));
    calibjet_phi.push_back(jet_phi->at(ij));
  }
}

void get_truthjet(std::vector<float>& goodtruthjet_pt, std::vector<float>& goodtruthjet_eta, std::vector<float>& goodtruthjet_phi, std::vector<bool>& truthjet_filter, std::vector<float>* jet_pt, std::vector<float>* jet_eta, std::vector<float>* jet_phi) {
  goodtruthjet_pt.clear();
  goodtruthjet_eta.clear();
  goodtruthjet_phi.clear();
  for (int ij = 0; ij < truthjet_filter.size(); ++ij) {
    if (truthjet_filter.at(ij)) continue;
    if (jet_pt->at(ij) < truthptbins[0] || jet_pt->at(ij) > truthptbins[truthnpt]) continue;
    goodtruthjet_pt.push_back(jet_pt->at(ij));
    goodtruthjet_eta.push_back(jet_eta->at(ij));
    goodtruthjet_phi.push_back(jet_phi->at(ij));
  }
}

void match_meas_truth(std::vector<float>& meas_eta, std::vector<float>& meas_phi, std::vector<float>& meas_matched, std::vector<float>& truth_eta, std::vector<float>& truth_phi, std::vector<float>& truth_matched, float jet_radius, bool has_zvertex) {
  meas_matched.assign(meas_eta.size(), -1);
  truth_matched.assign(truth_eta.size(), -1);
  float max_match_dR = jet_radius * 0.75;
  
  // Build list of all valid pairs with their dR
  std::vector<std::tuple<float, int, int>> pairs;  // (dR, meas_idx, truth_idx)
  
  for (int im = 0; im < meas_eta.size(); ++im) {
    for (int it = 0; it < truth_eta.size(); ++it) {
      float dR;
      if (has_zvertex) dR = get_dR(meas_eta[im], meas_phi[im], truth_eta[it], truth_phi[it]);
      else dR = get_dR_ellipsoid(meas_eta[im], meas_phi[im], truth_eta[it], truth_phi[it]);
      
      if (dR < max_match_dR) {
        pairs.push_back(std::make_tuple(dR, im, it));
      }
    }
  }
  
  // Sort by dR (smallest first)
  std::sort(pairs.begin(), pairs.end());
  
  // Greedy matching: assign pairs in order of increasing dR
  for (auto& p : pairs) {
    int im = std::get<1>(p);
    int it = std::get<2>(p);
    
    // Only match if both are still unmatched
    if (meas_matched[im] == -1 && truth_matched[it] == -1) {
      meas_matched[im] = it;
      truth_matched[it] = im;
    }
  }
}

void fill_response_matrix(TH1D*& h_truth, TH1D*& h_meas, TH2D*& h_resp, TH1D*& h_fake, TH1D*& h_miss, double scale, TF1* f_reweight, std::vector<float>& meas_pt, std::vector<float>& meas_matched, std::vector<float>& truth_pt, std::vector<float>& truth_matched) {
  for (int it = 0; it < truth_pt.size(); ++it) {
    double reweight_factor = f_reweight->Eval(truth_pt[it]);
    if (truth_matched[it] < 0) {
      h_truth->Fill(truth_pt[it], scale*reweight_factor);
      h_miss->Fill(truth_pt[it], scale*reweight_factor);
    } else {
      h_truth->Fill(truth_pt[it], scale*reweight_factor);
      h_resp->Fill(meas_pt[truth_matched[it]], truth_pt[it], scale*reweight_factor);
    }
  }
  for (int im = 0; im < meas_pt.size(); ++im) {
    if (meas_matched[im] < 0) {
      h_meas->Fill(meas_pt[im], scale);
      h_fake->Fill(meas_pt[im], scale);
    } else {
      h_meas->Fill(meas_pt[im], scale);
    }
  }
}

void fill_response_matrix_altz(TH1D*& h_truth, TH1D*& h_meas, TH2D*& h_resp, TH1D*& h_fake, TH1D*& h_miss, double scale, TF1* f_reweight, std::vector<float>& meas_pt, std::vector<float>& meas_matched, std::vector<float>& truth_pt, std::vector<float>& truth_matched, bool is_mbd_zvertex60, bool is_truth_zvertex60, bool has_zvertex) {
  if (is_truth_zvertex60 && is_mbd_zvertex60) {
    for (int it = 0; it < truth_pt.size(); ++it) {
      double reweight_factor = f_reweight->Eval(truth_pt[it]);
      if (truth_matched[it] < 0) {
        h_truth->Fill(truth_pt[it], scale*reweight_factor);
        h_miss->Fill(truth_pt[it], scale*reweight_factor);
      } else {
        h_truth->Fill(truth_pt[it], scale*reweight_factor);
        h_resp->Fill(meas_pt[truth_matched[it]], truth_pt[it], scale*reweight_factor);
      }
    }
    for (int im = 0; im < meas_pt.size(); ++im) {
      if (meas_matched[im] < 0) {
        h_meas->Fill(meas_pt[im], scale);
        h_fake->Fill(meas_pt[im], scale);
      } else {
        h_meas->Fill(meas_pt[im], scale);
      }
    }
  } else if (is_truth_zvertex60 && !is_mbd_zvertex60 && has_zvertex) {
    for (int it = 0; it < truth_pt.size(); ++it) {
      double reweight_factor = f_reweight->Eval(truth_pt[it]);
      h_truth->Fill(truth_pt[it], scale*reweight_factor);
      h_miss->Fill(truth_pt[it], scale*reweight_factor);
    }
  } else if (is_truth_zvertex60 && !has_zvertex) {
    for (int it = 0; it < truth_pt.size(); ++it) {
      double reweight_factor = f_reweight->Eval(truth_pt[it]);
      h_truth->Fill(truth_pt[it], scale*reweight_factor);
      h_miss->Fill(truth_pt[it], scale*reweight_factor);
    }
  } else if (!is_truth_zvertex60 && is_mbd_zvertex60) {
    for (int im = 0; im < meas_pt.size(); ++im) {
      h_meas->Fill(meas_pt[im], scale);
      h_fake->Fill(meas_pt[im], scale);
    }
  }
}

void fill_response_matrix_stat(TH1DBootstrap*& h_truth, TH1DBootstrap*& h_meas, TH2DBootstrap*& h_resp, TH1DBootstrap*& h_fake, TH1DBootstrap*& h_miss, double scale, TF1* f_reweight, std::vector<float>& meas_pt, std::vector<float>& meas_matched, std::vector<float>& truth_pt, std::vector<float>& truth_matched) {
  for (int it = 0; it < truth_pt.size(); ++it) {
    double reweight_factor = f_reweight->Eval(truth_pt[it]);
    if (truth_matched[it] < 0) {
      h_truth->Fill(truth_pt[it], scale*reweight_factor);
      h_miss->Fill(truth_pt[it], scale*reweight_factor);
    } else {
      h_truth->Fill(truth_pt[it], scale*reweight_factor);
      h_resp->Fill(meas_pt[truth_matched[it]], truth_pt[it], scale*reweight_factor);
    }
  }
  for (int im = 0; im < meas_pt.size(); ++im) {
    if (meas_matched[im] < 0) {
      h_meas->Fill(meas_pt[im], scale);
      h_fake->Fill(meas_pt[im], scale);
    } else {
      h_meas->Fill(meas_pt[im], scale);
    }
  }
}

void fill_response_matrix_altz_stat(TH1DBootstrap*& h_truth, TH1DBootstrap*& h_meas, TH2DBootstrap*& h_resp, TH1DBootstrap*& h_fake, TH1DBootstrap*& h_miss, double scale, TF1* f_reweight, std::vector<float>& meas_pt, std::vector<float>& meas_matched, std::vector<float>& truth_pt, std::vector<float>& truth_matched, bool is_mbd_zvertex60, bool is_truth_zvertex60, bool has_zvertex) {
  if (is_truth_zvertex60 && is_mbd_zvertex60) {
    for (int it = 0; it < truth_pt.size(); ++it) {
      double reweight_factor = f_reweight->Eval(truth_pt[it]);
      if (truth_matched[it] < 0) {
        h_truth->Fill(truth_pt[it], scale*reweight_factor);
        h_miss->Fill(truth_pt[it], scale*reweight_factor);
      } else {
        h_truth->Fill(truth_pt[it], scale*reweight_factor);
        h_resp->Fill(meas_pt[truth_matched[it]], truth_pt[it], scale*reweight_factor);
      }
    }
    for (int im = 0; im < meas_pt.size(); ++im) {
      if (meas_matched[im] < 0) {
        h_meas->Fill(meas_pt[im], scale);
        h_fake->Fill(meas_pt[im], scale);
      } else {
        h_meas->Fill(meas_pt[im], scale);
      }
    }
  } else if (is_truth_zvertex60 && !is_mbd_zvertex60 && has_zvertex) {
    for (int it = 0; it < truth_pt.size(); ++it) {
      double reweight_factor = f_reweight->Eval(truth_pt[it]);
      h_truth->Fill(truth_pt[it], scale*reweight_factor);
      h_miss->Fill(truth_pt[it], scale*reweight_factor);
    }
  } else if (is_truth_zvertex60 && !has_zvertex) {
    for (int it = 0; it < truth_pt.size(); ++it) {
      double reweight_factor = f_reweight->Eval(truth_pt[it]);
      h_truth->Fill(truth_pt[it], scale*reweight_factor);
      h_miss->Fill(truth_pt[it], scale*reweight_factor);
    }
  } else if (!is_truth_zvertex60 && is_mbd_zvertex60) {
    for (int im = 0; im < meas_pt.size(); ++im) {
      h_meas->Fill(meas_pt[im], scale);
      h_fake->Fill(meas_pt[im], scale);
    }
  }
}

void build_hist(std::string tag, TH1D*& h_truth, TH1D*& h_measure, TH2D*& h_respmatrix, TH1D*& h_fake, TH1D*& h_miss) {
  h_truth = new TH1D(Form("h_truth_%s", tag.c_str()), ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  h_measure = new TH1D(Form("h_measure_%s", tag.c_str()), ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  h_respmatrix = new TH2D(Form("h_respmatrix_%s", tag.c_str()), ";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", calibnpt, calibptbins, truthnpt, truthptbins);
  h_fake = new TH1D(Form("h_fake_%s", tag.c_str()), ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  h_miss = new TH1D(Form("h_miss_%s", tag.c_str()), ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
}

void build_hist_stat(std::string tag, TH1DBootstrap*& h_truth, TH1DBootstrap*& h_measure, TH2DBootstrap*& h_respmatrix, TH1DBootstrap*& h_fake, TH1DBootstrap*& h_miss, int nrep, BootstrapGenerator* bsgen) {
  h_truth = new TH1DBootstrap(Form("h_truth_%s", tag.c_str()), ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins, nrep, bsgen);
  h_measure = new TH1DBootstrap(Form("h_measure_%s", tag.c_str()), ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins, nrep, bsgen);
  h_respmatrix = new TH2DBootstrap(Form("h_respmatrix_%s", tag.c_str()), ";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", calibnpt, calibptbins, truthnpt, truthptbins, nrep, bsgen);
  h_fake = new TH1DBootstrap(Form("h_fake_%s", tag.c_str()), ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins, nrep, bsgen);
  h_miss = new TH1DBootstrap(Form("h_miss_%s", tag.c_str()), ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins, nrep, bsgen);
}
