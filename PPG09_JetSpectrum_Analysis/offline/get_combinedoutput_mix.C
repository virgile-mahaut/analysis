#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include "BootstrapGenerator/BootstrapGenerator.h"
#include "BootstrapGenerator/TH1DBootstrap.h"
#include "BootstrapGenerator/TH2DBootstrap.h"

R__LOAD_LIBRARY(libBootstrapGenerator.so)

void combine_truthjet(/*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out);
void combine_1Dhist(string histname, string surfix, std::string datatype, /*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out);
void combine_2Dhist(string histname, string surfix, /*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out);
void combine_1Dhist_closure(string surfix, string tag, /*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out);
void combine_2Dhist_closure(string surfix, string tag, /*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out);
void combine_1Dhist_stat(string histname, string surfix, std::string datatype, /*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out);
void combine_2Dhist_stat(string histname, string surfix, /*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out);

const double mb_cross_section = 4.1970e+10;
const double jet5_cross_section = 1.3878e+8;
const double jet12_cross_section = 1.4903e+6; 
const double jet20_cross_section = 6.2623e+4; 
const double jet30_cross_section = 2.5298e+3; 
const double jet40_cross_section = 1.3553e+2;
const double jet50_cross_section = 7.3113;
const double jet60_cross_section = 3.3261e-01;

void get_combinedoutput_mix(int radius_index = 4) {

  //TFile* f_MB = new TFile(Form("output_sim_hadd/output_sim_r0%d_MB.root", radius_index), "READ");
  //TFile* f_Jet5GeV = new TFile(Form("output_sim_hadd/output_sim_r0%d_Jet5GeV.root", radius_index), "READ");
  TFile* f_Jet12GeV = new TFile(Form("output_sim_hadd/output_sim_r0%d_Jet12GeV.root", radius_index), "READ");
  TFile* f_Jet20GeV = new TFile(Form("output_sim_hadd/output_sim_r0%d_Jet20GeV.root", radius_index), "READ");
  TFile* f_Jet30GeV = new TFile(Form("output_sim_hadd/output_sim_r0%d_Jet30GeV.root", radius_index), "READ");
  TFile* f_Jet40GeV = new TFile(Form("output_sim_hadd/output_sim_r0%d_Jet40GeV.root", radius_index), "READ");
  TFile* f_Jet50GeV = new TFile(Form("output_sim_hadd/output_sim_r0%d_Jet50GeV.root", radius_index), "READ");
  //TFile* f_Jet60GeV = new TFile(Form("output_sim_hadd/output_sim_r0%d_Jet60GeV.root", radius_index), "READ");
  if (/*!f_MB || !f_Jet5GeV ||*/ !f_Jet12GeV || !f_Jet20GeV || !f_Jet30GeV || !f_Jet40GeV || !f_Jet50GeV /*|| !f_Jet60GeV*/) {
    std::cout << "Error: cannot open one or more input files." << std::endl;
    return;
  }
  TFile* f_combined = new TFile(Form("output_sim_r0%d.root", radius_index), "RECREATE");

  combine_truthjet(/*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);

  combine_1Dhist("h_recojet_pt_record_nocut", "", "reco", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_recojet_pt_record", "", "reco", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);

  // Combine nominal histograms
  combine_1Dhist("h_truth", "", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_measure", "", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_2Dhist("h_respmatrix", "", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_fake", "", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_miss", "", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);

  combine_1Dhist_stat("h_truth", "_stat", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist_stat("h_measure", "_stat", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_2Dhist_stat("h_respmatrix", "_stat", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist_stat("h_fake", "_stat", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist_stat("h_miss", "_stat", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);

  // Combine unreweighted histograms
  combine_1Dhist("h_truth", "_unreweighted", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_measure", "_unreweighted", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_2Dhist("h_respmatrix", "_unreweighted", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_fake", "_unreweighted", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_miss", "_unreweighted", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);

  // Combine JES up variation histograms
  combine_1Dhist("h_truth", "_jesup", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_measure", "_jesup", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_2Dhist("h_respmatrix", "_jesup", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_fake", "_jesup", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_miss", "_jesup", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);

  // Combine JES down variation histograms
  combine_1Dhist("h_truth", "_jesdown", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_measure", "_jesdown", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_2Dhist("h_respmatrix", "_jesdown", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_fake", "_jesdown", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_miss", "_jesdown", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);

  // Combine JER up variation histograms
  combine_1Dhist("h_truth", "_jerup", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_measure", "_jerup", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_2Dhist("h_respmatrix", "_jerup", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_fake", "_jerup", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_miss", "_jerup", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);

  // Combine JER down variation histograms
  combine_1Dhist("h_truth", "_jerdown", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_measure", "_jerdown", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_2Dhist("h_respmatrix", "_jerdown", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_fake", "_jerdown", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_miss", "_jerdown", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);

  // Combine jet trigger up variation histograms
  combine_1Dhist("h_truth", "_jetup", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_measure", "_jetup", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_2Dhist("h_respmatrix", "_jetup", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_fake", "_jetup", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_miss", "_jetup", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);

  // Combine jet trigger down variation histograms
  combine_1Dhist("h_truth", "_jetdown", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_measure", "_jetdown", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_2Dhist("h_respmatrix", "_jetdown", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_fake", "_jetdown", "meas", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist("h_miss", "_jetdown", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);

  // full closure test
  combine_1Dhist_closure("h_fullclosure_", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist_closure("h_fullclosure_", "measure", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_2Dhist_closure("h_fullclosure_", "respmatrix", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist_closure("h_fullclosure_", "fake", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist_closure("h_fullclosure_", "miss", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);

  combine_1Dhist_closure("h_halfclosure_", "inputmeasure", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist_closure("h_halfclosure_", "truth", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist_closure("h_halfclosure_", "measure", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_2Dhist_closure("h_halfclosure_", "respmatrix", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist_closure("h_halfclosure_", "fake", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
  combine_1Dhist_closure("h_halfclosure_", "miss", /*f_MB, f_Jet5GeV, */f_Jet12GeV, f_Jet20GeV, f_Jet30GeV, f_Jet40GeV, f_Jet50GeV, /*f_Jet60GeV, */f_combined);
}

void combine_truthjet(/*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out) {
  //TH1D *h_MB_event_all = (TH1D*)f_MB->Get("h_event_all"); int mb_event_all = h_MB_event_all->GetBinContent(1); double mb_scale_all = mb_cross_section / (double)mb_event_all;
  //TH1D *h_Jet5GeV_event_all = (TH1D*)f_Jet5GeV->Get("h_event_all"); int jet5_event_all = h_Jet5GeV_event_all->GetBinContent(1); double jet5_scale_all = jet5_cross_section / (double)jet5_event_all;
  TH1D *h_Jet12GeV_event_all = (TH1D*)f_Jet12GeV->Get("h_event_all"); int jet12_event_all = h_Jet12GeV_event_all->GetBinContent(1); double jet12_scale_all = jet12_cross_section / (double)jet12_event_all;
  TH1D *h_Jet20GeV_event_all = (TH1D*)f_Jet20GeV->Get("h_event_all"); int jet20_event_all = h_Jet20GeV_event_all->GetBinContent(1); double jet20_scale_all = jet20_cross_section / (double)jet20_event_all;
  TH1D *h_Jet30GeV_event_all = (TH1D*)f_Jet30GeV->Get("h_event_all"); int jet30_event_all = h_Jet30GeV_event_all->GetBinContent(1); double jet30_scale_all = jet30_cross_section / (double)jet30_event_all;
  TH1D *h_Jet40GeV_event_all = (TH1D*)f_Jet40GeV->Get("h_event_all"); int jet40_event_all = h_Jet40GeV_event_all->GetBinContent(1); double jet40_scale_all = jet40_cross_section / (double)jet40_event_all;
  TH1D *h_Jet50GeV_event_all = (TH1D*)f_Jet50GeV->Get("h_event_all"); int jet50_event_all = h_Jet50GeV_event_all->GetBinContent(1); double jet50_scale_all = jet50_cross_section / (double)jet50_event_all;
  //TH1D *h_Jet60GeV_event_all = (TH1D*)f_Jet60GeV->Get("h_event_all"); int jet60_event_all = h_Jet60GeV_event_all->GetBinContent(1); double jet60_scale_all = jet60_cross_section / (double)jet60_event_all;

  string histname_all = "h_truthjet_pt_record_all";

  //TH1D* h_MB_all_forcombine = (TH1D*)f_MB->Get(histname_all.c_str());
  //TH1D* h_Jet5GeV_all_forcombine = (TH1D*)f_Jet5GeV->Get(histname_all.c_str());
  TH1D* h_Jet12GeV_all_forcombine = (TH1D*)f_Jet12GeV->Get(histname_all.c_str());
  TH1D* h_Jet20GeV_all_forcombine = (TH1D*)f_Jet20GeV->Get(histname_all.c_str());
  TH1D* h_Jet30GeV_all_forcombine = (TH1D*)f_Jet30GeV->Get(histname_all.c_str());
  TH1D* h_Jet40GeV_all_forcombine = (TH1D*)f_Jet40GeV->Get(histname_all.c_str());
  TH1D* h_Jet50GeV_all_forcombine = (TH1D*)f_Jet50GeV->Get(histname_all.c_str());
  //TH1D* h_Jet60GeV_all_forcombine = (TH1D*)f_Jet60GeV->Get(histname_all.c_str());

  //TH1D* h_all_combined = (TH1D*)h_MB_all_forcombine->Clone(histname_all.c_str());
  //h_all_combined->Scale(mb_scale_all);
  TH1D* h_all_combined = (TH1D*)h_Jet12GeV_all_forcombine->Clone(histname_all.c_str());
  h_all_combined->Scale(jet12_scale_all);
  //h_all_combined->Add(h_Jet5GeV_all_forcombine, jet5_scale_all);
  //h_all_combined->Add(h_Jet12GeV_all_forcombine, jet12_scale_all);
  h_all_combined->Add(h_Jet20GeV_all_forcombine, jet20_scale_all);
  h_all_combined->Add(h_Jet30GeV_all_forcombine, jet30_scale_all);
  h_all_combined->Add(h_Jet40GeV_all_forcombine, jet40_scale_all);
  h_all_combined->Add(h_Jet50GeV_all_forcombine, jet50_scale_all);
  //h_all_combined->Add(h_Jet60GeV_all_forcombine, jet60_scale_all);

  f_out->cd();
  h_all_combined->Write();
}

void combine_1Dhist(string histname, string surfix, std::string datatype, /*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out) {
  //TH1D *h_MB_event_all = (TH1D*)f_MB->Get("h_event_all"); int mb_event_all = h_MB_event_all->GetBinContent(1); double mb_scale_all = mb_cross_section / (double)mb_event_all;
  //TH1D *h_Jet5GeV_event_all = (TH1D*)f_Jet5GeV->Get("h_event_all"); int jet5_event_all = h_Jet5GeV_event_all->GetBinContent(1); double jet5_scale_all = jet5_cross_section / (double)jet5_event_all;
  TH1D *h_Jet12GeV_event_all = (TH1D*)f_Jet12GeV->Get("h_event_all"); int jet12_event_all = h_Jet12GeV_event_all->GetBinContent(1); double jet12_scale_all = jet12_cross_section / (double)jet12_event_all;
  TH1D *h_Jet20GeV_event_all = (TH1D*)f_Jet20GeV->Get("h_event_all"); int jet20_event_all = h_Jet20GeV_event_all->GetBinContent(1); double jet20_scale_all = jet20_cross_section / (double)jet20_event_all;
  TH1D *h_Jet30GeV_event_all = (TH1D*)f_Jet30GeV->Get("h_event_all"); int jet30_event_all = h_Jet30GeV_event_all->GetBinContent(1); double jet30_scale_all = jet30_cross_section / (double)jet30_event_all;
  TH1D *h_Jet40GeV_event_all = (TH1D*)f_Jet40GeV->Get("h_event_all"); int jet40_event_all = h_Jet40GeV_event_all->GetBinContent(1); double jet40_scale_all = jet40_cross_section / (double)jet40_event_all;
  TH1D *h_Jet50GeV_event_all = (TH1D*)f_Jet50GeV->Get("h_event_all"); int jet50_event_all = h_Jet50GeV_event_all->GetBinContent(1); double jet50_scale_all = jet50_cross_section / (double)jet50_event_all;
  //TH1D *h_Jet60GeV_event_all = (TH1D*)f_Jet60GeV->Get("h_event_all"); int jet60_event_all = h_Jet60GeV_event_all->GetBinContent(1); double jet60_scale_all = jet60_cross_section / (double)jet60_event_all;

  string histname_all = histname + "_all" + surfix;
  string histname_zvertex60 = histname + "_zvertex60" + surfix;

  //TH1D* h_MB_all_forcombine = (TH1D*)f_MB->Get(histname_all.c_str());
  //TH1D* h_Jet5GeV_all_forcombine = (TH1D*)f_Jet5GeV->Get(histname_all.c_str());
  TH1D* h_Jet12GeV_all_forcombine = (TH1D*)f_Jet12GeV->Get(histname_all.c_str());
  TH1D* h_Jet20GeV_all_forcombine = (TH1D*)f_Jet20GeV->Get(histname_all.c_str());
  TH1D* h_Jet30GeV_all_forcombine = (TH1D*)f_Jet30GeV->Get(histname_all.c_str());
  TH1D* h_Jet40GeV_all_forcombine = (TH1D*)f_Jet40GeV->Get(histname_all.c_str());
  TH1D* h_Jet50GeV_all_forcombine = (TH1D*)f_Jet50GeV->Get(histname_all.c_str());
  //TH1D* h_Jet60GeV_all_forcombine = (TH1D*)f_Jet60GeV->Get(histname_all.c_str());

  //TH1D* h_all_combined = (TH1D*)h_MB_all_forcombine->Clone(histname_all.c_str());
  //h_all_combined->Scale(mb_scale_all);
  TH1D* h_all_combined = (TH1D*)h_Jet12GeV_all_forcombine->Clone(histname_all.c_str());
  h_all_combined->Scale(jet12_scale_all);
  //h_all_combined->Add(h_Jet5GeV_all_forcombine, jet5_scale_all);
  //h_all_combined->Add(h_Jet12GeV_all_forcombine, jet12_scale_all);
  h_all_combined->Add(h_Jet20GeV_all_forcombine, jet20_scale_all);
  h_all_combined->Add(h_Jet30GeV_all_forcombine, jet30_scale_all);
  h_all_combined->Add(h_Jet40GeV_all_forcombine, jet40_scale_all);
  h_all_combined->Add(h_Jet50GeV_all_forcombine, jet50_scale_all);
  //h_all_combined->Add(h_Jet60GeV_all_forcombine, jet60_scale_all);

  //TH1D* h_MB_zvertex60_forcombine = (TH1D*)f_MB->Get(histname_zvertex60.c_str());
  //TH1D* h_Jet5GeV_zvertex60_forcombine = (TH1D*)f_Jet5GeV->Get(histname_zvertex60.c_str());
  TH1D* h_Jet12GeV_zvertex60_forcombine = (TH1D*)f_Jet12GeV->Get(histname_zvertex60.c_str());
  TH1D* h_Jet20GeV_zvertex60_forcombine = (TH1D*)f_Jet20GeV->Get(histname_zvertex60.c_str());
  TH1D* h_Jet30GeV_zvertex60_forcombine = (TH1D*)f_Jet30GeV->Get(histname_zvertex60.c_str());
  TH1D* h_Jet40GeV_zvertex60_forcombine = (TH1D*)f_Jet40GeV->Get(histname_zvertex60.c_str());
  TH1D* h_Jet50GeV_zvertex60_forcombine = (TH1D*)f_Jet50GeV->Get(histname_zvertex60.c_str());
  //TH1D* h_Jet60GeV_zvertex60_forcombine = (TH1D*)f_Jet60GeV->Get(histname_zvertex60.c_str());

  //TH1D* h_zvertex60_combined = (TH1D*)h_MB_zvertex60_forcombine->Clone(histname_zvertex60.c_str());
  //h_zvertex60_combined->Scale(mb_scale_all);
  TH1D* h_zvertex60_combined = (TH1D*)h_Jet12GeV_zvertex60_forcombine->Clone(histname_zvertex60.c_str());
  h_zvertex60_combined->Scale(jet12_scale_all);
  //h_zvertex60_combined->Add(h_Jet5GeV_zvertex60_forcombine, jet5_scale_all);
  //h_zvertex60_combined->Add(h_Jet12GeV_zvertex60_forcombine, jet12_scale_all);
  h_zvertex60_combined->Add(h_Jet20GeV_zvertex60_forcombine, jet20_scale_all);
  h_zvertex60_combined->Add(h_Jet30GeV_zvertex60_forcombine, jet30_scale_all);
  h_zvertex60_combined->Add(h_Jet40GeV_zvertex60_forcombine, jet40_scale_all);
  h_zvertex60_combined->Add(h_Jet50GeV_zvertex60_forcombine, jet50_scale_all);
  //h_zvertex60_combined->Add(h_Jet60GeV_zvertex60_forcombine, jet60_scale_all);

  f_out->cd();
  h_all_combined->Write();
  h_zvertex60_combined->Write();
}

void combine_2Dhist(string histname, string surfix, /*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out) {
  //TH1D *h_MB_event_all = (TH1D*)f_MB->Get("h_event_all"); int mb_event_all = h_MB_event_all->GetBinContent(1); double mb_scale_all = mb_cross_section / (double)mb_event_all;
  //TH1D *h_Jet5GeV_event_all = (TH1D*)f_Jet5GeV->Get("h_event_all"); int jet5_event_all = h_Jet5GeV_event_all->GetBinContent(1); double jet5_scale_all = jet5_cross_section / (double)jet5_event_all;
  TH1D *h_Jet12GeV_event_all = (TH1D*)f_Jet12GeV->Get("h_event_all"); int jet12_event_all = h_Jet12GeV_event_all->GetBinContent(1); double jet12_scale_all = jet12_cross_section / (double)jet12_event_all;
  TH1D *h_Jet20GeV_event_all = (TH1D*)f_Jet20GeV->Get("h_event_all"); int jet20_event_all = h_Jet20GeV_event_all->GetBinContent(1); double jet20_scale_all = jet20_cross_section / (double)jet20_event_all;
  TH1D *h_Jet30GeV_event_all = (TH1D*)f_Jet30GeV->Get("h_event_all"); int jet30_event_all = h_Jet30GeV_event_all->GetBinContent(1); double jet30_scale_all = jet30_cross_section / (double)jet30_event_all;
  TH1D *h_Jet40GeV_event_all = (TH1D*)f_Jet40GeV->Get("h_event_all"); int jet40_event_all = h_Jet40GeV_event_all->GetBinContent(1); double jet40_scale_all = jet40_cross_section / (double)jet40_event_all;
  TH1D *h_Jet50GeV_event_all = (TH1D*)f_Jet50GeV->Get("h_event_all"); int jet50_event_all = h_Jet50GeV_event_all->GetBinContent(1); double jet50_scale_all = jet50_cross_section / (double)jet50_event_all;
  //TH1D *h_Jet60GeV_event_all = (TH1D*)f_Jet60GeV->Get("h_event_all"); int jet60_event_all = h_Jet60GeV_event_all->GetBinContent(1); double jet60_scale_all = jet60_cross_section / (double)jet60_event_all;

  string histname_all = histname + "_all" + surfix;
  string histname_zvertex60 = histname + "_zvertex60" + surfix;

  //TH2D* h_MB_all_forcombine = (TH2D*)f_MB->Get(histname_all.c_str());
  //TH2D* h_Jet5GeV_all_forcombine = (TH2D*)f_Jet5GeV->Get(histname_all.c_str());
  TH2D* h_Jet12GeV_all_forcombine = (TH2D*)f_Jet12GeV->Get(histname_all.c_str());
  TH2D* h_Jet20GeV_all_forcombine = (TH2D*)f_Jet20GeV->Get(histname_all.c_str());
  TH2D* h_Jet30GeV_all_forcombine = (TH2D*)f_Jet30GeV->Get(histname_all.c_str());
  TH2D* h_Jet40GeV_all_forcombine = (TH2D*)f_Jet40GeV->Get(histname_all.c_str());
  TH2D* h_Jet50GeV_all_forcombine = (TH2D*)f_Jet50GeV->Get(histname_all.c_str());
  //TH2D* h_Jet60GeV_all_forcombine = (TH2D*)f_Jet60GeV->Get(histname_all.c_str());

  //TH2D* h_all_combined = (TH2D*)h_MB_all_forcombine->Clone(histname_all.c_str());
  //h_all_combined->Scale(mb_scale_all);
  TH2D* h_all_combined = (TH2D*)h_Jet12GeV_all_forcombine->Clone(histname_all.c_str());
  h_all_combined->Scale(jet12_scale_all);
  //h_all_combined->Add(h_Jet5GeV_all_forcombine, jet5_scale_all);
  //h_all_combined->Add(h_Jet12GeV_all_forcombine, jet12_scale_all);
  h_all_combined->Add(h_Jet20GeV_all_forcombine, jet20_scale_all);
  h_all_combined->Add(h_Jet30GeV_all_forcombine, jet30_scale_all);
  h_all_combined->Add(h_Jet40GeV_all_forcombine, jet40_scale_all);
  h_all_combined->Add(h_Jet50GeV_all_forcombine, jet50_scale_all);
  //h_all_combined->Add(h_Jet60GeV_all_forcombine, jet60_scale_all);

  //TH2D* h_MB_zvertex60_forcombine = (TH2D*)f_MB->Get(histname_zvertex60.c_str());
  //TH2D* h_Jet5GeV_zvertex60_forcombine = (TH2D*)f_Jet5GeV->Get(histname_zvertex60.c_str());
  TH2D* h_Jet12GeV_zvertex60_forcombine = (TH2D*)f_Jet12GeV->Get(histname_zvertex60.c_str());
  TH2D* h_Jet20GeV_zvertex60_forcombine = (TH2D*)f_Jet20GeV->Get(histname_zvertex60.c_str());
  TH2D* h_Jet30GeV_zvertex60_forcombine = (TH2D*)f_Jet30GeV->Get(histname_zvertex60.c_str());
  TH2D* h_Jet40GeV_zvertex60_forcombine = (TH2D*)f_Jet40GeV->Get(histname_zvertex60.c_str());
  TH2D* h_Jet50GeV_zvertex60_forcombine = (TH2D*)f_Jet50GeV->Get(histname_zvertex60.c_str());
  //TH2D* h_Jet60GeV_zvertex60_forcombine = (TH2D*)f_Jet60GeV->Get(histname_zvertex60.c_str());

  //TH2D* h_zvertex60_combined = (TH2D*)h_MB_zvertex60_forcombine->Clone(histname_zvertex60.c_str());
  //h_zvertex60_combined->Scale(mb_scale_all);
  TH2D* h_zvertex60_combined = (TH2D*)h_Jet12GeV_zvertex60_forcombine->Clone(histname_zvertex60.c_str());
  h_zvertex60_combined->Scale(jet12_scale_all);
  //h_zvertex60_combined->Add(h_Jet5GeV_zvertex60_forcombine, jet5_scale_all);
  //h_zvertex60_combined->Add(h_Jet12GeV_zvertex60_forcombine, jet12_scale_all);
  h_zvertex60_combined->Add(h_Jet20GeV_zvertex60_forcombine, jet20_scale_all);
  h_zvertex60_combined->Add(h_Jet30GeV_zvertex60_forcombine, jet30_scale_all);
  h_zvertex60_combined->Add(h_Jet40GeV_zvertex60_forcombine, jet40_scale_all);
  h_zvertex60_combined->Add(h_Jet50GeV_zvertex60_forcombine, jet50_scale_all);
  //h_zvertex60_combined->Add(h_Jet60GeV_zvertex60_forcombine, jet60_scale_all);

  f_out->cd();
  h_all_combined->Write();
  h_zvertex60_combined->Write();
}

void combine_1Dhist_stat(string histname, string surfix, std::string datatype, /*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out) {
  //TH1D *h_MB_event_all = (TH1D*)f_MB->Get("h_event_all"); int mb_event_all = h_MB_event_all->GetBinContent(1); double mb_scale_all = mb_cross_section / (double)mb_event_all;
  //TH1D *h_Jet5GeV_event_all = (TH1D*)f_Jet5GeV->Get("h_event_all"); int jet5_event_all = h_Jet5GeV_event_all->GetBinContent(1); double jet5_scale_all = jet5_cross_section / (double)jet5_event_all;
  TH1D *h_Jet12GeV_event_all = (TH1D*)f_Jet12GeV->Get("h_event_all"); int jet12_event_all = h_Jet12GeV_event_all->GetBinContent(1); double jet12_scale_all = jet12_cross_section / (double)jet12_event_all;
  TH1D *h_Jet20GeV_event_all = (TH1D*)f_Jet20GeV->Get("h_event_all"); int jet20_event_all = h_Jet20GeV_event_all->GetBinContent(1); double jet20_scale_all = jet20_cross_section / (double)jet20_event_all;
  TH1D *h_Jet30GeV_event_all = (TH1D*)f_Jet30GeV->Get("h_event_all"); int jet30_event_all = h_Jet30GeV_event_all->GetBinContent(1); double jet30_scale_all = jet30_cross_section / (double)jet30_event_all;
  TH1D *h_Jet40GeV_event_all = (TH1D*)f_Jet40GeV->Get("h_event_all"); int jet40_event_all = h_Jet40GeV_event_all->GetBinContent(1); double jet40_scale_all = jet40_cross_section / (double)jet40_event_all;
  TH1D *h_Jet50GeV_event_all = (TH1D*)f_Jet50GeV->Get("h_event_all"); int jet50_event_all = h_Jet50GeV_event_all->GetBinContent(1); double jet50_scale_all = jet50_cross_section / (double)jet50_event_all;
  //TH1D *h_Jet60GeV_event_all = (TH1D*)f_Jet60GeV->Get("h_event_all"); int jet60_event_all = h_Jet60GeV_event_all->GetBinContent(1); double jet60_scale_all = jet60_cross_section / (double)jet60_event_all;

  string histname_all = histname + "_all" + surfix;
  string histname_zvertex60 = histname + "_zvertex60" + surfix;

  //TH1DBootstrap* h_all_combined = (TH1DBootstrap*)f_MB->Get(histname_all.c_str());
  //TH1DBootstrap* h_Jet5GeV_all_forcombine = (TH1DBootstrap*)f_Jet5GeV->Get(histname_all.c_str());
  //TH1DBootstrap* h_Jet12GeV_all_forcombine = (TH1DBootstrap*)f_Jet12GeV->Get(histname_all.c_str());
  TH1DBootstrap* h_all_combined = (TH1DBootstrap*)f_Jet12GeV->Get(histname_all.c_str());
  TH1DBootstrap* h_Jet20GeV_all_forcombine = (TH1DBootstrap*)f_Jet20GeV->Get(histname_all.c_str());
  TH1DBootstrap* h_Jet30GeV_all_forcombine = (TH1DBootstrap*)f_Jet30GeV->Get(histname_all.c_str());
  TH1DBootstrap* h_Jet40GeV_all_forcombine = (TH1DBootstrap*)f_Jet40GeV->Get(histname_all.c_str());
  TH1DBootstrap* h_Jet50GeV_all_forcombine = (TH1DBootstrap*)f_Jet50GeV->Get(histname_all.c_str());
  //TH1DBootstrap* h_Jet60GeV_all_forcombine = (TH1DBootstrap*)f_Jet60GeV->Get(histname_all.c_str());
  //h_all_combined->Scale(mb_scale_all);
  //h_all_combined->Add(h_Jet5GeV_all_forcombine, jet5_scale_all);
  //h_all_combined->Add(h_Jet12GeV_all_forcombine, jet12_scale_all);
  h_all_combined->Scale(jet12_scale_all);
  h_all_combined->Add(h_Jet20GeV_all_forcombine, jet20_scale_all);
  h_all_combined->Add(h_Jet30GeV_all_forcombine, jet30_scale_all);
  h_all_combined->Add(h_Jet40GeV_all_forcombine, jet40_scale_all);
  h_all_combined->Add(h_Jet50GeV_all_forcombine, jet50_scale_all);
  //h_all_combined->Add(h_Jet60GeV_all_forcombine, jet60_scale_all);

  //TH1DBootstrap* h_zvertex60_combined = (TH1DBootstrap*)f_MB->Get(histname_zvertex60.c_str());
  //TH1DBootstrap* h_Jet5GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet5GeV->Get(histname_zvertex60.c_str());
  //TH1DBootstrap* h_Jet12GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet12GeV->Get(histname_zvertex60.c_str());
  TH1DBootstrap* h_zvertex60_combined = (TH1DBootstrap*)f_Jet12GeV->Get(histname_zvertex60.c_str());
  TH1DBootstrap* h_Jet20GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet20GeV->Get(histname_zvertex60.c_str());
  TH1DBootstrap* h_Jet30GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet30GeV->Get(histname_zvertex60.c_str());
  TH1DBootstrap* h_Jet40GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet40GeV->Get(histname_zvertex60.c_str());
  TH1DBootstrap* h_Jet50GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet50GeV->Get(histname_zvertex60.c_str());
  //TH1DBootstrap* h_Jet60GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet60GeV->Get(histname_zvertex60.c_str());
  //h_zvertex60_combined->Scale(mb_scale_all);
  //h_zvertex60_combined->Add(h_Jet5GeV_zvertex60_forcombine, jet5_scale_all);
  //h_zvertex60_combined->Add(h_Jet12GeV_zvertex60_forcombine, jet12_scale_all);
  h_zvertex60_combined->Scale(jet12_scale_all);
  h_zvertex60_combined->Add(h_Jet20GeV_zvertex60_forcombine, jet20_scale_all);
  h_zvertex60_combined->Add(h_Jet30GeV_zvertex60_forcombine, jet30_scale_all);
  h_zvertex60_combined->Add(h_Jet40GeV_zvertex60_forcombine, jet40_scale_all);
  h_zvertex60_combined->Add(h_Jet50GeV_zvertex60_forcombine, jet50_scale_all);
  //h_zvertex60_combined->Add(h_Jet60GeV_zvertex60_forcombine, jet60_scale_all);

  f_out->cd();
  h_all_combined->Write();
  h_zvertex60_combined->Write();
}

void combine_2Dhist_stat(string histname, string surfix, /*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out) {
  //TH1D *h_MB_event_all = (TH1D*)f_MB->Get("h_event_all"); int mb_event_all = h_MB_event_all->GetBinContent(1); double mb_scale_all = mb_cross_section / (double)mb_event_all;
  //TH1D *h_Jet5GeV_event_all = (TH1D*)f_Jet5GeV->Get("h_event_all"); int jet5_event_all = h_Jet5GeV_event_all->GetBinContent(1); double jet5_scale_all = jet5_cross_section / (double)jet5_event_all;
  TH1D *h_Jet12GeV_event_all = (TH1D*)f_Jet12GeV->Get("h_event_all"); int jet12_event_all = h_Jet12GeV_event_all->GetBinContent(1); double jet12_scale_all = jet12_cross_section / (double)jet12_event_all;
  TH1D *h_Jet20GeV_event_all = (TH1D*)f_Jet20GeV->Get("h_event_all"); int jet20_event_all = h_Jet20GeV_event_all->GetBinContent(1); double jet20_scale_all = jet20_cross_section / (double)jet20_event_all;
  TH1D *h_Jet30GeV_event_all = (TH1D*)f_Jet30GeV->Get("h_event_all"); int jet30_event_all = h_Jet30GeV_event_all->GetBinContent(1); double jet30_scale_all = jet30_cross_section / (double)jet30_event_all;
  TH1D *h_Jet40GeV_event_all = (TH1D*)f_Jet40GeV->Get("h_event_all"); int jet40_event_all = h_Jet40GeV_event_all->GetBinContent(1); double jet40_scale_all = jet40_cross_section / (double)jet40_event_all;
  TH1D *h_Jet50GeV_event_all = (TH1D*)f_Jet50GeV->Get("h_event_all"); int jet50_event_all = h_Jet50GeV_event_all->GetBinContent(1); double jet50_scale_all = jet50_cross_section / (double)jet50_event_all;
  //TH1D *h_Jet60GeV_event_all = (TH1D*)f_Jet60GeV->Get("h_event_all"); int jet60_event_all = h_Jet60GeV_event_all->GetBinContent(1); double jet60_scale_all = jet60_cross_section / (double)jet60_event_all;

  string histname_all = histname + "_all" + surfix;
  string histname_zvertex60 = histname + "_zvertex60" + surfix;

  //TH2DBootstrap* h_all_combined = (TH2DBootstrap*)f_MB->Get(histname_all.c_str());
  //TH2DBootstrap* h_Jet5GeV_all_forcombine = (TH2DBootstrap*)f_Jet5GeV->Get(histname_all.c_str());
  //TH2DBootstrap* h_Jet12GeV_all_forcombine = (TH2DBootstrap*)f_Jet12GeV->Get(histname_all.c_str());
  TH2DBootstrap* h_all_combined = (TH2DBootstrap*)f_Jet12GeV->Get(histname_all.c_str());
  TH2DBootstrap* h_Jet20GeV_all_forcombine = (TH2DBootstrap*)f_Jet20GeV->Get(histname_all.c_str());
  TH2DBootstrap* h_Jet30GeV_all_forcombine = (TH2DBootstrap*)f_Jet30GeV->Get(histname_all.c_str());
  TH2DBootstrap* h_Jet40GeV_all_forcombine = (TH2DBootstrap*)f_Jet40GeV->Get(histname_all.c_str());
  TH2DBootstrap* h_Jet50GeV_all_forcombine = (TH2DBootstrap*)f_Jet50GeV->Get(histname_all.c_str());
  //TH2DBootstrap* h_Jet60GeV_all_forcombine = (TH2DBootstrap*)f_Jet60GeV->Get(histname_all.c_str());
  //h_all_combined->Scale(mb_scale_all);
  //h_all_combined->Add(h_Jet5GeV_all_forcombine, jet5_scale_all);
  //h_all_combined->Add(h_Jet12GeV_all_forcombine, jet12_scale_all);
  h_all_combined->Scale(jet12_scale_all);
  h_all_combined->Add(h_Jet20GeV_all_forcombine, jet20_scale_all);
  h_all_combined->Add(h_Jet30GeV_all_forcombine, jet30_scale_all);
  h_all_combined->Add(h_Jet40GeV_all_forcombine, jet40_scale_all);
  h_all_combined->Add(h_Jet50GeV_all_forcombine, jet50_scale_all);
  //h_all_combined->Add(h_Jet60GeV_all_forcombine, jet60_scale_all);

  //TH2DBootstrap* h_zvertex60_combined = (TH2DBootstrap*)f_MB->Get(histname_zvertex60.c_str());
  //TH2DBootstrap* h_Jet5GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet5GeV->Get(histname_zvertex60.c_str());
  //TH2DBootstrap* h_Jet12GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet12GeV->Get(histname_zvertex60.c_str());
  TH2DBootstrap* h_zvertex60_combined = (TH2DBootstrap*)f_Jet12GeV->Get(histname_zvertex60.c_str());
  TH2DBootstrap* h_Jet20GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet20GeV->Get(histname_zvertex60.c_str());
  TH2DBootstrap* h_Jet30GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet30GeV->Get(histname_zvertex60.c_str());
  TH2DBootstrap* h_Jet40GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet40GeV->Get(histname_zvertex60.c_str());
  TH2DBootstrap* h_Jet50GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet50GeV->Get(histname_zvertex60.c_str());
  //TH2DBootstrap* h_Jet60GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet60GeV->Get(histname_zvertex60.c_str());
  //h_zvertex60_combined->Scale(mb_scale_all);
  //h_zvertex60_combined->Add(h_Jet5GeV_zvertex60_forcombine, jet5_scale_all);
  //h_zvertex60_combined->Add(h_Jet12GeV_zvertex60_forcombine, jet12_scale_all);
  h_zvertex60_combined->Scale(jet12_scale_all);
  h_zvertex60_combined->Add(h_Jet20GeV_zvertex60_forcombine, jet20_scale_all);
  h_zvertex60_combined->Add(h_Jet30GeV_zvertex60_forcombine, jet30_scale_all);
  h_zvertex60_combined->Add(h_Jet40GeV_zvertex60_forcombine, jet40_scale_all);
  h_zvertex60_combined->Add(h_Jet50GeV_zvertex60_forcombine, jet50_scale_all);
  //h_zvertex60_combined->Add(h_Jet60GeV_zvertex60_forcombine, jet60_scale_all);

  f_out->cd();
  h_all_combined->Write();
  h_zvertex60_combined->Write();
}

void combine_1Dhist_closure(string surfix, string tag, /*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out) {
  //TH1D *h_MB_event_all = (TH1D*)f_MB->Get("h_event_all"); int mb_event_all = h_MB_event_all->GetBinContent(1); double mb_scale_all = mb_cross_section / (double)mb_event_all;
  //TH1D *h_Jet5GeV_event_all = (TH1D*)f_Jet5GeV->Get("h_event_all"); int jet5_event_all = h_Jet5GeV_event_all->GetBinContent(1); double jet5_scale_all = jet5_cross_section / (double)jet5_event_all;
  TH1D *h_Jet12GeV_event_all = (TH1D*)f_Jet12GeV->Get("h_event_all"); int jet12_event_all = h_Jet12GeV_event_all->GetBinContent(1); double jet12_scale_all = jet12_cross_section / (double)jet12_event_all;
  TH1D *h_Jet20GeV_event_all = (TH1D*)f_Jet20GeV->Get("h_event_all"); int jet20_event_all = h_Jet20GeV_event_all->GetBinContent(1); double jet20_scale_all = jet20_cross_section / (double)jet20_event_all;
  TH1D *h_Jet30GeV_event_all = (TH1D*)f_Jet30GeV->Get("h_event_all"); int jet30_event_all = h_Jet30GeV_event_all->GetBinContent(1); double jet30_scale_all = jet30_cross_section / (double)jet30_event_all;
  TH1D *h_Jet40GeV_event_all = (TH1D*)f_Jet40GeV->Get("h_event_all"); int jet40_event_all = h_Jet40GeV_event_all->GetBinContent(1); double jet40_scale_all = jet40_cross_section / (double)jet40_event_all;
  TH1D *h_Jet50GeV_event_all = (TH1D*)f_Jet50GeV->Get("h_event_all"); int jet50_event_all = h_Jet50GeV_event_all->GetBinContent(1); double jet50_scale_all = jet50_cross_section / (double)jet50_event_all;
  //TH1D *h_Jet60GeV_event_all = (TH1D*)f_Jet60GeV->Get("h_event_all"); int jet60_event_all = h_Jet60GeV_event_all->GetBinContent(1); double jet60_scale_all = jet60_cross_section / (double)jet60_event_all;

  string histname = surfix + tag + "_zvertex60";

  //TH1D* h_MB_zvertex60_forcombine = (TH1D*)f_MB->Get(histname.c_str());
  //TH1D* h_Jet5GeV_zvertex60_forcombine = (TH1D*)f_Jet5GeV->Get(histname.c_str());
  TH1D* h_Jet12GeV_zvertex60_forcombine = (TH1D*)f_Jet12GeV->Get(histname.c_str());
  TH1D* h_Jet20GeV_zvertex60_forcombine = (TH1D*)f_Jet20GeV->Get(histname.c_str());
  TH1D* h_Jet30GeV_zvertex60_forcombine = (TH1D*)f_Jet30GeV->Get(histname.c_str());
  TH1D* h_Jet40GeV_zvertex60_forcombine = (TH1D*)f_Jet40GeV->Get(histname.c_str());
  TH1D* h_Jet50GeV_zvertex60_forcombine = (TH1D*)f_Jet50GeV->Get(histname.c_str());
  //TH1D* h_Jet60GeV_zvertex60_forcombine = (TH1D*)f_Jet60GeV->Get(histname.c_str());

  //TH1D* h_zvertex60_combined = (TH1D*)h_MB_zvertex60_forcombine->Clone(histname.c_str());
  //h_zvertex60_combined->Scale(mb_scale_all);
  TH1D* h_zvertex60_combined = (TH1D*)h_Jet12GeV_zvertex60_forcombine->Clone(histname.c_str());
  h_zvertex60_combined->Scale(jet12_scale_all);
  //h_zvertex60_combined->Add(h_Jet5GeV_zvertex60_forcombine, jet5_scale_all);
  //h_zvertex60_combined->Add(h_Jet12GeV_zvertex60_forcombine, jet12_scale_all);
  h_zvertex60_combined->Add(h_Jet20GeV_zvertex60_forcombine, jet20_scale_all);
  h_zvertex60_combined->Add(h_Jet30GeV_zvertex60_forcombine, jet30_scale_all);
  h_zvertex60_combined->Add(h_Jet40GeV_zvertex60_forcombine, jet40_scale_all);
  h_zvertex60_combined->Add(h_Jet50GeV_zvertex60_forcombine, jet50_scale_all);
  //h_zvertex60_combined->Add(h_Jet60GeV_zvertex60_forcombine, jet60_scale_all);

  f_out->cd();
  h_zvertex60_combined->Write();
}

void combine_2Dhist_closure(string surfix, string tag, /*TFile* f_MB, TFile* f_Jet5GeV, */TFile* f_Jet12GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet40GeV, TFile* f_Jet50GeV, /*TFile* f_Jet60GeV, */TFile* f_out) {
  //TH1D *h_MB_event_all = (TH1D*)f_MB->Get("h_event_all"); int mb_event_all = h_MB_event_all->GetBinContent(1); double mb_scale_all = mb_cross_section / (double)mb_event_all;
  //TH1D *h_Jet5GeV_event_all = (TH1D*)f_Jet5GeV->Get("h_event_all"); int jet5_event_all = h_Jet5GeV_event_all->GetBinContent(1); double jet5_scale_all = jet5_cross_section / (double)jet5_event_all;
  TH1D *h_Jet12GeV_event_all = (TH1D*)f_Jet12GeV->Get("h_event_all"); int jet12_event_all = h_Jet12GeV_event_all->GetBinContent(1); double jet12_scale_all = jet12_cross_section / (double)jet12_event_all;
  TH1D *h_Jet20GeV_event_all = (TH1D*)f_Jet20GeV->Get("h_event_all"); int jet20_event_all = h_Jet20GeV_event_all->GetBinContent(1); double jet20_scale_all = jet20_cross_section / (double)jet20_event_all;
  TH1D *h_Jet30GeV_event_all = (TH1D*)f_Jet30GeV->Get("h_event_all"); int jet30_event_all = h_Jet30GeV_event_all->GetBinContent(1); double jet30_scale_all = jet30_cross_section / (double)jet30_event_all;
  TH1D *h_Jet40GeV_event_all = (TH1D*)f_Jet40GeV->Get("h_event_all"); int jet40_event_all = h_Jet40GeV_event_all->GetBinContent(1); double jet40_scale_all = jet40_cross_section / (double)jet40_event_all;
  TH1D *h_Jet50GeV_event_all = (TH1D*)f_Jet50GeV->Get("h_event_all"); int jet50_event_all = h_Jet50GeV_event_all->GetBinContent(1); double jet50_scale_all = jet50_cross_section / (double)jet50_event_all;
  //TH1D *h_Jet60GeV_event_all = (TH1D*)f_Jet60GeV->Get("h_event_all"); int jet60_event_all = h_Jet60GeV_event_all->GetBinContent(1); double jet60_scale_all = jet60_cross_section / (double)jet60_event_all;

  string histname = surfix + tag + "_zvertex60";

  //TH2D* h_MB_zvertex60_forcombine = (TH2D*)f_MB->Get(histname.c_str());
  //TH2D* h_Jet5GeV_zvertex60_forcombine = (TH2D*)f_Jet5GeV->Get(histname.c_str());
  TH2D* h_Jet12GeV_zvertex60_forcombine = (TH2D*)f_Jet12GeV->Get(histname.c_str());
  TH2D* h_Jet20GeV_zvertex60_forcombine = (TH2D*)f_Jet20GeV->Get(histname.c_str());
  TH2D* h_Jet30GeV_zvertex60_forcombine = (TH2D*)f_Jet30GeV->Get(histname.c_str());
  TH2D* h_Jet40GeV_zvertex60_forcombine = (TH2D*)f_Jet40GeV->Get(histname.c_str());
  TH2D* h_Jet50GeV_zvertex60_forcombine = (TH2D*)f_Jet50GeV->Get(histname.c_str());
  //TH2D* h_Jet60GeV_zvertex60_forcombine = (TH2D*)f_Jet60GeV->Get(histname.c_str());

  //TH2D* h_zvertex60_combined = (TH2D*)h_MB_zvertex60_forcombine->Clone(histname.c_str());
  //h_zvertex60_combined->Scale(mb_scale_all);
  TH2D* h_zvertex60_combined = (TH2D*)h_Jet12GeV_zvertex60_forcombine->Clone(histname.c_str());
  h_zvertex60_combined->Scale(jet12_scale_all);
  //h_zvertex60_combined->Add(h_Jet5GeV_zvertex60_forcombine, jet5_scale_all);
  //h_zvertex60_combined->Add(h_Jet12GeV_zvertex60_forcombine, jet12_scale_all);
  h_zvertex60_combined->Add(h_Jet20GeV_zvertex60_forcombine, jet20_scale_all);
  h_zvertex60_combined->Add(h_Jet30GeV_zvertex60_forcombine, jet30_scale_all);
  h_zvertex60_combined->Add(h_Jet40GeV_zvertex60_forcombine, jet40_scale_all);
  h_zvertex60_combined->Add(h_Jet50GeV_zvertex60_forcombine, jet50_scale_all);
  //h_zvertex60_combined->Add(h_Jet60GeV_zvertex60_forcombine, jet60_scale_all);

  f_out->cd();
  h_zvertex60_combined->Write();
}