#include "AnNeutralMeson_micro_dst.h"
#include "ClusterSmallInfoContainer.h"
#include "ClusterSmallInfo.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <globalvertex/GlobalVertex.h>

// Spin DB
#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/resultsetmetadata.h>
#include <uspin/SpinDBContent.h>
#include <uspin/SpinDBContentv1.h>
#include <uspin/SpinDBOutput.h>

#include <jetbase/JetContainer.h>
#include <jetbase/Jet.h>

#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <numeric> // for std::iota
#include <cstdlib> // for rand
#include <cmath> // for fmod

#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <Math/GenVector/VectorUtil.h> // For DeltaR
#include <Math/LorentzRotation.h>
#include <Math/Rotation3D.h>
#include <Math/AxisAngle.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <bitset>

// To check virtual and resident memory usage
#include <Riostream.h>
#include <TSystem.h>
#include <malloc.h>

void monitorMemoryUsage(const std::string &label)
{
  ProcInfo_t procInfo;
  gSystem->GetProcInfo(&procInfo);
  std::cout << label << " - Memory usage: "
            << "Resident = " << procInfo.fMemResident << " kB, "
            << "Virtual = " << procInfo.fMemVirtual << " kB" << std::endl;
}

AnNeutralMeson_micro_dst::AnNeutralMeson_micro_dst(const std::string &name, const int runnb,
                                                   const std::string &outputfilename,
                                                   const std::string &outputfiletreename,
                                                   const int seednb)
  : SubsysReco(name)
  , runnumber(runnb)
  , outfilename(outputfilename)
  , outtreename(outputfiletreename)
  , seednumber(seednb)
{
}

AnNeutralMeson_micro_dst::~AnNeutralMeson_micro_dst()
{
  delete outfile;
  delete outfile_tree;
}

int AnNeutralMeson_micro_dst::Init(PHCompositeNode *)
{
  if (do_smearing)
  {
    rnd = new TRandom3(seednumber * 1e5 + runnumber);
    const std::string kSmearParFile =
      "/sphenix/u/virgilemahaut/fromBlair/function_compare_wide.root";
    TFile* smear_fin = TFile::Open(kSmearParFile.c_str(), "READ");
    if (!smear_fin || smear_fin->IsZombie())
    {
      std::cout << Name() << "::Init - WARNING: could not open photon smearing parameter file "
                << kSmearParFile << std::endl;
    }
    else
  {
      static const std::array<const char*, 3> kEnergyFuncNames = {
          "f_energy_difference_global_minimum", "f_energy_difference_cE_0p08", "f_energy_difference_cE_0p05"};
      static const std::array<const char*, 3> kPositionFuncNames = {"f_position_difference_global_minimum",
                                                                     "f_position_difference_cE_0p08",
                                                                     "f_position_difference_cE_0p05"};
      static const std::array<const char*, 3> kEnergyFuncNamesLegacy = {
          "g_energy_difference_global_minimum", "g_energy_difference_cE_0p08", "g_energy_difference_cE_0p05"};
      static const std::array<const char*, 3> kPositionFuncNamesLegacy = {
          "g_position_difference_global_minimum", "g_position_difference_cE_0p08", "g_position_difference_cE_0p05"};
      for (int k = 0; k < 3; ++k)
      {
        TF1* fe = dynamic_cast<TF1*>(smear_fin->Get(kEnergyFuncNames[k]));
        if (!fe)
        {
          fe = dynamic_cast<TF1*>(smear_fin->Get(kEnergyFuncNamesLegacy[k]));
        }
        if (!fe && k == 0)
        {
          fe = dynamic_cast<TF1*>(smear_fin->Get("f_energy_difference"));
        }
        if (!fe && k == 0)
        {
          fe = dynamic_cast<TF1*>(smear_fin->Get("g_energy_difference"));
        }
        if (!fe)
        {
          std::cout << Name() << "::Init - WARNING: energy smear function " << kEnergyFuncNames[k]
                    << " not found in " << kSmearParFile << std::endl;
        }
        else
        {
          f_energy_smear[k] =
              static_cast<TF1*>(fe->Clone(("f_energy_smear" + std::to_string(k) + "_clone").c_str()));
        }

        TF1* fp = dynamic_cast<TF1*>(smear_fin->Get(kPositionFuncNames[k]));
        if (!fp)
        {
          fp = dynamic_cast<TF1*>(smear_fin->Get(kPositionFuncNamesLegacy[k]));
        }
        if (!fp && k == 0)
        {
          fp = dynamic_cast<TF1*>(smear_fin->Get("f_position_difference"));
        }
        if (!fp && k == 0)
        {
          fp = dynamic_cast<TF1*>(smear_fin->Get("g_position_difference"));
        }
        if (!fp)
        {
          std::cout << Name() << "::Init - WARNING: position smear function " << kPositionFuncNames[k]
                    << " not found in " << kSmearParFile << std::endl;
        }
        else
        {
          f_position_smear[k] =
              static_cast<TF1*>(fp->Clone(("f_position_smear" + std::to_string(k) + "_clone").c_str()));
        }
      }
      smear_fin->Close();
      delete smear_fin;
    }
  }
  
  // Book tree (used for nano analysis)
  if (store_tree)
  {
    outfile_tree = new TFile(outtreename.c_str(),
                             "RECREATE");
    output_tree = new TTree(
      "tree_diphoton_compact",
      "Minimal diphoton info");
    output_tree->Branch(
      "diphoton_vertex_z", &diphoton_vertex_z,
      "diphoton_vertex_z/F");
    output_tree->Branch(
      "diphoton_bunchnumber", &diphoton_bunchnumber,
      "diphoton_bunchnumber/I");
    output_tree->Branch(
      "diphoton_mass", &diphoton_mass,
      "diphoton_mass/F");
    output_tree->Branch(
      "diphoton_eta", &diphoton_eta,
      "diphoton_eta/F");
    output_tree->Branch(
      "diphoton_phi", &diphoton_phi,
      "diphoton_phi/F");
    output_tree->Branch(
      "diphoton_pt", &diphoton_pt,
      "diphoton_pt/F");
    output_tree->Branch(
      "diphoton_xf", &diphoton_xf,
      "diphoton_xf/F");
  }
  
  // Book histograms
  outfile = new TFile(outfilename.c_str(),
                      "RECREATE");

  if (require_emulator_matching) {
    h_emulator_selection = new TH1I(
      "h_emulator_selection",
      ";Trigger Tile Count; Count",
      10,0,10);
  }
  
  if (require_efficiency_matching && require_emulator_matching)
  {
    h_efficiency_x_matching = new TH1I(
      "h_efficiency_x_matching",
      ";case;Counts",
      4, 1, 4);

    for (int iPt = 0; iPt < nPtBins; iPt++) {
      std::stringstream h_efficiency_x_matching_name;
      h_efficiency_x_matching_name << "h_efficiency_x_matching_pt_" << iPt;

      h_efficiency_x_matching_pt[iPt] = new TH1I(h_efficiency_x_matching_name.str().c_str(),
                                                 ";case;Counts",
                                                 4, 1, 4);
    }
  }

  if (do_smearing) {
    h_smear_E = new TH1F(
      "h_smear_E",
      "; E_{Smeared}/E_{Nominal}; Counts",
      100, 0.5, 1.5);

    h_smear_eta = new TH1F(
      "h_smear_eta",
      "; #eta_{Smeared}-#eta_{Nominal}; Counts",
      100, -0.025, 0.025);

    h_smear_phi = new TH1F(
      "h_smear_phi",
      "; #phi_{Smeared}-#phi_{Nominal}; Counts",
      100, -0.025, 0.025);

    h_smear_E_dE = new TH2F(
      "h_smear_E_dE",
      "; E [GeV]; E_{Smeared}/E_{Nominal}",
      100, 0, 10,
      100, 0.5, 1.5);

    h_smear_E_deta = new TH2F(
      "h_smear_E_deta",
      "; E [GeV]; #eta_{Smeared}-#eta_{Nominal}",
      100, 0, 10,
      100, -0.025, 0.025);

    h_smear_E_dphi = new TH2F(
      "h_smear_E_dphi",
      "; E [GeV]; #phi_{Smeared}-#phi_{Nominal}",
      100, 0, 10,
      100, -0.025, 0.025);
  }
  
  if (store_qa)
  {
    // Event QA
    h_multiplicity_efficiency = new TH1F(
      "h_multiplicity_efficiency",
      ";multiplicity; counts",
      10, 0, 10);

    h_multiplicity_emulator = new TH1F(
      "h_multiplicity_emulator",
      ";multiplicity; counts",
      10, 0, 10);

    h_event_zvtx = new TH1F(
      "h_event_zvtx",
      ";vertex z [cm]; counts",
      200, -200, 200);
  
    h_triggered_event_zvtx = new TH1F(
      "h_triggered_event_zvtx",
      ";vertex z [cm]; counts",
      200, -200, 200);

    h_triggered_event_photons = new TH1F(
      "h_triggered_event_photons",
      ";# photons; counts",
      20, 0, 10);

    h_analysis_event_zvtx = new TH1F(
      "h_analysis_event_zvtx",
      ";vertex z [cm]; counts",
      200, -200, 200);
  
    // event trigger QA
    h_trigger_live =
      new TH1F(
          "h_trigger_live",
          "Trigger Live;;Counts", 32, -0.5, 32 - 0.5);
    h_trigger_scaled = new TH1F(
      "h_trigger_scaled",
      "Trigger Scaled;;Counts", 32,
      -0.5, 32 - 0.5);

    for (const auto &[trigger_number,
                    trigger_name] : trigger_names)
    {
      h_trigger_live->GetXaxis()->SetBinLabel(trigger_number + 1,
                                              trigger_name.c_str());
      h_trigger_scaled->GetXaxis()->SetBinLabel(trigger_number + 1,
                                                trigger_name.c_str());
    }

    h_count_diphoton_nomatch = new TH1I(
      "h_count_diphoton_nomatch",
      ";>= 1 diphoton; Count",
      1,0,1);

    // Event mixing QA
    h_current_E = new TH1F(
      "h_current_E",
      ";E [GeV]; counts",
      200, 0, 20);
    h_current_eta = new TH1F(
      "h_current_eta",
      ";#eta [rad];counts",
      200, -2.0, 2.0);
    h_current_phi = new TH1F(
      "h_current_phi",
      ";#phi [rad];counts",
      128, -M_PI, M_PI);
    h_pool_E = new TH1F(
      "h_pool_E",
      ";E [GeV]; counts",
      200, 0, 20);
    h_pool_eta = new TH1F(
      "h_pool_eta",
      ";#eta [rad];counts",
      200, -2.0, 2.0);
    h_pool_phi = new TH1F(
      "h_pool_phi",
      ";#phi [rad];counts",
      128, -M_PI, M_PI);

    h_same_delta_eta = new TH1F(
      "h_same_delta_eta",
      ";#Delta #eta [rad];counts",
      500, 0, 10);
    h_same_delta_phi = new TH1F(
      "h_same_delta_phi",
      ";#Delta #phi [rad];counts",
      500, 0, 10);
    h_same_delta_R = new TH1F(
      "h_same_delta_R",
      ";#Delta R [rad];counts",
      500, 0, 10);
    h_same_alpha = new TH1F(
      "h_same_alpha",
      ";#alpha [rad];counts",
      500, 0, 1);
    h_mixed_delta_eta = new TH1F(
      "h_mixed_delta_eta",
      ";#Delta #eta [rad];counts",
      500, 0, 10);
    h_mixed_delta_phi = new TH1F(
      "h_mixed_delta_phi",
      ";#Delta #phi [rad];counts",
      500, 0, 10);
    h_mixed_delta_R = new TH1F(
      "h_mixed_delta_R",
      ";#Delta R [rad];counts",
      500, 0, 10);
    h_mixed_alpha = new TH1F(
      "h_mixed_alpha",
      ";#alpha [rad];counts",
      500, 0, 1);

    // Single photon QA
    h_photon_eta = new TH1F(
      "h_photon_eta",
      ";#eta [rad];counts",
      200, -2.0, 2.0);
    h_photon_phi = new TH1F(
      "h_photon_phi",
      ";#phi [rad];counts",
      128, -M_PI, M_PI);
    h_photon_pt = new TH1F(
      "h_photon_pt",
      ";p_{T} [GeV];counts",
      200, 0, 20);
    h_photon_zvtx = new TH1F(
      "h_photon_zvtx",
      ";z vertex [cm];counts",
      200, -200, 200);
  
    h_photon_eta_phi = new TH2F(
      "h_photon_eta_phi",
      ";#eta [rad];#phi [rad]",
      200, -2.0, 2.0, 128, -M_PI, M_PI);
    h_photon_eta_pt = new TH2F(
      "h_photon_eta_pt",
      ";#eta [rad];p_{T} [GeV]",
      200, -2.0, 2.0, 200, 0, 20);
    h_photon_eta_zvtx = new TH2F(
      "h_photon_eta_zvtx",
      ";#eta [rad];vertex z [cm]; counts",
      200, -2, 2, 200, -200, 200);
    h_photon_phi_pt = new TH2F(
      "h_photon_phi_pt",
      ";#phi [rad];p_{T} [GeV]; counts",
      128, -M_PI, M_PI, 200, 0, 20);
    h_photon_phi_zvtx = new TH2F(
      "h_photon_phi_zvtx",
      ";#phi [rad];vertex z [cm]; counts",
      128, -M_PI, M_PI, 200, -200, 200);
    h_photon_pt_zvtx = new TH2F(
      "h_photon_pt_zvtx",
      ";p_{T} [GeV];vertex z [cm]; counts",
      200, 0, 20, 200, -200, 200);

    for (int iPt = 0; iPt < nPtBins; iPt++) {
      h_photon_eta_phi_bin[iPt] = new TH2F(
       ("h_photon_eta_phi_bin_" + std::to_string(iPt)).c_str(),
       ";#eta; #phi [rad]",
       200, -2.0, 2.0, 128, -M_PI, M_PI);
      h_photon_phi_pt = new TH2F(
        ("h_photon_phi_pt_bin_" + std::to_string(iPt)).c_str(),
        ";#phi [rad];p_{T} [GeV]; counts",
        128, -M_PI, M_PI, 200, 0, 20);
    }
    
    h_selected_photon_eta = new TH1F(
      "h_selected_photon_eta",
      ";#eta [rad];counts",
      200, -2.0, 2.0);
    h_selected_photon_phi = new TH1F(
      "h_selected_photon_phi",
      ";#phi [rad];counts",
      128, -M_PI, M_PI);
    h_selected_photon_pt = new TH1F(
      "h_selected_photon_pt",
      ";p_{T} [GeV];counts",
      200, 0, 20);
    h_selected_photon_zvtx = new TH1F(
      "h_selected_photon_zvtx",
      ";z vertex [cm];counts",
      200, -200, 200);
  
    h_selected_photon_eta_phi = new TH2F(
      "h_selected_photon_eta_phi",
      ";#eta [rad];#phi [rad]",
      200, -2.0, 2.0, 128, -M_PI, M_PI);
    h_selected_photon_eta_pt = new TH2F(
      "h_selected_photon_eta_pt",
      ";#eta [rad];p_{T} [GeV]",
      200, -2.0, 2.0, 200, 0, 20);
    h_selected_photon_eta_zvtx = new TH2F(
      "h_selected_photon_eta_zvtx",
      ";#eta [rad];vertex z [cm]; counts",
      200, -2, 2, 200, -200, 200);
    h_selected_photon_phi_pt = new TH2F(
      "h_selected_photon_phi_pt",
      ";#phi [rad];p_{T} [GeV]; counts",
      128, -M_PI, M_PI, 200, 0, 20);
    h_selected_photon_phi_zvtx = new TH2F(
      "h_selected_photon_phi_zvtx",
      ";#phi [rad];vertex z [cm]; counts",
      128, -M_PI, M_PI, 200, -200, 200);
    h_selected_photon_pt_zvtx = new TH2F(
      "h_selected_photon_pt_zvtx",
      ";p_{T} [GeV];vertex z [cm]; counts",
      200, 0, 20, 200, -200, 200);

    for (int iPt = 0; iPt < nPtBins; iPt++) {
      h_selected_photon_eta_phi_bin[iPt] = new TH2F(
       ("h_selected_photon_eta_phi_bin_" + std::to_string(iPt)).c_str(),
       ";#eta; #phi [rad]",
       200, -2.0, 2.0, 128, -M_PI, M_PI);
      h_selected_photon_phi_pt = new TH2F(
        ("h_selected_photon_phi_pt_bin_" + std::to_string(iPt)).c_str(),
        ";#phi [rad];p_{T} [GeV]; counts",
        128, -M_PI, M_PI, 200, 0, 20);
    }

    // diphoton QA
    h_pair_eta = new TH1F(
      "h_pair_eta",
      ";#eta [rad];counts",
      200, -2.0, 2.0);
    h_pair_phi = new TH1F(
      "h_pair_phi",
      ";#phi [rad];counts",
      128, -M_PI, M_PI);
    h_pair_pt = new TH1F(
      "h_pair_pt",
      ";p_{T} [GeV];counts",
      200, 0, 20);
    h_pair_zvtx = new TH1F(
      "h_pair_zvtx",
      ";z vertex [cm];counts",
      200, -200, 200);
    h_pair_DeltaR = new TH1F(
      "h_pair_DeltaR",
      ";#Delta R [rad];counts",
      500, 0, 5);
    h_pair_alpha = new TH1F(
      "h_pair_alpha",
      ";#alpha [rad];counts",
      200, 0, 2);
    h_pair_alpha_pt = new TH2F(
      "h_pair_alpha_pt",
      ";#alpha; p_{T}; counts",
      200, 0, 2, 200, 0, 20);
    h_pair_eta_phi = new TH2F(
      "h_pair_eta_phi",
      ";#eta [rad];#phi [rad]",
      200, -2.0, 2.0, 128, -M_PI, M_PI);
    h_pair_eta_pt = new TH2F(
      "h_pair_eta_pt",
      ";#eta [rad];p_{T} [GeV]",
      200, -2.0, 2.0, 200, 0, 20);
    h_pair_eta_zvtx = new TH2F(
      "h_pair_eta_zvtx",
      ";#eta [rad];vertex z [cm]; counts",
      200, -2, 2, 200, -200, 200);
    h_pair_phi_pt = new TH2F(
      "h_pair_phi_pt",
      ";#phi [rad];p_{T} [GeV]; counts",
      128, -M_PI, M_PI, 200, 0, 20);
    h_pair_phi_zvtx = new TH2F(
      "h_pair_phi_zvtx",
      ";#phi [rad];vertex z [cm]; counts",
      128, -M_PI, M_PI, 200, -200, 200);
    h_pair_pt_zvtx = new TH2F(
      "h_pair_pt_zvtx",
      ";p_{T} [GeV];vertex z [cm]; counts",
      200, 0, 20, 200, -200, 200);

    // diphoton eta pT-bin by pT-bin
    h_pair_pi0_eta_pt_1 = new TH1F( // 1-2 GeV
      "h_pair_pi0_eta_pt_1",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_pi0_eta_pt_2 = new TH1F( // 2-3 GeV
      "h_pair_pi0_eta_pt_2",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_pi0_eta_pt_3 = new TH1F( // 3-4 GeV
      "h_pair_pi0_eta_pt_3",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_pi0_eta_pt_4 = new TH1F( // 4-5 GeV
      "h_pair_pi0_eta_pt_4",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_pi0_eta_pt_5 = new TH1F( // 5-6 GeV
      "h_pair_pi0_eta_pt_5",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_pi0_eta_pt_6 = new TH1F( // 6-7 GeV
      "h_pair_pi0_eta_pt_6",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_pi0_eta_pt_7 = new TH1F( // 7-8 GeV
      "h_pair_pi0_eta_pt_7",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_pi0_eta_pt_8 = new TH1F( // 8-10 GeV
      "h_pair_pi0_eta_pt_8",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_pi0_eta_pt_9 = new TH1F( // 10-20 GeV
      "h_pair_pi0_eta_pt_9",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_eta_eta_pt_1 = new TH1F( // 1-2 GeV
      "h_pair_eta_eta_pt_1",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_eta_eta_pt_2 = new TH1F( // 2-3 GeV
      "h_pair_eta_eta_pt_2",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_eta_eta_pt_3 = new TH1F( // 3-4 GeV
      "h_pair_eta_eta_pt_3",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_eta_eta_pt_4 = new TH1F( // 4-5 GeV
      "h_pair_eta_eta_pt_4",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_eta_eta_pt_5 = new TH1F( // 5-6 GeV
      "h_pair_eta_eta_pt_5",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_eta_eta_pt_6 = new TH1F( // 6-7 GeV
      "h_pair_eta_eta_pt_6",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_eta_eta_pt_7 = new TH1F( // 7-8 GeV
      "h_pair_eta_eta_pt_7",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_eta_eta_pt_8 = new TH1F( // 8-10 GeV
      "h_pair_eta_eta_pt_8",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    h_pair_eta_eta_pt_9 = new TH1F( // 10-20 GeV
      "h_pair_eta_eta_pt_9",
      ";#eta [rad]; Counts / [20 mrad]",
      200,-2.0,2.0);

    // diphoton xF pT-bin by pT-bin

    h_pair_pi0_xf_pt_1 = new TH1F( // 1-2 GeV
      "h_pair_pi0_xf_pt_1",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_pi0_xf_pt_2 = new TH1F( // 2-3 GeV
      "h_pair_pi0_xf_pt_2",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_pi0_xf_pt_3 = new TH1F( // 3-4 GeV
      "h_pair_pi0_xf_pt_3",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_pi0_xf_pt_4 = new TH1F( // 4-5 GeV
      "h_pair_pi0_xf_pt_4",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_pi0_xf_pt_5 = new TH1F( // 5-6 GeV
      "h_pair_pi0_xf_pt_5",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_pi0_xf_pt_6 = new TH1F( // 6-7 GeV
      "h_pair_pi0_xf_pt_6",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_pi0_xf_pt_7 = new TH1F( // 7-8 GeV
      "h_pair_pi0_xf_pt_7",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_pi0_xf_pt_8 = new TH1F( // 8-10 GeV
      "h_pair_pi0_xf_pt_8",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_pi0_xf_pt_9 = new TH1F( // 10-20 GeV
      "h_pair_pi0_xf_pt_9",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_eta_xf_pt_1 = new TH1F( // 1-2 GeV
      "h_pair_eta_xf_pt_1",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_eta_xf_pt_2 = new TH1F( // 2-3 GeV
      "h_pair_eta_xf_pt_2",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_eta_xf_pt_3 = new TH1F( // 3-4 GeV
      "h_pair_eta_xf_pt_3",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_eta_xf_pt_4 = new TH1F( // 4-5 GeV
      "h_pair_eta_xf_pt_4",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_eta_xf_pt_5 = new TH1F( // 5-6 GeV
      "h_pair_eta_xf_pt_5",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_eta_xf_pt_6 = new TH1F( // 6-7 GeV
      "h_pair_eta_xf_pt_6",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_eta_xf_pt_7 = new TH1F( // 7-8 GeV
      "h_pair_eta_xf_pt_7",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_eta_xf_pt_8 = new TH1F( // 8-10 GeV
      "h_pair_eta_xf_pt_8",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_pair_eta_xf_pt_9 = new TH1F( // 10-20 GeV
      "h_pair_eta_xf_pt_9",
      ";x_{F} [rad]; Counts / [2 mrad]",
      200,-0.2,0.2);

    h_meson_pi0_pt = new TH1F(
      "h_meson_pi0_pt",
      ";p_{T} [GeV];counts",
      200, 0, 20);

    h_meson_pi0_E = new TH1F(
      "h_meson_pi0_E",
      ";E [GeV];counts",
      200, 0, 20);
  
    h_meson_eta_pt = new TH1F(
      "h_meson_eta_pt",
      ";p_{T} [GeV];counts",
      200, 0, 20);

    h_meson_eta_E = new TH1F(
      "h_meson_eta_E",
      ";E [GeV];counts",
      200, 0, 20);
  }
  
  // diphoton invariant mass
  h_pair_mass = new TH1F(
      "h_pair_mass",
      ";M_{#gamma} [GeV];counts", 500, 0, 1);

  h_pair_mass_mixing = new TH1F(
      "h_pair_mass_mixing",
      ";M_{#gamma #gamma} [GeV];counts", 500, 0, 1);

  for (int iPt = 0; iPt < nPtBins; iPt++)
  {
    std::stringstream h_pair_mass_name;
    h_pair_mass_name << std::fixed << std::setprecision(0) << "h_pair_mass_pt_"
                     << iPt;
    std::stringstream h_pair_mass_title;
    h_pair_mass_title << std::fixed << std::setprecision(1) << "diphoton mass ["
                      << pTBins[iPt] << " < p_{T} [GeV/c] < "
                      << pTBins[iPt + 1]
                      << "];M_{#gamma#gamma} [GeV/c^{2}];counts";
    h_pair_mass_pt[iPt] = new TH1F(h_pair_mass_name.str().c_str(),
                                   h_pair_mass_title.str().c_str(), 500, 0, 1);
    h_pair_mass_name.str("");
    h_pair_mass_name << std::fixed << std::setprecision(0) << "h_pair_mass_mixing_pt_"
                     << iPt;
    h_pair_mass_mixing_pt[iPt] = new TH1F(h_pair_mass_name.str().c_str(),
                                            h_pair_mass_title.str().c_str(), 500, 0, 1);
  }

  for (int izvtx = 0; izvtx < nZvtxBins; izvtx++)
  {
    std::stringstream h_pair_mass_name;
    h_pair_mass_name << std::fixed << std::setprecision(0) << "h_pair_mass_zvtx_"
                     << izvtx;
    std::stringstream h_pair_mass_title;
    h_pair_mass_title << std::fixed << std::setprecision(1) << "diphoton mass ["
                      << - zvtxBins[izvtx + 1] << " < z_{vtx} [cm] < "
                      << zvtxBins[izvtx + 1]
                      << "];M_{#gamma#gamma} [GeV/c^{2}];counts";
    h_pair_mass_zvtx[izvtx] = new TH1F(h_pair_mass_name.str().c_str(),
                                       h_pair_mass_title.str().c_str(), 500, 0, 1);
  }

  for (int ietaBin = 0; ietaBin < nEtaBins; ietaBin++)
  {
    std::stringstream h_pair_mass_name;
    h_pair_mass_name << std::fixed << std::setprecision(0) << "h_pair_mass_eta_"
                     << ietaBin;
    std::stringstream h_pair_mass_title;
    h_pair_mass_title << std::fixed << std::setprecision(2) << "diphoton mass ["
                      << etaBins[ietaBin] << " < #eta < "
                      << etaBins[ietaBin + 1] << "];M_{#gamma#gamma};counts";
    h_pair_mass_eta[ietaBin] =
        new TH1F(h_pair_mass_name.str().c_str(),
                 h_pair_mass_title.str().c_str(), 500, 0, 1);
  }

  for (int ixfBin = 0; ixfBin < nXfBins; ixfBin++)
  {
    std::stringstream h_pair_mass_name;
    h_pair_mass_name << std::fixed << std::setprecision(0) << "h_pair_mass_xf_"
                     << ixfBin;
    std::stringstream h_pair_mass_title;
    h_pair_mass_title << std::fixed << std::setprecision(2) << "diphoton mass ["
                      << xfBins[ixfBin] << " < x_{F} < " << xfBins[ixfBin + 1]
                      << "];M_{#gamma#gamma};counts";
    h_pair_mass_xf[ixfBin] =
        new TH1F(h_pair_mass_name.str().c_str(),
                 h_pair_mass_title.str().c_str(), 500, 0, 1);
  }

  // Histogram for the average bin value
  h_average_pt[0] = new TH1F(
      "h_average_pt_pi0",
      ";iPt; pT average [GeV/c]", nPtBins,
      0, nPtBins);
  h_average_eta[0] = new TH1F(
      "h_average_eta_pi0",
      ";ieta; #eta average", nEtaBins, 0,
      nEtaBins);
  h_average_xf[0] =
      new TH1F(
          "h_average_xf_pi0",
          ";ixf; pT average", nXfBins, 0, nXfBins);
  h_norm_pt[0] = new TH1F(
      "h_norm_pt_pi0",
      ";iPt; pT norm [GeV/c]", nPtBins,
      0, nPtBins);
  h_norm_eta[0] = new TH1F(
      "h_norm_eta_pi0",
      ";ieta; #eta norm", nEtaBins, 0,
      nEtaBins);
  h_norm_xf[0] =
      new TH1F(
          "h_norm_xf_pi0",
          ";ixf; pT norm", nXfBins, 0, nXfBins);
  h_average_pt[1] = new TH1F(
      "h_average_pt_eta",
      ";iPt; pT average [GeV/c]", nPtBins,
      0, nPtBins);
  h_average_eta[1] = new TH1F(
      "h_average_eta_eta",
      ";ieta; #eta average", nEtaBins, 0,
      nEtaBins);
  h_average_xf[1] =
      new TH1F(
          "h_average_xf_eta",
          ";ixf; pT average", nXfBins, 0, nXfBins);
  h_norm_pt[1] = new TH1F(
      "h_norm_pt_eta",
      ";iPt; pT norm [GeV/c]", nPtBins,
      0, nPtBins);
  h_norm_eta[1] = new TH1F(
      "h_norm_eta_eta",
      ";ieta; #eta norm", nEtaBins, 0,
      nEtaBins);
  h_norm_xf[1] =
      new TH1F(
          "h_norm_xf_eta",
          ";ixf; pT norm", nXfBins, 0, nXfBins);

  // Beam- spin- and kinematic-dependent yields -> pT dependent
  for (int iB = 0; iB < nBeams; iB++)
  {
    for (int iP = 0; iP < nParticles; iP++)
    {
      for (int iR = 0; iR < nRegions; iR++)
      {
        for (int iS = 0; iS < nSpins; iS++)
        {
          for (int iPt = 0; iPt < nPtBins; iPt++)
          {
            for (int ieta = 0; ieta < nEtaRegions; ieta++)
            {
              // low |eta| vs high |eta|
              std::stringstream h_yield_name;
              h_yield_name << "h_yield_" << beams[iB] << "_" << particle[iP]
                           << "_" << regions[iR] << "_"
                           << "pT_" << iPt << "_"
                           << etaRegions[ieta]
                           << "_" << spins[iS];
              std::stringstream h_yield_title;
              h_yield_title << std::fixed << std::setprecision(1)
                            << "P_{T} #in [" << pTBins[iPt] << ", "
                            << pTBins[iPt + 1] << "] "
                            << (ieta == 0 ? "|#eta| < 0.35"
                                          : "|#eta| > 0.35")
                            << " " << (iP == 0 ? "(#pi^{0} " : "(#eta ")
                            << regions[iR] << " range);"
                            << "#phi_{" << (iB == 0 ? "yellow" : "blue")
                            << "};counts";
              h_yield_1[iB][iP][iR][iPt][ieta][iS] =
                  new TH1F(h_yield_name.str().c_str(),
                           h_yield_title.str().c_str(), 12, -M_PI, M_PI);
              for (int iDir = 0; iDir < nDirections; iDir++)
              {
                if (ieta == 0)
                {
                  // forward-going vs backward-going
                  h_yield_name.str(
                      "");
                  h_yield_name << "h_yield_" << beams[iB] << "_" << particle[iP]
                               << "_" << regions[iR] << "_"
                               << "pT_" << iPt << "_"
                               << directions[iDir]
                               << "_" << spins[iS];
                  h_yield_title.str(
                      "");
                  h_yield_title << std::fixed << std::setprecision(1)
                                << "P_{T} #in [" << pTBins[iPt] << ", "
                                << pTBins[iPt + 1] << "] " << directions[iDir]
                                << " " << (iP == 0 ? "(#pi^{0} " : "(#eta ")
                                << regions[iR] << " range);"
                                << "#phi_{" << (iB == 0 ? "yellow" : "blue")
                                << "};counts";
                  h_yield_2[iB][iP][iR][iPt][iDir][iS] =
                      new TH1F(h_yield_name.str().c_str(),
                               h_yield_title.str().c_str(), 12, -M_PI, M_PI);
                }
                // low |eta| vs high |eta| AND forward-going vs backward-going
                h_yield_name.str(
                    "");
                h_yield_name << "h_yield_" << beams[iB] << "_" << particle[iP]
                             << "_" << regions[iR] << "_"
                             << "pT_" << iPt << "_"
                             << etaRegions[ieta]
                             << "_" << directions[iDir] << "_" << spins[iS];
                h_yield_title.str(
                    "");
                h_yield_title << std::fixed << std::setprecision(1)
                              << "P_{T} #in [" << pTBins[iPt] << ", "
                              << pTBins[iPt + 1] << "] "
                              << (ieta == 0 ? "|#eta| < 0.35"
                                            : "|#eta| > 0.35")
                              << " " << directions[iDir] << " "
                              << (iP == 0 ? "(#pi^{0} "
                                          : "(#eta ")
                              << regions[iR] << " range);"
                              << "#phi_{" << (iB == 0 ? "yellow" : "blue")
                              << "};counts";
                h_yield_3[iB][iP][iR][iPt][ieta][iDir][iS] =
                    new TH1F(h_yield_name.str().c_str(),
                             h_yield_title.str().c_str(), 12, -M_PI, M_PI);
              }
            }
          }
        }
      }
    }
  }

  // Beam- spin- and kinematic-dependent yields -> eta-dependent
  for (int iB = 0; iB < nBeams; iB++)
  {
    for (int iP = 0; iP < nParticles; iP++)
    {
      for (int iR = 0; iR < nRegions; iR++)
      {
        for (int iS = 0; iS < nSpins; iS++)
        {
          for (int ietaBin = 0; ietaBin < nEtaBins; ietaBin++)
          {
            std::stringstream h_yield_name;
            h_yield_name << "h_yield_" << beams[iB] << "_" << particle[iP]
                         << "_" << regions[iR] << "_"
                         << "eta_" << (int) ietaBin << "_" << spins[iS];
            std::stringstream h_yield_title;
            h_yield_title << std::fixed << std::setprecision(1) << "#eta #in ["
                          << etaBins[ietaBin] << ", " << etaBins[ietaBin]
                          << "] " << (iP == 0 ? "(#pi^{0} " : "(#eta ")
                          << regions[iR] << " range);"
                          << "#phi_{" << (iB == 0 ? "yellow" : "blue")
                          << "};counts";
            h_yield_eta[iB][iP][iR][ietaBin][iS] =
                new TH1F(h_yield_name.str().c_str(),
                         h_yield_title.str().c_str(), 12, -M_PI, M_PI);
          }
        }
      }
    }
  }

  // Beam- spin- and kinematic-dependent yields -> xf-dependent
  for (int iB = 0; iB < nBeams; iB++)
  {
    for (int iP = 0; iP < nParticles; iP++)
    {
      for (int iR = 0; iR < nRegions; iR++)
      {
        for (int iS = 0; iS < nSpins; iS++)
        {
          for (int ixfBin = 0; ixfBin < nXfBins; ixfBin++)
          {
            std::stringstream h_yield_name;
            h_yield_name << "h_yield_" << beams[iB] << "_" << particle[iP]
                         << "_" << regions[iR] << "_"
                         << "xf_" << (int) ixfBin << "_" << spins[iS];
            std::stringstream h_yield_title;
            h_yield_title << std::fixed << std::setprecision(1) << "x_{F} #in ["
                          << xfBins[ixfBin] << ", " << xfBins[ixfBin + 1]
                          << "] " << (iP == 0 ? "(#pi^{0} " : "(#eta ")
                          << regions[iR] << " range);"
                          << "#phi_{" << (iB == 0 ? "yellow" : "blue")
                          << "};counts";
            h_yield_xf[iB][iP][iR][ixfBin][iS] =
                new TH1F(h_yield_name.str().c_str(),
                         h_yield_title.str().c_str(), 12, -M_PI, M_PI);
          }
        }
      }
    }
  }

  // Histogram to illustrate cluster level cuts
  chi2_cuts_labels.resize(n_chi2_cuts);
  ecore_cuts_labels.resize(n_ecore_cuts);

  for (int i = 0; i < n_chi2_cuts; i++)
  {
    std::stringstream label;
    label << std::fixed << std::setprecision(1) << chi2_cuts[i];
    chi2_cuts_labels[i] = label.str();
  }

  for (int i = 0; i < n_ecore_cuts; i++)
  {
    std::stringstream label;
    label << std::fixed << std::setprecision(2) << ecore_cuts[i] << " GeV";
    ecore_cuts_labels[i] = label.str();
  }

  // Only store cluster level cuts if several cuts are given in (E,chi2):

  if (n_chi2_cuts > 1 || n_ecore_cuts > 1)
  {
    h_cluster_level_cuts_total = new TH2F(
      "h_cluster_level_cuts_total",
      "Total number of Diphotons;Cluster #chi^{2} Cut;Min Cluster Ecore Cut",
      n_chi2_cuts, 0, (float) n_chi2_cuts, n_ecore_cuts, 0, (float) n_ecore_cuts);
    for (int il = 0; il < n_chi2_cuts; il++)
      h_cluster_level_cuts_total->GetXaxis()->SetBinLabel(
        il + 1, chi2_cuts_labels[il].c_str());
    for (int il = 0; il < n_ecore_cuts; il++)
      h_cluster_level_cuts_total->GetYaxis()->SetBinLabel(
        il + 1, ecore_cuts_labels[il].c_str());

    for (int iP = 0; iP < nParticles; iP++)
    {
      for (int iR = 0; iR < nRegions; iR++)
      {
        std::stringstream h_name;
        h_name << "h_cluster_level_cuts_" << particle[iP] << "_" << regions[iR];
        std::stringstream h_title;
        h_title << "Total Number of Diphotons (" << (iP == 0 ? "#pi^{0}" : "#eta")
                << " " << regions[iR]
                << ");Cluster #chi^{2} Cut;Min Cluster Energy Cut";
        h_cluster_level_cuts_particle[iP][iR] =
          new TH2F(h_name.str().c_str(), h_title.str().c_str(), n_chi2_cuts, 0,
                   (float) n_chi2_cuts, n_ecore_cuts, 0, (float) n_ecore_cuts);
        for (int il = 0; il < n_chi2_cuts; il++)
          h_cluster_level_cuts_particle[iP][iR]->GetXaxis()->SetBinLabel(
            il + 1, chi2_cuts_labels[il].c_str());
        for (int il = 0; il < n_ecore_cuts; il++)
          h_cluster_level_cuts_particle[iP][iR]->GetYaxis()->SetBinLabel(
            il + 1, ecore_cuts_labels[il].c_str());

        for (int iPt = 0; iPt < nPtBins; iPt++)
        {
          if ((iP == 0) && (iR == 0))
          {
            h_name.str("");
            h_name << "h_cluster_level_cuts_total_pt_" << iPt;
            h_title.str("");
            h_title << "Total Number of Diphotons (p_{T} #in [" << pTBins[iPt]
                    << "," << pTBins[iPt + 1]
                    << "] GeV/c);Cluster #chi^{2} Cut;Min Cluster Energy Cut";
            h_cluster_level_cuts_total_pt[iPt] = new TH2F(
              h_name.str().c_str(), h_title.str().c_str(), n_chi2_cuts, 0,
              (float) n_chi2_cuts, n_ecore_cuts, 0, (float) n_ecore_cuts);
            for (int il = 0; il < n_chi2_cuts; il++)
              h_cluster_level_cuts_total_pt[iPt]->GetXaxis()->SetBinLabel(
                il + 1, chi2_cuts_labels[il].c_str());
            for (int il = 0; il < n_ecore_cuts; il++)
              h_cluster_level_cuts_total_pt[iPt]->GetYaxis()->SetBinLabel(
                il + 1, ecore_cuts_labels[il].c_str());
          }
          h_name.str(
                     "");
          h_name << "h_cluster_level_cuts_" << particle[iP] << "_" << regions[iR]
                 << "_pt_" << iPt;
          h_title.str(
                      "");
          h_title << "Total Number of Diphotons ("
                  << (iP == 0 ? "#pi^{0}" : "#eta")
                  << " " << regions[iR]
                  << ", p_{T} #in [" << pTBins[iPt] << "," << pTBins[iPt + 1]
                  << "] GeV/c);Cluster #chi^{2} Cut;Min Cluster Energy Cut";
          h_cluster_level_cuts_particle_pt[iP][iR][iPt] = new TH2F(
            h_name.str().c_str(), h_title.str().c_str(), n_chi2_cuts, 0,
            (float) n_chi2_cuts, n_ecore_cuts, 0, (float) n_ecore_cuts);
          for (int il = 0; il < n_chi2_cuts; il++)
            h_cluster_level_cuts_particle_pt[iP][iR][iPt]
              ->GetXaxis()
              ->SetBinLabel(il + 1, chi2_cuts_labels[il].c_str());
          for (int il = 0; il < n_ecore_cuts; il++)
            h_cluster_level_cuts_particle_pt[iP][iR][iPt]
              ->GetYaxis()
              ->SetBinLabel(il + 1, ecore_cuts_labels[il].c_str());
        }
      }
    }
  }

  return 0;
}

int AnNeutralMeson_micro_dst::InitRun(PHCompositeNode *topNode)
{
  syncobject = findNode::getClass<SyncObject>(topNode, syncdefs::SYNCNODENAME);
  if (!syncobject)
  {
    std::cout << "syncobject not found" << std::endl;
  }
  
  towergeom =
    findNode::getClass<RawTowerGeomContainer>(topNode,
                                              "TOWERGEOM_CEMC");

  if (require_emulator_matching)
  {
    emcaltiles =
      findNode::getClass<TriggerTile>(topNode,
                                    "TILES_EMCAL");
    if (!emcaltiles)
    {
      std::cerr << "AnNeutralMeson_micro_dst TILES_EMCAL node is missing for emulator matching" << "\n";
      std::cerr << "Either load DST_TRIGGER_EMULATOR files with Fun4AllDstInputManager" << "\n";
      std::cerr << "Or call method set_photon_trigger_emulator_matching_requirement(false)" << std::endl;
    }
  }
  
  vertexmap =
    findNode::getClass<GlobalVertexMap>(topNode,
                                        "GlobalVertexMap");
  if (!vertexmap)
  {
      std::cerr << "AnNeutralMeson_micro_dst GlobalVertexMap node is missing"
                << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
  }
  
  _smallclusters = findNode::getClass<ClusterSmallInfoContainer>(topNode, "CLUSTER_SMALLINFO_CEMC");
  if (!_smallclusters)
  {
    std::cerr << "AnNeutralMeson_micro_dst CLUSTER_SMALLINFO_CEMC node is missing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Check chi2 and ecore cuts are non-empty
  if (n_chi2_cuts == 0 || n_ecore_cuts == 0)
  {
    std::cerr << "No chi2 or no ecore cuts defined. Abort run" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // See how the trigger logic is defined for this run in the CDB.
  odbc::Connection *con = nullptr;
  try
  {
    con = odbc::DriverManager::getConnection(db_name, user_name, "");
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  std::stringstream cmd;
  cmd << "SELECT * FROM " << table_name << " WHERE runnumber=" << runnumber << ";";
  odbc::Statement *stmt = con->createStatement();
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (rs->next() == 0)
  {
    std::cout << std::format("Error: Can't find GL1 scaledown data for run {}", runnumber) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  int ncols = rs->getMetaData()->getColumnCount();
  for (int icol = 0; icol < ncols; icol++)
  {
    std::string col_name = rs->getMetaData()->getColumnName(icol + 1);
    // i.e. scaledownXX where XX is the trigger number (00 to 63)
    if (col_name.starts_with("scaledown")) 
    {
      int scaledown_index = std::stoi(col_name.substr(9,2));
      std::string cont = rs->getString(col_name);
      scaledown[scaledown_index] = std::stoi(cont);
    }
  }


  // Get livetime fraction
  for (int index = 0; index < 64; index++)
  {
    cmd.str("");
    cmd << "SELECT raw, live FROM gl1_scalers WHERE runnumber=" << runnumber << " AND index=" << index << ";";
    stmt = nullptr; stmt = con->createStatement();
    rs = nullptr;
    try
    {
      rs = stmt->executeQuery(cmd.str());
    }
    catch (odbc::SQLException &e)
    {
      std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
      delete rs;
      delete stmt;
      delete con;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    if (rs->next() == 0)
    {
      delete rs;
      delete stmt;
      livetime[index] = 0;
      continue;
    }

    double total_live = std::stod(rs->getString("live"));
    double total_raw = std::stod(rs->getString("raw"));
    livetime[index] = total_live / total_raw;
  }

  // Spin pattern from SpinDB
  // Get spin pattern from SpinDB
  // Use the default QA level
  //unsigned int qa_level = 0xffff;  // or equivalently 65535 in decimal notation -> Default
  SpinDBOutput spin_out("phnxrc");
  SpinDBContent *spin_cont = new SpinDBContentv1();
  spin_out.StoreDBContent(runnumber, runnumber);
  spin_out.GetDBContentStore(spin_cont, runnumber);

  crossingshift = spin_cont->GetCrossingShift();

  for (int i = 0; i < nBunches; i++)
  {
    beamspinpat[0][i] = -1 * spin_cont->GetSpinPatternYellow(i);
    beamspinpat[1][i] = -1 * spin_cont->GetSpinPatternBlue(i);
    // spin pattern at sPHENIX is -1*(CDEV pattern)
    // The spin pattern corresponds to the expected spin at IP12, before the
    // siberian snake
  }

  _eventcounter = -1;
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnNeutralMeson_micro_dst::process_event(PHCompositeNode *topNode)
{
  _eventcounter++;

  if (_eventcounter % 100'000 == 0)
  {
    monitorMemoryUsage("check");
    std::cout << "event: " << _eventcounter << std::endl;
  }

  // Store event-level entries
  live_trigger = _smallclusters->get_live_trigger();
  scaled_trigger = _smallclusters->get_scaled_trigger();
  diphoton_bunchnumber = _smallclusters->get_bunch_number();
  cluster_number = _smallclusters->size();

  if (vertexmap && !vertexmap->empty())
  {
    GlobalVertex *vtx = vertexmap->begin()->second;
    if (vtx)
    {
      // vtx_x and vtx_y always estimated to 0
      vertex_z = vtx->get_z();
    }
    else
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  else
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (store_qa)
  {
  // Trigger count per bit (QA)
    for (const auto &[trigger_number,
                        trigger_name] : trigger_names)
    {
      if ((live_trigger >> trigger_number & 0x1U) == 0x1U)
      {
        h_trigger_live->Fill(trigger_number + 1);
      }
      if ((scaled_trigger >> trigger_number & 0x1U) == 0x1U)
      {
        h_trigger_scaled->Fill(trigger_number + 1);
      }
    }
  }

  // Vertex cut
  if (std::abs(vertex_z) > vertex_max)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  diphoton_vertex_z = vertex_z;

  if (store_qa) h_event_zvtx->Fill(vertex_z);

  // Additional trigger selection
  if (require_mbd_trigger_bit && !require_photon_trigger_bit)
  {
    // Minimum Bias Event
    if (!mbd_trigger_bit())
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    mbd_trigger_bit_event = true;
  }
  else if (require_photon_trigger_bit && !require_mbd_trigger_bit)
  {
    // Photon Trigger Event
    if (!photon_trigger_bit())
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    photon_trigger_bit_event = true;
  }
  else { // If both triggers are used, apply a pT threshold between them
    mbd_trigger_bit_event = mbd_trigger_bit();
    photon_trigger_bit_event = photon_trigger_bit();
  }

  if (store_qa) {
    h_triggered_event_photons->Fill(cluster_number);
    
    // Fill vertex info, for QA
    h_triggered_event_zvtx->Fill(vertex_z);
  }

  // If more than one cut for ecore and chi2 is given, apply the following method for cluster number comparison (QA only)
  if (n_chi2_cuts > 1 || n_ecore_cuts > 1)
  {
    cluster_cuts();
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // Select the "good" photon clusters
  good_photons.clear();
  num_photons = 0;
  for (int iclus = 0; iclus < cluster_number; iclus++)
  {
    ClusterSmallInfo *smallcluster = _smallclusters->get_cluster_at(iclus);
    float eta = smallcluster->get_eta();
    float phi = smallcluster->get_phi();
    float chi2 = smallcluster->get_chi2();
    float ecore = smallcluster->get_ecore();

    // Smear the photon's kinematics based on the EMCal given resolution
    if (do_smearing)
    {
      float eta_nominal = eta;
      float phi_nominal = phi;
      float E_nominal = ecore;
      smear_photon(eta, phi, ecore);
      
      h_smear_E->Fill(ecore/E_nominal);
      h_smear_eta->Fill(eta - eta_nominal);
      h_smear_phi->Fill(phi - phi_nominal);

      h_smear_E_dE->Fill(E_nominal, ecore/E_nominal);
      h_smear_E_deta->Fill(E_nominal, eta - eta_nominal);
      h_smear_E_dphi->Fill(E_nominal, phi - phi_nominal);
    }
    
    if (chi2 <= chi2_cuts[0] &&
        ecore >= ecore_cuts[0])
    {
      ROOT::Math::PtEtaPhiMVector photon(ecore / std::cosh(eta), eta, phi, 0);
      if (store_qa)
      {
        h_photon_eta->Fill(photon.Eta());
        h_photon_phi->Fill(photon.Phi());
        h_photon_pt->Fill(photon.Pt());
        h_photon_zvtx->Fill(vertex_z);
        h_photon_eta_phi->Fill(photon.Eta(), photon.Phi());
        int iPt = FindBinBinary(diphoton_pt, pTBins, nPtBins + 1);
        if (!(iPt < 0 || iPt >= nPtBins)) {
          h_photon_eta_phi_bin[iPt]->Fill(photon.Eta(), photon.Phi());
        }
        h_photon_eta_pt->Fill(photon.Eta(), photon.Pt());
        h_photon_eta_zvtx->Fill(photon.Eta(), vertex_z);
        h_photon_phi_pt->Fill(photon.Phi(), photon.Pt());
        h_photon_phi_zvtx->Fill(photon.Phi(), vertex_z);
        h_photon_pt_zvtx->Fill(photon.Pt(), vertex_z);
      }
      Cluster cluster(photon, false); // by default, a cluster does not activate the trigger
      good_photons.push_back(cluster);
      num_photons++;
    }
  }

  // Trigger emulator -> determine what cluster(s) fired the trigger (cluster.isTrigger = true)
  if (require_emulator_matching)
  {
    trigger_emulator_input();
  }

  // Select the photon pairs
  int multiplicity_efficiency = 0;
  int multiplicity_emulator = 0;
  bool first_diphoton_nomatch = true;
  for (int iclus1 = 0; iclus1 < num_photons - 1; iclus1++)
  {
    ROOT::Math::PtEtaPhiMVector photon1 = good_photons[iclus1].p4;
    for (int iclus2 = iclus1 + 1; iclus2 < num_photons; iclus2++)
    {
      ROOT::Math::PtEtaPhiMVector photon2 = good_photons[iclus2].p4;
      
      ROOT::Math::PtEtaPhiMVector diphoton = photon1 + photon2;

      diphoton_pt = diphoton.Pt();

      float eta_diff = photon1.Eta() - photon2.Eta();
      float phi_diff = WrapAngleDifference(photon1.Phi(), photon2.Phi());
      if (store_qa) h_pair_DeltaR->Fill(std::sqrt(std::pow(eta_diff, 2) + std::pow(phi_diff, 2)));

      // General diphoton cut
      if (diphoton_cut(photon1, photon2, diphoton)) continue;

      if (store_qa && diphoton.M() > 0.25 && diphoton.M() < 0.38)
      {
        h_same_delta_eta->Fill(std::abs(photon1.Eta() - photon2.Eta()));
        h_same_delta_phi->Fill(std::abs(ROOT::Math::VectorUtil::DeltaPhi(photon1, photon2)));
        h_same_delta_R->Fill(ROOT::Math::VectorUtil::DeltaR(photon1, photon2));
        h_same_alpha->Fill(std::abs(photon1.E() - photon2.E())/(photon1.E() + photon2.E()));
      }

      if (store_qa && first_diphoton_nomatch)
      {
        h_count_diphoton_nomatch->Fill(1);
        first_diphoton_nomatch = false;
      }

      // Distinct diphoton selection between MBD and Photon Trigger

      // For a photon-triggered event
      // Keep candidates at high-pT which satisfy the trigger matching criterion

      // Trigger emulator matching cut
      if (require_emulator_matching && !require_efficiency_matching)
      {
        emulator_match = (good_photons[iclus1].isTrigger || good_photons[iclus2].isTrigger);
        if (!(emulator_match))
        {
          continue;
        }
        multiplicity_emulator++;
      }

      // Trigger efficiency matching cut ("Proxy method")
      if (require_efficiency_matching && !require_emulator_matching)
      {
        efficiency_match = trigger_efficiency_matching(photon1, photon2, diphoton);
        if (!(efficiency_match))
        {
          continue;
        }
      }

      // Emulator vs proxy comparison
      if (require_efficiency_matching && require_emulator_matching)
      {
        int iPt = FindBinBinary(diphoton_pt, pTBins, nPtBins + 1);
        if (iPt < 0 || iPt >= nPtBins) continue;
        emulator_match = (good_photons[iclus1].isTrigger || good_photons[iclus2].isTrigger);
        efficiency_match = trigger_efficiency_matching(photon1, photon2, diphoton);
        if (emulator_match && efficiency_match)
        {
          h_efficiency_x_matching_pt[iPt]->AddBinContent(1);
          h_efficiency_x_matching->AddBinContent(1);
        }
        if (!emulator_match && !efficiency_match)
        {
          h_efficiency_x_matching_pt[iPt]->AddBinContent(2);
          h_efficiency_x_matching->AddBinContent(2);
        }
        if (emulator_match && !efficiency_match)
        {
          h_efficiency_x_matching_pt[iPt]->AddBinContent(3);
          h_efficiency_x_matching->AddBinContent(3);
        }
        if (!emulator_match && efficiency_match) {
          h_efficiency_x_matching_pt[iPt]->AddBinContent(4);
          h_efficiency_x_matching->AddBinContent(4);
        }
        continue;
      }

      if (store_qa)
      {
        h_pair_alpha->Fill(std::abs(photon1.E() - photon2.E())/(photon1.E() + photon2.E()));
        h_pair_alpha_pt->Fill(std::abs(photon1.E() - photon2.E())/(photon1.E() + photon2.E()), diphoton_pt);
      }
        
      if (multiplicity_efficiency == 1) {
        h_analysis_event_zvtx->Fill(vertex_z);
      }

      // Store diphoton properties
      diphoton_mass = diphoton.mag();
      diphoton_eta = diphoton.Eta();
      diphoton_phi = diphoton.Phi();
      diphoton_xf = 2 * diphoton.Pz() / Vs;

      if (store_qa)
      {
        h_selected_photon_eta->Fill(photon1.Eta());
        h_selected_photon_phi->Fill(photon1.Phi());
        h_selected_photon_pt->Fill(photon1.Pt());
        h_selected_photon_zvtx->Fill(vertex_z);
        h_selected_photon_eta_phi->Fill(photon1.Eta(), photon1.Phi());
        h_selected_photon_eta_pt->Fill(photon1.Eta(), photon1.Pt());
        h_selected_photon_eta_zvtx->Fill(photon1.Eta(), vertex_z);
        h_selected_photon_phi_pt->Fill(photon1.Phi(), photon1.Pt());
        h_selected_photon_phi_zvtx->Fill(photon1.Phi(), vertex_z);
        h_selected_photon_pt_zvtx->Fill(photon1.Pt(), vertex_z);
        h_selected_photon_eta->Fill(photon2.Eta());
        h_selected_photon_phi->Fill(photon2.Phi());
        h_selected_photon_pt->Fill(photon2.Pt());
        h_selected_photon_zvtx->Fill(vertex_z);
        h_selected_photon_eta_phi->Fill(photon2.Eta(), photon2.Phi());
        h_selected_photon_eta_pt->Fill(photon2.Eta(), photon2.Pt());
        h_selected_photon_eta_zvtx->Fill(photon2.Eta(), vertex_z);
        h_selected_photon_phi_pt->Fill(photon2.Phi(), photon2.Pt());
        h_selected_photon_phi_zvtx->Fill(photon2.Phi(), vertex_z);
        h_selected_photon_pt_zvtx->Fill(photon2.Pt(), vertex_z);

        int iPt = FindBinBinary(diphoton_pt, pTBins, nPtBins + 1);
        if (!(iPt < 0 || iPt >= nPtBins)) {
          h_selected_photon_eta_phi_bin[iPt]->Fill(photon1.Eta(), photon1.Phi());
          h_selected_photon_eta_phi_bin[iPt]->Fill(photon2.Eta(), photon2.Phi());
        }
        
        h_pair_eta->Fill(diphoton_eta);
        h_pair_phi->Fill(diphoton_phi);
        h_pair_pt->Fill(diphoton_pt);
        h_pair_zvtx->Fill(vertex_z);
        h_pair_eta_phi->Fill(diphoton_eta, diphoton_phi);
        h_pair_eta_pt->Fill(diphoton_eta, diphoton_pt);
        h_pair_eta_zvtx->Fill(diphoton_eta, vertex_z);
        h_pair_phi_pt->Fill(diphoton_phi, diphoton_pt);
        h_pair_phi_zvtx->Fill(diphoton_phi, vertex_z);
        h_pair_pt_zvtx->Fill(diphoton_pt, vertex_z);
      }

      // Select the pt bin:
      int iPt = FindBinBinary(diphoton_pt, pTBins, nPtBins + 1);

      // Select the zvtx bin:
      int izvtx = FindBinBinary(std::abs(vertex_z), zvtxBins, nZvtxBins + 1);

      // Select the eta bin:
      int ietaBin = FindBinBinary(diphoton_eta, etaBins, nEtaBins + 1);

      // Select the xf bin:
      int ixfBin = FindBinBinary(diphoton_xf, xfBins, nXfBins + 1);

      // Select the eta region (small eta vs high eta)
      int ietaRegion = (std::abs(diphoton_eta) < etaThreshold ? 0 : 1);

      // Invariant mass distributions
      h_pair_mass->Fill(diphoton_mass);
      if (izvtx >= 0 && izvtx < nZvtxBins)
        for (int izvtxtmp = izvtx; izvtxtmp < nZvtxBins; izvtxtmp++)
          h_pair_mass_zvtx[izvtxtmp]->Fill(diphoton_mass);
      if (iPt >= 0 && iPt < nPtBins)
        h_pair_mass_pt[iPt]->Fill(diphoton_mass);
      if (ietaBin >= 0 && ietaBin < nEtaBins)
        h_pair_mass_eta[ietaBin]->Fill(diphoton_mass);
      if (ixfBin >= 0 && ixfBin < nXfBins)
        h_pair_mass_xf[ixfBin]->Fill(diphoton_mass);

      // Fill minimal output ttree
      // Can be used for future short analysis (see AnNeutralMeson_nano)
      if (store_tree)
      {
        output_tree->Fill();
      }

      // Select particle and region index;
      int iP = -1;
      int iR = -1;
      int ival = FindBinBinary(diphoton_mass, band_limits,
                               nParticles * (nRegions + 1) * 2);
      // Ignore any diphoton which is not within the side band or the peak band
      if (ival >= 0 && ival % 2 == 0)
      {
        iP = ival / 2 / (nRegions + 1);            // particle index
        iR = (ival / 2 % (nRegions + 1) + 1) % 2;  // region index
      }

      // Store pT, eta and xF values in order to compute the average value of
      // each bin
      if (iP != -1) 
      {
        if (iPt >= 0 && iPt < nPtBins) 
        {
          h_average_pt[iP]->Fill(iPt, diphoton_pt);
          h_norm_pt[iP]->Fill(iPt, 1);
        }
        if (ietaBin >= 0 && ietaBin < nEtaBins)
        {
          h_average_eta[iP]->Fill(ietaBin, diphoton_eta);
          h_norm_eta[iP]->Fill(ietaBin, 1);
        }
        if (ixfBin >= 0 && ixfBin < nXfBins)
        {
          h_average_xf[iP]->Fill(ixfBin, diphoton_xf);
          h_norm_xf[iP]->Fill(ixfBin, 1);
        }
      }

      // Get Kinematic distributions:
      if (store_qa)
      {
        if (iP == 0 && iR == 0)
        {
          h_meson_pi0_pt->Fill(diphoton_pt);
          h_meson_pi0_E->Fill(diphoton.E());
          if (1 < diphoton_pt && diphoton_pt < 2)
          {
            h_pair_pi0_eta_pt_1->Fill(diphoton_eta);
            h_pair_pi0_xf_pt_1->Fill(diphoton_xf);
          }
          else if (2 < diphoton_pt && diphoton_pt < 3)
          {
            h_pair_pi0_eta_pt_2->Fill(diphoton_eta);
            h_pair_pi0_xf_pt_2->Fill(diphoton_xf);
          }
          else if (3 < diphoton_pt && diphoton_pt < 4)
          {
            h_pair_pi0_eta_pt_3->Fill(diphoton_eta);
            h_pair_pi0_xf_pt_3->Fill(diphoton_xf);
          }
          else if (4 < diphoton_pt && diphoton_pt < 5)
          {
            h_pair_pi0_eta_pt_4->Fill(diphoton_eta);
            h_pair_pi0_xf_pt_4->Fill(diphoton_xf);
          }
          else if (5 < diphoton_pt && diphoton_pt < 6)
          {
            h_pair_pi0_eta_pt_5->Fill(diphoton_eta);
            h_pair_pi0_xf_pt_5->Fill(diphoton_xf);
          }
          else if (6 < diphoton_pt && diphoton_pt < 7)
          {
            h_pair_pi0_eta_pt_6->Fill(diphoton_eta);
            h_pair_pi0_xf_pt_6->Fill(diphoton_xf);
          }
          else if (7 < diphoton_pt && diphoton_pt < 8)
          {
            h_pair_pi0_eta_pt_7->Fill(diphoton_eta);
            h_pair_pi0_xf_pt_7->Fill(diphoton_xf);
          }
          else if (8 < diphoton_pt && diphoton_pt < 10)
          {
            h_pair_pi0_eta_pt_8->Fill(diphoton_eta);
            h_pair_pi0_xf_pt_8->Fill(diphoton_xf);
          }
          else if (10 < diphoton_pt && diphoton_pt < 20)
          {
            h_pair_pi0_eta_pt_9->Fill(diphoton_eta);
            h_pair_pi0_xf_pt_9->Fill(diphoton_xf);
          }
        }

        if (iP == 1 && iR == 0)
        {
          h_meson_eta_pt->Fill(diphoton_pt);
          h_meson_eta_E->Fill(diphoton.E());
          if (1 < diphoton_pt && diphoton_pt < 2)
          {
            h_pair_eta_eta_pt_1->Fill(diphoton_eta);
            h_pair_eta_xf_pt_1->Fill(diphoton_xf);
          }
          else if (2 < diphoton_pt && diphoton_pt < 3)
          {
            h_pair_eta_eta_pt_2->Fill(diphoton_eta);
            h_pair_eta_xf_pt_2->Fill(diphoton_xf);
          }
          else if (3 < diphoton_pt && diphoton_pt < 4)
          {
            h_pair_eta_eta_pt_3->Fill(diphoton_eta);
            h_pair_eta_xf_pt_3->Fill(diphoton_xf);
          }
          else if (4 < diphoton_pt && diphoton_pt < 5)
          {
            h_pair_eta_eta_pt_4->Fill(diphoton_eta);
            h_pair_eta_xf_pt_4->Fill(diphoton_xf);
          }
          else if (5 < diphoton_pt && diphoton_pt < 6)
          {
            h_pair_eta_eta_pt_5->Fill(diphoton_eta);
            h_pair_eta_xf_pt_5->Fill(diphoton_xf);
          }
          else if (6 < diphoton_pt && diphoton_pt < 7)
          {
            h_pair_eta_eta_pt_6->Fill(diphoton_eta);
            h_pair_eta_xf_pt_6->Fill(diphoton_xf);
          }
          else if (7 < diphoton_pt && diphoton_pt < 8)
          {
            h_pair_eta_eta_pt_7->Fill(diphoton_eta);
            h_pair_eta_xf_pt_7->Fill(diphoton_xf);
          }
          else if (8 < diphoton_pt && diphoton_pt < 10)
          {
            h_pair_eta_eta_pt_8->Fill(diphoton_eta);
            h_pair_eta_xf_pt_8->Fill(diphoton_xf);
          }
          else if (10 < diphoton_pt && diphoton_pt < 20)
          {
            h_pair_eta_eta_pt_9->Fill(diphoton_eta);
            h_pair_eta_xf_pt_9->Fill(diphoton_xf);
          }
        }
      }

      // Compute yield
      float phi_shift = 0;
      float phi_beam = 0;
      for (int iB = 0; iB < nBeams; iB++)
      {
        int iDir = (diphoton_eta * beamDirection[iB] > 0 ? 0 : 1);
        int ietaBinRelative =
            (beamDirection[iB] > 0 ? ietaBin : nEtaBins - 1 - ietaBin);
        int ixfBinRelative =
            (beamDirection[iB] > 0 ? ixfBin : nXfBins - 1 - ixfBin);

        // sanity check (eta, xf > 0 -> forward)
        if (((ietaBin >= 0 && ietaBin < nEtaBins) && (ixfBin >= 0 && ixfBin < nXfBins)) &&
            (((iDir == 0) && ((etaBins[ietaBinRelative + 1] <= 0) ||
                              (xfBins[ixfBinRelative + 1] <= 0))) ||
             ((iDir == 1) && ((etaBins[ietaBinRelative] >= 0) ||
                              (xfBins[ixfBinRelative] >= 0)))))
        {
          std::cerr << "Error: Bad definition of relative eta, xf bins.\n";
          std::cerr << (iDir == 0 ? "forward"
                                  : "backward")
                    << " -> ixf = " << ixfBinRelative << " ("
                    << xfBins[ixfBinRelative] << ", "
                    << xfBins[ixfBinRelative + 1]
                    << "), ieta = " << ietaBinRelative << " ("
                    << etaBins[ietaBinRelative] << ", "
                    << etaBins[ietaBinRelative + 1] << ")\n";
        }

        phi_shift = (iB == 0 ? M_PI / 2 : -M_PI / 2);
        phi_beam = diphoton_phi + phi_shift;  // phi global to phi yellow/blue
        WrapAngle(phi_beam);
        if (iR != -1 && iP != -1)
        {
          if (beamspinpat[iB][(crossingshift + diphoton_bunchnumber) % nBunches] == 1)
          {
            if (iPt >= 0 && iPt < nPtBins)
            {
              h_yield_1[iB][iP][iR][iPt][ietaRegion][0]->Fill(phi_beam);
              h_yield_2[iB][iP][iR][iPt][iDir][0]->Fill(phi_beam);
              h_yield_3[iB][iP][iR][iPt][ietaRegion][iDir][0]->Fill(phi_beam);
            }
            if (ietaBin >= 0 && ietaBin < nEtaBins)
            {
              h_yield_eta[iB][iP][iR < 2 ? iR : 1][ietaBinRelative][0]->Fill(phi_beam);
            }
            if (ixfBin >= 0 && ixfBin < nXfBins)
            {
              h_yield_xf[iB][iP][iR < 2 ? iR : 1][ixfBinRelative][0]->Fill(phi_beam);
            }
          }
          else if (beamspinpat[iB][(crossingshift + diphoton_bunchnumber) % nBunches] == -1)
          {
            if (iPt >= 0 && iPt < nPtBins)
            {
              h_yield_1[iB][iP][iR < 2 ? iR : 1][iPt][ietaRegion][1]->Fill(phi_beam);
              h_yield_2[iB][iP][iR < 2 ? iR : 1][iPt][iDir][1]->Fill(phi_beam);
              h_yield_3[iB][iP][iR < 2 ? iR : 1][iPt][ietaRegion][iDir][1]->Fill(phi_beam);
            }
            if (ietaBin >= 0 && ietaBin < nEtaBins)
            {
              h_yield_eta[iB][iP][iR < 2 ? iR : 1][ietaBinRelative][1]->Fill(phi_beam);
            }
            if (ixfBin >= 0 && ixfBin < nXfBins)
            {
              h_yield_xf[iB][iP][iR < 2 ? iR : 1][ixfBinRelative][1]->Fill(phi_beam);
            }
          }
        }
      }
    }
  }

  if (store_qa)
  {
    h_multiplicity_efficiency->Fill(multiplicity_efficiency);
    h_multiplicity_emulator->Fill(multiplicity_emulator);
  }

  // Event Mixing: If there is at least one good pair, store the event in the pool
  if (use_event_mixing)
  {
    if (require_mbd_trigger_bit)
    {
      event_mixing_mbd(topNode);
    }
    else if (require_photon_trigger_bit)
    {
      event_mixing_photon();
    }
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnNeutralMeson_micro_dst::End(PHCompositeNode *)
{
  outfile->cd();
  outfile->Write();
  outfile->Close();
  delete outfile;
  outfile = nullptr; // To prevent double deletion with the destructor
  if (store_tree)
  {
    outfile_tree->cd();
    outfile_tree->Write();
    outfile_tree->Close();
    delete outfile_tree;
    outfile_tree = nullptr;
  }
  return 0;
}

void AnNeutralMeson_micro_dst::cluster_cuts()
{
  // Define the good clusters with respect to the chi2 and energy (corrected) cuts
  // set before the execution
  std::vector<std::vector<ROOT::Math::PtEtaPhiMVector>> good_photons_by_cut(n_chi2_cuts * n_ecore_cuts);
  for (int iclus = 0; iclus < cluster_number; iclus++)
  {
    const ClusterSmallInfo *smallcluster = _smallclusters->get_cluster_at(iclus);
    float eta = smallcluster->get_eta();
    float phi = smallcluster->get_phi();
    float chi2 = smallcluster->get_chi2();
    float ecore = smallcluster->get_ecore();
    int ichi2cut_max = FindBinBinary(chi2, chi2_cuts.data(), n_chi2_cuts);
    int iecorecut_max = FindBinBinary(ecore, ecore_cuts.data(), n_ecore_cuts);
    ROOT::Math::PtEtaPhiMVector photon(ecore / eta, eta, phi, 0);
    for (int ichi2cut = 0; ichi2cut <= ichi2cut_max; ichi2cut++)
    {
      for (int iecorecut = 0; iecorecut <= iecorecut_max; iecorecut++)
      {
        good_photons_by_cut[ichi2cut * n_ecore_cuts + iecorecut].push_back(photon);
      }
    }
  }

  for (int ichi2cut = 0; ichi2cut < n_chi2_cuts; ichi2cut++)
  {
    for (int iecorecut = 0; iecorecut < n_ecore_cuts; iecorecut++)
    {
      num_photons = good_photons_by_cut[ichi2cut * n_ecore_cuts + iecorecut].size();
      for (int iclus1 = 0; iclus1 < num_photons; iclus1++)
      {
        ROOT::Math::PtEtaPhiMVector photon1 = good_photons_by_cut[ichi2cut * n_ecore_cuts + iecorecut][iclus1];
        for (int iclus2 = iclus1 + 1; iclus2 < num_photons; iclus2++)
        {
          ROOT::Math::PtEtaPhiMVector photon2 = good_photons_by_cut[ichi2cut * n_ecore_cuts + iecorecut][iclus2];
          ROOT::Math::PtEtaPhiMVector diphoton = photon1 + photon2;

          if (diphoton_cut(photon1, photon2, diphoton)) continue;

          diphoton_mass = diphoton.mag();

          // Select particle and region index;
          int iP = -1;
          int iR = -1;
          int ival = FindBinBinary(diphoton_mass, band_limits,
                                   nParticles * (nRegions + 1) * 2);
          // Ignore any diphoton which is not within the side band or the peak band
          if (ival >= 0 && ival % 2 == 0)
          {
            iP = ival / 2 / (nRegions + 1);            // particle index
            iR = (ival / 2 % (nRegions + 1) + 1) % 2;  // region index
          }

          int iPt = FindBinBinary(diphoton_pt, pTBins, nPtBins + 1);

          if (iPt >= 0 && iPt < nPtBins)
            h_cluster_level_cuts_total_pt[iPt]
                ->Fill((float) ichi2cut + 0.5, (float) iecorecut + 0.5);
          // Check within pi0 and eta mass range
          if (iP != -1 && iR != -1)
          {
            // OK
            h_cluster_level_cuts_particle[iP][iR]
              ->Fill((float) ichi2cut + 0.5, (float) iecorecut + 0.5);
            if (iPt >= 0 && iPt < nPtBins)
              h_cluster_level_cuts_particle_pt[iP][iR][iPt]
                ->Fill((float) ichi2cut + 0.5, (float) iecorecut + 0.5);
          }
        }
      }
    }
  }
}

bool AnNeutralMeson_micro_dst::diphoton_cut(ROOT::Math::PtEtaPhiMVector p1,
                                            ROOT::Math::PtEtaPhiMVector p2,
                                            ROOT::Math::PtEtaPhiMVector ppair)
{
  float alpha = std::abs(p1.E() - p2.E()) / (p1.E() + p2.E());
  float pt = ppair.Pt();
  return (alpha > alphaCut || pt < pTCutMin || pt > pTCutMax);
}

bool AnNeutralMeson_micro_dst::mbd_trigger_bit()
{
  trigger_mbd_any_vtx = ((scaled_trigger >> 10 & 0x1U) == 0x1U);
  trigger_mbd_vtx_10 = ((scaled_trigger >> 12 & 0x1U) == 0x1U);
  trigger_mbd = (trigger_mbd_any_vtx || trigger_mbd_vtx_10);

  return trigger_mbd;
}

bool AnNeutralMeson_micro_dst::photon_trigger_bit()
{
  trigger_mbd_photon_3 = false;
  trigger_mbd_photon_4 = false;
  
  // If (photon 3 GeV + MBD NS >= 1) is enabled for recording
  if ((scaledown[25] != -1) &&
      ((scaled_trigger >> 25 & 0x1U) == 0x1U))
  {
    trigger_mbd_photon_3 = true;
  }
  // Else if (photon 3 GeV + MBD NS >=1, vtx < 10) is enabled for recording
  if ((scaledown[36] != -1) &&
      ((scaled_trigger >> 36 & 0x1U) == 0x1U))
  {
    trigger_mbd_photon_3 = true;
  }
  // Else if (photon 3 GeV) is enabled for recording with no prescale
  if ((scaledown[29] == 0) &&
      (((live_trigger >> 25 & 0x1U) == 0x1U) ||
       ((live_trigger >> 36 & 0x1U) == 0x1U)))
  {
    if ((live_trigger >> 25 & 0x1U) == 0x1U)
    {
      trigger_mbd_photon_3 = true;
    }
    if ((live_trigger >> 36 & 0x1U) == 0x1U)
    {
      trigger_mbd_photon_3 = true;
    }
  }

  // If (photon 4 GeV + MBD NS >= 1) is enabled for recording
  if ((scaledown[26] != -1) &&
      ((scaled_trigger >> 26 & 0x1U) == 0x1U))
  {
    trigger_mbd_photon_4 = true;
  }
  // Else if (photon 4 GeV + MBD NS >=1, vtx < 10) is enabled for recording
  if ((scaledown[37] != -1) &&
      ((scaled_trigger >> 37 & 0x1U) == 0x1U))
  {
    trigger_mbd_photon_4 = true;
  }
  // Else if (photon 4 GeV) is enabled for recording with no prescale
  if ((scaledown[30] == 0) &&
      (((live_trigger >> 26 & 0x1U) == 0x1U) ||
       ((live_trigger >> 37 & 0x1U) == 0x1U)))
  {
    if ((live_trigger >> 26 & 0x1U) == 0x1U)
    {
      trigger_mbd_photon_4 = true;
    }
    if ((live_trigger >> 37 & 0x1U) == 0x1U)
    {
      trigger_mbd_photon_4 = true;
    }
  }
  return (trigger_mbd_photon_3 || trigger_mbd_photon_4);
}

void AnNeutralMeson_micro_dst::trigger_emulator_input()
{
  // 1) Identify the tile(s) (8x8 tower windows) which fired the trigger
  // And if so, to what energy: 3 GeV or 4 GeV
  // In most cases, only one tile fired
  emulator_selection = 0;
  fired_indices.clear();
  
  tileNumber = emcaltiles->get_tile_number();
  if (trigger_mbd_photon_3)
  {
    for (int itile = 0; itile < tileNumber; itile++)
    {
      if (emcaltiles->get_tile_energy_adc(itile) >= adc_threshold_3)
      {
        fired_indices.push_back(itile);
        emulator_selection++;
      }
    }
  }
  else if (trigger_mbd_photon_4)
  {
    for (int itile = 0; itile < tileNumber; itile++)
    {
      if (emcaltiles->get_tile_energy_adc(itile) >= adc_threshold_4)
      {
        fired_indices.push_back(itile);
        emulator_selection++;
      }
    }
  }
  h_emulator_selection->Fill(emulator_selection);

  // Identify all the good photon clusters
  // located in the firing tile
  // And store them by decreasing energy
  emulator_energies.clear();
  emulator_clusters.clear();
  if (emulator_selection) // i.e. emulator detects at least one firing tile
  {
    for (int i = 0; i < emulator_selection; i++)
    {
      // Tile dimensions
      int itile = fired_indices[i];
      int ieta_min = emcaltiles->get_tile_eta(itile) * 8;
      int ieta_max = (emcaltiles->get_tile_eta(itile) + 1) * 8 - 1;
      int iphi_min = emcaltiles->get_tile_phi(itile) * 8;
      int iphi_max = (emcaltiles->get_tile_phi(itile) + 1) * 8 - 1;
      double eta_min = towergeom->get_etabounds(std::max(0,ieta_min)).first;
      double eta_max = towergeom->get_etabounds(std::min(ieta_max,95)).second;
      double z_min = radius * std::sinh(eta_min);
      double z_max = radius * std::sinh(eta_max);
      double phi_min = towergeom->get_phibounds(std::max(0,iphi_min)).first;
      double phi_max = towergeom->get_phibounds(std::min(255,iphi_max)).second;
      phi_min = (phi_min >= 0 ? phi_min : phi_min + 2 * M_PI);
      phi_max = (phi_max >= 0 ? phi_max : phi_max + 2 * M_PI);
      for (int iclus1 = 0; iclus1 < num_photons; iclus1++)
      {
        ROOT::Math::PtEtaPhiMVector photon1 = good_photons[iclus1].p4;
        double photon1_z = vertex_z + radius * std::sinh(photon1.Eta());
        double photon1_phi = (photon1.Phi() >= 0 ? photon1.Phi() : photon1.Phi() + 2 * M_PI);
        
        if ((photon1_z >= z_min && photon1_z <= z_max) &&
            ((phi_max >= phi_min && photon1_phi >= phi_min && photon1_phi <= phi_max) ||
             (phi_max <= phi_min && (photon1_phi <= phi_max || photon1_phi >= phi_min)))
            )
        {
          emulator_energies[itile].push_back(photon1.E());
          emulator_clusters[itile].push_back(iclus1);
        }
      }
      if (auto it = emulator_energies.find(itile); it != emulator_energies.end() && !it->second.empty())
      {
        emulator_multiplicities[itile] = emulator_energies[itile].size();

        // Sort the clusters by decreasing energy
        std::vector<int> indices(emulator_energies[itile].size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](int a, int b) -> bool { return emulator_energies[itile][a] > emulator_energies[itile][b]; });
        std::vector<int> vec_temp_clusters(emulator_clusters[itile]);
        std::vector<float> vec_temp_energies(emulator_energies[itile]);
        for (int iclus = 0; iclus < (int) emulator_energies[itile].size(); iclus++)
        {
          emulator_clusters[itile][iclus] = vec_temp_clusters[indices[iclus]];
          emulator_energies[itile][iclus] = vec_temp_energies[indices[iclus]];
        }
        vec_temp_clusters.clear();
        vec_temp_energies.clear();
      }
    }
  }

  // Identify the trigger firing clusters:
  // Clusters which account for more than 20% of the tile energy
  // Usually, there is only 1 per event
  std::map<int, std::vector<int>>::const_iterator it_trigger;
  std::map<int, std::vector<float>>::const_iterator it_energy;
  it_trigger = emulator_clusters.begin();
  it_energy = emulator_energies.begin();
  for (; it_trigger != emulator_clusters.end(); it_trigger++, it_energy++)
  {
    std::vector<float>::const_iterator it_local_energy;
    float trigger_window_sum = std::accumulate(it_energy->second.begin(), it_energy->second.end(), 0.0);
    float emulator_iclus = 0;
    float clus_energy = trigger_window_sum;
    for (size_t local_index = 0; local_index < it_trigger->second.size(); local_index++)
    {
      emulator_iclus = it_trigger->second[local_index];
      clus_energy = it_energy->second[local_index];
      if (clus_energy < 0.2 * trigger_window_sum)
      {
        break;
      }
      good_photons[emulator_iclus].isTrigger = true;
    }
  }
}

bool AnNeutralMeson_micro_dst::trigger_efficiency_matching(const ROOT::Math::PtEtaPhiMVector& photon1,
                                                           const ROOT::Math::PtEtaPhiMVector& photon2,
                                                           const ROOT::Math::PtEtaPhiMVector& diphoton)
                                            
{
  float energy_threshold = 0;

  if (trigger_mbd_photon_3) energy_threshold = energy_threshold_3; // Photon 3 GeV trigger
  else if (trigger_mbd_photon_4) energy_threshold = energy_threshold_4; // Photon 4 GeV trigger
  
  else 
  {
    std::cerr << "Error: either photon 3 GeV or photon 4 GeV scaled trigger bit should be fired at this point.\n";
    exit(1);
  }
  
  // If one of the two clusters has sufficient energy to fire the specific trigger
  // that led to the event written on disk, keep the event
  if (!(photon1.E() >= energy_threshold || photon2.E() >= energy_threshold))
  {
    // If no, check if:
    // - distance between two clusters is small enough for them to have been in the same trigger tile
    // - if their energy sum is sufficient to fire the specific trigger.
    float delta_eta = std::abs(photon1.Eta() - photon2.Eta());
    float delta_phi = WrapAngleDifference(photon1.Phi(), photon2.Phi());
    if (!((diphoton.E() > energy_threshold) &&
          (delta_eta < delta_eta_threshold) &&
          (delta_phi < delta_phi_threshold)))
    {
      // Otherwise, discard the diphoton
      return false;
    }
  }
  return true;
}

void AnNeutralMeson_micro_dst::smear_photon(float& eta,
                                            float& phi,
                                            float& ecore)
{
  // Defaults if parameter file missing or Eval invalid
  float E_res = 0.08f;
  float P_res = 0.0001f;

  std::vector<float> E_vector;
  std::vector<float> P_vector;
  E_vector.reserve(3);
  P_vector.reserve(3);
  
  for (int smear_idx = 0; smear_idx < nDifferences; smear_idx++)
  {
    const double se = f_energy_smear[smear_idx]->Eval(ecore);
    if (std::isfinite(se) && se > 0.0)
    {
      E_vector.push_back(static_cast<float>(se));
    }

    const double sp = f_position_smear[smear_idx]->Eval(ecore);
    if (std::isfinite(sp) && sp > 0.0)
    {
      P_vector.push_back(static_cast<float>(sp));
    }
  }
  if (!E_vector.empty())
  {
    E_res = *std::max_element(E_vector.begin(), E_vector.end());
  }
  if (!P_vector.empty())
  {
    P_res = *std::max_element(P_vector.begin(), P_vector.end());
  }

  float smeared_energy = std::max(0.0, rnd->Gaus(ecore, ecore * E_res));
  float smeared_phi = rnd->Gaus(phi, P_res);
  WrapAngle(smeared_phi);
  float smeared_eta = rnd->Gaus(eta, P_res);

  if (do_smear_eta) 
    eta = smeared_eta;
  if (do_smear_phi)
    phi = smeared_phi;
  if (do_smear_ecore)
    ecore = smeared_energy;
}

void AnNeutralMeson_micro_dst::event_mixing_mbd(PHCompositeNode *topNode)
{
  
  EventInPool current_event;
  current_event.zvtx = vertex_z;
  current_event.clusters = good_photons;
  
  // Associate cluster from current event to
  // clusters in pooled events
  int iPoolZvtxBin = FindBinBinary(vertex_z, poolZvtxBins, nPoolZvtxBins + 1);
  int iMultiplicityBin = FindBinBinary(num_photons, multiplicityBins, nMultiplicityBins + 1);

  // Determine the angle of the leading jet
  std::vector<ROOT::Math::PtEtaPhiMVector> photon_vectors;
  for (int iclus = 0; iclus < cluster_number; iclus++)
  {
    ClusterSmallInfo *smallcluster = _smallclusters->get_cluster_at(iclus);
    float eta = smallcluster->get_eta();
    float phi = smallcluster->get_phi();
    float ecore = smallcluster->get_ecore();
    ROOT::Math::PtEtaPhiMVector photon(ecore / std::cosh(eta), eta, phi, 0);
    photon_vectors.push_back(photon);
  }

  // Determine  the leading jet position
  JetContainer *jet_cont = findNode::getClass<JetContainer>(topNode, "AntiKt_Tower_r04");

  if (!jet_cont)
  {
    std::cout << "Jet Container is empty" << std::endl;
    return;
  }

  Jet *leading_jet = nullptr;
  float leading_pt = 0;
  float sum_relative_pt = 0;
  float sum_absolute_pt = 0;
  for (Jet *jet : *jet_cont)
  {
    std::cout << "(" << jet->get_pt() << ", " << jet->get_eta() << ", " << jet->get_phi() << ")" << std::endl; 
    sum_relative_pt += jet->get_pt();
    sum_absolute_pt += std::abs(jet->get_pt());
    if (jet->get_pt() > leading_pt)
    {
      leading_jet = jet;
    }
  }
  std::cout << std::endl;

  // Sharpness (to define "jetty" events):
  float S = std::abs(sum_relative_pt) / sum_absolute_pt;
  bool is_jet_event = false;
  if (S > 0.3) is_jet_event = true;

  std::cout << "Sharpness = " << S << std::endl;

  if (is_jet_event)
  {
    std::cout << "Jetty event (" << leading_jet->get_eta() << ", " << leading_jet->get_phi() << ")" << std::endl; 
    current_event.eta_lead = leading_jet->get_eta();
    current_event.phi_lead = leading_jet->get_phi();
  }
  else
  {
    std::cout << "Not a jetty event" << std::endl;
  }

  int iPoolEtaBin = FindBinBinary(current_event.eta_lead, poolEtaBins, nPoolEtaBins + 1);
  
  if ((iPoolZvtxBin < 0 || iPoolZvtxBin >= nPoolZvtxBins) ||
      (iMultiplicityBin < 0 || iMultiplicityBin >= nMultiplicityBins) ||
      (iPoolEtaBin < 0 || iPoolEtaBin >= nPoolEtaBins))
  {
    return;
  }
  
  int CurrentNmix = std::min(int(pool[iPoolZvtxBin][iMultiplicityBin][iPoolEtaBin].size()), Nmix);
  float sum_multiplicities_pool = 0;
  for (int ipool = 0; ipool < CurrentNmix; ipool++)
  {
    EventInPool pool_event = pool[iPoolZvtxBin][iMultiplicityBin][iPoolEtaBin][ipool];
    sum_multiplicities_pool += (float) pool_event.clusters.size();
  }

  float correction_weight = (sum_multiplicities_pool > 0.5 ? ((float) current_event.clusters.size() - 1)/ sum_multiplicities_pool : 0);

  for (auto current_cluster: current_event.clusters)
  {
    for (int ipool = 0; ipool < CurrentNmix; ipool++)
    {
      EventInPool pool_event = pool[iPoolZvtxBin][iMultiplicityBin][iPoolEtaBin][ipool];
      for (auto pool_cluster: pool_event.clusters)
      {
        // Adapt pseudorapidity direction to the current vertex position 
        float z_pool = pool_event.zvtx + radius * std::sinh(pool_cluster.p4.Eta());
        float eta_pool = std::asinh((z_pool - current_event.zvtx) / radius);
          
        ROOT::Math::PtEtaPhiMVector pool_photon(0, 0, 0, 0);
        if (is_jet_event) // rotation
        {
          // Adapt pseudorapidity direction to the current vertex position 
          float z_lead_pool = pool_event.zvtx + radius * std::sinh(pool_event.eta_lead);
          float eta_lead_pool = std::asinh((z_lead_pool - current_event.zvtx) / radius);
          
          // Rotate the event to align the jet axes
          float eta_rot = eta_pool + current_event.eta_lead - eta_lead_pool;
          float phi_rot = pool_cluster.p4.Phi() + current_event.phi_lead - pool_event.phi_lead;
          pool_photon.SetCoordinates(pool_cluster.p4.Pt(), eta_rot, phi_rot, 0);
        }
        else // Naïve event mixing
        {
          pool_photon.SetCoordinates(pool_cluster.p4.Pt(), eta_pool, pool_cluster.p4.Phi(), 0);
        }
          
        ROOT::Math::PtEtaPhiMVector mixed_pair = current_cluster.p4 + pool_photon;
        if (false)
          mixed_pair = current_cluster.p4 + pool_photon;

        if (store_qa)
        {
          h_current_E->Fill(current_cluster.p4.E());
          h_current_eta->Fill(current_cluster.p4.Eta());
          h_current_phi->Fill(current_cluster.p4.Phi());
          h_pool_E->Fill(pool_photon.E());
          h_pool_eta->Fill(pool_photon.Eta());
          h_pool_phi->Fill(pool_photon.Phi());
        }
        
        // The mixed event pair should satisfy the main cut as main analysis
        if (diphoton_cut(current_cluster.p4, pool_photon, mixed_pair)) continue;
        // In addition, require a minimum separation between two clusters
        // In reconstruction, two clusters cannot be too close
        float dR = ROOT::Math::VectorUtil::DeltaR(current_cluster.p4, pool_photon);
        if (dR < dR_min) continue;
        if (dR < 0.05) continue;

        if (store_qa && mixed_pair.M() > 0.25 && mixed_pair.M() < 0.38)
        {
          h_mixed_delta_eta->Fill(std::abs(current_cluster.p4.Eta() - pool_photon.Eta()), correction_weight);
          h_mixed_delta_phi->Fill(std::abs(ROOT::Math::VectorUtil::DeltaPhi(current_cluster.p4, pool_photon)), correction_weight);
          h_mixed_delta_R->Fill(dR, correction_weight);
          h_mixed_alpha->Fill(std::abs(current_cluster.p4.E() - pool_photon.E()) / (current_cluster.p4.E() + pool_photon.E()), correction_weight);
        }
        
        int iPt = FindBinBinary(mixed_pair.Pt(), pTBins, nPtBins + 1);
        
        float mixed_mass = mixed_pair.M();
        
        h_pair_mass_mixing->Fill(mixed_mass, correction_weight);
        h_pair_mass_mixing_pt[iPt]->Fill(mixed_mass, correction_weight);
      }
    }
  }
  pool[iPoolZvtxBin][iMultiplicityBin][iPoolEtaBin].push_front(current_event);
  if (CurrentNmix == Nmix)
  {
    pool[iPoolZvtxBin][iMultiplicityBin][iPoolEtaBin].pop_back();
  }
}

void AnNeutralMeson_micro_dst::event_mixing_photon()
{
  // if ((std::abs(vertex_z) >= 50) || (cluster_number > 11))
  //   return;
  
  // EventInPool current_event;
  // current_event.zvtx = vertex_z;
  // current_event.clusters = good_photons;
  
  // // Associate trigger cluster from current event to
  // // non-trigger clusters in pooled events
  // int iPoolZvtxBin = FindBinBinary(vertex_z, poolZvtxBins, nPoolZvtxBins + 1);
  // int iMultiplicityBin = FindBinBinary(cluster_number, multiplicityBins, nMultiplicityBins + 1);
  // int CurrentNmix = std::min(int(pool[iPoolZvtxBin][iMultiplicityBin].size()), Nmix);
  // for (auto current_cluster: current_event.clusters)
  // {
  //   if (! current_cluster.isTrigger) continue;
  //   for (int ipool = 0; ipool < CurrentNmix; ipool++)
  //   {
  //     EventInPool pool_event = pool[iPoolZvtxBin][iMultiplicityBin][ipool];
  //     for (auto pool_cluster: pool_event.clusters)
  //     {
  //       if (pool_cluster.isTrigger) continue;

  //       ROOT::Math::PtEtaPhiMVector mixed_pair = current_cluster.p4 + pool_cluster.p4;
  //       if (diphoton_cut(current_cluster.p4, pool_cluster.p4, mixed_pair)) continue;
  //       int iPt = FindBinBinary(mixed_pair.Pt(), pTBins, nPtBins + 1);
        
  //       float mixed_mass = mixed_pair.M();
        
  //       h_pair_mass_mixing->Fill(mixed_mass);
  //       h_pair_mass_mixing_pt[iPt]->Fill(mixed_mass);
  //     }
  //   }
  // }
  // pool[iPoolZvtxBin][iMultiplicityBin].push_front(current_event);
  // if (CurrentNmix == Nmix)
  // {
  //   pool[iPoolZvtxBin][iMultiplicityBin].pop_back();
  // }
} 

bool AnNeutralMeson_micro_dst::startswith(const std::string& str, const std::string& cmp)
{
  return str.compare(0, cmp.length(), cmp) == 0;
}

void AnNeutralMeson_micro_dst::WrapAngle(float& phi)
{
  while (phi > static_cast<float>(M_PI)) {
    phi -= 2.0 * M_PI;
  }
  while (phi < - static_cast<float>(M_PI)) {
    phi += 2.0 * M_PI;
  }
}

float AnNeutralMeson_micro_dst::WrapAngleDifference(const float& phi1, const float& phi2)
{
  float difference = std::fmod(phi1 - phi2, 2 * M_PI);
  if (difference < 0)
  {
    difference += 2 * M_PI;
  }
  return std::min(difference, (float)(2 * M_PI - difference));
}

ROOT::Math::XYZVector AnNeutralMeson_micro_dst::EnergyWeightedAverageP3(
    const std::vector<ROOT::Math::PtEtaPhiMVector>& v4s)
{
  ROOT::Math::XYZVector sum_vector(0.,0.,0.);
  double sumE = 0.0;

  for (const auto& v : v4s) {
    const double E = v.E();
    sum_vector  += E * ROOT::Math::XYZVector(v.Px(), v.Py(), v.Pz());
    sumE += E;
  }

  return (sumE > 0.0) ? (sum_vector * (1.0 / sumE)) : ROOT::Math::XYZVector(0.,0.,0.);
}
