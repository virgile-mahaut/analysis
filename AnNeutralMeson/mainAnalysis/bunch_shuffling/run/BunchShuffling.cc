#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "AsymmetryCalc/WrapperFinalAsymArrays.h"
#include "AsymmetryCalc/Constants.h"
#include "AsymmetryCalc/ShuffleBunches.h"

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

std::unordered_map<std::string, std::string> read_config(const std::string& filename)
{
  std::unordered_map<std::string, std::string> config;
  std::ifstream file(filename);
  if (!file.is_open())
  {
    std::cerr << "Warning: Config file " << filename << " not found. All config paths will be set to the default." << std::endl;
    return config;
  }

  std::string line;
  while (std::getline(file, line))
  {
    std::istringstream is_line(line);
    std::string key;
    if (std::getline(is_line, key, '='))
    {
      std::string value;
      if (std::getline(is_line, value, '='))
      {
        config[key] = value;
      }
    }
  }
  return config;
}

void book_shuffle_file(TFile* (&output_file_shuffle),
                       TH1F* (&h_AN_pt_shuffle)[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nDirections],
                       TH1F* (&h_AN_pt_nodir_shuffle)[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins],
                       TH1F* (&h_AN_eta_shuffle)[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nEtaBins],
                       TH1F* (&h_AN_xf_shuffle)[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nXfBins],
                       TH1F* (&h_AN_pt_shuffle_details)[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nDirections][ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nConfigs],
                       TH1F* (&h_AN_pt_nodir_shuffle_details)[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nConfigs],
                       TH1F* (&h_AN_eta_shuffle_details)[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nEtaBins][ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nConfigs],
                       TH1F* (&h_AN_xf_shuffle_details)[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nXfBins][ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nConfigs],
                       const std::string &outputfilename)
{
  output_file_shuffle = new TFile(outputfilename.c_str(), "RECREATE");
  output_file_shuffle->cd();
  for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++) {
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      {
        std::stringstream hname;
        hname << "h_AN_pt_nodir_shuffle_"
            << ASYM_CONSTANTS::particle[iP] << "_pt_"
            << iPt << "_sqrt";
        std::stringstream htitle;
        htitle << ";A_{N}/#delta A_{N}^{stat}; Counts";
        h_AN_pt_nodir_shuffle[iP][iPt] = new TH1F(
          hname.str().c_str(),
          htitle.str().c_str(),
          100, -4, 4);
        // Beam- and config- specific values
        for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++) {
          for (int iC = 0; iC < ASYM_CONSTANTS::nConfigs; iC++) {
            hname.str("");
            hname << "h_AN_pt_nodir_shuffle_"
                  << ASYM_CONSTANTS::beams[iB] << "_"
                  << ASYM_CONSTANTS::particle[iP] << "_pt_"
                  << iPt << "_"
                  << ASYM_CONSTANTS::configuration[iC] << "_sqrt";
            htitle.str("");
            htitle << ";A_{N}/#delta A_{N}^{stat}; Counts";
            h_AN_pt_nodir_shuffle_details[iP][iPt][iB][iC] = new TH1F(
              hname.str().c_str(),
              htitle.str().c_str(),
              100, -4, 4);
          }
        }
      }
      for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++) {
        std::stringstream hname;
        hname << "h_AN_pt_shuffle_"
            << ASYM_CONSTANTS::particle[iP] << "_"
            << ASYM_CONSTANTS::directions[iDir] << "_pt_"
            << iPt << "_sqrt";
        std::stringstream htitle;
        htitle << ";A_{N}/#delta A_{N}^{stat}; Counts";
        h_AN_pt_shuffle[iP][iPt][iDir] = new TH1F(
          hname.str().c_str(),
          htitle.str().c_str(),
          100, -4, 4);
        // Beam- and config- specific values
        for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++) {
          for (int iC = 0; iC < ASYM_CONSTANTS::nConfigs; iC++) {
            hname.str("");
            hname << "h_AN_pt_shuffle_"
                  << ASYM_CONSTANTS::beams[iB] << "_"
                  << ASYM_CONSTANTS::particle[iP] << "_"
                  << ASYM_CONSTANTS::directions[iDir] << "_pt_"
                  << iPt << "_"
                  << ASYM_CONSTANTS::configuration[iC] << "_sqrt";
            htitle.str("");
            htitle << ";A_{N}/#delta A_{N}^{stat}; Counts";
            h_AN_pt_shuffle_details[iP][iPt][iDir][iB][iC] = new TH1F(
              hname.str().c_str(),
              htitle.str().c_str(),
              100, -4, 4);
          }
        }
      }
    }
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      std::stringstream hname;
      hname << "h_AN_eta_shuffle_"
            << ASYM_CONSTANTS::particle[iP] << "_eta_"
            << iEta << "_sqrt";
      std::stringstream htitle;
      htitle << ";A_{N}/#delta A_{N}^{stat}; Counts";
      h_AN_eta_shuffle[iP][iEta] = new TH1F(
        hname.str().c_str(),
        htitle.str().c_str(),
        100, -4, 4);
      // Beam- and config- specific values
      for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++) {
        for (int iC = 0; iC < ASYM_CONSTANTS::nConfigs; iC++) {
          hname.str("");
          hname << "h_AN_eta_shuffle_"
                << ASYM_CONSTANTS::beams[iB] << "_"
                << ASYM_CONSTANTS::particle[iP] << "_eta_"
                << iEta << "_"
                << ASYM_CONSTANTS::configuration[iC] << "_sqrt";
          htitle.str("");
          htitle << ";A_{N}/#delta A_{N}^{stat}; Counts";
          h_AN_eta_shuffle_details[iP][iEta][iB][iC] = new TH1F(
            hname.str().c_str(),
            htitle.str().c_str(),
            100, -4, 4);
        }
      }
    }
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      std::stringstream hname;
      hname << "h_AN_xf_shuffle_"
            << ASYM_CONSTANTS::particle[iP] << "_xf_"
            << iXf << "_sqrt";
      std::stringstream htitle;
      htitle << ";A_{N}/#delta A_{N}^{stat}; Counts";
      h_AN_xf_shuffle[iP][iXf] = new TH1F(
        hname.str().c_str(),
        htitle.str().c_str(),
        100, -4, 4);
      // Beam- and config- specific values
      for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++) {
        for (int iC = 0; iC < ASYM_CONSTANTS::nConfigs; iC++) {
          hname.str("");
          hname << "h_AN_xf_shuffle_"
                << ASYM_CONSTANTS::beams[iB] << "_"
                << ASYM_CONSTANTS::particle[iP] << "_xf_"
                << iXf << "_"
                << ASYM_CONSTANTS::configuration[iC] << "_sqrt";
          htitle.str("");
          htitle << ";A_{N}/#delta A_{N}^{stat}; Counts";
          h_AN_xf_shuffle_details[iP][iXf][iB][iC] = new TH1F(
            hname.str().c_str(),
            htitle.str().c_str(),
            100, -4, 4);
        }
      }
    }
  }
}

void fill_shuffle_detailed_histos(
       AsymmetryCalc::ShuffleBunches *shuffle_bunches_0mrad_mbd,
       AsymmetryCalc::ShuffleBunches *shuffle_bunches_0mrad_photon,
       AsymmetryCalc::ShuffleBunches *shuffle_bunches_15mrad_mbd,
       AsymmetryCalc::ShuffleBunches *shuffle_bunches_15mrad_photon,
       TH1F* h_AN_pt_shuffle_details[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nDirections][ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nConfigs],
       TH1F* h_AN_pt_nodir_shuffle_details[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nConfigs],
       TH1F* h_AN_eta_shuffle_details[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nEtaBins][ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nConfigs],
       TH1F* h_AN_xf_shuffle_details[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nXfBins][ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nConfigs],
       const int iIter)
{
  // Total of 4 configurations
  const int nConfigs = ASYM_CONSTANTS::nConfigs;
  AsymmetryCalc::ShuffleBunches *shuffle_bunches[nConfigs] = {
    shuffle_bunches_0mrad_mbd,
    shuffle_bunches_0mrad_photon,
    shuffle_bunches_15mrad_mbd,
    shuffle_bunches_15mrad_photon
  };

  // Load corrected asymmetries from every configuration
  std::vector<const AsymmetryCalc::pt_corr_array*> corrected_mean_pt(nConfigs);
  std::vector<const AsymmetryCalc::pt_nodir_corr_array*> corrected_mean_pt_nodir(nConfigs);
  std::vector<const AsymmetryCalc::eta_corr_array*> corrected_mean_eta(nConfigs);
  std::vector<const AsymmetryCalc::xf_corr_array*> corrected_mean_xf(nConfigs);
  std::vector<const AsymmetryCalc::pt_corr_array*> corrected_unc_pt(nConfigs);
  std::vector<const AsymmetryCalc::pt_nodir_corr_array*> corrected_unc_pt_nodir(nConfigs);
  std::vector<const AsymmetryCalc::eta_corr_array*> corrected_unc_eta(nConfigs);
  std::vector<const AsymmetryCalc::xf_corr_array*> corrected_unc_xf(nConfigs);

  for (int iConfig = 0; iConfig < nConfigs; iConfig++) {
    corrected_mean_pt[iConfig] = &shuffle_bunches[iConfig]->get_corrected_mean_pt()[iIter];
    corrected_mean_pt_nodir[iConfig] = &shuffle_bunches[iConfig]->get_corrected_mean_pt_nodir()[iIter];
    corrected_mean_eta[iConfig] = &shuffle_bunches[iConfig]->get_corrected_mean_eta()[iIter];
    corrected_mean_xf[iConfig] = &shuffle_bunches[iConfig]->get_corrected_mean_xf()[iIter];
    corrected_unc_pt[iConfig] = &shuffle_bunches[iConfig]->get_corrected_unc_pt()[iIter];
    corrected_unc_pt_nodir[iConfig] = &shuffle_bunches[iConfig]->get_corrected_unc_pt_nodir()[iIter];
    corrected_unc_eta[iConfig] = &shuffle_bunches[iConfig]->get_corrected_unc_eta()[iIter];
    corrected_unc_xf[iConfig] = &shuffle_bunches[iConfig]->get_corrected_unc_xf()[iIter];
  }

  for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++) {
    for (int iC = 0; iC < ASYM_CONSTANTS::nConfigs; iC++) {
      for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++) {
        for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
          {
            // Beam- and config- specific values
            h_AN_pt_nodir_shuffle_details[iP][iPt][iB][iC]->Fill(
              (*corrected_mean_pt_nodir[iC])(iB, iP, iPt) /
              (*corrected_unc_pt_nodir[iC])(iB, iP, iPt));
          }
          for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++) {
            // Beam- and config- specific values
            h_AN_pt_shuffle_details[iP][iPt][iDir][iB][iC]->Fill(
              (*corrected_mean_pt[iC])(iB, iP, iPt, iDir) /
              (*corrected_unc_pt[iC])(iB, iP, iPt, iDir));
          }
        }
        for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
          h_AN_eta_shuffle_details[iP][iEta][iB][iC]->Fill(
            (*corrected_mean_eta[iC])(iB, iP, iEta) /
            (*corrected_unc_eta[iC])(iB, iP, iEta));
                                                       
        }
        for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
          h_AN_xf_shuffle_details[iP][iXf][iB][iC]->Fill(
            (*corrected_mean_xf[iC])(iB, iP, iXf) /
            (*corrected_unc_xf[iC])(iB, iP, iXf));
        }
      }
    }
  }
}

void fill_shuffle_histos(AsymmetryCalc::pt_final_array& merged_mean_pt,
                         AsymmetryCalc::pt_nodir_final_array& merged_mean_pt_nodir,
                         AsymmetryCalc::eta_final_array& merged_mean_eta,
                         AsymmetryCalc::xf_final_array& merged_mean_xf,
                         AsymmetryCalc::pt_final_array& merged_unc_pt,
                         AsymmetryCalc::pt_nodir_final_array& merged_unc_pt_nodir,
                         AsymmetryCalc::eta_final_array& merged_unc_eta,
                         AsymmetryCalc::xf_final_array& merged_unc_xf,
                         TH1F* h_AN_pt_shuffle[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nDirections],
                         TH1F* h_AN_pt_nodir_shuffle[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins],
                         TH1F* h_AN_eta_shuffle[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nEtaBins],
                         TH1F* h_AN_xf_shuffle[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nXfBins])
{
  for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++) {
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      {
        h_AN_pt_nodir_shuffle[iP][iPt]->Fill(merged_mean_pt_nodir(iP, iPt) / merged_unc_pt_nodir(iP, iPt));
      }
      for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++) {
        h_AN_pt_shuffle[iP][iPt][iDir]->Fill(merged_mean_pt(iP, iPt, iDir) / merged_unc_pt(iP, iPt, iDir));
      }
    }
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      h_AN_eta_shuffle[iP][iEta]->Fill(merged_mean_eta(iP, iEta) / merged_unc_eta(iP, iEta));
    }
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      h_AN_xf_shuffle[iP][iXf]->Fill(merged_mean_xf(iP, iXf) / merged_unc_xf(iP, iXf));
    }
  }
}

void save_shuffle_histos(TFile *output_file_shuffle)
{
  output_file_shuffle->cd();
  output_file_shuffle->Write();
  output_file_shuffle->Close();
  delete output_file_shuffle;
  output_file_shuffle = nullptr;
}

static bool valid_measurement(double val, double unc)
{
  return std::isfinite(val) && std::isfinite(unc) && unc > 0.0 && val != 0;
}

void compute_average_2(double& res_val, double& res_unc,
                       double in_val_1, double in_unc_1,
                       double in_val_2, double in_unc_2)
{
    const bool ok1 = valid_measurement(in_val_1, in_unc_1);
    const bool ok2 = valid_measurement(in_val_2, in_unc_2);

    if (!ok1 && !ok2) {
        res_val = 0.0;
        res_unc = 0.0;
        return;
    }
    if (!ok1) {
        res_val = in_val_2;
        res_unc = in_unc_2;
        return;
    }
    if (!ok2) {
        res_val = in_val_1;
        res_unc = in_unc_1;
        return;
    }

    const double w1 = 1.0 / (in_unc_1 * in_unc_1);
    const double w2 = 1.0 / (in_unc_2 * in_unc_2);

    res_val = (w1 * in_val_1 + w2 * in_val_2) / (w1 + w2);
    res_unc = 1.0 / std::sqrt(w1 + w2);
}                      

void save_final_histograms(const std::string &name,
                           AsymmetryCalc::pt_final_array& merged_mean_pt,
                           AsymmetryCalc::pt_nodir_final_array& merged_mean_pt_nodir,
                           AsymmetryCalc::eta_final_array& merged_mean_eta,
                           AsymmetryCalc::xf_final_array& merged_mean_xf,
                           AsymmetryCalc::pt_final_array& merged_unc_pt,
                           AsymmetryCalc::pt_nodir_final_array& merged_unc_pt_nodir,
                           AsymmetryCalc::eta_final_array& merged_unc_eta,
                           AsymmetryCalc::xf_final_array& merged_unc_xf)
{
  TFile *outputfile = new TFile(name.c_str(), "RECREATE");
  outputfile->cd();
  for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++) {
    for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++) {
      TGraphErrors *graph_AN = new TGraphErrors();
      std::stringstream graph_name;
      graph_name << "graph_AN_pt_final_" << ASYM_CONSTANTS::particle[iP] << "_"
                 << ASYM_CONSTANTS::directions[iDir] << "_sqrt";
      std::stringstream graph_title;
      graph_title << "; iPt; A_N";
      graph_AN->SetName(graph_name.str().c_str());
      graph_AN->SetTitle(graph_title.str().c_str());
      for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
        graph_AN->SetPoint(iPt, (double)(iPt+1), merged_mean_pt(iP, iPt, iDir));
        graph_AN->SetPointError(iPt, 0, merged_unc_pt(iP, iPt, iDir)); 
      }
      graph_AN->Write();
    }
    {
      TGraphErrors *graph_AN = new TGraphErrors();
      std::stringstream graph_name;
      graph_name << "graph_AN_eta_final_" << ASYM_CONSTANTS::particle[iP] << "_sqrt";
      std::stringstream graph_title;
      graph_title << "; iEta; A_N";
      graph_AN->SetName(graph_name.str().c_str());
      graph_AN->SetTitle(graph_title.str().c_str());
      for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
        graph_AN->SetPoint(iEta, (double)(iEta+1), merged_mean_eta(iP, iEta));
        graph_AN->SetPointError(iEta, 0, merged_unc_eta(iP, iEta)); 
      }
      graph_AN->Write();
    }
    {
      TGraphErrors *graph_AN = new TGraphErrors();
      std::stringstream graph_name;
      graph_name << "graph_AN_xf_final_" << ASYM_CONSTANTS::particle[iP] << "_sqrt";
      std::stringstream graph_title;
      graph_title << "; iXf; A_N";
      graph_AN->SetName(graph_name.str().c_str());
      graph_AN->SetTitle(graph_title.str().c_str());
      for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
        graph_AN->SetPoint(iXf, (double)(iXf+1), merged_mean_xf(iP, iXf));
        graph_AN->SetPointError(iXf, 0, merged_unc_xf(iP, iXf)); 
      }
      graph_AN->Write();
    }
  }
  outputfile->Close();
  delete outputfile;
}

void merge_configurations(AsymmetryCalc::ShuffleBunches *shuffle_bunches_0mrad_mbd,
                          AsymmetryCalc::ShuffleBunches *shuffle_bunches_0mrad_photon,
                          AsymmetryCalc::ShuffleBunches *shuffle_bunches_15mrad_mbd,
                          AsymmetryCalc::ShuffleBunches *shuffle_bunches_15mrad_photon,
                          AsymmetryCalc::pt_final_array& merged_mean_pt,
                          AsymmetryCalc::pt_nodir_final_array& merged_mean_pt_nodir,
                          AsymmetryCalc::eta_final_array& merged_mean_eta,
                          AsymmetryCalc::xf_final_array& merged_mean_xf,
                          AsymmetryCalc::pt_final_array& merged_unc_pt,
                          AsymmetryCalc::pt_nodir_final_array& merged_unc_pt_nodir,
                          AsymmetryCalc::eta_final_array& merged_unc_eta,
                          AsymmetryCalc::xf_final_array& merged_unc_xf,
                          const int iIter)
{
  // Total of 4 configurations
  const int nConfigs = ASYM_CONSTANTS::nConfigs;
  AsymmetryCalc::ShuffleBunches *shuffle_bunches[nConfigs] = {
    shuffle_bunches_0mrad_mbd,
    shuffle_bunches_0mrad_photon,
    shuffle_bunches_15mrad_mbd,
    shuffle_bunches_15mrad_photon
  };

  // Load corrected asymmetries from every configuration
  std::vector<const AsymmetryCalc::pt_corr_array*> corrected_mean_pt(nConfigs);
  std::vector<const AsymmetryCalc::pt_nodir_corr_array*> corrected_mean_pt_nodir(nConfigs);
  std::vector<const AsymmetryCalc::eta_corr_array*> corrected_mean_eta(nConfigs);
  std::vector<const AsymmetryCalc::xf_corr_array*> corrected_mean_xf(nConfigs);
  std::vector<const AsymmetryCalc::pt_corr_array*> corrected_unc_pt(nConfigs);
  std::vector<const AsymmetryCalc::pt_nodir_corr_array*> corrected_unc_pt_nodir(nConfigs);
  std::vector<const AsymmetryCalc::eta_corr_array*> corrected_unc_eta(nConfigs);
  std::vector<const AsymmetryCalc::xf_corr_array*> corrected_unc_xf(nConfigs);

  for (int iConfig = 0; iConfig < nConfigs; iConfig++) {
    corrected_mean_pt[iConfig] = &shuffle_bunches[iConfig]->get_corrected_mean_pt()[iIter];
    corrected_mean_pt_nodir[iConfig] = &shuffle_bunches[iConfig]->get_corrected_mean_pt_nodir()[iIter];
    corrected_mean_eta[iConfig] = &shuffle_bunches[iConfig]->get_corrected_mean_eta()[iIter];
    corrected_mean_xf[iConfig] = &shuffle_bunches[iConfig]->get_corrected_mean_xf()[iIter];
    corrected_unc_pt[iConfig] = &shuffle_bunches[iConfig]->get_corrected_unc_pt()[iIter];
    corrected_unc_pt_nodir[iConfig] = &shuffle_bunches[iConfig]->get_corrected_unc_pt_nodir()[iIter];
    corrected_unc_eta[iConfig] = &shuffle_bunches[iConfig]->get_corrected_unc_eta()[iIter];
    corrected_unc_xf[iConfig] = &shuffle_bunches[iConfig]->get_corrected_unc_xf()[iIter];
  }

  // Average the 4 configurations in each kinematic bin
  for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++) {
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      {
        // First average over trigger
        double mean_trigger[4] = {0};
        double unc_trigger[4] = {0};
        compute_average_2(
             mean_trigger[0], unc_trigger[0],
             (*corrected_mean_pt_nodir[0])(0, iP, iPt),
             (*corrected_unc_pt_nodir[0])(0, iP, iPt),
             (*corrected_mean_pt_nodir[1])(0, iP, iPt),
             (*corrected_unc_pt_nodir[1])(0, iP, iPt)); // Yellow 0 mrad

        compute_average_2(
             mean_trigger[1], unc_trigger[1],
             (*corrected_mean_pt_nodir[2])(0, iP, iPt),
             (*corrected_unc_pt_nodir[2])(0, iP, iPt),
             (*corrected_mean_pt_nodir[3])(0, iP, iPt),
             (*corrected_unc_pt_nodir[3])(0, iP, iPt)); // Yellow 1.5 mrad

        compute_average_2(
             mean_trigger[2], unc_trigger[2],
             (*corrected_mean_pt_nodir[0])(1, iP, iPt),
             (*corrected_unc_pt_nodir[0])(1, iP, iPt),
             (*corrected_mean_pt_nodir[1])(1, iP, iPt),
             (*corrected_unc_pt_nodir[1])(1, iP, iPt)); // Blue 0 mrad

        compute_average_2(
             mean_trigger[3], unc_trigger[3],
             (*corrected_mean_pt_nodir[2])(1, iP, iPt),
             (*corrected_unc_pt_nodir[2])(1, iP, iPt),
             (*corrected_mean_pt_nodir[3])(1, iP, iPt),
             (*corrected_unc_pt_nodir[3])(1, iP, iPt)); // Blue 1.5 mrad
        
        // Then average over crossing angle
        double mean_trigger_angle[2] = {0};
        double unc_trigger_angle[2] = {0};
        compute_average_2(
             mean_trigger_angle[0], unc_trigger_angle[0],
             mean_trigger[0], unc_trigger[0],
             mean_trigger[1], unc_trigger[1]); // Yellow
        compute_average_2(
             mean_trigger_angle[1], unc_trigger_angle[1],
             mean_trigger[2], unc_trigger[2],
             mean_trigger[3], unc_trigger[3]); // Blue

        // And finally average over beam
        double mean_final = 0;
        double unc_final = 0;
        compute_average_2(
             mean_final, unc_final,
             mean_trigger_angle[0], unc_trigger_angle[0],
             mean_trigger_angle[1], unc_trigger_angle[1]);
        merged_mean_pt_nodir(iP, iPt) = mean_final;
        merged_unc_pt_nodir(iP, iPt) = unc_final;
      }
      for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++) {
        
        // First average over trigger
        double mean_trigger[4] = {0};
        double unc_trigger[4] = {0};
        compute_average_2(
             mean_trigger[0], unc_trigger[0],
             (*corrected_mean_pt[0])(0, iP, iPt, iDir),
             (*corrected_unc_pt[0])(0, iP, iPt, iDir),
             (*corrected_mean_pt[1])(0, iP, iPt, iDir),
             (*corrected_unc_pt[1])(0, iP, iPt, iDir)); // Yellow 0 mrad

        compute_average_2(
             mean_trigger[1], unc_trigger[1],
             (*corrected_mean_pt[2])(0, iP, iPt, iDir),
             (*corrected_unc_pt[2])(0, iP, iPt, iDir),
             (*corrected_mean_pt[3])(0, iP, iPt, iDir),
             (*corrected_unc_pt[3])(0, iP, iPt, iDir)); // Yellow 1.5 mrad

        compute_average_2(
             mean_trigger[2], unc_trigger[2],
             (*corrected_mean_pt[0])(1, iP, iPt, iDir),
             (*corrected_unc_pt[0])(1, iP, iPt, iDir),
             (*corrected_mean_pt[1])(1, iP, iPt, iDir),
             (*corrected_unc_pt[1])(1, iP, iPt, iDir)); // Blue 0 mrad

        compute_average_2(
             mean_trigger[3], unc_trigger[3],
             (*corrected_mean_pt[2])(1, iP, iPt, iDir),
             (*corrected_unc_pt[2])(1, iP, iPt, iDir),
             (*corrected_mean_pt[3])(1, iP, iPt, iDir),
             (*corrected_unc_pt[3])(1, iP, iPt, iDir)); // Blue 1.5 mrad
        
        // Then average over crossing angle
        double mean_trigger_angle[2] = {0};
        double unc_trigger_angle[2] = {0};
        compute_average_2(
             mean_trigger_angle[0], unc_trigger_angle[0],
             mean_trigger[0], unc_trigger[0],
             mean_trigger[1], unc_trigger[1]); // Yellow
        compute_average_2(
             mean_trigger_angle[1], unc_trigger_angle[1],
             mean_trigger[2], unc_trigger[2],
             mean_trigger[3], unc_trigger[3]); // Blue

        // And finally average over beam
        double mean_final = 0;
        double unc_final = 0;
        compute_average_2(
             mean_final, unc_final,
             mean_trigger_angle[0], unc_trigger_angle[0],
             mean_trigger_angle[1], unc_trigger_angle[1]);
        merged_mean_pt(iP, iPt, iDir) = mean_final;
        merged_unc_pt(iP, iPt, iDir) = unc_final;
      }
    }
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      // First average over trigger
      double mean_trigger[4] = {0};
      double unc_trigger[4] = {0};
      
      compute_average_2(
           mean_trigger[0], unc_trigger[0],
           (*corrected_mean_eta[0])(0, iP, iEta),
           (*corrected_unc_eta[0])(0, iP, iEta),
           (*corrected_mean_eta[1])(0, iP, iEta),
           (*corrected_unc_eta[1])(0, iP, iEta)); // Yellow 0 mrad
      compute_average_2(
           mean_trigger[1], unc_trigger[1],
           (*corrected_mean_eta[2])(0, iP, iEta),
           (*corrected_unc_eta[2])(0, iP, iEta),
           (*corrected_mean_eta[3])(0, iP, iEta),
           (*corrected_unc_eta[3])(0, iP, iEta)); // Yellow 1.5 mrad
      compute_average_2(
           mean_trigger[2], unc_trigger[2],
           (*corrected_mean_eta[0])(1, iP, iEta),
           (*corrected_unc_eta[0])(1, iP, iEta),
           (*corrected_mean_eta[1])(1, iP, iEta),
           (*corrected_unc_eta[1])(1, iP, iEta)); // Blue 0 mrad
      compute_average_2(
           mean_trigger[3], unc_trigger[3],
           (*corrected_mean_eta[2])(1, iP, iEta),
           (*corrected_unc_eta[2])(1, iP, iEta),
           (*corrected_mean_eta[3])(1, iP, iEta),
           (*corrected_unc_eta[3])(1, iP, iEta)); // Blue 1.5 mrad

      // Then average over crossing angle
      double mean_trigger_angle[2] = {0};
      double unc_trigger_angle[2] = {0};
      compute_average_2(
           mean_trigger_angle[0], unc_trigger_angle[0],
           mean_trigger[0], unc_trigger[0],
           mean_trigger[1], unc_trigger[1]); // Yellow
      compute_average_2(
           mean_trigger_angle[1], unc_trigger_angle[1],
           mean_trigger[2], unc_trigger[2],
           mean_trigger[3], unc_trigger[3]); // Blue

      // And finally average over beam
      double mean_final = 0;
      double unc_final = 0;
      compute_average_2(
           mean_final, unc_final,
           mean_trigger_angle[0], unc_trigger_angle[0],
           mean_trigger_angle[1], unc_trigger_angle[1]);
      //std::cout << "final (" << mean_final << ", " << unc_final << ")\n\n" << std::endl;

      merged_mean_eta(iP, iEta) = mean_final;
      merged_unc_eta(iP, iEta) = unc_final;
    }
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      
      // First average over trigger
      double mean_trigger[4] = {0};
      double unc_trigger[4] = {0};

      compute_average_2(
           mean_trigger[0], unc_trigger[0],
           (*corrected_mean_xf[0])(0, iP, iXf),
           (*corrected_unc_xf[0])(0, iP, iXf),
           (*corrected_mean_xf[1])(0, iP, iXf),
           (*corrected_unc_xf[1])(0, iP, iXf)); // Yellow 0 mrad
      compute_average_2(
           mean_trigger[1], unc_trigger[1],
           (*corrected_mean_xf[2])(0, iP, iXf),
           (*corrected_unc_xf[2])(0, iP, iXf),
           (*corrected_mean_xf[3])(0, iP, iXf),
           (*corrected_unc_xf[3])(0, iP, iXf)); // Yellow 1.5 mrad
      compute_average_2(
           mean_trigger[2], unc_trigger[2],
           (*corrected_mean_xf[0])(1, iP, iXf),
           (*corrected_unc_xf[0])(1, iP, iXf),
           (*corrected_mean_xf[1])(1, iP, iXf),
           (*corrected_unc_xf[1])(1, iP, iXf)); // Blue 0 mrad
      compute_average_2(
           mean_trigger[3], unc_trigger[3],
           (*corrected_mean_xf[2])(1, iP, iXf),
           (*corrected_unc_xf[2])(1, iP, iXf),
           (*corrected_mean_xf[3])(1, iP, iXf),
           (*corrected_unc_xf[3])(1, iP, iXf)); // Blue 1.5 mrad

      // Then average over crossing angle
      double mean_trigger_angle[2] = {0};
      double unc_trigger_angle[2] = {0};
      compute_average_2(
           mean_trigger_angle[0], unc_trigger_angle[0],
           mean_trigger[0], unc_trigger[0],
           mean_trigger[1], unc_trigger[1]); // Yellow
      compute_average_2(
           mean_trigger_angle[1], unc_trigger_angle[1],
           mean_trigger[2], unc_trigger[2],
           mean_trigger[3], unc_trigger[3]); // Blue

      // And finally average over beam
      double mean_final = 0;
      double unc_final = 0;
      compute_average_2(
           mean_final, unc_final,
           mean_trigger_angle[0], unc_trigger_angle[0],
           mean_trigger_angle[1], unc_trigger_angle[1]);

      merged_mean_xf(iP, iXf) = mean_final;
      merged_unc_xf(iP, iXf) = unc_final;
    }
  }
}

int main(int argc, char *argv[])
{
  if (argc < 3) {
    std::cout << "Usage: ./BunchShuffling <iterMin> <iterMax> [<config file> (default: config.txt)]" << std::endl;
    return 1;
  }

  int shuffleIterMin = std::stoi(argv[1]);
  int shuffleIterMax = std::stoi(argv[2]);

  std::string config_filename = "config.txt";
  if (argc > 3)
  {
    config_filename = argv[3];
  }

  const int shuffleMaxStorage = 10; // To avoid having too much resident memory

  bool perform_shuffle = true; // Set to false to deactivate the shuffle (debug only)
  bool store_raw = false;
  bool store_final = false;

  // Default path
  std::unordered_map<std::string, std::string> config = read_config(config_filename);
  
  std::string path_to_bunch_dependent_yields_mbd = config["path_to_bunch_dependent_yields_mbd"];
  std::string path_to_bunch_dependent_yields_photon = config["path_to_bunch_dependent_yields_photon"];
  std::string path_to_fill_list_0mrad = config["path_to_fill_list_0mrad"];
  std::string path_to_fill_list_15mrad = config["path_to_fill_list_15mrad"];
  std::string path_to_spin_pattern_0mrad = config["path_to_spin_pattern_0mrad"];
  std::string path_to_spin_pattern_15mrad = config["path_to_spin_pattern_15mrad"];
  std::string path_to_csv_mean = config["path_to_csv_mean"];
  std::string path_to_csv_ratios = config["path_to_csv_ratios"];
  std::string suffix = config["suffix"]; // This parameter is allowed to be empty
  std::string production_date = config["production_date"]; // This parameter is allowed to be empty
  std::string stored_iter_asym_histos = config["stored_iter_asym_histos"];
  std::string stored_final_asym_histos = config["stored_final_asym_histos"]; 
  
  if (path_to_bunch_dependent_yields_mbd.empty())
  {
    path_to_bunch_dependent_yields_mbd = "/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/macros/comparison_micro_DST_analysis/nano_analysis/analysis_complete_ana509_01312026_MBD_pt_l4_0mrad/fast_analysis_";
  }
  if (path_to_bunch_dependent_yields_photon.empty())
  {
    path_to_bunch_dependent_yields_photon = "/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/macros/comparison_micro_DST_analysis/nano_analysis/analysis_complete_ana509_01312026_efficiency_30_pt_g3_0mrad/fast_analysis_";
  }
  if (path_to_fill_list_0mrad.empty())
  {
    path_to_fill_list_0mrad = "/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/macros/comparison_micro_DST_analysis/nano_analysis/run_to_fill_0mrad.txt";
  }
  if (path_to_fill_list_15mrad.empty())
  {
    path_to_fill_list_15mrad = "/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/macros/comparison_micro_DST_analysis/nano_analysis/run_to_fill_15mrad.txt";
  }
  if (path_to_spin_pattern_0mrad.empty())
  {
    path_to_spin_pattern_0mrad = "/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/macros/comparison_micro_DST_analysis/nano_analysis/spin_pattern_fill_0mrad.root";
  }
  if (path_to_spin_pattern_15mrad.empty())
  {
    path_to_spin_pattern_15mrad = "/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/macros/comparison_micro_DST_analysis/nano_analysis/spin_pattern_fill_15mrad.root";
  }
  if (path_to_csv_mean.empty())
  {
    path_to_csv_mean = "/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/macros/comparison_micro_DST_analysis/nano_analysis/csv_inputs";
  }
  if (path_to_csv_ratios.empty())
  {
    path_to_csv_ratios = "/sphenix/u/virgilemahaut/work/analysis/AnNeutralMeson/macros/comparison_micro_DST_analysis/nano_analysis/csv_ratios";
    std::cout << "default path_to_csv_ratios = " << path_to_csv_ratios << std::endl;
  }
  else {
    std::cout << "custom path_to_csv_ratios = " << path_to_csv_ratios << std::endl;
  }
  if (production_date.empty()) {
    production_date = "06052026";
  }

  if (!stored_iter_asym_histos.empty()) {
    store_raw = true;
  }
  if (!stored_final_asym_histos.empty()) {
    store_final = true;
  }

  // All pT/eta/xF means and ratios
  std::string input_csv_mbd_0mrad_pt_pi0_name = path_to_csv_mean + "/input_pt_mean_" + production_date + suffix + "_MBD_pt_l3_0mrad.csv";
  std::string input_csv_mbd_0mrad_pt_eta_name = path_to_csv_mean + "/input_pt_mean_" + production_date + suffix + "_MBD_pt_l4_0mrad.csv";
  std::string input_csv_mbd_0mrad_eta_pi0_name = path_to_csv_mean + "/input_eta_mean_" + production_date + suffix + "_MBD_pt_l3_0mrad.csv";
  std::string input_csv_mbd_0mrad_eta_eta_name = path_to_csv_mean + "/input_eta_mean_" + production_date + suffix + "_MBD_pt_l4_0mrad.csv";
  std::string input_csv_mbd_0mrad_xf_pi0_name = path_to_csv_mean + "/input_xf_mean_" + production_date + suffix + "_MBD_pt_l3_0mrad.csv";
  std::string input_csv_mbd_0mrad_xf_eta_name = path_to_csv_mean + "/input_xf_mean_" + production_date + suffix + "_MBD_pt_l4_0mrad.csv";
  std::string input_csv_mbd_0mrad_pt_ratios_pi0_name = path_to_csv_ratios + "/pt_ratios_" + production_date + suffix + "_MBD_pt_l3_0mrad.csv";
  std::string input_csv_mbd_0mrad_pt_ratios_eta_name = path_to_csv_ratios + "/pt_ratios_" + production_date + suffix + "_MBD_pt_l4_0mrad.csv";
  std::string input_csv_mbd_0mrad_eta_ratios_pi0_name = path_to_csv_ratios + "/eta_ratios_" + production_date + suffix + "_MBD_pt_l3_0mrad.csv";
  std::string input_csv_mbd_0mrad_eta_ratios_eta_name = path_to_csv_ratios + "/eta_ratios_" + production_date + suffix + "_MBD_pt_l4_0mrad.csv";
  std::string input_csv_mbd_0mrad_xf_ratios_pi0_name = path_to_csv_ratios + "/xf_ratios_" + production_date + suffix + "_MBD_pt_l3_0mrad.csv";
  std::string input_csv_mbd_0mrad_xf_ratios_eta_name = path_to_csv_ratios + "/xf_ratios_" + production_date + suffix + "_MBD_pt_l4_0mrad.csv";

  std::string input_csv_photon_0mrad_pt_pi0_name = path_to_csv_mean + "/input_pt_mean_" + production_date + suffix + "_efficiency_30_pt_g3_0mrad.csv";
  std::string input_csv_photon_0mrad_pt_eta_name = path_to_csv_mean + "/input_pt_mean_" + production_date + suffix + "_efficiency_30_pt_g4_0mrad.csv";
  std::string input_csv_photon_0mrad_eta_pi0_name = path_to_csv_mean + "/input_eta_mean_" + production_date + suffix + "_efficiency_30_pt_g3_0mrad.csv";
  std::string input_csv_photon_0mrad_eta_eta_name = path_to_csv_mean + "/input_eta_mean_" + production_date + suffix + "_efficiency_30_pt_g4_0mrad.csv";
  std::string input_csv_photon_0mrad_xf_pi0_name = path_to_csv_mean + "/input_xf_mean_" + production_date + suffix + "_efficiency_30_pt_g3_0mrad.csv";
  std::string input_csv_photon_0mrad_xf_eta_name = path_to_csv_mean + "/input_xf_mean_" + production_date + suffix + "_efficiency_30_pt_g4_0mrad.csv";
  std::string input_csv_photon_0mrad_pt_ratios_pi0_name = path_to_csv_ratios + "/pt_ratios_" + production_date + suffix + "_efficiency_30_pt_g3_0mrad.csv";
  std::string input_csv_photon_0mrad_pt_ratios_eta_name = path_to_csv_ratios + "/pt_ratios_" + production_date + suffix + "_efficiency_30_pt_g4_0mrad.csv";
  std::string input_csv_photon_0mrad_eta_ratios_pi0_name = path_to_csv_ratios + "/eta_ratios_" + production_date + suffix + "_efficiency_30_pt_g3_0mrad.csv";
  std::string input_csv_photon_0mrad_eta_ratios_eta_name = path_to_csv_ratios + "/eta_ratios_" + production_date + suffix + "_efficiency_30_pt_g4_0mrad.csv";
  std::string input_csv_photon_0mrad_xf_ratios_pi0_name = path_to_csv_ratios + "/xf_ratios_" + production_date + suffix + "_efficiency_30_pt_g3_0mrad.csv";
  std::string input_csv_photon_0mrad_xf_ratios_eta_name = path_to_csv_ratios + "/xf_ratios_" + production_date + suffix + "_efficiency_30_pt_g4_0mrad.csv";

  std::string input_csv_mbd_15mrad_pt_pi0_name = path_to_csv_mean + "/input_pt_mean_" + production_date + suffix + "_MBD_pt_l3_15mrad.csv";
  std::string input_csv_mbd_15mrad_pt_eta_name = path_to_csv_mean + "/input_pt_mean_" + production_date + suffix + "_MBD_pt_l4_15mrad.csv";
  std::string input_csv_mbd_15mrad_eta_pi0_name = path_to_csv_mean + "/input_eta_mean_" + production_date + suffix + "_MBD_pt_l3_15mrad.csv";
  std::string input_csv_mbd_15mrad_eta_eta_name = path_to_csv_mean + "/input_eta_mean_" + production_date + suffix + "_MBD_pt_l4_15mrad.csv";
  std::string input_csv_mbd_15mrad_xf_pi0_name = path_to_csv_mean + "/input_xf_mean_" + production_date + suffix + "_MBD_pt_l3_15mrad.csv";
  std::string input_csv_mbd_15mrad_xf_eta_name = path_to_csv_mean + "/input_xf_mean_" + production_date + suffix + "_MBD_pt_l4_15mrad.csv";
  std::string input_csv_mbd_15mrad_pt_ratios_pi0_name = path_to_csv_ratios + "/pt_ratios_" + production_date + suffix + "_MBD_pt_l3_15mrad.csv";
  std::string input_csv_mbd_15mrad_pt_ratios_eta_name = path_to_csv_ratios + "/pt_ratios_" + production_date + suffix + "_MBD_pt_l4_15mrad.csv";
  std::string input_csv_mbd_15mrad_eta_ratios_pi0_name = path_to_csv_ratios + "/eta_ratios_" + production_date + suffix + "_MBD_pt_l3_15mrad.csv";
  std::string input_csv_mbd_15mrad_eta_ratios_eta_name = path_to_csv_ratios + "/eta_ratios_" + production_date + suffix + "_MBD_pt_l4_15mrad.csv";
  std::string input_csv_mbd_15mrad_xf_ratios_pi0_name = path_to_csv_ratios + "/xf_ratios_" + production_date + suffix + "_MBD_pt_l3_15mrad.csv";
  std::string input_csv_mbd_15mrad_xf_ratios_eta_name = path_to_csv_ratios + "/xf_ratios_" + production_date + suffix + "_MBD_pt_l4_15mrad.csv";

  std::string input_csv_photon_15mrad_pt_pi0_name = path_to_csv_mean + "/input_pt_mean_" + production_date + suffix + "_efficiency_30_pt_g3_15mrad.csv";
  std::string input_csv_photon_15mrad_pt_eta_name = path_to_csv_mean + "/input_pt_mean_" + production_date + suffix + "_efficiency_30_pt_g4_15mrad.csv";
  std::string input_csv_photon_15mrad_eta_pi0_name = path_to_csv_mean + "/input_eta_mean_" + production_date + suffix + "_efficiency_30_pt_g3_15mrad.csv";
  std::string input_csv_photon_15mrad_eta_eta_name = path_to_csv_mean + "/input_eta_mean_" + production_date + suffix + "_efficiency_30_pt_g4_15mrad.csv";
  std::string input_csv_photon_15mrad_xf_pi0_name = path_to_csv_mean + "/input_xf_mean_" + production_date + suffix + "_efficiency_30_pt_g3_15mrad.csv";
  std::string input_csv_photon_15mrad_xf_eta_name = path_to_csv_mean + "/input_xf_mean_" + production_date + suffix + "_efficiency_30_pt_g4_15mrad.csv";
  std::string input_csv_photon_15mrad_pt_ratios_pi0_name = path_to_csv_ratios + "/pt_ratios_" + production_date + suffix + "_efficiency_30_pt_g3_15mrad.csv";
  std::string input_csv_photon_15mrad_pt_ratios_eta_name = path_to_csv_ratios + "/pt_ratios_" + production_date + suffix + "_efficiency_30_pt_g4_15mrad.csv";
  std::string input_csv_photon_15mrad_eta_ratios_pi0_name = path_to_csv_ratios + "/eta_ratios_" + production_date + suffix + "_efficiency_30_pt_g3_15mrad.csv";
  std::string input_csv_photon_15mrad_eta_ratios_eta_name = path_to_csv_ratios + "/eta_ratios_" + production_date + suffix + "_efficiency_30_pt_g4_15mrad.csv";
  std::string input_csv_photon_15mrad_xf_ratios_pi0_name = path_to_csv_ratios + "/xf_ratios_" + production_date + suffix + "_efficiency_30_pt_g3_15mrad.csv";
  std::string input_csv_photon_15mrad_xf_ratios_eta_name = path_to_csv_ratios + "/xf_ratios_" + production_date + suffix + "_efficiency_30_pt_g4_15mrad.csv";

  int nIterations = (shuffleIterMax - shuffleIterMin + 1);
  int nStore = nIterations / shuffleMaxStorage;
  if (nIterations % shuffleMaxStorage != 0) {
    nStore++;
  }

  std::stringstream shuffle_file_name;
  shuffle_file_name << "shuffle_" << shuffleIterMin << "_" << shuffleIterMax << ".root";
  TFile *output_file_shuffle = nullptr;
  TH1F *h_AN_pt_shuffle[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nDirections] = {nullptr};
  TH1F *h_AN_pt_nodir_shuffle[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins] = {nullptr};
  TH1F *h_AN_eta_shuffle[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nEtaBins] = {nullptr};
  TH1F *h_AN_xf_shuffle[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nXfBins] = {nullptr};
  TH1F *h_AN_pt_shuffle_details[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nDirections][ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nConfigs] {nullptr};
  TH1F *h_AN_pt_nodir_shuffle_details[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nConfigs] {nullptr};
  TH1F *h_AN_eta_shuffle_details[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nEtaBins][ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nConfigs] {nullptr};
  TH1F *h_AN_xf_shuffle_details[ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nXfBins][ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nConfigs] {nullptr};
  book_shuffle_file(output_file_shuffle,
                    h_AN_pt_shuffle,
                    h_AN_pt_nodir_shuffle,
                    h_AN_eta_shuffle,
                    h_AN_xf_shuffle,
                    h_AN_pt_shuffle_details,
                    h_AN_pt_nodir_shuffle_details,
                    h_AN_eta_shuffle_details,
                    h_AN_xf_shuffle_details,
                    shuffle_file_name.str());
  
  int shuffleMin = 0;
  int shuffleMax = 0;
  for (int iStore = 0; iStore < nStore; iStore++) {
    
    shuffleMin = iStore * shuffleMaxStorage + shuffleIterMin;
    shuffleMax = std::min((iStore + 1) * shuffleMaxStorage - 1 + shuffleIterMin,
                          shuffleIterMax);
    std::cout << "shuffleMin, shuffleMax = " << shuffleMin << ", " << shuffleMax << std::endl;
    
    AsymmetryCalc::ShuffleBunches *sb_mbd_0mrad = new AsymmetryCalc::ShuffleBunches(
        shuffleMin,
        shuffleMax,
        path_to_fill_list_0mrad,
        path_to_spin_pattern_0mrad,
        path_to_bunch_dependent_yields_mbd,
        input_csv_mbd_0mrad_pt_pi0_name,
        input_csv_mbd_0mrad_eta_pi0_name,
        input_csv_mbd_0mrad_xf_pi0_name,
        input_csv_mbd_0mrad_pt_eta_name,
        input_csv_mbd_0mrad_eta_eta_name,
        input_csv_mbd_0mrad_xf_eta_name,
        input_csv_mbd_0mrad_pt_ratios_pi0_name,
        input_csv_mbd_0mrad_eta_ratios_pi0_name,
        input_csv_mbd_0mrad_xf_ratios_pi0_name,
        input_csv_mbd_0mrad_pt_ratios_eta_name,
        input_csv_mbd_0mrad_eta_ratios_eta_name,
        input_csv_mbd_0mrad_xf_ratios_eta_name);
    sb_mbd_0mrad->set_shuffle(perform_shuffle);
    sb_mbd_0mrad->set_store_fill_histos(store_raw, (stored_iter_asym_histos + "/mbd_0mrad_"));
    sb_mbd_0mrad->run_1();
    sb_mbd_0mrad->run_2();

    AsymmetryCalc::ShuffleBunches *sb_photon_0mrad = new AsymmetryCalc::ShuffleBunches(
        shuffleMin,
        shuffleMax,
        path_to_fill_list_0mrad,
        path_to_spin_pattern_0mrad,
        path_to_bunch_dependent_yields_photon,
        input_csv_photon_0mrad_pt_pi0_name,
        input_csv_photon_0mrad_eta_pi0_name,
        input_csv_photon_0mrad_xf_pi0_name,
        input_csv_photon_0mrad_pt_eta_name,
        input_csv_photon_0mrad_eta_eta_name,
        input_csv_photon_0mrad_xf_eta_name,
        input_csv_photon_0mrad_pt_ratios_pi0_name,
        input_csv_photon_0mrad_eta_ratios_pi0_name,
        input_csv_photon_0mrad_xf_ratios_pi0_name,
        input_csv_photon_0mrad_pt_ratios_eta_name,
        input_csv_photon_0mrad_eta_ratios_eta_name,
        input_csv_photon_0mrad_xf_ratios_eta_name);
    sb_photon_0mrad->set_shuffle(perform_shuffle);
    sb_photon_0mrad->set_store_fill_histos(store_raw, (stored_iter_asym_histos + "/photon_0mrad_"));
    sb_photon_0mrad->run_1();
    sb_photon_0mrad->run_2();

    AsymmetryCalc::ShuffleBunches *sb_mbd_15mrad = new AsymmetryCalc::ShuffleBunches(
        shuffleMin,
        shuffleMax,
        path_to_fill_list_15mrad,
        path_to_spin_pattern_15mrad,
        path_to_bunch_dependent_yields_mbd,
        input_csv_mbd_15mrad_pt_pi0_name,
        input_csv_mbd_15mrad_eta_pi0_name,
        input_csv_mbd_15mrad_xf_pi0_name,
        input_csv_mbd_15mrad_pt_eta_name,
        input_csv_mbd_15mrad_eta_eta_name,
        input_csv_mbd_15mrad_xf_eta_name,
        input_csv_mbd_15mrad_pt_ratios_pi0_name,
        input_csv_mbd_15mrad_eta_ratios_pi0_name,
        input_csv_mbd_15mrad_xf_ratios_pi0_name,
        input_csv_mbd_15mrad_pt_ratios_eta_name,
        input_csv_mbd_15mrad_eta_ratios_eta_name,
        input_csv_mbd_15mrad_xf_ratios_eta_name);
    sb_mbd_15mrad->set_shuffle(perform_shuffle);
    sb_mbd_15mrad->set_store_fill_histos(store_raw, (stored_iter_asym_histos + "/mbd_15mrad_"));
    sb_mbd_15mrad->run_1();
    sb_mbd_15mrad->run_2();

    AsymmetryCalc::ShuffleBunches *sb_photon_15mrad = new AsymmetryCalc::ShuffleBunches(
        shuffleMin,
        shuffleMax,
        path_to_fill_list_15mrad,
        path_to_spin_pattern_15mrad,
        path_to_bunch_dependent_yields_photon,
        input_csv_photon_15mrad_pt_pi0_name,
        input_csv_photon_15mrad_eta_pi0_name,
        input_csv_photon_15mrad_xf_pi0_name,
        input_csv_photon_15mrad_pt_eta_name,
        input_csv_photon_15mrad_eta_eta_name,
        input_csv_photon_15mrad_xf_eta_name,
        input_csv_photon_15mrad_pt_ratios_pi0_name,
        input_csv_photon_15mrad_eta_ratios_pi0_name,
        input_csv_photon_15mrad_xf_ratios_pi0_name,
        input_csv_photon_15mrad_pt_ratios_eta_name,
        input_csv_photon_15mrad_eta_ratios_eta_name,
        input_csv_photon_15mrad_xf_ratios_eta_name);
    sb_photon_15mrad->set_shuffle(perform_shuffle);
    sb_photon_15mrad->set_store_fill_histos(store_raw, (stored_iter_asym_histos + "/photon_15mrad_"));
    sb_photon_15mrad->run_1();
    sb_photon_15mrad->run_2();

    int local_number_shuffles = (shuffleMax - shuffleMin + 1);
    for (int iIter = 0; iIter < local_number_shuffles; iIter++) {
      AsymmetryCalc::pt_final_array merged_mean_pt;
      AsymmetryCalc::pt_nodir_final_array merged_mean_pt_nodir;
      AsymmetryCalc::eta_final_array merged_mean_eta;
      AsymmetryCalc::xf_final_array merged_mean_xf;
      AsymmetryCalc::pt_final_array merged_unc_pt;
      AsymmetryCalc::pt_nodir_final_array merged_unc_pt_nodir;
      AsymmetryCalc::eta_final_array merged_unc_eta;
      AsymmetryCalc::xf_final_array merged_unc_xf;
      merge_configurations(
          sb_mbd_0mrad,
          sb_photon_0mrad,
          sb_mbd_15mrad,
          sb_photon_15mrad,
          merged_mean_pt,
          merged_mean_pt_nodir,
          merged_mean_eta,
          merged_mean_xf,
          merged_unc_pt,
          merged_unc_pt_nodir,
          merged_unc_eta,
          merged_unc_xf,
          iIter);

      fill_shuffle_detailed_histos(
         sb_mbd_0mrad,
         sb_photon_0mrad,
         sb_mbd_15mrad,
         sb_photon_15mrad,
         h_AN_pt_shuffle_details,
         h_AN_pt_nodir_shuffle_details,
         h_AN_eta_shuffle_details,
         h_AN_xf_shuffle_details,
         iIter);

      // Fill shuffle histogram

      if (store_final) {
        std::string output_final_name = stored_final_asym_histos;
        save_final_histograms(
         output_final_name,
         merged_mean_pt,
         merged_mean_pt_nodir,
         merged_mean_eta,
         merged_mean_xf,
         merged_unc_pt,
         merged_unc_pt_nodir,
         merged_unc_eta,
         merged_unc_xf);
      }

      // After filling shuffle histos, "merged" vectors are cleared automatically
      fill_shuffle_histos(merged_mean_pt,
                          merged_mean_pt_nodir,
                          merged_mean_eta,
                          merged_mean_xf,
                          merged_unc_pt,
                          merged_unc_pt_nodir,
                          merged_unc_eta,
                          merged_unc_xf,
                          h_AN_pt_shuffle,
                          h_AN_pt_nodir_shuffle,
                          h_AN_eta_shuffle,
                          h_AN_xf_shuffle);
    }

    delete sb_mbd_0mrad;
    delete sb_photon_0mrad;
    delete sb_mbd_15mrad;
    delete sb_photon_15mrad;
  }

  save_shuffle_histos(output_file_shuffle);
  
  return 0;
}
