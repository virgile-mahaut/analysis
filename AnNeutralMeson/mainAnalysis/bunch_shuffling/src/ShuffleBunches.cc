#include "AsymmetryCalc/ShuffleBunches.h"
#include "AsymmetryCalc/Constants.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <cstdlib>
#include <numeric> // for iota

namespace AsymmetryCalc {
  
  ShuffleBunches::ShuffleBunches(const int iMin,
                                 const int iMax,
                                 const std::string& inputfilename_fill_run,
                                 const std::string& inputfilename_spin_pattern,
                                 const std::string& inputfile_template,
                                 const std::string& input_csv_pt_pi0_name,
                                 const std::string& input_csv_eta_pi0_name,
                                 const std::string& input_csv_xf_pi0_name,
                                 const std::string& input_csv_pt_eta_name,
                                 const std::string& input_csv_eta_eta_name,
                                 const std::string& input_csv_xf_eta_name,
                                 const std::string& input_csv_pt_ratios_pi0_name,
                                 const std::string& input_csv_eta_ratios_pi0_name,
                                 const std::string& input_csv_xf_ratios_pi0_name,
                                 const std::string& input_csv_pt_ratios_eta_name,
                                 const std::string& input_csv_eta_ratios_eta_name,
                                 const std::string& input_csv_xf_ratios_eta_name)
    : iterMin(iMin),
      iterMax(iMax),
      infile_template(inputfile_template)
  {
    nIterations = iterMax - iterMin + 1;
    std::ifstream infile_fill_run(inputfilename_fill_run);
    if (!infile_fill_run) {
      std::cerr << "Could not open file [" << inputfilename_fill_run << "]" << std::endl;
      std::exit(1);
    }

    int fill_number, run_number;
    while (infile_fill_run >> fill_number >> run_number) {
      if (fill_to_runs.empty() || fill_to_runs.back().first != fill_number)
      {
        fill_to_runs.push_back({fill_number, {}});
      }
      fill_to_runs.back().second.push_back(run_number);
    }
    nFills = fill_to_runs.size();


    // Read average bin value
    get_average_bin_pt_pi0(input_csv_pt_pi0_name);
    get_average_bin_eta_pi0(input_csv_eta_pi0_name);
    get_average_bin_xf_pi0(input_csv_xf_pi0_name);
    get_average_bin_pt_eta(input_csv_pt_eta_name);
    get_average_bin_eta_eta(input_csv_eta_eta_name);
    get_average_bin_xf_eta(input_csv_xf_eta_name);

    // Read background ratio
    get_pt_ratios_pi0(input_csv_pt_ratios_pi0_name);
    get_eta_ratios_pi0(input_csv_eta_ratios_pi0_name);
    get_xf_ratios_pi0(input_csv_xf_ratios_pi0_name);
    get_pt_ratios_eta(input_csv_pt_ratios_eta_name);
    get_eta_ratios_eta(input_csv_eta_ratios_eta_name);
    get_xf_ratios_eta(input_csv_xf_ratios_eta_name);

    // Read Spin Pattern
    inputfile_spin = TFile::Open(inputfilename_spin_pattern.c_str());
    spin_patterns = (TTree*)inputfile_spin->Get("spin_patterns");
    spin_patterns->SetBranchAddress("fill_number", &fill_spin);
    spin_patterns->SetBranchAddress("yellow_polarization", &yellow_polarization);
    spin_patterns->SetBranchAddress("blue_polarization", &blue_polarization);
    spin_patterns->SetBranchAddress("yellow_spin_pattern", &yellow_spin_pattern);
    spin_patterns->SetBranchAddress("blue_spin_pattern", &blue_spin_pattern);
  }

  void ShuffleBunches::get_average_bin_pt_pi0(const std::string& input_csv_name)
  {
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      std::getline(line_pi0, entry, ',');
      float pTMean = std::stof(entry);
      pTMeans[0][iPt] = pTMean;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      std::getline(line_eta, entry, ','); // skip
    }
  }

  void ShuffleBunches::get_average_bin_pt_eta(const std::string& input_csv_name)
  {
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      std::getline(line_pi0, entry, ','); // skip
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      std::getline(line_eta, entry, ',');
      float pTMean = std::stof(entry);
      pTMeans[1][iPt] = pTMean;
    }
  }

  void ShuffleBunches::get_average_bin_eta_pi0(const std::string& input_csv_name)
  {
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      std::getline(line_pi0, entry, ',');
      float etaMean = std::stof(entry);
      etaMeans[0][iEta] = etaMean;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      std::getline(line_eta, entry, ','); // skip
    }
  }

  void ShuffleBunches::get_average_bin_eta_eta(const std::string& input_csv_name)
  {
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      std::getline(line_pi0, entry, ','); // skip
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      std::getline(line_eta, entry, ',');
      float etaMean = std::stof(entry);
      etaMeans[1][iEta] = etaMean;
    }
  }

  void ShuffleBunches::get_average_bin_xf_pi0(const std::string& input_csv_name)
  {
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      std::getline(line_pi0, entry, ',');
      float xfMean = std::stof(entry);
      xfMeans[0][iXf] = xfMean;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      std::getline(line_eta, entry, ','); // skip
    }
  }

  void ShuffleBunches::get_average_bin_xf_eta(const std::string& input_csv_name)
  {
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      std::getline(line_pi0, entry, ','); // skip
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      std::getline(line_eta, entry, ',');
      float xfMean = std::stof(entry);
      xfMeans[1][iXf] = xfMean;
    }
  }
  
  void ShuffleBunches::get_pt_ratios_pi0(const std::string& input_csv_name)
  {
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      std::getline(line_pi0, entry, ',');
      float pTRatio = std::stof(entry);
      bkg_ratio_pt[0][iPt] = pTRatio;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      std::getline(line_eta, entry, ','); // skip
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0_err(line);
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      std::getline(line_pi0_err, entry, ',');
      float pTRatioErr = std::stof(entry);
      bkg_ratio_err_pt[0][iPt] = pTRatioErr;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta_err(line);
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      std::getline(line_eta_err, entry, ','); // skip
    }
  }

  void ShuffleBunches::get_pt_ratios_eta(const std::string& input_csv_name)
  {
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      std::getline(line_pi0, entry, ','); // skip
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      std::getline(line_eta, entry, ',');
      float pTRatio = std::stof(entry);
      bkg_ratio_pt[1][iPt] = pTRatio;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0_err(line);
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      std::getline(line_pi0_err, entry, ','); // skip
    }
    std::getline(input_csv, line);
    std::stringstream line_eta_err(line);
    for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
      std::getline(line_eta_err, entry, ',');
      float pTRatioErr = std::stof(entry);
      bkg_ratio_err_pt[1][iPt] = pTRatioErr;
    }
  }

  void ShuffleBunches::get_eta_ratios_pi0(const std::string& input_csv_name)
  {
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      std::getline(line_pi0, entry, ',');
      float etaRatio = std::stof(entry);
      bkg_ratio_eta[0][iEta] = etaRatio;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      std::getline(line_eta, entry, ','); // skip
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0_err(line);
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      std::getline(line_pi0_err, entry, ',');
      float etaRatioErr = std::stof(entry);
      bkg_ratio_err_eta[0][iEta] = etaRatioErr;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta_err(line);
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      std::getline(line_eta_err, entry, ','); // skip
    }
  }

  void ShuffleBunches::get_eta_ratios_eta(const std::string& input_csv_name)
  {
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      std::getline(line_pi0, entry, ','); // skip
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      std::getline(line_eta, entry, ',');
      float etaRatio = std::stof(entry);
      bkg_ratio_eta[1][iEta] = etaRatio;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0_err(line);
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      std::getline(line_pi0_err, entry, ','); // skip
    }
    std::getline(input_csv, line);
    std::stringstream line_eta_err(line);
    for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
      std::getline(line_eta_err, entry, ',');
      float etaRatioErr = std::stof(entry);
      bkg_ratio_err_eta[1][iEta] = etaRatioErr;
    }
  }

  void ShuffleBunches::get_xf_ratios_pi0(const std::string& input_csv_name)
  {
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      std::getline(line_pi0, entry, ',');
      float xfRatio = std::stof(entry);
      bkg_ratio_xf[0][iXf] = xfRatio;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      std::getline(line_eta, entry, ','); // skip
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0_err(line);
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      std::getline(line_pi0_err, entry, ',');
      float xfRatioErr = std::stof(entry);
      bkg_ratio_err_xf[0][iXf] = xfRatioErr;
    }
    std::getline(input_csv, line);
    std::stringstream line_eta_err(line);
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      std::getline(line_eta_err, entry, ','); // skip
    }
  }

  void ShuffleBunches::get_xf_ratios_eta(const std::string& input_csv_name)
  {
    std::string line;
    std::stringstream linestream;
    std::string entry;
    std::ifstream input_csv;
    input_csv.open(input_csv_name);
    if (!input_csv) {
      std::cout << "File " << input_csv_name << " could not be opened." << std::endl;
      return;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0(line);
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      std::getline(line_pi0, entry, ','); // skip
    }
    std::getline(input_csv, line);
    std::stringstream line_eta(line);
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      std::getline(line_eta, entry, ',');
      float xfRatio = std::stof(entry);
      bkg_ratio_xf[1][iXf] = xfRatio;
    }
    std::getline(input_csv, line);
    std::stringstream line_pi0_err(line);
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      std::getline(line_pi0_err, entry, ','); // skip
    }
    std::getline(input_csv, line);
    std::stringstream line_eta_err(line);
    for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
      std::getline(line_eta_err, entry, ',');
      float xfRatioErr = std::stof(entry);
      bkg_ratio_err_xf[1][iXf] = xfRatioErr;
    }
  }

  ShuffleBunches::~ShuffleBunches()
  {
    clear_fill_yields();
    clear_fill_raw_asyms();
    clear_average_raw_asyms();
    clear_corr_asyms();
    clear_fit_asyms();
  }

  void ShuffleBunches::clear_fill_yields()
  {
    vec_yield_pt.clear();
    vec_yield_eta.clear();
    vec_yield_xf.clear();
  }

  void ShuffleBunches::clear_fill_raw_asyms()
  {
    fill_pt_asyms.clear();
    fill_eta_asyms.clear();
    fill_xf_asyms.clear();
    fill_pt_asym_errs.clear();
    fill_eta_asym_errs.clear();
    fill_xf_asym_errs.clear();
  }

  void ShuffleBunches::clear_average_raw_asyms()
  {
    mean_pt.clear();
    mean_eta.clear();
    mean_xf.clear();
    unc_pt.clear();
    unc_eta.clear();
    unc_xf.clear();
  }

  void ShuffleBunches::clear_corr_asyms()
  {
    corrected_mean_pt.clear();
    corrected_mean_eta.clear();
    corrected_mean_xf.clear();
    corrected_unc_pt.clear();
    corrected_unc_eta.clear();
    corrected_unc_xf.clear();
  }

  void ShuffleBunches::clear_fit_asyms()
  {
    fit_mean_pt.clear();
    fit_mean_eta.clear();
    fit_mean_xf.clear();
    fit_unc_pt.clear();
    fit_unc_eta.clear();
    fit_unc_xf.clear();
  }

  void ShuffleBunches::reserve_fill_yields()
  {
    vec_yield_pt.reserve(nIterations);
    vec_yield_eta.reserve(nIterations);
    vec_yield_xf.reserve(nIterations);
  }


  void ShuffleBunches::reserve_fill_raw_asyms()
  {
    fill_pt_asyms.reserve(nFills * nIterations);
    fill_eta_asyms.reserve(nFills * nIterations);
    fill_xf_asyms.reserve(nFills * nIterations);
    fill_pt_asym_errs.reserve(nFills * nIterations);
    fill_eta_asym_errs.reserve(nFills * nIterations);
    fill_xf_asym_errs.reserve(nFills * nIterations);
  }
  
  void ShuffleBunches::reserve_average_raw_asyms()
  {
    mean_pt.reserve(nIterations);
    mean_eta.reserve(nIterations);
    mean_xf.reserve(nIterations);
    unc_pt.reserve(nIterations);
    unc_eta.reserve(nIterations);
    unc_xf.reserve(nIterations);
  }

  void ShuffleBunches::reserve_corr_asyms()
  {
    corrected_mean_pt.reserve(nIterations);
    corrected_mean_eta.reserve(nIterations);
    corrected_mean_xf.reserve(nIterations);
    corrected_unc_pt.reserve(nIterations);
    corrected_unc_eta.reserve(nIterations);
    corrected_unc_xf.reserve(nIterations);
  }

  void ShuffleBunches::reserve_fit_asyms()
  {
    fit_mean_pt.reserve(nIterations);
    fit_mean_eta.reserve(nIterations);
    fit_mean_xf.reserve(nIterations);
    fit_unc_pt.reserve(nIterations);
    fit_unc_eta.reserve(nIterations);
    fit_unc_xf.reserve(nIterations);
  }

  void ShuffleBunches::reset()
  {
    clear_fill_yields();
    clear_fill_raw_asyms();
    clear_average_raw_asyms();
    clear_corr_asyms();
    clear_fit_asyms();
  }

  void ShuffleBunches::run_1()
  {
    reserve_fill_yields();
    reserve_fill_raw_asyms();
    global_irun = 0;
    std::cout << "Compute " << nFills << " fills" << std::endl;
    for (int iFill = 0; iFill < nFills; iFill++)
    {
      compute_fill(iFill);
    }
  }

  void ShuffleBunches::run_2()
  {
    reserve_average_raw_asyms();
    for (int iIter = 0; iIter < nIterations; iIter++) {
      std::cout << "average iter " << (iterMin + iIter) << std::endl;
      if (store_iter_histos) {
        std::stringstream outputfilename;
        outputfilename << outfilename_iter_template << "raw_" << (iterMin + iIter) << ".root";
        book_outfile_iter_histograms(outputfilename.str().c_str());
      }
      average_asyms(iIter);
      if (store_iter_histos) {
        get_iter_raw_asym_histograms();
        save_outfile_iter_histograms();
      }
    }
    clear_fill_raw_asyms();
    reserve_fit_asyms();
    for (int iIter = 0; iIter < nIterations; iIter++) {
      std::cout << "fit iter " << (iterMin + iIter) << std::endl;
      if (store_iter_histos) {
        std::stringstream outputfilename;
        outputfilename << outfilename_iter_template << "fit_" << (iterMin + iIter) << ".root";
        book_outfile_iter_histograms(outputfilename.str().c_str());
      }
      fit_raw_asyms(iIter);
      if (store_iter_histos) {
        get_iter_fit_asym_histograms();
        save_outfile_iter_histograms();
      }
    }
    clear_average_raw_asyms();
    reserve_corr_asyms();
    for (int iIter = 0; iIter < nIterations; iIter++) {
      std::cout << "corr iter " << (iterMin + iIter) << std::endl;
      if (store_iter_histos) {
        std::stringstream outputfilename;
        outputfilename << outfilename_iter_template << "corr_" << (iterMin + iIter) << ".root";
        book_outfile_iter_histograms(outputfilename.str().c_str());
      }
      compute_corrected_asyms(iIter);
      if (store_iter_histos) {
        get_iter_corr_asym_histograms();
        save_outfile_iter_histograms();
      }
    }
  }

  void ShuffleBunches::set_store_fill_histos(bool val, const std::string& outputfilename_template)
  {
    if (val && nIterations > 1)
    {
      std::cerr << "Error. Can't store histograms for more than one shuffle iteration" << std::endl;
      store_fill_histos = false;
    }
    else if (val)
    {
      store_fill_histos = true;
      outfilename_fill_template = outputfilename_template;
    }
    else
    {
      store_fill_histos = false;
    }
  }

  void ShuffleBunches::set_store_iter_histos(bool val, const std::string& outputfilename_template)
  {
    if (val && nIterations > 1)
    {
      std::cerr << "Error. Can't store histograms for more than one shuffle iteration" << std::endl;
      store_iter_histos = false;
    }
    else if (val)
    {
      store_iter_histos = true;
      outfilename_iter_template = outputfilename_template;
    }
    else
    {
      store_iter_histos = false;
    }
  }

  void ShuffleBunches::book_outfile_fill_histograms(const std::string& outputfilename)
  {
    // Book histograms
    outfile_fill_histograms = new TFile(outputfilename.c_str(), "RECREATE");
  }

  void ShuffleBunches::book_outfile_iter_histograms(const std::string& outputfilename)
  {
    // Book histograms
    outfile_iter_histograms = new TFile(outputfilename.c_str(), "RECREATE");
  }

  void ShuffleBunches::save_outfile_fill_histograms()
  {
    outfile_fill_histograms->cd();
    outfile_fill_histograms->Write();
    outfile_fill_histograms->Close();
    delete outfile_fill_histograms;
    outfile_fill_histograms = nullptr;
  }

  void ShuffleBunches::save_outfile_iter_histograms()
  {
    outfile_iter_histograms->cd();
    outfile_iter_histograms->Write();
    outfile_iter_histograms->Close();
    delete outfile_iter_histograms;
    outfile_iter_histograms = nullptr;
  }
  
  void ShuffleBunches::get_yield_histograms()
  {
    pt_yield_array &arr_pt_yield = vec_yield_pt.back();
    eta_yield_array &arr_eta_yield = vec_yield_eta.back();
    xf_yield_array &arr_xf_yield = vec_yield_xf.back();
    
    outfile_fill_histograms->cd();
    for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++)
    {
      for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++)
      {
        for (int iR = 0; iR < ASYM_CONSTANTS::nRegions; iR++)
        {
          for (int iS = 0; iS < ASYM_CONSTANTS::nSpins; iS++)
          {
            for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++)
            {
              for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++)
              {
                // Book histogram
                std::stringstream h_yield_name;
                h_yield_name << "h_yield_" << ASYM_CONSTANTS::beams[iB] << "_" << ASYM_CONSTANTS::particle[iP]
                             << "_" << ASYM_CONSTANTS::regions[iR] << "_"
                             << "pT_" << iPt << "_"
                             << ASYM_CONSTANTS::directions[iDir]
                             << "_" << ASYM_CONSTANTS::spins[iS];
                std::stringstream h_yield_title;
                h_yield_title << "#phi_{" << (iB == 0 ? "yellow" : "blue")
                              << "};counts";
                h_yield_pt[iB][iP][iR][iPt][iDir][iS] =
                  new TH1F(h_yield_name.str().c_str(),
                           h_yield_title.str().c_str(), 12, -M_PI, M_PI);

                // Fill histogram
                for (int iphi = 0; iphi < ASYM_CONSTANTS::nPhiBins; iphi++)
                {
                  h_yield_pt[iB][iP][iR][iPt][iDir][iS]->SetBinContent(
                    iphi + 1,
                    arr_pt_yield(iB, iP, iR, iPt, iDir, iS, iphi));
                }
              }
            }
            for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++)
            {
              // Book histogram
              std::stringstream h_yield_name;
              h_yield_name << "h_yield_" << ASYM_CONSTANTS::beams[iB] << "_" << ASYM_CONSTANTS::particle[iP]
                           << "_" << ASYM_CONSTANTS::regions[iR] << "_"
                           << "eta_" << (int) iEta << "_" << ASYM_CONSTANTS::spins[iS];
              std::stringstream h_yield_title;
              h_yield_title << "#phi_{" << (iB == 0 ? "yellow" : "blue")
                            << "};counts";
              h_yield_eta[iB][iP][iR][iEta][iS] =
                  new TH1F(h_yield_name.str().c_str(),
                           h_yield_title.str().c_str(), 12, -M_PI, M_PI);

              // Fill histogram
              for (int iphi = 0; iphi < ASYM_CONSTANTS::nPhiBins; iphi++)
              {
                h_yield_eta[iB][iP][iR][iEta][iS]->SetBinContent(
                  iphi + 1,
                  arr_eta_yield(iB, iP, iR, iEta, iS, iphi));
              }
            }
            for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++)
            {
              // Book histogram
              std::stringstream h_yield_name;
              h_yield_name << "h_yield_" << ASYM_CONSTANTS::beams[iB] << "_" << ASYM_CONSTANTS::particle[iP]
                           << "_" << ASYM_CONSTANTS::regions[iR] << "_"
                           << "xf_" << (int) iXf << "_" << ASYM_CONSTANTS::spins[iS];
              std::stringstream h_yield_title;
              h_yield_title << "#phi_{" << (iB == 0 ? "yellow" : "blue")
                            << "};counts";
              h_yield_xf[iB][iP][iR][iXf][iS] =
                  new TH1F(h_yield_name.str().c_str(),
                           h_yield_title.str().c_str(), 12, -M_PI, M_PI);
              
              // Fill histogram
              for (int iphi = 0; iphi < ASYM_CONSTANTS::nPhiBins; iphi++)
              {
                h_yield_xf[iB][iP][iR][iXf][iS]->SetBinContent(
                  iphi + 1,
                  arr_xf_yield(iB, iP, iR, iXf, iS, iphi));
              }
            }
          }
        }
      }
    }
  }

  void ShuffleBunches::get_fill_raw_asym_histograms()
  {
    pt_asym_array& arr_pt_asym = fill_pt_asyms.back();
    eta_asym_array& arr_eta_asym = fill_eta_asyms.back();
    xf_asym_array& arr_xf_asym = fill_xf_asyms.back();
    pt_asym_array& arr_pt_asym_err = fill_pt_asym_errs.back();
    eta_asym_array& arr_eta_asym_err = fill_eta_asym_errs.back();
    xf_asym_array& arr_xf_asym_err = fill_xf_asym_errs.back();
    
    outfile_fill_histograms->cd();
    for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++)
    {
      for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++)
      {
        for (int iR = 0; iR < ASYM_CONSTANTS::nRegions; iR++)
        {
          for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++)
          {
            for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++)
            {
              TGraphErrors *graph_AN = new TGraphErrors();
              std::stringstream graph_AN_name;
              graph_AN_name << "graph_asym_pt_"
                            << ASYM_CONSTANTS::beams[iB] << "_" << ASYM_CONSTANTS::particle[iP] << "_"
                            << ASYM_CONSTANTS::regions[iR] << "_"
                            << "pT_" << iPt << "_"
                            << ASYM_CONSTANTS::directions[iDir] << "_"
                            << "sqrt";
              graph_AN->SetName(graph_AN_name.str().c_str());
              std::stringstream graph_title;
              graph_title << "; #phi; " << (iR == 0 ? "Raw" : "Bkg") << " Asymmetry";
              graph_AN->SetTitle(graph_title.str().c_str());
              for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++)
              {
                float phi = -M_PI + M_PI / 12 + iPhi * M_PI / 6.;
                graph_AN->SetPoint(iPhi, phi, arr_pt_asym(iB, iP, iR, iPt, iDir, iPhi));
                graph_AN->SetPointError(iPhi, 0, arr_pt_asym_err(iB, iP, iR, iPt, iDir, iPhi));
              }
              graph_AN->Write();
            }
          }
          for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++)
          {
            TGraphErrors *graph_AN = new TGraphErrors();
            std::stringstream graph_AN_name;
            graph_AN_name << "graph_asym_eta_"
                          << ASYM_CONSTANTS::beams[iB] << "_" << ASYM_CONSTANTS::particle[iP] << "_"
                          << ASYM_CONSTANTS::regions[iR] << "_"
                          << "eta_" << iEta << "_sqrt";
            graph_AN->SetName(graph_AN_name.str().c_str());
            std::stringstream graph_title;
            graph_title << "; #phi; " << (iR == 0 ? "Raw" : "Bkg") << " Asymmetry";
            graph_AN->SetTitle(graph_title.str().c_str());
            for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++)
            {
              float phi = -M_PI + M_PI / 12 + iPhi * M_PI / 6.;
              graph_AN->SetPoint(iPhi, phi, arr_eta_asym(iB, iP, iR, iEta, iPhi));
              graph_AN->SetPointError(iPhi, 0, arr_eta_asym_err(iB, iP, iR, iEta, iPhi));
            }
            graph_AN->Write();
          }
          for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++)
          {
            TGraphErrors *graph_AN = new TGraphErrors();
            std::stringstream graph_AN_name;
            graph_AN_name << "graph_asym_xf_"
                          << ASYM_CONSTANTS::beams[iB] << "_" << ASYM_CONSTANTS::particle[iP] << "_"
                          << ASYM_CONSTANTS::regions[iR] << "_"
                          << "xf_" << iXf << "_sqrt";
            graph_AN->SetName(graph_AN_name.str().c_str());
            std::stringstream graph_title;
            graph_title << "; #phi; " << (iR == 0 ? "Raw" : "Bkg") << " Asymmetry";
            graph_AN->SetTitle(graph_title.str().c_str());

            // Fill
            for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++)
            {
              float phi = -M_PI + M_PI / 12 + iPhi * M_PI / 6.;
              graph_AN->SetPoint(iPhi, phi, arr_xf_asym(iB, iP, iR, iXf, iPhi));
              graph_AN->SetPointError(iPhi, 0, arr_xf_asym_err(iB, iP, iR, iXf, iPhi));
            }
            graph_AN->Write();
          }
        }
      }
    }
  }
  
  void ShuffleBunches::compute_fill(const int ifill)
  {
    //std::cout << "Compute fill " << ifill << std::endl;
    //std::cout << "yields.size = " << vec_yield_xf.size() << std::endl;
    const int fill_number = fill_to_runs[ifill].first;
    const std::vector<int>& runs = fill_to_runs[ifill].second;
    const int nRuns = runs.size();

    // Get fill-dependent spin pattern
    spin_patterns->GetEntry((unsigned long)ifill);
    if (fill_number != fill_spin) {
      std::cerr << "Mismatch in fill number: " << fill_number << " != " << fill_spin << std::endl;
      return;
    }
    for (int i = 0; i < ASYM_CONSTANTS::nBunches; i++) {
      spin_pattern[0][i] = yellow_spin_pattern[i];
      spin_pattern[1][i] = blue_spin_pattern[i];
    }
    polarization[0] = yellow_polarization;
    polarization[1] = blue_polarization;

    if (store_fill_histos)
    {
      std::stringstream outputfilename;
      outputfilename << outfilename_fill_template << fill_number << ".root";
      book_outfile_fill_histograms(outputfilename.str());
    }

    for (int iRun = 0; iRun < nRuns; iRun++)
    {
      std::string infilename = infile_template + std::to_string(runs[iRun]) + ".bin";
      read_binary_array(infilename);
      for (int iIter = 0; iIter < nIterations; iIter++)
      {
        int iter = iterMin + iIter;
        // Shuffle bunch indices and flip spin
        int seednumber = (do_shuffle ? (iter + 100000 * (global_irun + iRun)) : 0);
        shuffle_bunch_indices(seednumber);
        
        // Sum shuffled bunches
        sum_yields(iIter);
      }
    }

    //std::cout << "yields.size 2 = " << vec_yield_xf.size() << std::endl;

    if (store_fill_histos) {
      get_yield_histograms();
    }

    for (int iIter = 0; iIter < nIterations; iIter++)
    {
      // Get (Square Root) Raw Asymmetries
      compute_raw_asyms(iIter);
    }
    clear_fill_yields();

    //std::cout << "yields.size 3 = " << vec_yield_xf.size() << std::endl;

    global_irun += nRuns;

    if (store_fill_histos)
    {
      get_fill_raw_asym_histograms();
      save_outfile_fill_histograms(); // includes deletion
    }
  }
  
  void ShuffleBunches::read_binary_array(const std::string& inputfilename)
  {
    std::ifstream inbinary(inputfilename, std::ios::binary);

    if (!inbinary.is_open()) {
      std::cerr << "Error: Cannot open input file " << inputfilename << std::endl;
      exit(1);
    }
    
    //std::cout << "size of pt array = " << sizeof(array_yield_pt_bunches) << std::endl;
    inbinary.read(reinterpret_cast<char*>(array_yield_pt_bunches), sizeof(array_yield_pt_bunches));
    inbinary.read(reinterpret_cast<char*>(array_yield_eta_bunches), sizeof(array_yield_eta_bunches));
    inbinary.read(reinterpret_cast<char*>(array_yield_xf_bunches), sizeof(array_yield_xf_bunches));
    
    uint64_t sum_yield_pt = 0;
    uint64_t sum_yield_eta = 0;
    uint64_t sum_yield_xf = 0;
    for (int iS = 0; iS < 2; iS++) {
      for (int iBunch = 0; iBunch < ASYM_CONSTANTS::nBunches; ++iBunch) {
        for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; ++iB) {
          for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; ++iP) {
            for (int iR = 0; iR < ASYM_CONSTANTS::nRegions; ++iR) {
              for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; ++iPhi) {
                for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; ++iPt) {
                  for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; ++iDir) {
                    sum_yield_pt += array_yield_pt_bunches[iB][iP][iR][iPt][iDir][iS][iPhi][iBunch];
                  }
                }
                for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; ++iEta) {
                  sum_yield_eta += array_yield_eta_bunches[iB][iP][iR][iEta][iS][iPhi][iBunch];
                }
                for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; ++iXf) {
                  sum_yield_xf += array_yield_xf_bunches[iB][iP][iR][iXf][iS][iPhi][iBunch];
                }
              }
            }
          }
        }
      }
    }

    // std::cout << "input file = " << inputfilename << std::endl;
    // std::cout << "yields = (" << sum_yield_pt << "," << sum_yield_eta << "," << sum_yield_xf << ")" << std::endl;
    
  }
  
  void ShuffleBunches::shuffle_bunch_indices(const int seednumber)
  {
    // Shuffle the spin pattern if the seed number is non zero
    // Only consider the first 111/(nBunches = 120) bunches, the last 9 are
    // always empty.

    // Reset indices before the shuffle
    std::iota(shuffled_indices, shuffled_indices + ASYM_CONSTANTS::nBunches, 0);
    if (seednumber != 0)
    {
      // Pseudo-random generator: Mersenne twister with a period of 2^19937 - 1
      //std::cout << "seed = " << seednumber << std::endl;
      std::mt19937 rng(seednumber);
      std::shuffle(shuffled_indices, shuffled_indices + ASYM_CONSTANTS::nBunches, rng);

      // With a probability 1/2, flip all the bunch spins
      std::uniform_int_distribution<std::mt19937::result_type> dist2(0, 1);
      int factor = 2 * (int) dist2(rng) - 1;
      spin_flip = (factor == -1 ? true : false);
    }
    
    //std::cout << "shuffled_indices = ";
    // for (int i = 0; i < ASYM_CONSTANTS::nBunches; i++)
    //   std::cout << shuffled_indices[i] << " ";
    // std::cout << std::endl;
  }

  void ShuffleBunches::sum_yields(const int iIter)
  {
    pt_yield_array arr_yield_pt;
    eta_yield_array arr_yield_eta;
    xf_yield_array arr_yield_xf;
    
    // Use the shuffled spins
    for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++) {
      for (int iBunch = 0; iBunch < ASYM_CONSTANTS::nBunches; iBunch++) {
        shuffled_spin_pattern[iB][iBunch] = spin_pattern[iB][shuffled_indices[iBunch]];
      }
    }
    
    auto spin_to_index = [](int8_t s) -> int {
      if (s == 1)  return 0;
      if (s == -1) return 1;
      std::cerr << "Error for spin value (expected +1 or -1)" << std::endl;
      exit(1);
    };

    for (int iBunch = 0; iBunch < ASYM_CONSTANTS::nBunches; ++iBunch) {
      for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; ++iB) {
        const int8_t true_spin = spin_pattern[iB][iBunch];

        int8_t shuffled_spin = shuffled_spin_pattern[iB][iBunch];
        if (spin_flip) shuffled_spin = -shuffled_spin;

        const int iS_true    = spin_to_index(true_spin);
        const int iS_shuffle = spin_to_index(shuffled_spin);

        for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; ++iP) {
          for (int iR = 0; iR < ASYM_CONSTANTS::nRegions; ++iR) {
            for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; ++iPhi) {
              for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; ++iPt) {
                for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; ++iDir) {
                  arr_yield_pt(iB, iP, iR, iPt, iDir, iS_shuffle, iPhi) +=
                    array_yield_pt_bunches[iB][iP][iR][iPt][iDir][iS_true][iPhi][iBunch];
                  // std::cout << "array_yield_pt_bunches[" << iB << "][" << iP << "][" << iR << "][" << iPt << "][" << iDir << "][" << iS_true << "][" << iPhi << "][" << iBunch << "] = "
                  //           << array_yield_pt_bunches[iB][iP][iR][iPt][iDir][iS_true][iPhi][iBunch]
                  //           << std::endl;
                }
              }

              for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; ++iEta) {
                arr_yield_eta(iB, iP, iR, iEta, iS_shuffle, iPhi) +=
                  array_yield_eta_bunches[iB][iP][iR][iEta][iS_true][iPhi][iBunch];
                // std::cout << "array_yield_eta_bunches[" << iB << "][" << iP << "][" << iR << "][" << iEta << "][" << iS_true << "][" << iPhi << "][" << iBunch << "] = "
                //           << array_yield_eta_bunches[iB][iP][iR][iEta][iS_true][iPhi][iBunch]
                //           << std::endl;
              }

              for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; ++iXf) {
                arr_yield_xf(iB, iP, iR, iXf, iS_shuffle, iPhi) +=
                  array_yield_xf_bunches[iB][iP][iR][iXf][iS_true][iPhi][iBunch];
                // std::cout << "array_yield_xf_bunches[" << iB << "][" << iP << "][" << iR << "][" << iXf << "][" << iS_true << "][" << iPhi << "][" << iBunch << "] = "
                //           << array_yield_xf_bunches[iB][iP][iR][iXf][iS_true][iPhi][iBunch]
                //           << std::endl;
              }
            }
          }
        }
      }
    }
    
    if (vec_yield_pt.size() <= iIter) {
      vec_yield_pt.resize(iIter+1);
      vec_yield_eta.resize(iIter+1);
      vec_yield_xf.resize(iIter+1);
    }
    vec_yield_pt[iIter] += arr_yield_pt;
    vec_yield_eta[iIter] += arr_yield_eta;
    vec_yield_xf[iIter] += arr_yield_xf;
  }

  void ShuffleBunches::compute_raw_asyms(const int iIter)
  {
    pt_yield_array &arr_pt_yield = vec_yield_pt[iIter];
    eta_yield_array &arr_eta_yield = vec_yield_eta[iIter];
    xf_yield_array &arr_xf_yield = vec_yield_xf[iIter];
    
    pt_asym_array arr_pt_asym;
    pt_nodir_asym_array arr_pt_nodir_asym;
    eta_asym_array arr_eta_asym;
    xf_asym_array arr_xf_asym;
    pt_asym_array arr_pt_asym_err;
    pt_nodir_asym_array arr_pt_nodir_asym_err;
    eta_asym_array arr_eta_asym_err;
    xf_asym_array arr_xf_asym_err;
    for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++) {
      for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++) {
        for (int iR = 0; iR < ASYM_CONSTANTS::nRegions; iR++) {
          for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
            {
              double array_asym[ASYM_CONSTANTS::nPhiBins];
              double array_asym_err[ASYM_CONSTANTS::nPhiBins];
              uint32_t array_yield[ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins];
              for (int iS = 0; iS < ASYM_CONSTANTS::nSpins; iS++) {
                for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
                  array_yield[iS][iPhi] =
                    arr_pt_yield(iB, iP, iR, iPt, 0, iS, iPhi) +
                    arr_pt_yield(iB, iP, iR, iPt, 1, iS, iPhi);
                }
              }
              compute_sqrt_asym(array_yield, polarization[iB], array_asym, array_asym_err);
              for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
                arr_pt_nodir_asym(iB, iP, iR, iPt, iPhi) = array_asym[iPhi];
                arr_pt_nodir_asym_err(iB, iP, iR, iPt, iPhi) = array_asym_err[iPhi];
              }
            }
            for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++) {
              double array_asym[ASYM_CONSTANTS::nPhiBins];
              double array_asym_err[ASYM_CONSTANTS::nPhiBins];
              uint32_t array_yield[ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins];
              for (int iS = 0; iS < ASYM_CONSTANTS::nSpins; iS++) {
                for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
                  array_yield[iS][iPhi] = arr_pt_yield(iB, iP, iR, iPt, iDir, iS, iPhi);
                  // std::cout << "arr_pt_yield(" << iB << ", " << iP << ", " << iR << ", " << iPt << ", " << iDir << ", " << iS << ", " << iPhi << ") = "
                  //           << arr_pt_yield(iB, iP, iR, iPt, iDir, iS, iPhi)
                  //           << std::endl;
                }
              }
              compute_sqrt_asym(array_yield, polarization[iB], array_asym, array_asym_err);
              for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
                arr_pt_asym(iB, iP, iR, iPt, iDir, iPhi) = array_asym[iPhi];
                arr_pt_asym_err(iB, iP, iR, iPt, iDir, iPhi) = array_asym_err[iPhi];
                // std::cout << "arr_pt_asym(" << iB << ", " << iP << ", " << iR << ", " << iPt << ", " << iDir << ", " << iPhi << ") = " 
                //           << arr_pt_asym(iB, iP, iR, iPt, iDir, iPhi) << std::endl;
              }
            }
          }
          for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
            double array_asym[ASYM_CONSTANTS::nPhiBins];
            double array_asym_err[ASYM_CONSTANTS::nPhiBins];
            uint32_t array_yield[ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins];
            for (int iS = 0; iS < ASYM_CONSTANTS::nSpins; iS++) {
              for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
                array_yield[iS][iPhi] = arr_eta_yield(iB, iP, iR, iEta, iS, iPhi);
              }
            }
            compute_sqrt_asym(array_yield, polarization[iB], array_asym, array_asym_err);
            for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
              arr_eta_asym(iB, iP, iR, iEta, iPhi) = array_asym[iPhi];
              arr_eta_asym_err(iB, iP, iR, iEta, iPhi) = array_asym_err[iPhi];
            }
          }
          for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
            double array_asym[ASYM_CONSTANTS::nPhiBins];
            double array_asym_err[ASYM_CONSTANTS::nPhiBins];
            uint32_t array_yield[ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins];
            for (int iS = 0; iS < ASYM_CONSTANTS::nSpins; iS++) {
              for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
                array_yield[iS][iPhi] = arr_xf_yield(iB, iP, iR, iXf, iS, iPhi);
              }
            }
            compute_sqrt_asym(array_yield, polarization[iB], array_asym, array_asym_err);
            for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
              arr_xf_asym(iB, iP, iR, iXf, iPhi) = array_asym[iPhi];
              arr_xf_asym_err(iB, iP, iR, iXf, iPhi) = array_asym_err[iPhi];
            }
          }
        }
      }
    }
    fill_pt_asyms.push_back(arr_pt_asym);
    fill_pt_nodir_asyms.push_back(arr_pt_nodir_asym);
    fill_eta_asyms.push_back(arr_eta_asym);
    fill_xf_asyms.push_back(arr_xf_asym);
    fill_pt_asym_errs.push_back(arr_pt_asym_err);
    fill_pt_nodir_asym_errs.push_back(arr_pt_nodir_asym_err);
    fill_eta_asym_errs.push_back(arr_eta_asym_err);
    fill_xf_asym_errs.push_back(arr_xf_asym_err);
  }

  void ShuffleBunches::average_asyms(const int iIter)
  {
    double sum_weights_pt[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nDirections][ASYM_CONSTANTS::nPhiBins] = {0};
    double sum_weights_pt_nodir[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nPhiBins] = {0};
    double sum_weights_eta[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nEtaBins][ASYM_CONSTANTS::nPhiBins] = {0};
    double sum_weights_xf[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nXfBins][ASYM_CONSTANTS::nPhiBins] = {0};

    double sum_weighted_asym_pt[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nDirections][ASYM_CONSTANTS::nPhiBins] = {0};
    double sum_weighted_asym_pt_nodir[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nPtBins][ASYM_CONSTANTS::nPhiBins] = {0};
    double sum_weighted_asym_eta[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nEtaBins][ASYM_CONSTANTS::nPhiBins] = {0};
    double sum_weighted_asym_xf[ASYM_CONSTANTS::nBeams][ASYM_CONSTANTS::nParticles][ASYM_CONSTANTS::nRegions][ASYM_CONSTANTS::nXfBins][ASYM_CONSTANTS::nPhiBins] = {0};
    
    for (size_t iFill = 0; iFill < nFills; iFill++) {
      for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++) {
        for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++) {
          for (int iR = 0; iR < ASYM_CONSTANTS::nRegions; iR++) {
            for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
              for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
                {
                  const double& asym = fill_pt_nodir_asyms[nIterations * iFill + iIter](iB, iP, iR, iPt, iPhi);
                  const double& asym_err = fill_pt_nodir_asym_errs[nIterations * iFill + iIter](iB, iP, iR, iPt, iPhi);
                  if (asym == asym && asym_err == asym_err && asym_err > 0)
                  {
                    sum_weights_pt_nodir[iB][iP][iR][iPt][iPhi] += std::pow(asym_err, -2);
                    sum_weighted_asym_pt_nodir[iB][iP][iR][iPt][iPhi] += asym * std::pow(asym_err, -2);
                  }
                }
                for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++) {
                  const double& asym = fill_pt_asyms[nIterations * iFill + iIter](iB, iP, iR, iPt, iDir, iPhi);
                  const double& asym_err = fill_pt_asym_errs[nIterations * iFill + iIter](iB, iP, iR, iPt, iDir, iPhi);
                  if (asym == asym && asym_err == asym_err && asym_err > 0)
                  {
                    sum_weights_pt[iB][iP][iR][iPt][iDir][iPhi] += std::pow(asym_err, -2);
                    sum_weighted_asym_pt[iB][iP][iR][iPt][iDir][iPhi] += asym * std::pow(asym_err, -2);
                  }
                }
              }
              for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
                const double& asym = fill_eta_asyms[nIterations * iFill + iIter](iB, iP, iR, iEta, iPhi);
                const double& asym_err = fill_eta_asym_errs[nIterations * iFill + iIter](iB, iP, iR, iEta, iPhi);
                if (asym == asym && asym_err == asym_err && asym_err > 0)
                {
                  sum_weights_eta[iB][iP][iR][iEta][iPhi] += std::pow(asym_err, -2);
                  sum_weighted_asym_eta[iB][iP][iR][iEta][iPhi] += asym * std::pow(asym_err, -2);
                }
              }
              for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
                const double& asym = fill_xf_asyms[nIterations * iFill + iIter](iB, iP, iR, iXf, iPhi);
                const double& asym_err = fill_xf_asym_errs[nIterations * iFill + iIter](iB, iP, iR, iXf, iPhi);
                if (asym == asym && asym_err == asym_err && asym_err > 0)
                {
                  
                  sum_weights_xf[iB][iP][iR][iXf][iPhi] += std::pow(asym_err, -2);
                  sum_weighted_asym_xf[iB][iP][iR][iXf][iPhi] += asym * std::pow(asym_err, -2);
                  // if (iB == 1 && iP == 1 && iR == 1 && iXf == 4 && iPhi == 0) {
                  //   std::cout << "iFill, asym, asym_err, sum_weights, sum_weighted_asyms =\n"
                  //             << iFill << ", "
                  //             << asym << ", "
                  //             << asym_err << ", "
                  //             << sum_weights_xf[iB][iP][iR][iXf][iPhi] << ", "
                  //             << sum_weighted_asym_xf[iB][iP][iR][iXf][iPhi] << std::endl;
                  // }
                }
              }
            }
          }
        }
      }
    }

    pt_asym_array average_mean_pt;
    pt_nodir_asym_array average_mean_pt_nodir;
    eta_asym_array average_mean_eta;
    xf_asym_array average_mean_xf;

    pt_asym_array average_unc_pt;
    pt_nodir_asym_array average_unc_pt_nodir;
    eta_asym_array average_unc_eta;
    xf_asym_array average_unc_xf;
    
    for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++) {
      for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++) {
        for (int iR = 0; iR < ASYM_CONSTANTS::nRegions; iR++) {
          for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
            for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
              {
                if (sum_weights_pt_nodir[iB][iP][iR][iPt][iPhi] > 0) {
                  average_mean_pt_nodir(iB, iP, iR, iPt, iPhi) =
                    sum_weighted_asym_pt_nodir[iB][iP][iR][iPt][iPhi] /
                    sum_weights_pt_nodir[iB][iP][iR][iPt][iPhi];
                  average_unc_pt_nodir(iB, iP, iR, iPt, iPhi) = 1. /
                    std::sqrt(sum_weights_pt_nodir[iB][iP][iR][iPt][iPhi]);
                }
                else {
                  average_mean_pt_nodir(iB, iP, iR, iPt, iPhi) = 0;
                  average_unc_pt_nodir(iB, iP, iR, iPt, iPhi) = 0;
                }
              }
              for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++) {
                if (sum_weights_pt[iB][iP][iR][iPt][iDir][iPhi] > 0) {
                  average_mean_pt(iB, iP, iR, iPt, iDir, iPhi) =
                    sum_weighted_asym_pt[iB][iP][iR][iPt][iDir][iPhi] /
                    sum_weights_pt[iB][iP][iR][iPt][iDir][iPhi];
                  average_unc_pt(iB, iP, iR, iPt, iDir, iPhi) = 1. /
                    std::sqrt(sum_weights_pt[iB][iP][iR][iPt][iDir][iPhi]);
                }
                else {
                  average_mean_pt(iB, iP, iR, iPt, iDir, iPhi) = 0;
                  average_unc_pt(iB, iP, iR, iPt, iDir, iPhi) = 0;
                }
              }
            }
            for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
              if (sum_weights_eta[iB][iP][iR][iEta][iPhi] > 0) {
                average_mean_eta(iB, iP, iR, iEta, iPhi) =
                  sum_weighted_asym_eta[iB][iP][iR][iEta][iPhi] /
                  sum_weights_eta[iB][iP][iR][iEta][iPhi];
                average_unc_eta(iB, iP, iR, iEta, iPhi) = 1. /
                  std::sqrt(sum_weights_eta[iB][iP][iR][iEta][iPhi]);
              } else {
                average_mean_eta(iB, iP, iR, iEta, iPhi) = 0;
                average_unc_eta(iB, iP, iR, iEta, iPhi) = 0;
              }
            }
            for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
              if (sum_weights_xf[iB][iP][iR][iXf][iPhi] > 0) {
                average_mean_xf(iB, iP, iR, iXf, iPhi) =
                  sum_weighted_asym_xf[iB][iP][iR][iXf][iPhi] /
                  sum_weights_xf[iB][iP][iR][iXf][iPhi];
                average_unc_xf(iB, iP, iR, iXf, iPhi) = 1. /
                  std::sqrt(sum_weights_xf[iB][iP][iR][iXf][iPhi]);
              } else {
                average_mean_xf(iB, iP, iR, iXf, iPhi) = 0;
                average_unc_xf(iB, iP, iR, iXf, iPhi) = 0;
              }
            }
          }
        }
      }
    }

    // // Fill-dependent raw asymmetries are no longer needed
    // fill_pt_asyms.clear();
    // fill_eta_asyms.clear();
    // fill_xf_asyms.clear();
    // fill_pt_asym_errs.clear();
    // fill_eta_asym_errs.clear();
    // fill_xf_asym_errs.clear();

    // Now fill the iteration-dependent vector
    mean_pt.push_back(average_mean_pt);
    mean_pt_nodir.push_back(average_mean_pt_nodir);
    mean_eta.push_back(average_mean_eta);
    mean_xf.push_back(average_mean_xf);
    unc_pt.push_back(average_unc_pt);
    unc_pt_nodir.push_back(average_unc_pt_nodir);
    unc_eta.push_back(average_unc_eta);
    unc_xf.push_back(average_unc_xf);
  }

  void ShuffleBunches::fit_raw_asyms(const int iIter)
  {
    const pt_asym_array& average_mean_pt = mean_pt[iIter];
    const pt_nodir_asym_array& average_mean_pt_nodir = mean_pt_nodir[iIter];
    const eta_asym_array& average_mean_eta = mean_eta[iIter];
    const xf_asym_array& average_mean_xf = mean_xf[iIter];
    const pt_asym_array& average_unc_pt = unc_pt[iIter];
    const pt_nodir_asym_array& average_unc_pt_nodir = unc_pt_nodir[iIter];
    const eta_asym_array& average_unc_eta = unc_eta[iIter];
    const xf_asym_array& average_unc_xf = unc_xf[iIter];

    pt_fit_array fit_average_mean_pt;
    pt_nodir_fit_array fit_average_mean_pt_nodir;
    eta_fit_array fit_average_mean_eta;
    xf_fit_array fit_average_mean_xf;
    pt_fit_array fit_average_unc_pt;
    pt_nodir_fit_array fit_average_unc_pt_nodir;
    eta_fit_array fit_average_unc_eta;
    xf_fit_array fit_average_unc_xf;

    for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++) {
      for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++) {
        for (int iR = 0; iR < ASYM_CONSTANTS::nRegions; iR++) {
          for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
            {
              double asym_azimuth[ASYM_CONSTANTS::nPhiBins] = {0};
              double asym_err_azimuth[ASYM_CONSTANTS::nPhiBins] = {0};
              for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
                asym_azimuth[iPhi] = average_mean_pt_nodir(iB, iP, iR, iPt, iPhi);
                asym_err_azimuth[iPhi] = average_unc_pt_nodir(iB, iP, iR, iPt, iPhi);
              }
              fit_asym(fit_average_mean_pt_nodir(iB, iP, iR, iPt),
                       fit_average_unc_pt_nodir(iB, iP, iR, iPt),
                       asym_azimuth,
                       asym_err_azimuth);
            }
            for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++) {
              double asym_azimuth[ASYM_CONSTANTS::nPhiBins] = {0};
              double asym_err_azimuth[ASYM_CONSTANTS::nPhiBins] = {0};
              for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
                asym_azimuth[iPhi] = average_mean_pt(iB, iP, iR, iPt, iDir, iPhi);
                asym_err_azimuth[iPhi] = average_unc_pt(iB, iP, iR, iPt, iDir, iPhi);
              }
              fit_asym(fit_average_mean_pt(iB, iP, iR, iPt, iDir),
                       fit_average_unc_pt(iB, iP, iR, iPt, iDir),
                       asym_azimuth,
                       asym_err_azimuth);
            }
          }
          for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
            double asym_azimuth[ASYM_CONSTANTS::nPhiBins] = {0};
            double asym_err_azimuth[ASYM_CONSTANTS::nPhiBins] = {0};
            for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
              asym_azimuth[iPhi] = average_mean_eta(iB, iP, iR, iEta, iPhi);
              asym_err_azimuth[iPhi] = average_unc_eta(iB, iP, iR, iEta, iPhi);
            }
            fit_asym(fit_average_mean_eta(iB, iP, iR, iEta),
                     fit_average_unc_eta(iB, iP, iR, iEta),
                     asym_azimuth,
                     asym_err_azimuth);
          }
          for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
            double asym_azimuth[ASYM_CONSTANTS::nPhiBins] = {0};
            double asym_err_azimuth[ASYM_CONSTANTS::nPhiBins] = {0};
            for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
              asym_azimuth[iPhi] = average_mean_xf(iB, iP, iR, iXf, iPhi);
              asym_err_azimuth[iPhi] = average_unc_xf(iB, iP, iR, iXf, iPhi);
            }
            fit_asym(fit_average_mean_xf(iB, iP, iR, iXf),
                     fit_average_unc_xf(iB, iP, iR, iXf),
                     asym_azimuth,
                     asym_err_azimuth);
          }
        }
      }
    }
    fit_mean_pt.push_back(fit_average_mean_pt);
    fit_mean_pt_nodir.push_back(fit_average_mean_pt_nodir);
    fit_mean_eta.push_back(fit_average_mean_eta);
    fit_mean_xf.push_back(fit_average_mean_xf);
    fit_unc_pt.push_back(fit_average_unc_pt);
    fit_unc_pt_nodir.push_back(fit_average_unc_pt_nodir);
    fit_unc_eta.push_back(fit_average_unc_eta);
    fit_unc_xf.push_back(fit_average_unc_xf);
  }

  void ShuffleBunches::fit_asym(double& asym_amplitude,
                                double& asym_amplitude_err,
                                const double (&asym_azimuth_values)[ASYM_CONSTANTS::nPhiBins],
                                const double (&asym_azimuth_uncertainties)[ASYM_CONSTANTS::nPhiBins],
                                bool fit_geom)
  {
    int N = (fit_geom ? ASYM_CONSTANTS::nPhiBins / 2 : ASYM_CONSTANTS::nPhiBins);

    // Output defaults
    asym_amplitude     = 0.0;
    asym_amplitude_err = 0.0;

    // 12 equal bins from -pi to +pi  => use bin centers
    double phi[N];
    double y[N];
    double ex[N];
    double ey[N];

    double phiMin = (fit_geom ? - M_PI / 2 : -M_PI);
    double phiMax = (fit_geom ? M_PI / 2 : M_PI);
    const double dphi = (phiMax - phiMin) / N;

    for (int i = 0; i < N; ++i) {
      int iphi = (fit_geom ? i + 3 : i);
      phi[i] = phiMin + (i + 0.5) * dphi;
      y[i]   = asym_azimuth_values[iphi];
      ex[i]  = 0.0;                      
      
      // Protect against zero/negative uncertainties (ROOT chi2 fit needs positive errors)
      ey[i] = (asym_azimuth_uncertainties[iphi] > 0.0)
        ? asym_azimuth_uncertainties[iphi]
        : 1.0;
    }

    TGraphErrors gr(N, phi, y, ex, ey);
    gr.SetName("gr_asym_phi");
    gr.SetTitle("Asymmetry vs #phi;#phi;Asymmetry");

    // Model: f(phi) = -[0] * sin(phi), where [0] is the amplitude A
    TF1 f_asym("f_asym", "-[0]*sin(x)", phiMin, phiMax);
    f_asym.SetParName(0, "A");

    // Weighted fit, quiet mode, return fit result
    // "Q" quiet, "S" return TFitResultPtr
    TFitResultPtr fitResult = gr.Fit(&f_asym, "QS");

    asym_amplitude     = f_asym.GetParameter(0);
    asym_amplitude_err = f_asym.GetParError(0);
  }

    void ShuffleBunches::compute_corrected_asyms(const int iIter)
  {
    const pt_fit_array& fit_average_mean_pt = fit_mean_pt[iIter];
    const pt_nodir_fit_array& fit_average_mean_pt_nodir = fit_mean_pt_nodir[iIter];
    const eta_fit_array& fit_average_mean_eta = fit_mean_eta[iIter];
    const xf_fit_array& fit_average_mean_xf = fit_mean_xf[iIter];
    const pt_fit_array& fit_average_unc_pt = fit_unc_pt[iIter];
    const pt_nodir_fit_array& fit_average_unc_pt_nodir = fit_unc_pt_nodir[iIter];
    const eta_fit_array& fit_average_unc_eta = fit_unc_eta[iIter];
    const xf_fit_array& fit_average_unc_xf = fit_unc_xf[iIter];

    pt_corr_array corrected_average_mean_pt;
    pt_nodir_corr_array corrected_average_mean_pt_nodir;
    eta_corr_array corrected_average_mean_eta;
    xf_corr_array corrected_average_mean_xf;
    pt_corr_array corrected_average_unc_pt;
    pt_nodir_corr_array corrected_average_unc_pt_nodir;
    eta_corr_array corrected_average_unc_eta;
    xf_corr_array corrected_average_unc_xf;

    for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++) {
      for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++) {
        {
          for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
            compute_corrected_asym(corrected_average_mean_pt_nodir(iB, iP, iPt),
                                   corrected_average_unc_pt_nodir(iB, iP, iPt),
                                   fit_average_mean_pt_nodir(iB, iP, 0, iPt),
                                   fit_average_mean_pt_nodir(iB, iP, 1, iPt),
                                   fit_average_unc_pt_nodir(iB, iP, 0, iPt),
                                   fit_average_unc_pt_nodir(iB, iP, 1, iPt),
                                   bkg_ratio_pt[iP][iPt],
                                   bkg_ratio_err_pt[iP][iPt]);
          }
        }
        for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++) {
          for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
            compute_corrected_asym(corrected_average_mean_pt(iB, iP, iPt, iDir),
                                   corrected_average_unc_pt(iB, iP, iPt, iDir),
                                   fit_average_mean_pt(iB, iP, 0, iPt, iDir),
                                   fit_average_mean_pt(iB, iP, 1, iPt, iDir),
                                   fit_average_unc_pt(iB, iP, 0, iPt, iDir),
                                   fit_average_unc_pt(iB, iP, 1, iPt, iDir),
                                   bkg_ratio_pt[iP][iPt],
                                   bkg_ratio_err_pt[iP][iPt]);
          }
        }
        for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++) {
          compute_corrected_asym(corrected_average_mean_eta(iB, iP, iEta),
                                 corrected_average_unc_eta(iB, iP, iEta),
                                 fit_average_mean_eta(iB, iP, 0, iEta),
                                 fit_average_mean_eta(iB, iP, 1, iEta),
                                 fit_average_unc_eta(iB, iP, 0, iEta),
                                 fit_average_unc_eta(iB, iP, 1, iEta),
                                 bkg_ratio_eta[iP][iEta],
                                 bkg_ratio_err_eta[iP][iEta]);
        }
        for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++) {
          compute_corrected_asym(corrected_average_mean_xf(iB, iP, iXf),
                                 corrected_average_unc_xf(iB, iP, iXf),
                                 fit_average_mean_xf(iB, iP, 0, iXf),
                                 fit_average_mean_xf(iB, iP, 1, iXf),
                                 fit_average_unc_xf(iB, iP, 0, iXf),
                                 fit_average_unc_xf(iB, iP, 1, iXf),
                                 bkg_ratio_xf[iP][iXf],
                                 bkg_ratio_err_xf[iP][iXf]);
        }
      }
    }
    corrected_mean_pt.push_back(corrected_average_mean_pt);
    corrected_mean_pt_nodir.push_back(corrected_average_mean_pt_nodir);
    corrected_mean_eta.push_back(corrected_average_mean_eta);
    corrected_mean_xf.push_back(corrected_average_mean_xf);
    corrected_unc_pt.push_back(corrected_average_unc_pt);
    corrected_unc_pt_nodir.push_back(corrected_average_unc_pt_nodir);
    corrected_unc_eta.push_back(corrected_average_unc_eta);
    corrected_unc_xf.push_back(corrected_average_unc_xf);
  }

  void ShuffleBunches::compute_corrected_asym(double &corrected_asym,
                                              double &corrected_asym_err,
                                              const double &raw_asym,
                                              const double &bkg_asym,
                                              const double &raw_asym_err,
                                              const double &bkg_asym_err,
                                              const double &R,
                                              const double &R_err)
  {
    const double Delta_R = 0;//(R > 1 ? 0 : (raw_asym - bkg_asym) / (1. - R) * R_err);
    const double Delta_raw = (R > 1 ? 0 : raw_asym_err / (1 - R));
    const double Delta_bkg = (R > 1 ? 0 : -R * bkg_asym_err / (1 - R));
    corrected_asym = (R > 1 ? 0 : (raw_asym - R * bkg_asym) / (1 - R));
    corrected_asym_err = (R > 1 ? 0 : std::sqrt(std::pow(Delta_R, 2) + std::pow(Delta_raw, 2) + std::pow(Delta_bkg, 2)));
  }

  void ShuffleBunches::get_iter_raw_asym_histograms()
  {
    pt_asym_array& arr_pt_asym = mean_pt.back();
    eta_asym_array& arr_eta_asym = mean_eta.back();
    xf_asym_array& arr_xf_asym = mean_xf.back();
    pt_asym_array& arr_pt_asym_err = unc_pt.back();
    eta_asym_array& arr_eta_asym_err = unc_eta.back();
    xf_asym_array& arr_xf_asym_err = unc_xf.back();
    
    outfile_iter_histograms->cd();
    for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++)
    {
      for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++)
      {
        for (int iR = 0; iR < ASYM_CONSTANTS::nRegions; iR++)
        {
          for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++)
          {
            for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++)
            {
              TGraphErrors *graph_AN = new TGraphErrors();
              std::stringstream graph_AN_name;
              graph_AN_name << "graph_asym_pt_"
                            << ASYM_CONSTANTS::beams[iB] << "_" << ASYM_CONSTANTS::particle[iP] << "_"
                            << ASYM_CONSTANTS::regions[iR] << "_"
                            << "pT_" << iPt << "_"
                            << ASYM_CONSTANTS::directions[iDir] << "_"
                            << "sqrt";
              graph_AN->SetName(graph_AN_name.str().c_str());
              std::stringstream graph_title;
              graph_title << "; #phi; " << (iR == 0 ? "Raw" : "Bkg") << " Asymmetry";
              graph_AN->SetTitle(graph_title.str().c_str());
              for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++)
              {
                float phi = -M_PI + M_PI / 12 + iPhi * M_PI / 6.;
                graph_AN->SetPoint(iPhi, phi, arr_pt_asym(iB, iP, iR, iPt, iDir, iPhi));
                graph_AN->SetPointError(iPhi, 0, arr_pt_asym_err(iB, iP, iR, iPt, iDir, iPhi));
              }
              graph_AN->Write();
            }
          }
          for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++)
          {
            TGraphErrors *graph_AN = new TGraphErrors();
            std::stringstream graph_AN_name;
            graph_AN_name << "graph_asym_eta_"
                          << ASYM_CONSTANTS::beams[iB] << "_" << ASYM_CONSTANTS::particle[iP] << "_"
                          << ASYM_CONSTANTS::regions[iR] << "_"
                          << "eta_" << iEta << "_sqrt";
            graph_AN->SetName(graph_AN_name.str().c_str());
            std::stringstream graph_title;
            graph_title << "; #phi; " << (iR == 0 ? "Raw" : "Bkg") << " Asymmetry";
            graph_AN->SetTitle(graph_title.str().c_str());
            for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++)
            {
              float phi = -M_PI + M_PI / 12 + iPhi * M_PI / 6.;
              graph_AN->SetPoint(iPhi, phi, arr_eta_asym(iB, iP, iR, iEta, iPhi));
              graph_AN->SetPointError(iPhi, 0, arr_eta_asym_err(iB, iP, iR, iEta, iPhi));
            }
            graph_AN->Write();
          }
          for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++)
          {
            TGraphErrors *graph_AN = new TGraphErrors();
            std::stringstream graph_AN_name;
            graph_AN_name << "graph_asym_xf_"
                          << ASYM_CONSTANTS::beams[iB] << "_" << ASYM_CONSTANTS::particle[iP] << "_"
                          << ASYM_CONSTANTS::regions[iR] << "_"
                          << "xf_" << iXf << "_sqrt";
            graph_AN->SetName(graph_AN_name.str().c_str());
            std::stringstream graph_title;
            graph_title << "; #phi; " << (iR == 0 ? "Raw" : "Bkg") << " Asymmetry";
            graph_AN->SetTitle(graph_title.str().c_str());

            // Iter
            for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++)
            {
              float phi = -M_PI + M_PI / 12 + iPhi * M_PI / 6.;
              graph_AN->SetPoint(iPhi, phi, arr_xf_asym(iB, iP, iR, iXf, iPhi));
              graph_AN->SetPointError(iPhi, 0, arr_xf_asym_err(iB, iP, iR, iXf, iPhi));
            }
            graph_AN->Write();
          }
        }
      }
    }
  }

  void ShuffleBunches::get_iter_fit_asym_histograms()
  {
    pt_fit_array& arr_pt_asym = fit_mean_pt.back();
    eta_fit_array& arr_eta_asym = fit_mean_eta.back();
    xf_fit_array& arr_xf_asym = fit_mean_xf.back();
    pt_fit_array& arr_pt_asym_err = fit_unc_pt.back();
    eta_fit_array& arr_eta_asym_err = fit_unc_eta.back();
    xf_fit_array& arr_xf_asym_err = fit_unc_xf.back();
    
    outfile_iter_histograms->cd();
    for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++)
    {
      for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++)
      {
        for (int iR = 0; iR < ASYM_CONSTANTS::nRegions; iR++)
        {
          for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++)
          {
            TGraphErrors *graph_AN = new TGraphErrors();
            std::stringstream graph_AN_name;
            graph_AN_name << "graph_asym_pt_"
                          << ASYM_CONSTANTS::beams[iB] << "_"
                          << ASYM_CONSTANTS::particle[iP] << "_"
                          << ASYM_CONSTANTS::regions[iR] << "_"
                          << ASYM_CONSTANTS::directions[iDir] << "_"
                          << "sqrt";
            graph_AN->SetName(graph_AN_name.str().c_str());
            std::stringstream graph_title;
            graph_title << "; p_{T} [GeV]; A_{N}^{" << (iR == 0 ? "Raw" : "Bkg") << "}";
            graph_AN->SetTitle(graph_title.str().c_str());
            for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++)
            {
              float pT = pTMeans[iP][iPt];
              graph_AN->SetPoint(iPt, pT, arr_pt_asym(iB, iP, iR, iPt, iDir));
              graph_AN->SetPointError(iPt, 0, arr_pt_asym_err(iB, iP, iR, iPt, iDir));
              graph_AN->Write();
            }
          }
          {
            TGraphErrors *graph_AN = new TGraphErrors();
            std::stringstream graph_AN_name;
            graph_AN_name << "graph_asym_eta_"
                          << ASYM_CONSTANTS::beams[iB] << "_"
                          << ASYM_CONSTANTS::particle[iP] << "_"
                          << ASYM_CONSTANTS::regions[iR] << "_"
                          << "sqrt";
            graph_AN->SetName(graph_AN_name.str().c_str());
            std::stringstream graph_title;
            graph_title << "; #eta; A_{N}^{" << (iR == 0 ? "Raw" : "Bkg") << "}";
            graph_AN->SetTitle(graph_title.str().c_str());
            for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++)
            {
              float eta = etaMeans[iP][iEta];
              graph_AN->SetPoint(iEta, eta, arr_eta_asym(iB, iP, iR, iEta));
              graph_AN->SetPointError(iEta, 0, arr_eta_asym_err(iB, iP, iR, iEta));

            }
            graph_AN->Write();
          }
          {
            TGraphErrors *graph_AN = new TGraphErrors();
            std::stringstream graph_AN_name;
            graph_AN_name << "graph_asym_xf_"
                          << ASYM_CONSTANTS::beams[iB] << "_"
                          << ASYM_CONSTANTS::particle[iP] << "_"
                          << ASYM_CONSTANTS::regions[iR] << "_"
                          << "sqrt";
            graph_AN->SetName(graph_AN_name.str().c_str());
            std::stringstream graph_title;
            graph_title << "; x_{F}; A_{N}^{" << (iR == 0 ? "Raw" : "Bkg") << "}";
            graph_AN->SetTitle(graph_title.str().c_str());
            for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++)
            {
              float xf = xfMeans[iP][iXf];
              graph_AN->SetPoint(iXf, xf, arr_xf_asym(iB, iP, iR, iXf));
              graph_AN->SetPointError(iXf, 0, arr_xf_asym_err(iB, iP, iR, iXf));

            }
            graph_AN->Write();
          }
        }
      }
    }
  }

  void ShuffleBunches::get_iter_corr_asym_histograms()
  {
    pt_corr_array& arr_pt_asym = corrected_mean_pt.back();
    eta_corr_array& arr_eta_asym = corrected_mean_eta.back();
    xf_corr_array& arr_xf_asym = corrected_mean_xf.back();
    pt_corr_array& arr_pt_asym_err = corrected_unc_pt.back();
    eta_corr_array& arr_eta_asym_err = corrected_unc_eta.back();
    xf_corr_array& arr_xf_asym_err = corrected_unc_xf.back();
    
    outfile_iter_histograms->cd();
    for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++)
    {
      for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++)
      {
        for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++)
        {
          TGraphErrors *graph_AN = new TGraphErrors();
          std::stringstream graph_AN_name;
          graph_AN_name << "graph_asym_pt_"
                        << ASYM_CONSTANTS::beams[iB] << "_" << ASYM_CONSTANTS::particle[iP] << "_"
                        << ASYM_CONSTANTS::directions[iDir] << "_"
                        << "sqrt";
          graph_AN->SetName(graph_AN_name.str().c_str());
          std::stringstream graph_title;
          graph_title << "; p_{T} [GeV]; A_{N}";
          graph_AN->SetTitle(graph_title.str().c_str());
          for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++)
          {
            float pT = pTMeans[iP][iPt];
            graph_AN->SetPoint(iPt, pT, arr_pt_asym(iB, iP, iPt, iDir));
            graph_AN->SetPointError(iPt, 0, arr_pt_asym_err(iB, iP, iPt, iDir));
            graph_AN->Write();
          }
        }
        {
          TGraphErrors *graph_AN = new TGraphErrors();
          std::stringstream graph_AN_name;
          graph_AN_name << "graph_asym_eta_"
                        << ASYM_CONSTANTS::beams[iB] << "_" << ASYM_CONSTANTS::particle[iP] << "_sqrt";
          graph_AN->SetName(graph_AN_name.str().c_str());
          std::stringstream graph_title;
          graph_title << "; #eta; A_{N}";
          graph_AN->SetTitle(graph_title.str().c_str());
          for (int iEta = 0; iEta < ASYM_CONSTANTS::nEtaBins; iEta++)
          {
            float eta = etaMeans[iP][iEta];
            graph_AN->SetPoint(iEta, eta, arr_eta_asym(iB, iP, iEta));
            graph_AN->SetPointError(iEta, 0, arr_eta_asym_err(iB, iP, iEta));

          }
          graph_AN->Write();
        }
        {
          TGraphErrors *graph_AN = new TGraphErrors();
          std::stringstream graph_AN_name;
          graph_AN_name << "graph_asym_xf_"
                        << ASYM_CONSTANTS::beams[iB] << "_" << ASYM_CONSTANTS::particle[iP] << "_"
                        << "sqrt";
          graph_AN->SetName(graph_AN_name.str().c_str());
          std::stringstream graph_title;
          graph_title << "; x_{F}; A_{N}";
          graph_AN->SetTitle(graph_title.str().c_str());
          for (int iXf = 0; iXf < ASYM_CONSTANTS::nXfBins; iXf++)
          {
            float xf = xfMeans[iP][iXf];
            graph_AN->SetPoint(iXf, xf, arr_xf_asym(iB, iP, iXf));
            graph_AN->SetPointError(iXf, 0, arr_xf_asym_err(iB, iP, iXf));

          }
          graph_AN->Write();
        }
      }
    }
  }
  
  void ShuffleBunches::check_array()
  {
    std::memset(array_yield_pt, 0, sizeof(array_yield_pt));
    std::memset(array_yield_eta, 0, sizeof(array_yield_eta));
    std::memset(array_yield_xf, 0, sizeof(array_yield_xf));
    
    for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++) {
      for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++) {
        for (int iR = 0; iR < ASYM_CONSTANTS::nRegions; iR++) {
          for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
            for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++) {
              for (int iS = 0; iS < ASYM_CONSTANTS::nSpins; iS++) {
                for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
                  for(int iBunch = 0; iBunch < ASYM_CONSTANTS::nBunches; iBunch++) {
                    // Increment pT yield
                    array_yield_pt[iB][iP][iR][iPt][iDir][iS][iPhi] +=
                      array_yield_pt_bunches[iB][iP][iR][iPt][iDir][iS][iPhi][iBunch];
                  }
                }
              }
            }
          }
        }
      }
    }
    
    for (int iB = 0; iB < ASYM_CONSTANTS::nBeams; iB++) {
      for (int iP = 0; iP < ASYM_CONSTANTS::nParticles; iP++) {
        for (int iR = 0; iR < ASYM_CONSTANTS::nRegions; iR++) {
          for (int iPt = 0; iPt < ASYM_CONSTANTS::nPtBins; iPt++) {
            for (int iDir = 0; iDir < ASYM_CONSTANTS::nDirections; iDir++) {
              for (int iS = 0; iS < ASYM_CONSTANTS::nSpins; iS++) {
                for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
                  std::cout << "array_yield_pt[" << iB << "][" << iP << "][" << iR << "][" << iPt << "][" << iDir << "][" << iS << "][" << iPhi << "] = "
                            << array_yield_pt[iB][iP][iR][iPt][iDir][iS][iPhi] << std::endl;
                }
              }
            }
          }
        }
      }
    }
  }

  void ShuffleBunches::compute_rellum_asym(const uint32_t (&array_yield)[ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins],
                                           const double pol,
                                           double (&array_asym)[ASYM_CONSTANTS::nPhiBins],
                                           double (&array_asym_err)[ASYM_CONSTANTS::nPhiBins])
  {
    auto rellum_asym = [&](double Nup, double Ndown, double RL) -> double {
      double numerator = Nup - RL * Ndown;
      double denominator = Nup + RL * Ndown;
      double  asym = numerator / denominator / pol * 100;
      return asym;
    };

    auto rellum_asym_err = [&](double Nup, double Ndown, double RL) -> double {
      // hypothesis: negligible uncertainty on the relative luminosity compared to the yields
      double NupErr = std::sqrt(Nup);
      double NdownErr = std::sqrt(Ndown);
      double t1 = 2 * RL * Ndown * NupErr;
      double t2 = 2 * RL * Nup * NdownErr;
      double denominator = pow(Nup + RL * Ndown,2);
      double asym_err = sqrt(pow(t1, 2) + pow(t2, 2)) / denominator / pol * 100;
      return asym_err;
    };

    for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
      double Nup = (double) array_yield[0][iPhi];
      double Ndown = (double) array_yield[1][iPhi];
      
      array_asym[iPhi] = rellum_asym(Nup, Ndown, relative_luminosity);
      array_asym_err[iPhi] = rellum_asym_err(Nup, Ndown, relative_luminosity);
    }
  }

  void ShuffleBunches::compute_sqrt_asym(const uint32_t (&array_yield)[ASYM_CONSTANTS::nSpins][ASYM_CONSTANTS::nPhiBins],
                                         const double pol,
                                         double (&array_asym)[ASYM_CONSTANTS::nPhiBins],
                                         double (&array_asym_err)[ASYM_CONSTANTS::nPhiBins])
  {
    auto sqrt_asym = [&](double NLup, double NLdown, double NRup, double NRdown) -> double {
      double numerator = sqrt(NLup*NRdown)-sqrt(NRup*NLdown);
      double denominator = sqrt(NLup*NRdown)+sqrt(NRup*NLdown);
      double asym = numerator / denominator / pol * 100; // Polarization given in %
      return asym;
    };

    auto sqrt_asym_err = [&](double NLup, double NLdown, double NRup, double NRdown) -> double {
      double NLupErr = std::sqrt(NLup);
      double NLdownErr = std::sqrt(NLdown);
      double NRupErr = std::sqrt(NRup);
      double NRdownErr = std::sqrt(NRdown);
      double t1 = sqrt(NLdown*NRup*NRdown/NLup)*NLupErr;
      double t2 = sqrt(NLup*NRup*NRdown/NLdown)*NLdownErr;
      double t3 = sqrt(NLup*NLdown*NRdown/NRup)*NRupErr;
      double t4 = sqrt(NLup*NLdown*NRup/NRdown)*NRdownErr;
      double denominator = pow(sqrt(NLup*NRdown)+sqrt(NRup*NLdown),2);
      double asym_err = sqrt(pow(t1,2)+pow(t2,2)+pow(t3,2)+pow(t4,2)) / denominator / pol * 100;
      return asym_err;
    };

    for (int iPhi = 0; iPhi < ASYM_CONSTANTS::nPhiBins; iPhi++) {
      int iPhiL = iPhi;
      int iPhiR = (iPhi + (int)(ASYM_CONSTANTS::nPhiBins/2)) % ASYM_CONSTANTS::nPhiBins;
      double NLup = (double) array_yield[0][iPhiL];
      double NLdown = (double) array_yield[1][iPhiL];
      double NRup = (double) array_yield[0][iPhiR];
      double NRdown = (double) array_yield[1][iPhiR];
      
      array_asym[iPhi] = sqrt_asym(NLup, NLdown, NRup, NRdown);
      array_asym_err[iPhi] = sqrt_asym_err(NLup, NLdown, NRup, NRdown);
    }
  }
};
