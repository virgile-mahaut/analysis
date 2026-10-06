#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TList.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "BootstrapGenerator/TH1DBootstrap.h"
#include "BootstrapGenerator/TH2DBootstrap.h"

R__LOAD_LIBRARY(libBootstrapGenerator.so)

void read_luminosity(const std::string& filename, double& lumi_all, double& lumi_zvertex60) {
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cerr << "Error: Cannot open " << filename << std::endl;
    lumi_all = 0;
    lumi_zvertex60 = 0;
    return;
  }
  double dummy1;
  infile >> dummy1 >> lumi_zvertex60 >> lumi_all;
  infile.close();
  std::cout << "Read from " << filename << ": lumi_all=" << lumi_all << ", lumi_zvertex60=" << lumi_zvertex60 << std::endl;
}

void combine_one(const std::string& pattern, int mode, int jet_radius_index, double frac_0mrad_all, double frac_1p5mrad_all, double frac_0mrad_zvertex60, double frac_1p5mrad_zvertex60) {
  
  std::vector<std::string> filenames;
  std::vector<double> weights_all;
  std::vector<double> weights_zvertex60;
  std::string outname;
  
  if (mode == 1) {
    outname = Form("output_sim_hadd/output_sim_r0%d_%s.root", jet_radius_index, pattern.c_str());
    filenames.push_back(Form("output_sim_hadd/output_sim_0mrad_single_r0%d_%s.root", jet_radius_index, pattern.c_str()));
    filenames.push_back(Form("output_sim_hadd/output_sim_1p5mrad_single_r0%d_%s.root", jet_radius_index, pattern.c_str()));
    weights_all = {frac_0mrad_all, frac_1p5mrad_all};
    weights_zvertex60 = {frac_0mrad_zvertex60, frac_1p5mrad_zvertex60};
  }
  else if (mode == 2) {
    outname = Form("output_sim_hadd/output_sim_r0%d_%s.root", jet_radius_index, pattern.c_str());
    filenames.push_back(Form("output_sim_hadd/output_sim_0mrad_single_r0%d_%s.root", jet_radius_index, pattern.c_str()));
    filenames.push_back(Form("output_sim_hadd/output_sim_0mrad_double_r0%d_%s.root", jet_radius_index, pattern.c_str()));
    weights_all = {0.78, 0.22};
    weights_zvertex60 = {0.78, 0.22};
  }
  else if (mode == 3) {
    outname = Form("output_sim_hadd/output_sim_r0%d_%s.root", jet_radius_index, pattern.c_str());
    filenames.push_back(Form("output_sim_hadd/output_sim_1p5mrad_single_r0%d_%s.root", jet_radius_index, pattern.c_str()));
    filenames.push_back(Form("output_sim_hadd/output_sim_1p5mrad_double_r0%d_%s.root", jet_radius_index, pattern.c_str()));
    weights_all = {0.921, 0.079};
    weights_zvertex60 = {0.921, 0.079};
  }
  else if (mode == 4) {
    outname = Form("output_sim_hadd/output_sim_r0%d_%s.root", jet_radius_index, pattern.c_str());
    filenames.push_back(Form("output_sim_hadd/output_sim_0mrad_single_r0%d_%s.root", jet_radius_index, pattern.c_str()));
    filenames.push_back(Form("output_sim_hadd/output_sim_0mrad_double_r0%d_%s.root", jet_radius_index, pattern.c_str()));
    filenames.push_back(Form("output_sim_hadd/output_sim_1p5mrad_single_r0%d_%s.root", jet_radius_index, pattern.c_str()));
    filenames.push_back(Form("output_sim_hadd/output_sim_1p5mrad_double_r0%d_%s.root", jet_radius_index, pattern.c_str()));
    weights_all = {0.78*frac_0mrad_all, 0.22*frac_0mrad_all, 0.921*frac_1p5mrad_all, 0.079*frac_1p5mrad_all};
    weights_zvertex60 = {0.78*frac_0mrad_zvertex60, 0.22*frac_0mrad_zvertex60, 0.921*frac_1p5mrad_zvertex60, 0.079*frac_1p5mrad_zvertex60};
  }
  
  // Open files
  std::vector<TFile*> files;
  for (const auto& fname : filenames) {
    TFile* f = new TFile(fname.c_str(), "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "Error: Cannot open " << fname << std::endl;
      for (auto ff : files) { ff->Close(); delete ff; }
      return;
    }
    files.push_back(f);
  }
  
  TFile* f_out = new TFile(outname.c_str(), "RECREATE");
  
  // Loop over histograms
  TList* keys = files[0]->GetListOfKeys();
  TIter next(keys);
  TKey* key;
  
  while ((key = (TKey*)next())) {
    std::string histname = key->GetName();
    std::string classname = key->GetClassName();
    
    // Choose weights based on histogram name
    bool use_zvertex60 = (histname.find("zvertex60") != std::string::npos);
    const std::vector<double>& weights = use_zvertex60 ? weights_zvertex60 : weights_all;
    
    if (classname == "TH1D") {
      TH1D* h_combined = nullptr;
      for (int i = 0; i < files.size(); i++) {
        TH1D* h = (TH1D*)files[i]->Get(histname.c_str());
        if (!h) continue;
        if (!h_combined) {
          f_out->cd();
          h_combined = new TH1D(*h);
          h_combined->SetName(histname.c_str());
          h_combined->Scale(weights[i]);
        } else {
          h_combined->Add(h, weights[i]);
        }
      }
      if (h_combined) h_combined->Write();
    }
    else if (classname == "TH2D") {
      TH2D* h_combined = nullptr;
      for (int i = 0; i < files.size(); i++) {
        TH2D* h = (TH2D*)files[i]->Get(histname.c_str());
        if (!h) continue;
        if (!h_combined) {
          f_out->cd();
          h_combined = new TH2D(*h);
          h_combined->SetName(histname.c_str());
          h_combined->Scale(weights[i]);
        } else {
          h_combined->Add(h, weights[i]);
        }
      }
      if (h_combined) h_combined->Write();
    }
    else if (classname == "TH1DBootstrap") {
      TH1DBootstrap* h_combined = nullptr;
      for (int i = 0; i < files.size(); i++) {
        TH1DBootstrap* h = (TH1DBootstrap*)files[i]->Get(histname.c_str());
        if (!h) continue;
        if (!h_combined) {
          f_out->cd();
          h_combined = new TH1DBootstrap(*h);
          h_combined->SetName(histname.c_str());
          h_combined->Scale(weights[i]);
        } else {
          TH1DBootstrap h_temp(*h);
          h_temp.Scale(weights[i]);
          h_combined->Add(&h_temp);
        }
      }
      if (h_combined) h_combined->Write();
    }
    else if (classname == "TH2DBootstrap") {
      TH2DBootstrap* h_combined = nullptr;
      for (int i = 0; i < files.size(); i++) {
        TH2DBootstrap* h = (TH2DBootstrap*)files[i]->Get(histname.c_str());
        if (!h) continue;
        if (!h_combined) {
          f_out->cd();
          h_combined = new TH2DBootstrap(*h);
          h_combined->SetName(histname.c_str());
          h_combined->Scale(weights[i]);
        } else {
          TH2DBootstrap h_temp(*h);
          h_temp.Scale(weights[i]);
          h_combined->Add(&h_temp);
        }
      }
      if (h_combined) h_combined->Write();
    }
  }
  
  f_out->Close();
  for (auto f : files) { f->Close(); delete f; }
  delete f_out;
  
  std::cout << "Output: " << outname << std::endl;
}

void get_mixoutput(int mode, int jet_radius_index = 4) {
  // mode 1: 0mrad single + 1p5mrad single
  // mode 2: 0mrad single + 0mrad double
  // mode 3: 1p5mrad single + 1p5mrad double
  // mode 4: 0mrad single + 0mrad double + 1p5mrad single + 1p5mrad double
  
  double lumi_0mrad_all, lumi_1p5mrad_all;
  double lumi_0mrad_zvertex60, lumi_1p5mrad_zvertex60;
  read_luminosity("luminosity_0mrad.txt", lumi_0mrad_all, lumi_0mrad_zvertex60);
  read_luminosity("luminosity_1p5mrad.txt", lumi_1p5mrad_all, lumi_1p5mrad_zvertex60);
  
  double frac_0mrad_all = lumi_0mrad_all / (lumi_0mrad_all + lumi_1p5mrad_all);
  double frac_1p5mrad_all = lumi_1p5mrad_all / (lumi_0mrad_all + lumi_1p5mrad_all);

  double frac_0mrad_zvertex60 = lumi_0mrad_zvertex60 / (lumi_0mrad_zvertex60 + lumi_1p5mrad_zvertex60);
  double frac_1p5mrad_zvertex60 = lumi_1p5mrad_zvertex60 / (lumi_0mrad_zvertex60 + lumi_1p5mrad_zvertex60);
  
  std::cout << "Lumi fractions: 0mrad=" << frac_0mrad_all << ", 1p5mrad=" << frac_1p5mrad_all << std::endl;
  std::cout << "Lumi fractions (zvertex60): 0mrad=" << frac_0mrad_zvertex60 << ", 1p5mrad=" << frac_1p5mrad_zvertex60 << std::endl;
  
  std::vector<std::string> patterns = {"Jet12GeV", "Jet20GeV", "Jet30GeV", "Jet40GeV", "Jet50GeV"};
  
  for (const auto& pattern : patterns) {
    std::cout << "\nCombining " << pattern << "..." << std::endl;
    combine_one(pattern, mode, jet_radius_index, frac_0mrad_all, frac_1p5mrad_all, frac_0mrad_zvertex60, frac_1p5mrad_zvertex60);
  }
  
  std::cout << "\nAll done!" << std::endl;
}
