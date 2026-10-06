#include "../corrections/EfficiencyCorrection.h"
#include "../corrections/LambdaFeedDownCorrection.h"
#include "../corrections/GeoAcceptanceCorrection.h"
#include "../corrections/TrivialEfficiencyCorrection.h"

#include "ResonanceRatio.h"
#include "LambdaModel.h"
#include "KshortModel.h"
#include "../bco_correction/V0DuplicateReader.h"

void Lambda_Kshort_ratio()
{
  TFile* Ks_file = TFile::Open("/sphenix/user/aopatton/Ana538CombinedOutput/outputLFPPG/KShort6RunCombined.root");
  TFile* lambda_file = TFile::Open("/sphenix/user/aopatton/Ana538CombinedOutput/outputLFPPG/Lambda6RunCombined.root");

  TTree* Ks_tree = (TTree*)Ks_file->Get("DecayTree");
  TTree* lambda_tree = (TTree*)lambda_file->Get("DecayTree");

  std::vector<HistogramInfo> variables =
  {
    BinInfo::final_pt_bins,
    BinInfo::final_eta_bins,
    BinInfo::final_rapidity_bins,
    BinInfo::final_phi_bins,
  };

  RooArgList Ks_args;
  RooArgList lambda_args;

  RooRealVar m_ks("K_S0_mass","K_S0_mass",0.45,0.55);
  RooRealVar m_lambda("Lambda0_mass","Lambda0_mass",1.1,1.14);

  Ks_args.add(m_ks);
  lambda_args.add(m_lambda);

  std::vector<RooRealVar> Ks_vars;
  std::vector<RooRealVar> lambda_vars;

  for (HistogramInfo& hinfo : variables)
  {
    std::string Ks_branchname = "K_S0_" + hinfo.name;
    std::string lambda_branchname = "Lambda0_" + hinfo.name;
    std::cout << Ks_branchname << " " << lambda_branchname << std::endl;

    RooRealVar Ks_var(Ks_branchname.c_str(), Ks_branchname.c_str(), hinfo.bins.front(), hinfo.bins.back());
    RooRealVar lambda_var(lambda_branchname.c_str(), lambda_branchname.c_str(), hinfo.bins.front(), hinfo.bins.back());

    Ks_vars.push_back(Ks_var);
    lambda_vars.push_back(lambda_var);
  }

  for (int i = 0; i < (int)Ks_vars.size(); i++)
  {
    Ks_args.add(Ks_vars[i]);
    lambda_args.add(lambda_vars[i]);
  }

  Ks_args.Print();
  lambda_args.Print();

  RooDataSet Ks_ds("K_S0","K_S0",Ks_args);
  RooDataSet lambda_ds("Lambda0","Lambda0",lambda_args);

  V0DuplicateReader ks_reader(Ks_tree, V0DuplicateReader::ParticleType::K0s);
  V0DuplicateReader lambda_reader(lambda_tree, V0DuplicateReader::ParticleType::Lambda);

  ks_reader.enableDeltaBCOCut(0, 350);
  lambda_reader.enableDeltaBCOCut(0, 350);

  for (Long64_t i = 0; i < ks_reader.entries(); ++i)
  {
    ks_reader.loadEntry(i);

    if (!ks_reader.passesDeltaBCOCut()) continue;
    if (!ks_reader.isCurrentEntryUnique()) continue;

    m_ks.setVal(ks_reader.mass());

    for (int ivar = 0; ivar < (int)variables.size(); ++ivar)
    {
      const std::string& name = variables[ivar].name;

      if (name == "pT") Ks_vars[ivar].setVal(ks_reader.pt());
      else if (name == "pseudorapidity") Ks_vars[ivar].setVal(ks_reader.eta());
      else if (name == "rapidity") Ks_vars[ivar].setVal(ks_reader.rapidity());
      else if (name == "phi") Ks_vars[ivar].setVal(ks_reader.phi());
    }

    Ks_ds.add(Ks_args);
  }

  for (Long64_t i = 0; i < lambda_reader.entries(); ++i)
  {
    lambda_reader.loadEntry(i);

    if (!lambda_reader.passesDeltaBCOCut()) continue;
    if (!lambda_reader.isCurrentEntryUnique()) continue;

    m_lambda.setVal(lambda_reader.mass());

    for (int ivar = 0; ivar < (int)variables.size(); ++ivar)
    {
      const std::string& name = variables[ivar].name;

      if (name == "pT") lambda_vars[ivar].setVal(lambda_reader.pt());
      else if (name == "pseudorapidity") lambda_vars[ivar].setVal(lambda_reader.eta());
      else if (name == "rapidity") lambda_vars[ivar].setVal(lambda_reader.rapidity());
      else if (name == "phi") lambda_vars[ivar].setVal(lambda_reader.phi());
    }

    lambda_ds.add(lambda_args);
  }

  std::cout << "K0s unique candidates      = " << ks_reader.numberOfUniqueCandidates() << std::endl;
  std::cout << "K0s duplicate candidates   = " << ks_reader.numberOfDuplicateCandidates() << std::endl;
  std::cout << "Lambda unique candidates   = " << lambda_reader.numberOfUniqueCandidates() << std::endl;
  std::cout << "Lambda duplicate candidates= " << lambda_reader.numberOfDuplicateCandidates() << std::endl;

  KshortModel kshort_model;
  LambdaModel lambda_model;

  std::string fd_filename = "/sphenix/tg/tg01/hf/hjheng/HF-analysis/simulation/Pythia_ppMinBias/cascade_feeddown/Cascade_feeddown_fraction.root";
  std::vector<std::vector<std::shared_ptr<CorrectionHistogram1D>>> corrections(variables.size());

  corrections[0].push_back(std::make_shared<LambdaFeedDownCorrection>(fd_filename,"h_feeddown_frac_xi_all"));
  corrections[0].push_back(std::make_shared<EfficiencyCorrection>());
  corrections[0].push_back(std::make_shared<GeoAcceptanceCorrection>("/sphenix/u/cdean/analysis/LightFlavorRatios/geometric_acceptance/analysis/plots/Lambda0_to_KS0_geometric_acceptance_ratio.root","pT"));

  corrections[1].push_back(std::make_shared<LambdaFeedDownCorrection>(fd_filename,"h_feeddown_frac_xi_eta_all"));
  corrections[1].push_back(std::make_shared<TrivialEfficiencyCorrection>("awef"));
  corrections[1].push_back(std::make_shared<GeoAcceptanceCorrection>("/sphenix/u/cdean/analysis/LightFlavorRatios/geometric_acceptance/analysis/plots/Lambda0_to_KS0_geometric_acceptance_ratio_eta.root","Lambda0_inGeo_#eta"));

  corrections[2].push_back(std::make_shared<LambdaFeedDownCorrection>(fd_filename,"h_feeddown_frac_xi_rapidity_all"));
  corrections[2].push_back(std::make_shared<TrivialEfficiencyCorrection>("aweg"));
  corrections[2].push_back(std::make_shared<GeoAcceptanceCorrection>("/sphenix/u/cdean/analysis/LightFlavorRatios/geometric_acceptance/analysis/plots/Lambda0_to_KS0_geometric_acceptance_ratio_rap.root","Lambda0_inGeo_y"));

  corrections[3].push_back(std::make_shared<LambdaFeedDownCorrection>(fd_filename,"h_feeddown_frac_xi_phi_all"));
  corrections[3].push_back(std::make_shared<TrivialEfficiencyCorrection>("aweh"));
  corrections[3].push_back(std::make_shared<GeoAcceptanceCorrection>("/sphenix/u/cdean/analysis/LightFlavorRatios/geometric_acceptance/analysis/plots/Lambda0_to_KS0_geometric_acceptance_ratio_phi.root","Lambda0_inGeo_#phi"));

  TFile* fout = new TFile("fits.root","RECREATE");

  ResonanceRatio analyzer(lambda_model, kshort_model,
                          fout,"lambdaKsratio","#Lambda/2K_{S}^{0} ratio",1./2.,true,
                          variables,corrections);

  analyzer.calculate_ratios_unbinned(lambda_ds,Ks_ds);
}
