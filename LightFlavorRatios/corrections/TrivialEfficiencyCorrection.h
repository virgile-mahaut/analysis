#ifndef TRIVIAL_EFFICIENCY_CORRECTION_H
#define TRIVIAL_EFFICIENCY_CORRECTION_H

#include "CorrectionHistogram1D.h"

struct TrivialEfficiencyCorrection : CorrectionHistogram1D
{
  std::vector<float> hist_bins = {0.5,0.8,1.1,1.4,1.8,2.2,3.,4.};

  TrivialEfficiencyCorrection(std::string a)
  {
    h_corr = new TH1F(a.c_str(),"efficiency",1,-10.,10.);

    for(int i=1; i<=h_corr->GetNbinsX(); i++)
    {
      h_corr->SetBinContent(1,1.);
      h_corr->SetBinError(1,0.);
    }
    h_corr->SetDirectory(nullptr);
    name = "eff";
    title = "efficiency";
  }

  TrivialEfficiencyCorrection(const TrivialEfficiencyCorrection& e)
  {
    h_corr = e.h_corr;
    name = e.name;
    title = e.title;
  }

  void apply_correction(float xlow, float xhigh, TH1F* h, int bin) override
  {
  }
};

#endif
