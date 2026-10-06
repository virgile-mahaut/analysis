#ifndef EFFICIENCY_CORRECTION_H
#define EFFICIENCY_CORRECTION_H

#include "CorrectionHistogram1D.h"

struct EfficiencyCorrection : CorrectionHistogram1D
{
  std::vector<float> pt_center = {0.65,      0.95,     1.25,      1.6,       2.,        2.6,       3.5     };
  std::vector<float> pt_low    = {0.5,       0.8,      1.1,       1.4,       1.8,       2.2,       3.      };
  std::vector<float> pt_high   = {0.8,       1.1,      1.4,       1.8,       2.2,       3.,        4.      };
  std::vector<float> eff       = {0.975164,  1.02684,  1.09048,   1.11564,   1.05061,   1.02402,   1.03747 };
  std::vector<float> eff_err   = {0.0167643, 0.013714, 0.0173775, 0.0218095, 0.0325241, 0.0452745, 0.111922};

  std::vector<float> hist_bins = {0.5,0.8,1.1,1.4,1.8,2.2,3.,4.};

  EfficiencyCorrection()
  {
    h_corr = new TH1F("eff","efficiency",hist_bins.size()-1,hist_bins.data());

    for(int i=1; i<=h_corr->GetNbinsX(); i++)
    {
      h_corr->SetBinContent(i,eff[i-1]);
      h_corr->SetBinError(i,eff_err[i-1]);
    }
    h_corr->SetDirectory(nullptr);
    name = "eff";
    title = "efficiency";
  }

  EfficiencyCorrection(const EfficiencyCorrection& e)
  {
    h_corr = e.h_corr;
    name = e.name;
    title = e.title;
  }

  void apply_correction(float xlow, float xhigh, TH1F* h, int bin) override
  {
    std::pair<double,double> corr_and_err = get_val_and_error(xlow,xhigh);
    float new_val = h->GetBinContent(bin)/corr_and_err.first;
    float new_err = new_val * sqrt(pow(h->GetBinError(bin)/h->GetBinContent(bin),2.)+pow(corr_and_err.second/corr_and_err.first,2.));
    h->SetBinContent(bin,new_val);
    h->SetBinError(bin,new_err);
  }
};

/* from Tony:

pT 0.65 ptlow 0.5 pthigh 0.8 double ratio 0.975164 error 0.0167643
pT 0.95 ptlow 0.8 pthigh 1.1 double ratio 1.02684 error 0.013714
pT 1.25 ptlow 1.1 pthigh 1.4 double ratio 1.09048 error 0.0173775
pT 1.6 ptlow 1.4 pthigh 1.8 double ratio 1.11564 error 0.0218095
pT 2 ptlow 1.8 pthigh 2.2 double ratio 1.05061 error 0.0325241
pT 2.6 ptlow 2.2 pthigh 3 double ratio 1.02402 error 0.0452745
pT 3.5 ptlow 3 pthigh 4 double ratio 1.03747 error 0.111922

*/

#endif
