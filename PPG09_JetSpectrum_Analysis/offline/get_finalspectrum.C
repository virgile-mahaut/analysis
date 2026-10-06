#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <vector>
#include "TLatex.h"
#include <string>
#include "unfold_Def.h"
#include "/sphenix/user/hanpuj/CaloDataAna24_skimmed/src/draw_template.C"               

void do_normalization(TH1D* h_unfolded, double luminosity, float jet_radius) {
  double delta_eta = 2 * (1.1 - jet_radius);
  for (int i = 1; i <= h_unfolded->GetNbinsX(); ++i) {
    double binwidth = h_unfolded->GetBinWidth(i);
    h_unfolded->SetBinContent(i, h_unfolded->GetBinContent(i) / (double)(luminosity * binwidth * delta_eta));
    h_unfolded->SetBinError(i, h_unfolded->GetBinError(i) / (double)(luminosity * binwidth * delta_eta));
  }
}

TH1D* get_statistical_uncertainty(TH1D* h_unfold, string name) {
  TH1D* h_uncertainty = (TH1D*)h_unfold->Clone(name.c_str());
  h_uncertainty->Reset();
  for (int i = 1; i <= h_unfold->GetNbinsX(); ++i) {
    double content = h_unfold->GetBinContent(i);
    if (content == 0) continue;
    double error = h_unfold->GetBinError(i);
    double uncertainty = error / content;
    h_uncertainty->SetBinContent(i, uncertainty);
    h_uncertainty->SetBinError(i, 0);
  }
  return h_uncertainty;
}

void get_uncertainty(TH1D* &h_uncertainty_up, TH1D* &h_uncertainty_down, TH1D* h_nominal, TH1D* h_shiftup, TH1D* h_shiftdown, string name) {
  h_uncertainty_up = (TH1D*)h_nominal->Clone((name + "_up").c_str());
  h_uncertainty_down = (TH1D*)h_nominal->Clone((name + "_down").c_str());
  h_uncertainty_up->Reset();
  h_uncertainty_down->Reset();
  for (int i = 1; i <= h_nominal->GetNbinsX(); ++i) {
    if (h_nominal->GetBinContent(i) == 0) continue;
    double diff_up = (h_shiftup->GetBinContent(i) - h_nominal->GetBinContent(i)) / (double)h_nominal->GetBinContent(i);
    double diff_down = (h_shiftdown->GetBinContent(i) - h_nominal->GetBinContent(i)) / (double)h_nominal->GetBinContent(i);
    double up = 0, down = 0;
    if (diff_up >= 0 && diff_down < 0) {
      up = diff_up;
      down = TMath::Abs(diff_down);
    } else if (diff_up < 0 && diff_down >= 0) {
      up = diff_down;
      down = TMath::Abs(diff_up);
    } else if (diff_up >= 0 && diff_down >= 0) {
      up = TMath::Max(diff_up, diff_down);
      down = 0; 
    } else if (diff_up < 0 && diff_down < 0) {
      up = 0;
      down = TMath::Max(TMath::Abs(diff_up), TMath::Abs(diff_down));
    }

    h_uncertainty_up->SetBinContent(i, up);
    h_uncertainty_up->SetBinError(i, 0);
    h_uncertainty_down->SetBinContent(i, down);
    h_uncertainty_down->SetBinError(i, 0);
  }
}

void get_unfold_uncertainty(TH1D* &h_uncertainty_up, TH1D* &h_uncertainty_down, TH1D* h_nominal, TH1D* h_shift, string name) {
  h_uncertainty_up = (TH1D*)h_nominal->Clone((name + "_up").c_str());
  h_uncertainty_down = (TH1D*)h_nominal->Clone((name + "_down").c_str());
  h_uncertainty_up->Reset();
  h_uncertainty_down->Reset();
  for (int i = 1; i <= h_nominal->GetNbinsX(); ++i) {
    double up = 0, down = 0;
    if (h_nominal->GetBinContent(i) == 0) continue;

    double diff = (h_shift->GetBinContent(i) - h_nominal->GetBinContent(i)) / (double)h_nominal->GetBinContent(i);
    up = TMath::Abs(diff);
    down = TMath::Abs(diff);

    h_uncertainty_up->SetBinContent(i, up);
    h_uncertainty_up->SetBinError(i, 0);
    h_uncertainty_down->SetBinContent(i, down);
    h_uncertainty_down->SetBinError(i, 0);
  }
}

void get_priority_uncertainty(TH1D* &h_uncertainty_up, TH1D* &h_uncertainty_down, TH1D* h_nominal, TH1D* h_alt, string name) {
  h_uncertainty_up = (TH1D*)h_nominal->Clone((name + "_up").c_str());
  h_uncertainty_down = (TH1D*)h_nominal->Clone((name + "_down").c_str());
  h_uncertainty_up->Reset();
  h_uncertainty_down->Reset();
  for (int i = 1; i <= h_nominal->GetNbinsX(); ++i) {
    if (h_nominal->GetBinContent(i) == 0) continue;
    double diff = TMath::Abs((h_alt->GetBinContent(i) - h_nominal->GetBinContent(i)) / (double)h_nominal->GetBinContent(i));
    //std::cout << "Bin " << i << ": nominal = " << h_nominal->GetBinContent(i) << ", alt = " << h_alt->GetBinContent(i) << ", diff = " << diff << std::endl;
    h_uncertainty_up->SetBinContent(i, diff);
    h_uncertainty_up->SetBinError(i, 0);
    h_uncertainty_down->SetBinContent(i, diff);
    h_uncertainty_down->SetBinError(i, 0);
  }
}

void get_total_uncertainty(TH1D* &h_total_up, TH1D* &h_total_down,TH1D* h_statistical, std::vector<TH1D*> h_systematic_up, std::vector<TH1D*> h_systematic_down) {
  h_total_up = (TH1D*)h_statistical->Clone("h_total_up");
  h_total_down = (TH1D*)h_statistical->Clone("h_total_down");
  h_total_up->Reset();
  h_total_down->Reset();
  for (int i = 1; i <= h_statistical->GetNbinsX(); ++i) {
    double stat_error = h_statistical->GetBinContent(i);

    double sys_error_up = 0;
    for (auto& h_up : h_systematic_up) {
      if (h_up->GetBinContent(i) > 0) {
        sys_error_up += h_up->GetBinContent(i) * h_up->GetBinContent(i);
      }
    }
    double total_error_up = TMath::Sqrt(stat_error * stat_error + sys_error_up);
    h_total_up->SetBinContent(i, total_error_up);
    h_total_up->SetBinError(i, 0);

    double sys_error_down = 0;
    for (auto& h_down : h_systematic_down) {
      if (h_down->GetBinContent(i) > 0) {
        sys_error_down += h_down->GetBinContent(i) * h_down->GetBinContent(i);
      }
    }
    double total_error_down = TMath::Sqrt(stat_error * stat_error + sys_error_down);
    h_total_down->SetBinContent(i, total_error_down);
    h_total_down->SetBinError(i, 0);
  }
}

void draw_total_uncertainty(TH1D* h_total_up, TH1D* h_total_down, TH1D* h_statistical, std::vector<TH1D*> h_systematic_up, std::vector<TH1D*> h_systematic_down, string label, std::vector<int> color, std::vector<string> legend, float jet_radius) {
  TFile* f_uncout = new TFile(Form("output_uncertainty_%s.root", label.c_str()), "RECREATE");

  TH1D* h_statistical_up = (TH1D*)h_statistical->Clone("h_statistical_up");
  TH1D* h_statistical_down = (TH1D*)h_statistical->Clone("h_statistical_down");
  h_statistical_down->Scale(-1.0);
  h_statistical_up->SetLineColor(kViolet-1);
  h_statistical_down->SetLineColor(kViolet-1);
  for (int i = 0; i < h_systematic_down.size(); ++i) {
    h_systematic_up[i]->SetLineColor(color[i]);
    h_systematic_down[i]->SetLineColor(color[i]);
    h_systematic_down[i]->Scale(-1.0);
  }
  h_total_up->SetMarkerStyle(20);
  h_total_down->SetMarkerStyle(20);
  h_total_up->SetLineColor(kBlack);
  h_total_down->SetLineColor(kBlack);
  h_total_up->SetMarkerColor(kBlack);
  h_total_down->SetMarkerColor(kBlack);
  h_total_down->Scale(-1.0);

  TCanvas *can_unc = new TCanvas("can_unc", "", 850, 800);
  gStyle->SetPalette(57);
  can_unc->SetTopMargin(0.03);
  can_unc->SetLeftMargin(0.15);
  can_unc->SetBottomMargin(0.12);
  can_unc->SetRightMargin(0.05);
  TH2F *frame_unc = new TH2F("frame_unc", "", 10, calibptbins[0], calibptbins[calibnpt], 2, -1, 1);
  frame_unc->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_unc->GetYaxis()->SetTitle("Relative Uncertainty");
  frame_unc->GetXaxis()->SetTitleOffset(0.98);
  frame_unc->GetYaxis()->SetTitleOffset(1.17);
  frame_unc->GetXaxis()->SetLabelSize(0.045);
  frame_unc->GetYaxis()->SetLabelSize(0.045);
  frame_unc->GetXaxis()->CenterTitle();
  frame_unc->GetYaxis()->CenterTitle();
  frame_unc->Draw();
  h_statistical_up->Draw("hist same");
  h_statistical_down->Draw("hist same");
  h_statistical_up->Write();
  h_statistical_down->Write();
  for (int i = 0; i < h_systematic_up.size(); ++i) {
    h_systematic_up[i]->Draw("hist same");
    h_systematic_down[i]->Draw("hist same");
    h_systematic_up[i]->Write();
    h_systematic_down[i]->Write();
  }
  h_total_up->Draw("p same");
  h_total_down->Draw("p same");
  h_total_up->Write();
  h_total_down->Write();
  myText(0.17, 0.92, 1, "#bf{#it{sPHENIX}} Internal", 0.05);
  myText(0.17, 0.86, 1, "p+p #sqrt{s} = 200 GeV", 0.05);
  myText(0.17, 0.8, 1, Form("anti-k_{t} #kern[-0.3]{#it{R}} = %.1f", jet_radius), 0.05);
  myText(0.17, 0.74, 1, Form("|#eta^{jet}| < %.1f", (1.1-jet_radius)), 0.05);
  TLegend *leg_unc1 = new TLegend(0.15, 0.13, 0.75, 0.25);
  leg_unc1->SetNColumns(2);
  leg_unc1->SetBorderSize(0);
  leg_unc1->SetFillStyle(0);
  leg_unc1->SetTextFont(42);
  leg_unc1->AddEntry(h_statistical_up, "Statistical Uncertainty", "l");
  for (int i = 0; i < h_systematic_up.size(); ++i) {
    leg_unc1->AddEntry(h_systematic_up[i], legend[i].c_str(), "l");
  }
  leg_unc1->AddEntry(h_total_up, "Total Uncertainty", "p");
  leg_unc1->Draw();
  TLine *line = new TLine(calibptbins[0], 0, calibptbins[calibnpt], 0);
  line->SetLineStyle(9);
  line->Draw("same");
  can_unc->SaveAs(("figure/spectrum_uncertainty" + label + ".png").c_str());

  delete can_unc;
  delete frame_unc;
  f_uncout->Close();
}

float low_edge_yaxis(float low_point) {
  if (low_point <= 0) return 0.0001;
  else return std::pow(10.0, std::floor(std::log10(low_point)))/10.;
}

float high_edge_yaxis(float high_point) {
  if (high_point <= 0) return 0.0001;
  else return std::pow(10.0, std::ceil(std::log10(high_point)))*10.;
}

double IntegrateSpline(TSpline3* spline, double xlow, double xhigh, int nsteps = 100) {
    double dx = (xhigh - xlow) / nsteps;
    double integral = 0.0;
    
    for (int i = 0; i < nsteps; i++) {
        double x1 = xlow + i * dx;
        double x2 = xlow + (i + 1) * dx;
        integral += 0.5 * (spline->Eval(x1) + spline->Eval(x2)) * dx;
    }
    return integral;
}

void get_finalspectrum() {

  TFile *f_out = new TFile("output_final_spectrum.root", "RECREATE");
  TFile *f_out_forratio = new TFile("output_forratio.root", "RECREATE");

  ifstream f_luminosity("luminosity.txt");
  double luminosity22_30, luminosity22_60, luminosity22_all;
  f_luminosity >> luminosity22_30 >> luminosity22_60 >> luminosity22_all;
  std::cout << "Luminosity: " << luminosity22_30 << "  " << luminosity22_60 << "  " << luminosity22_all << std::endl;

  TFile* f_r02_spectrum = new TFile("output_unfolded_r02.root", "READ");
  TFile* f_r03_spectrum = new TFile("output_unfolded_r03.root", "READ");
  TFile* f_r04_spectrum = new TFile("output_unfolded_r04.root", "READ");
  TFile* f_r05_spectrum = new TFile("output_unfolded_r05.root", "READ");
  TFile* f_r06_spectrum = new TFile("output_unfolded_r06.root", "READ");
  TFile* f_r08_spectrum = new TFile("output_unfolded_r08.root", "READ");
  if (!f_r02_spectrum || f_r02_spectrum->IsZombie()) {
    std::cerr << "Error opening file: output_unfolded_r02.root" << std::endl;
    return;
  }
  if (!f_r03_spectrum || f_r03_spectrum->IsZombie()) {
    std::cerr << "Error opening file: output_unfolded_r03.root" << std::endl;
    return;
  }
  if (!f_r04_spectrum || f_r04_spectrum->IsZombie()) {
    std::cerr << "Error opening file: output_unfolded_r04.root" << std::endl;
    return;
  }
  if (!f_r05_spectrum || f_r05_spectrum->IsZombie()) {
    std::cerr << "Error opening file: output_unfolded_r05.root" << std::endl;
    return;
  }
  if (!f_r06_spectrum || f_r06_spectrum->IsZombie()) {
    std::cerr << "Error opening file: output_unfolded_r06.root" << std::endl;
    return;
  }
  if (!f_r08_spectrum || f_r08_spectrum->IsZombie()) {
    std::cerr << "Error opening file: output_unfolded_r08.root" << std::endl;
    return;
  }

  TH1D *h_r02_unfold_all, *h_r02_unfold_all_stat, *h_r02_unfold_all_unreweighted, *h_r02_unfold_all_jesup, *h_r02_unfold_all_jesdown, *h_r02_unfold_all_jerup, *h_r02_unfold_all_jerdown, *h_r02_unfold_all_jetup, *h_r02_unfold_all_jetdown, *h_r02_unfold_all_unfoldunc;
  h_r02_unfold_all = (TH1D*)f_r02_spectrum->Get("h_unfold_all"); h_r02_unfold_all_stat = (TH1D*)((TH1DBootstrap*)f_r02_spectrum->Get("h_unfold_all_stat"))->GetNominal(); h_r02_unfold_all_unreweighted = (TH1D*)f_r02_spectrum->Get("h_unfold_all_unreweighted"); h_r02_unfold_all_jesup = (TH1D*)f_r02_spectrum->Get("h_unfold_all_jesup"); h_r02_unfold_all_jesdown = (TH1D*)f_r02_spectrum->Get("h_unfold_all_jesdown"); h_r02_unfold_all_jerup = (TH1D*)f_r02_spectrum->Get("h_unfold_all_jerup"); h_r02_unfold_all_jerdown = (TH1D*)f_r02_spectrum->Get("h_unfold_all_jerdown"); h_r02_unfold_all_jetup = (TH1D*)f_r02_spectrum->Get("h_unfold_all_jetup"); h_r02_unfold_all_jetdown = (TH1D*)f_r02_spectrum->Get("h_unfold_all_jetdown"); h_r02_unfold_all_unfoldunc = (TH1D*)f_r02_spectrum->Get("h_unfold_all_unfoldunc");
  do_normalization(h_r02_unfold_all, luminosity22_all, 0.2); do_normalization(h_r02_unfold_all_stat, luminosity22_all, 0.2); do_normalization(h_r02_unfold_all_unreweighted, luminosity22_all, 0.2); do_normalization(h_r02_unfold_all_jesup, luminosity22_all, 0.2); do_normalization(h_r02_unfold_all_jesdown, luminosity22_all, 0.2); do_normalization(h_r02_unfold_all_jerup, luminosity22_all, 0.2); do_normalization(h_r02_unfold_all_jerdown, luminosity22_all, 0.2); do_normalization(h_r02_unfold_all_jetup, luminosity22_all, 0.2); do_normalization(h_r02_unfold_all_jetdown, luminosity22_all, 0.2); do_normalization(h_r02_unfold_all_unfoldunc, luminosity22_all, 0.2);
  TH1D* h_r02_uncertainty_all_stat = get_statistical_uncertainty(h_r02_unfold_all_stat, "h_r02_uncertainty_all_stat");
  TH1D* h_r02_uncertainty_all_jesup, *h_r02_uncertainty_all_jesdown; get_uncertainty(h_r02_uncertainty_all_jesup, h_r02_uncertainty_all_jesdown, h_r02_unfold_all, h_r02_unfold_all_jesup, h_r02_unfold_all_jesdown, "h_r02_uncertainty_all_jes");
  TH1D* h_r02_uncertainty_all_jerup, *h_r02_uncertainty_all_jerdown; get_uncertainty(h_r02_uncertainty_all_jerup, h_r02_uncertainty_all_jerdown, h_r02_unfold_all, h_r02_unfold_all_jerup, h_r02_unfold_all_jerdown, "h_r02_uncertainty_all_jer");
  TH1D* h_r02_uncertainty_all_jetup, *h_r02_uncertainty_all_jetdown; get_uncertainty(h_r02_uncertainty_all_jetup, h_r02_uncertainty_all_jetdown, h_r02_unfold_all, h_r02_unfold_all_jetup, h_r02_unfold_all_jetdown, "h_r02_uncertainty_all_jet");
  TH1D* h_r02_uncertainty_all_unfoldup, *h_r02_uncertainty_all_unfolddown; get_unfold_uncertainty(h_r02_uncertainty_all_unfoldup, h_r02_uncertainty_all_unfolddown, h_r02_unfold_all, h_r02_unfold_all_unfoldunc, "h_r02_uncertainty_all_unfold");
  TH1D* h_r02_uncertainty_all_priorityup, *h_r02_uncertainty_all_prioritydown; get_priority_uncertainty(h_r02_uncertainty_all_priorityup, h_r02_uncertainty_all_prioritydown, h_r02_unfold_all, h_r02_unfold_all_unreweighted, "h_r02_uncertainty_all_priority");
  std::vector<TH1D*> h_r02_uncertainty_all_up = {h_r02_uncertainty_all_jesup, h_r02_uncertainty_all_jerup, h_r02_uncertainty_all_jetup, h_r02_uncertainty_all_unfoldup, h_r02_uncertainty_all_priorityup};
  std::vector<TH1D*> h_r02_uncertainty_all_down = {h_r02_uncertainty_all_jesdown, h_r02_uncertainty_all_jerdown, h_r02_uncertainty_all_jetdown, h_r02_uncertainty_all_unfolddown, h_r02_uncertainty_all_prioritydown};
  TH1D* h_r02_uncertainty_all_total_up, *h_r02_uncertainty_all_total_down; get_total_uncertainty(h_r02_uncertainty_all_total_up, h_r02_uncertainty_all_total_down, h_r02_uncertainty_all_stat, h_r02_uncertainty_all_up, h_r02_uncertainty_all_down);

  TH1D *h_r03_unfold_all, *h_r03_unfold_all_stat, *h_r03_unfold_all_unreweighted, *h_r03_unfold_all_jesup, *h_r03_unfold_all_jesdown, *h_r03_unfold_all_jerup, *h_r03_unfold_all_jerdown, *h_r03_unfold_all_jetup, *h_r03_unfold_all_jetdown, *h_r03_unfold_all_unfoldunc;
  h_r03_unfold_all = (TH1D*)f_r03_spectrum->Get("h_unfold_all"); h_r03_unfold_all_stat = (TH1D*)((TH1DBootstrap*)f_r03_spectrum->Get("h_unfold_all_stat"))->GetNominal(); h_r03_unfold_all_unreweighted = (TH1D*)f_r03_spectrum->Get("h_unfold_all_unreweighted"); h_r03_unfold_all_jesup = (TH1D*)f_r03_spectrum->Get("h_unfold_all_jesup"); h_r03_unfold_all_jesdown = (TH1D*)f_r03_spectrum->Get("h_unfold_all_jesdown"); h_r03_unfold_all_jerup = (TH1D*)f_r03_spectrum->Get("h_unfold_all_jerup"); h_r03_unfold_all_jerdown = (TH1D*)f_r03_spectrum->Get("h_unfold_all_jerdown"); h_r03_unfold_all_jetup = (TH1D*)f_r03_spectrum->Get("h_unfold_all_jetup"); h_r03_unfold_all_jetdown = (TH1D*)f_r03_spectrum->Get("h_unfold_all_jetdown"); h_r03_unfold_all_unfoldunc = (TH1D*)f_r03_spectrum->Get("h_unfold_all_unfoldunc");
  do_normalization(h_r03_unfold_all, luminosity22_all, 0.3); do_normalization(h_r03_unfold_all_stat, luminosity22_all, 0.3); do_normalization(h_r03_unfold_all_unreweighted, luminosity22_all, 0.3); do_normalization(h_r03_unfold_all_jesup, luminosity22_all, 0.3); do_normalization(h_r03_unfold_all_jesdown, luminosity22_all, 0.3); do_normalization(h_r03_unfold_all_jerup, luminosity22_all, 0.3); do_normalization(h_r03_unfold_all_jerdown, luminosity22_all, 0.3); do_normalization(h_r03_unfold_all_jetup, luminosity22_all, 0.3); do_normalization(h_r03_unfold_all_jetdown, luminosity22_all, 0.3); do_normalization(h_r03_unfold_all_unfoldunc, luminosity22_all, 0.3);
  TH1D* h_r03_uncertainty_all_stat = get_statistical_uncertainty(h_r03_unfold_all_stat, "h_r03_uncertainty_all_stat");
  TH1D* h_r03_uncertainty_all_jesup, *h_r03_uncertainty_all_jesdown; get_uncertainty(h_r03_uncertainty_all_jesup, h_r03_uncertainty_all_jesdown, h_r03_unfold_all, h_r03_unfold_all_jesup, h_r03_unfold_all_jesdown, "h_r03_uncertainty_all_jes");
  TH1D* h_r03_uncertainty_all_jerup, *h_r03_uncertainty_all_jerdown; get_uncertainty(h_r03_uncertainty_all_jerup, h_r03_uncertainty_all_jerdown, h_r03_unfold_all, h_r03_unfold_all_jerup, h_r03_unfold_all_jerdown, "h_r03_uncertainty_all_jer");
  TH1D* h_r03_uncertainty_all_jetup, *h_r03_uncertainty_all_jetdown; get_uncertainty(h_r03_uncertainty_all_jetup, h_r03_uncertainty_all_jetdown, h_r03_unfold_all, h_r03_unfold_all_jetup, h_r03_unfold_all_jetdown, "h_r03_uncertainty_all_jet");
  TH1D* h_r03_uncertainty_all_unfoldup, *h_r03_uncertainty_all_unfolddown; get_unfold_uncertainty(h_r03_uncertainty_all_unfoldup, h_r03_uncertainty_all_unfolddown, h_r03_unfold_all, h_r03_unfold_all_unfoldunc, "h_r03_uncertainty_all_unfold");
  TH1D* h_r03_uncertainty_all_priorityup, *h_r03_uncertainty_all_prioritydown; get_priority_uncertainty(h_r03_uncertainty_all_priorityup, h_r03_uncertainty_all_prioritydown, h_r03_unfold_all, h_r03_unfold_all_unreweighted, "h_r03_uncertainty_all_priority");
  std::vector<TH1D*> h_r03_uncertainty_all_up = {h_r03_uncertainty_all_jesup, h_r03_uncertainty_all_jerup, h_r03_uncertainty_all_jetup, h_r03_uncertainty_all_unfoldup, h_r03_uncertainty_all_priorityup};
  std::vector<TH1D*> h_r03_uncertainty_all_down = {h_r03_uncertainty_all_jesdown, h_r03_uncertainty_all_jerdown, h_r03_uncertainty_all_jetdown, h_r03_uncertainty_all_unfolddown, h_r03_uncertainty_all_prioritydown};
  TH1D* h_r03_uncertainty_all_total_up, *h_r03_uncertainty_all_total_down; get_total_uncertainty(h_r03_uncertainty_all_total_up, h_r03_uncertainty_all_total_down, h_r03_uncertainty_all_stat, h_r03_uncertainty_all_up, h_r03_uncertainty_all_down);

  TH1D *h_r04_unfold_all, *h_r04_unfold_all_stat, *h_r04_unfold_all_unreweighted, *h_r04_unfold_all_jesup, *h_r04_unfold_all_jesdown, *h_r04_unfold_all_jerup, *h_r04_unfold_all_jerdown, *h_r04_unfold_all_jetup, *h_r04_unfold_all_jetdown, *h_r04_unfold_all_unfoldunc;
  h_r04_unfold_all = (TH1D*)f_r04_spectrum->Get("h_unfold_all"); h_r04_unfold_all_stat = (TH1D*)((TH1DBootstrap*)f_r04_spectrum->Get("h_unfold_all_stat"))->GetNominal(); h_r04_unfold_all_unreweighted = (TH1D*)f_r04_spectrum->Get("h_unfold_all_unreweighted"); h_r04_unfold_all_jesup = (TH1D*)f_r04_spectrum->Get("h_unfold_all_jesup"); h_r04_unfold_all_jesdown = (TH1D*)f_r04_spectrum->Get("h_unfold_all_jesdown"); h_r04_unfold_all_jerup = (TH1D*)f_r04_spectrum->Get("h_unfold_all_jerup"); h_r04_unfold_all_jerdown = (TH1D*)f_r04_spectrum->Get("h_unfold_all_jerdown"); h_r04_unfold_all_jetup = (TH1D*)f_r04_spectrum->Get("h_unfold_all_jetup"); h_r04_unfold_all_jetdown = (TH1D*)f_r04_spectrum->Get("h_unfold_all_jetdown"); h_r04_unfold_all_unfoldunc = (TH1D*)f_r04_spectrum->Get("h_unfold_all_unfoldunc");
  do_normalization(h_r04_unfold_all, luminosity22_all, 0.4); do_normalization(h_r04_unfold_all_stat, luminosity22_all, 0.4); do_normalization(h_r04_unfold_all_unreweighted, luminosity22_all, 0.4); do_normalization(h_r04_unfold_all_jesup, luminosity22_all, 0.4); do_normalization(h_r04_unfold_all_jesdown, luminosity22_all, 0.4); do_normalization(h_r04_unfold_all_jerup, luminosity22_all, 0.4); do_normalization(h_r04_unfold_all_jerdown, luminosity22_all, 0.4); do_normalization(h_r04_unfold_all_jetup, luminosity22_all, 0.4); do_normalization(h_r04_unfold_all_jetdown, luminosity22_all, 0.4); do_normalization(h_r04_unfold_all_unfoldunc, luminosity22_all, 0.4);
  TH1D* h_r04_uncertainty_all_stat = get_statistical_uncertainty(h_r04_unfold_all_stat, "h_r04_uncertainty_all_stat");
  TH1D* h_r04_uncertainty_all_jesup, *h_r04_uncertainty_all_jesdown; get_uncertainty(h_r04_uncertainty_all_jesup, h_r04_uncertainty_all_jesdown, h_r04_unfold_all, h_r04_unfold_all_jesup, h_r04_unfold_all_jesdown, "h_r04_uncertainty_all_jes");
  TH1D* h_r04_uncertainty_all_jerup, *h_r04_uncertainty_all_jerdown; get_uncertainty(h_r04_uncertainty_all_jerup, h_r04_uncertainty_all_jerdown, h_r04_unfold_all, h_r04_unfold_all_jerup, h_r04_unfold_all_jerdown, "h_r04_uncertainty_all_jer");
  TH1D* h_r04_uncertainty_all_jetup, *h_r04_uncertainty_all_jetdown; get_uncertainty(h_r04_uncertainty_all_jetup, h_r04_uncertainty_all_jetdown, h_r04_unfold_all, h_r04_unfold_all_jetup, h_r04_unfold_all_jetdown, "h_r04_uncertainty_all_jet");
  TH1D* h_r04_uncertainty_all_unfoldup, *h_r04_uncertainty_all_unfolddown; get_unfold_uncertainty(h_r04_uncertainty_all_unfoldup, h_r04_uncertainty_all_unfolddown, h_r04_unfold_all, h_r04_unfold_all_unfoldunc, "h_r04_uncertainty_all_unfold");
  TH1D* h_r04_uncertainty_all_priorityup, *h_r04_uncertainty_all_prioritydown; get_priority_uncertainty(h_r04_uncertainty_all_priorityup, h_r04_uncertainty_all_prioritydown, h_r04_unfold_all, h_r04_unfold_all_unreweighted, "h_r04_uncertainty_all_priority");
  std::vector<TH1D*> h_r04_uncertainty_all_up = {h_r04_uncertainty_all_jesup, h_r04_uncertainty_all_jerup, h_r04_uncertainty_all_jetup, h_r04_uncertainty_all_unfoldup, h_r04_uncertainty_all_priorityup};
  std::vector<TH1D*> h_r04_uncertainty_all_down = {h_r04_uncertainty_all_jesdown, h_r04_uncertainty_all_jerdown, h_r04_uncertainty_all_jetdown, h_r04_uncertainty_all_unfolddown, h_r04_uncertainty_all_prioritydown};
  TH1D* h_r04_uncertainty_all_total_up, *h_r04_uncertainty_all_total_down; get_total_uncertainty(h_r04_uncertainty_all_total_up, h_r04_uncertainty_all_total_down, h_r04_uncertainty_all_stat, h_r04_uncertainty_all_up, h_r04_uncertainty_all_down);
  std::vector<TH1D*> h_r04_uncertainty_all_up_cancelled = {h_r04_uncertainty_all_unfoldup, h_r04_uncertainty_all_priorityup};
  std::vector<TH1D*> h_r04_uncertainty_all_down_cancelled = {h_r04_uncertainty_all_unfolddown, h_r04_uncertainty_all_prioritydown};
  TH1D* h_r04_uncertainty_all_total_up_cancelled, *h_r04_uncertainty_all_total_down_cancelled; get_total_uncertainty(h_r04_uncertainty_all_total_up_cancelled, h_r04_uncertainty_all_total_down_cancelled, h_r04_uncertainty_all_stat, h_r04_uncertainty_all_up_cancelled, h_r04_uncertainty_all_down_cancelled);

  TH1D *h_r04_unfold_zvertex60, *h_r04_unfold_zvertex60_stat, *h_r04_unfold_zvertex60_unreweighted, *h_r04_unfold_zvertex60_jesup, *h_r04_unfold_zvertex60_jesdown, *h_r04_unfold_zvertex60_jerup, *h_r04_unfold_zvertex60_jerdown, *h_r04_unfold_zvertex60_jetup, *h_r04_unfold_zvertex60_jetdown, *h_r04_unfold_zvertex60_unfoldunc;
  h_r04_unfold_zvertex60 = (TH1D*)f_r04_spectrum->Get("h_unfold_zvertex60"); h_r04_unfold_zvertex60_stat = (TH1D*)((TH1DBootstrap*)f_r04_spectrum->Get("h_unfold_zvertex60_stat"))->GetNominal(); h_r04_unfold_zvertex60_unreweighted = (TH1D*)f_r04_spectrum->Get("h_unfold_zvertex60_unreweighted"); h_r04_unfold_zvertex60_jesup = (TH1D*)f_r04_spectrum->Get("h_unfold_zvertex60_jesup"); h_r04_unfold_zvertex60_jesdown = (TH1D*)f_r04_spectrum->Get("h_unfold_zvertex60_jesdown"); h_r04_unfold_zvertex60_jerup = (TH1D*)f_r04_spectrum->Get("h_unfold_zvertex60_jerup"); h_r04_unfold_zvertex60_jerdown = (TH1D*)f_r04_spectrum->Get("h_unfold_zvertex60_jerdown"); h_r04_unfold_zvertex60_jetup = (TH1D*)f_r04_spectrum->Get("h_unfold_zvertex60_jetup"); h_r04_unfold_zvertex60_jetdown = (TH1D*)f_r04_spectrum->Get("h_unfold_zvertex60_jetdown"); h_r04_unfold_zvertex60_unfoldunc = (TH1D*)f_r04_spectrum->Get("h_unfold_zvertex60_unfoldunc");
  do_normalization(h_r04_unfold_zvertex60, luminosity22_60, 0.4); do_normalization(h_r04_unfold_zvertex60_stat, luminosity22_60, 0.4); do_normalization(h_r04_unfold_zvertex60_unreweighted, luminosity22_60, 0.4); do_normalization(h_r04_unfold_zvertex60_jesup, luminosity22_60, 0.4); do_normalization(h_r04_unfold_zvertex60_jesdown, luminosity22_60, 0.4); do_normalization(h_r04_unfold_zvertex60_jerup, luminosity22_60, 0.4); do_normalization(h_r04_unfold_zvertex60_jerdown, luminosity22_60, 0.4); do_normalization(h_r04_unfold_zvertex60_jetup, luminosity22_60, 0.4); do_normalization(h_r04_unfold_zvertex60_jetdown, luminosity22_60, 0.4); do_normalization(h_r04_unfold_zvertex60_unfoldunc, luminosity22_60, 0.4);
  TH1D* h_r04_uncertainty_zvertex60_stat = get_statistical_uncertainty(h_r04_unfold_zvertex60_stat, "h_r04_uncertainty_zvertex60_stat");
  TH1D* h_r04_uncertainty_zvertex60_jesup,* h_r04_uncertainty_zvertex60_jesdown; get_uncertainty(h_r04_uncertainty_zvertex60_jesup, h_r04_uncertainty_zvertex60_jesdown, h_r04_unfold_zvertex60, h_r04_unfold_zvertex60_jesup, h_r04_unfold_zvertex60_jesdown, "h_r04_uncertainty_zvertex60_jes");
  TH1D* h_r04_uncertainty_zvertex60_jerup,* h_r04_uncertainty_zvertex60_jerdown; get_uncertainty(h_r04_uncertainty_zvertex60_jerup, h_r04_uncertainty_zvertex60_jerdown, h_r04_unfold_zvertex60, h_r04_unfold_zvertex60_jerup, h_r04_unfold_zvertex60_jerdown, "h_r04_uncertainty_zvertex60_jer");
  TH1D* h_r04_uncertainty_zvertex60_jetup,* h_r04_uncertainty_zvertex60_jetdown; get_uncertainty(h_r04_uncertainty_zvertex60_jetup, h_r04_uncertainty_zvertex60_jetdown, h_r04_unfold_zvertex60, h_r04_unfold_zvertex60_jetup, h_r04_unfold_zvertex60_jetdown, "h_r04_uncertainty_zvertex60_jet");
  TH1D* h_r04_uncertainty_zvertex60_unfoldup,* h_r04_uncertainty_zvertex60_unfolddown; get_unfold_uncertainty(h_r04_uncertainty_zvertex60_unfoldup, h_r04_uncertainty_zvertex60_unfolddown, h_r04_unfold_zvertex60, h_r04_unfold_zvertex60_unfoldunc, "h_r04_uncertainty_zvertex60_unfold");
  TH1D* h_r04_uncertainty_zvertex60_priorityup,* h_r04_uncertainty_zvertex60_prioritydown; get_priority_uncertainty(h_r04_uncertainty_zvertex60_priorityup, h_r04_uncertainty_zvertex60_prioritydown, h_r04_unfold_zvertex60, h_r04_unfold_zvertex60_unreweighted, "h_r04_uncertainty_zvertex60_priority");
  std::vector<TH1D*> h_r04_uncertainty_zvertex60_up = {h_r04_uncertainty_zvertex60_jesup, h_r04_uncertainty_zvertex60_jerup, h_r04_uncertainty_zvertex60_jetup, h_r04_uncertainty_zvertex60_unfoldup, h_r04_uncertainty_zvertex60_priorityup};
  std::vector<TH1D*> h_r04_uncertainty_zvertex60_down = {h_r04_uncertainty_zvertex60_jesdown, h_r04_uncertainty_zvertex60_jerdown, h_r04_uncertainty_zvertex60_jetdown, h_r04_uncertainty_zvertex60_unfolddown, h_r04_uncertainty_zvertex60_prioritydown};
  TH1D* h_r04_uncertainty_zvertex60_total_up, *h_r04_uncertainty_zvertex60_total_down; get_total_uncertainty(h_r04_uncertainty_zvertex60_total_up, h_r04_uncertainty_zvertex60_total_down, h_r04_uncertainty_zvertex60_stat, h_r04_uncertainty_zvertex60_up, h_r04_uncertainty_zvertex60_down);
  std::vector<TH1D*> h_r04_uncertainty_zvertex60_up_cancelled = {h_r04_uncertainty_zvertex60_unfoldup, h_r04_uncertainty_zvertex60_priorityup};
  std::vector<TH1D*> h_r04_uncertainty_zvertex60_down_cancelled = {h_r04_uncertainty_zvertex60_unfolddown, h_r04_uncertainty_zvertex60_prioritydown};
  TH1D* h_r04_uncertainty_zvertex60_total_up_cancelled, *h_r04_uncertainty_zvertex60_total_down_cancelled; get_total_uncertainty(h_r04_uncertainty_zvertex60_total_up_cancelled, h_r04_uncertainty_zvertex60_total_down_cancelled, h_r04_uncertainty_zvertex60_stat, h_r04_uncertainty_zvertex60_up_cancelled, h_r04_uncertainty_zvertex60_down_cancelled);

  TH1D *h_r05_unfold_all, *h_r05_unfold_all_stat, *h_r05_unfold_all_unreweighted, *h_r05_unfold_all_jesup, *h_r05_unfold_all_jesdown, *h_r05_unfold_all_jerup, *h_r05_unfold_all_jerdown, *h_r05_unfold_all_jetup, *h_r05_unfold_all_jetdown, *h_r05_unfold_all_unfoldunc;
  h_r05_unfold_all = (TH1D*)f_r05_spectrum->Get("h_unfold_all"); h_r05_unfold_all_stat = (TH1D*)((TH1DBootstrap*)f_r05_spectrum->Get("h_unfold_all_stat"))->GetNominal(); h_r05_unfold_all_unreweighted = (TH1D*)f_r05_spectrum->Get("h_unfold_all_unreweighted"); h_r05_unfold_all_jesup = (TH1D*)f_r05_spectrum->Get("h_unfold_all_jesup"); h_r05_unfold_all_jesdown = (TH1D*)f_r05_spectrum->Get("h_unfold_all_jesdown"); h_r05_unfold_all_jerup = (TH1D*)f_r05_spectrum->Get("h_unfold_all_jerup"); h_r05_unfold_all_jerdown = (TH1D*)f_r05_spectrum->Get("h_unfold_all_jerdown"); h_r05_unfold_all_jetup = (TH1D*)f_r05_spectrum->Get("h_unfold_all_jetup"); h_r05_unfold_all_jetdown = (TH1D*)f_r05_spectrum->Get("h_unfold_all_jetdown"); h_r05_unfold_all_unfoldunc = (TH1D*)f_r05_spectrum->Get("h_unfold_all_unfoldunc");
  do_normalization(h_r05_unfold_all, luminosity22_all, 0.5); do_normalization(h_r05_unfold_all_stat, luminosity22_all, 0.5); do_normalization(h_r05_unfold_all_unreweighted, luminosity22_all, 0.5); do_normalization(h_r05_unfold_all_jesup, luminosity22_all, 0.5); do_normalization(h_r05_unfold_all_jesdown, luminosity22_all, 0.5); do_normalization(h_r05_unfold_all_jerup, luminosity22_all, 0.5); do_normalization(h_r05_unfold_all_jerdown, luminosity22_all, 0.5); do_normalization(h_r05_unfold_all_jetup, luminosity22_all, 0.5); do_normalization(h_r05_unfold_all_jetdown, luminosity22_all, 0.5); do_normalization(h_r05_unfold_all_unfoldunc, luminosity22_all, 0.5);
  TH1D* h_r05_uncertainty_all_stat = get_statistical_uncertainty(h_r05_unfold_all_stat, "h_r05_uncertainty_all_stat");
  TH1D* h_r05_uncertainty_all_jesup, *h_r05_uncertainty_all_jesdown; get_uncertainty(h_r05_uncertainty_all_jesup, h_r05_uncertainty_all_jesdown, h_r05_unfold_all, h_r05_unfold_all_jesup, h_r05_unfold_all_jesdown, "h_r05_uncertainty_all_jes");
  TH1D* h_r05_uncertainty_all_jerup, *h_r05_uncertainty_all_jerdown; get_uncertainty(h_r05_uncertainty_all_jerup, h_r05_uncertainty_all_jerdown, h_r05_unfold_all, h_r05_unfold_all_jerup, h_r05_unfold_all_jerdown, "h_r05_uncertainty_all_jer");
  TH1D* h_r05_uncertainty_all_jetup, *h_r05_uncertainty_all_jetdown; get_uncertainty(h_r05_uncertainty_all_jetup, h_r05_uncertainty_all_jetdown, h_r05_unfold_all, h_r05_unfold_all_jetup, h_r05_unfold_all_jetdown, "h_r05_uncertainty_all_jet");
  TH1D* h_r05_uncertainty_all_unfoldup, *h_r05_uncertainty_all_unfolddown; get_unfold_uncertainty(h_r05_uncertainty_all_unfoldup, h_r05_uncertainty_all_unfolddown, h_r05_unfold_all, h_r05_unfold_all_unfoldunc, "h_r05_uncertainty_all_unfold");
  TH1D* h_r05_uncertainty_all_priorityup, *h_r05_uncertainty_all_prioritydown; get_priority_uncertainty(h_r05_uncertainty_all_priorityup, h_r05_uncertainty_all_prioritydown, h_r05_unfold_all, h_r05_unfold_all_unreweighted, "h_r05_uncertainty_all_priority");
  std::vector<TH1D*> h_r05_uncertainty_all_up = {h_r05_uncertainty_all_jesup, h_r05_uncertainty_all_jerup, h_r05_uncertainty_all_jetup, h_r05_uncertainty_all_unfoldup, h_r05_uncertainty_all_priorityup};
  std::vector<TH1D*> h_r05_uncertainty_all_down = {h_r05_uncertainty_all_jesdown, h_r05_uncertainty_all_jerdown, h_r05_uncertainty_all_jetdown, h_r05_uncertainty_all_unfolddown, h_r05_uncertainty_all_prioritydown};
  TH1D* h_r05_uncertainty_all_total_up, *h_r05_uncertainty_all_total_down; get_total_uncertainty(h_r05_uncertainty_all_total_up, h_r05_uncertainty_all_total_down, h_r05_uncertainty_all_stat, h_r05_uncertainty_all_up, h_r05_uncertainty_all_down);

  TH1D *h_r06_unfold_all, *h_r06_unfold_all_stat, *h_r06_unfold_all_unreweighted, *h_r06_unfold_all_jesup, *h_r06_unfold_all_jesdown, *h_r06_unfold_all_jerup, *h_r06_unfold_all_jerdown, *h_r06_unfold_all_jetup, *h_r06_unfold_all_jetdown, *h_r06_unfold_all_unfoldunc;
  h_r06_unfold_all = (TH1D*)f_r06_spectrum->Get("h_unfold_all"); h_r06_unfold_all_stat = (TH1D*)((TH1DBootstrap*)f_r06_spectrum->Get("h_unfold_all_stat"))->GetNominal(); h_r06_unfold_all_unreweighted = (TH1D*)f_r06_spectrum->Get("h_unfold_all_unreweighted"); h_r06_unfold_all_jesup = (TH1D*)f_r06_spectrum->Get("h_unfold_all_jesup"); h_r06_unfold_all_jesdown = (TH1D*)f_r06_spectrum->Get("h_unfold_all_jesdown"); h_r06_unfold_all_jerup = (TH1D*)f_r06_spectrum->Get("h_unfold_all_jerup"); h_r06_unfold_all_jerdown = (TH1D*)f_r06_spectrum->Get("h_unfold_all_jerdown"); h_r06_unfold_all_jetup = (TH1D*)f_r06_spectrum->Get("h_unfold_all_jetup"); h_r06_unfold_all_jetdown = (TH1D*)f_r06_spectrum->Get("h_unfold_all_jetdown"); h_r06_unfold_all_unfoldunc = (TH1D*)f_r06_spectrum->Get("h_unfold_all_unfoldunc");
  do_normalization(h_r06_unfold_all, luminosity22_all, 0.6); do_normalization(h_r06_unfold_all_stat, luminosity22_all, 0.6); do_normalization(h_r06_unfold_all_unreweighted, luminosity22_all, 0.6); do_normalization(h_r06_unfold_all_jesup, luminosity22_all, 0.6); do_normalization(h_r06_unfold_all_jesdown, luminosity22_all, 0.6); do_normalization(h_r06_unfold_all_jerup, luminosity22_all, 0.6); do_normalization(h_r06_unfold_all_jerdown, luminosity22_all, 0.6); do_normalization(h_r06_unfold_all_jetup, luminosity22_all, 0.6); do_normalization(h_r06_unfold_all_jetdown, luminosity22_all, 0.6); do_normalization(h_r06_unfold_all_unfoldunc, luminosity22_all, 0.6);
  TH1D* h_r06_uncertainty_all_stat = get_statistical_uncertainty(h_r06_unfold_all_stat, "h_r06_uncertainty_all_stat");
  TH1D* h_r06_uncertainty_all_jesup, *h_r06_uncertainty_all_jesdown; get_uncertainty(h_r06_uncertainty_all_jesup, h_r06_uncertainty_all_jesdown, h_r06_unfold_all, h_r06_unfold_all_jesup, h_r06_unfold_all_jesdown, "h_r06_uncertainty_all_jes");
  TH1D* h_r06_uncertainty_all_jerup, *h_r06_uncertainty_all_jerdown; get_uncertainty(h_r06_uncertainty_all_jerup, h_r06_uncertainty_all_jerdown, h_r06_unfold_all, h_r06_unfold_all_jerup, h_r06_unfold_all_jerdown, "h_r06_uncertainty_all_jer");
  TH1D* h_r06_uncertainty_all_jetup, *h_r06_uncertainty_all_jetdown; get_uncertainty(h_r06_uncertainty_all_jetup, h_r06_uncertainty_all_jetdown, h_r06_unfold_all, h_r06_unfold_all_jetup, h_r06_unfold_all_jetdown, "h_r06_uncertainty_all_jet");
  TH1D* h_r06_uncertainty_all_unfoldup, *h_r06_uncertainty_all_unfolddown; get_unfold_uncertainty(h_r06_uncertainty_all_unfoldup, h_r06_uncertainty_all_unfolddown, h_r06_unfold_all, h_r06_unfold_all_unfoldunc, "h_r06_uncertainty_all_unfold");
  TH1D* h_r06_uncertainty_all_priorityup, *h_r06_uncertainty_all_prioritydown; get_priority_uncertainty(h_r06_uncertainty_all_priorityup, h_r06_uncertainty_all_prioritydown, h_r06_unfold_all, h_r06_unfold_all_unreweighted, "h_r06_uncertainty_all_priority");
  std::vector<TH1D*> h_r06_uncertainty_all_up = {h_r06_uncertainty_all_jesup, h_r06_uncertainty_all_jerup, h_r06_uncertainty_all_jetup, h_r06_uncertainty_all_unfoldup, h_r06_uncertainty_all_priorityup};
  std::vector<TH1D*> h_r06_uncertainty_all_down = {h_r06_uncertainty_all_jesdown, h_r06_uncertainty_all_jerdown, h_r06_uncertainty_all_jetdown, h_r06_uncertainty_all_unfolddown, h_r06_uncertainty_all_prioritydown};
  TH1D* h_r06_uncertainty_all_total_up, *h_r06_uncertainty_all_total_down; get_total_uncertainty(h_r06_uncertainty_all_total_up, h_r06_uncertainty_all_total_down, h_r06_uncertainty_all_stat, h_r06_uncertainty_all_up, h_r06_uncertainty_all_down);

  TH1D *h_r08_unfold_all, *h_r08_unfold_all_stat, *h_r08_unfold_all_unreweighted, *h_r08_unfold_all_jesup, *h_r08_unfold_all_jesdown, *h_r08_unfold_all_jerup, *h_r08_unfold_all_jerdown, *h_r08_unfold_all_jetup, *h_r08_unfold_all_jetdown, *h_r08_unfold_all_unfoldunc;
  h_r08_unfold_all = (TH1D*)f_r08_spectrum->Get("h_unfold_all"); h_r08_unfold_all_stat = (TH1D*)((TH1DBootstrap*)f_r08_spectrum->Get("h_unfold_all_stat"))->GetNominal(); h_r08_unfold_all_unreweighted = (TH1D*)f_r08_spectrum->Get("h_unfold_all_unreweighted"); h_r08_unfold_all_jesup = (TH1D*)f_r08_spectrum->Get("h_unfold_all_jesup"); h_r08_unfold_all_jesdown = (TH1D*)f_r08_spectrum->Get("h_unfold_all_jesdown"); h_r08_unfold_all_jerup = (TH1D*)f_r08_spectrum->Get("h_unfold_all_jerup"); h_r08_unfold_all_jerdown = (TH1D*)f_r08_spectrum->Get("h_unfold_all_jerdown"); h_r08_unfold_all_jetup = (TH1D*)f_r08_spectrum->Get("h_unfold_all_jetup"); h_r08_unfold_all_jetdown = (TH1D*)f_r08_spectrum->Get("h_unfold_all_jetdown"); h_r08_unfold_all_unfoldunc = (TH1D*)f_r08_spectrum->Get("h_unfold_all_unfoldunc");
  do_normalization(h_r08_unfold_all, luminosity22_all, 0.8); do_normalization(h_r08_unfold_all_stat, luminosity22_all, 0.8); do_normalization(h_r08_unfold_all_unreweighted, luminosity22_all, 0.8); do_normalization(h_r08_unfold_all_jesup, luminosity22_all, 0.8); do_normalization(h_r08_unfold_all_jesdown, luminosity22_all, 0.8); do_normalization(h_r08_unfold_all_jerup, luminosity22_all, 0.8); do_normalization(h_r08_unfold_all_jerdown, luminosity22_all, 0.8); do_normalization(h_r08_unfold_all_jetup, luminosity22_all, 0.8); do_normalization(h_r08_unfold_all_jetdown, luminosity22_all, 0.8); do_normalization(h_r08_unfold_all_unfoldunc, luminosity22_all, 0.8);
  TH1D* h_r08_uncertainty_all_stat = get_statistical_uncertainty(h_r08_unfold_all_stat, "h_r08_uncertainty_all_stat");
  TH1D* h_r08_uncertainty_all_jesup, *h_r08_uncertainty_all_jesdown; get_uncertainty(h_r08_uncertainty_all_jesup, h_r08_uncertainty_all_jesdown, h_r08_unfold_all, h_r08_unfold_all_jesup, h_r08_unfold_all_jesdown, "h_r08_uncertainty_all_jes");
  TH1D* h_r08_uncertainty_all_jerup, *h_r08_uncertainty_all_jerdown; get_uncertainty(h_r08_uncertainty_all_jerup, h_r08_uncertainty_all_jerdown, h_r08_unfold_all, h_r08_unfold_all_jerup, h_r08_unfold_all_jerdown, "h_r08_uncertainty_all_jer");
  TH1D* h_r08_uncertainty_all_jetup, *h_r08_uncertainty_all_jetdown; get_uncertainty(h_r08_uncertainty_all_jetup, h_r08_uncertainty_all_jetdown, h_r08_unfold_all, h_r08_unfold_all_jetup, h_r08_unfold_all_jetdown, "h_r08_uncertainty_all_jet");
  TH1D* h_r08_uncertainty_all_unfoldup, *h_r08_uncertainty_all_unfolddown; get_unfold_uncertainty(h_r08_uncertainty_all_unfoldup, h_r08_uncertainty_all_unfolddown, h_r08_unfold_all, h_r08_unfold_all_unfoldunc, "h_r08_uncertainty_all_unfold");
  TH1D* h_r08_uncertainty_all_priorityup, *h_r08_uncertainty_all_prioritydown; get_priority_uncertainty(h_r08_uncertainty_all_priorityup, h_r08_uncertainty_all_prioritydown, h_r08_unfold_all, h_r08_unfold_all_unreweighted, "h_r08_uncertainty_all_priority");
  std::vector<TH1D*> h_r08_uncertainty_all_up = {h_r08_uncertainty_all_jesup, h_r08_uncertainty_all_jerup, h_r08_uncertainty_all_jetup, h_r08_uncertainty_all_unfoldup, h_r08_uncertainty_all_priorityup};
  std::vector<TH1D*> h_r08_uncertainty_all_down = {h_r08_uncertainty_all_jesdown, h_r08_uncertainty_all_jerdown, h_r08_uncertainty_all_jetdown, h_r08_uncertainty_all_unfolddown, h_r08_uncertainty_all_prioritydown};
  TH1D* h_r08_uncertainty_all_total_up, *h_r08_uncertainty_all_total_down; get_total_uncertainty(h_r08_uncertainty_all_total_up, h_r08_uncertainty_all_total_down, h_r08_uncertainty_all_stat, h_r08_uncertainty_all_up, h_r08_uncertainty_all_down);

  std::vector<int> color;
  std::vector<string> legend;

  color = {kRed, kYellow+1, kAzure+2, kGreen+2, kMagenta};
  legend = {"Jet Energy Scale", "Jet Energy Resolution", "Jet Trigger Efficiency", "Unfolding Iteration", "Unfolding Priority"};
  draw_total_uncertainty(h_r02_uncertainty_all_total_up, h_r02_uncertainty_all_total_down, h_r02_uncertainty_all_stat, h_r02_uncertainty_all_up, h_r02_uncertainty_all_down, "r02_all", color, legend, 0.2);
  color.clear(); legend.clear();

  color = {kRed, kYellow+1, kAzure+2, kGreen+2, kMagenta};
  legend = {"Jet Energy Scale", "Jet Energy Resolution", "Jet Trigger Efficiency", "Unfolding Iteration", "Unfolding Priority"};
  draw_total_uncertainty(h_r03_uncertainty_all_total_up, h_r03_uncertainty_all_total_down, h_r03_uncertainty_all_stat, h_r03_uncertainty_all_up, h_r03_uncertainty_all_down, "r03_all", color, legend, 0.3);
  color.clear(); legend.clear();

  color = {kRed, kYellow+1, kAzure+2, kGreen+2, kMagenta};
  legend = {"Jet Energy Scale", "Jet Energy Resolution", "Jet Trigger Efficiency", "Unfolding Iteration", "Unfolding Priority"};
  draw_total_uncertainty(h_r04_uncertainty_all_total_up, h_r04_uncertainty_all_total_down, h_r04_uncertainty_all_stat, h_r04_uncertainty_all_up, h_r04_uncertainty_all_down, "r04_all", color, legend, 0.4);
  color.clear(); legend.clear();
  color = {kRed, kYellow+1, kAzure+2, kGreen+2, kMagenta, kOrange-2};
  legend = {"Jet Energy Scale", "Jet Energy Resolution", "Jet Trigger Efficiency", "Unfolding Iteration", "Unfolding Priority"};
  draw_total_uncertainty(h_r04_uncertainty_zvertex60_total_up, h_r04_uncertainty_zvertex60_total_down, h_r04_uncertainty_zvertex60_stat, h_r04_uncertainty_zvertex60_up, h_r04_uncertainty_zvertex60_down, "r04_zvertex60", color, legend, 0.4);
  color.clear(); legend.clear();

  color = {kRed, kYellow+1, kAzure+2, kGreen+2, kMagenta};
  legend = {"Jet Energy Scale", "Jet Energy Resolution", "Jet Trigger Efficiency", "Unfolding Iteration", "Unfolding Priority"};
  draw_total_uncertainty(h_r05_uncertainty_all_total_up, h_r05_uncertainty_all_total_down, h_r05_uncertainty_all_stat, h_r05_uncertainty_all_up, h_r05_uncertainty_all_down, "r05_all", color, legend, 0.6);
  color.clear(); legend.clear();

  color = {kRed, kYellow+1, kAzure+2, kGreen+2, kMagenta};
  legend = {"Jet Energy Scale", "Jet Energy Resolution", "Jet Trigger Efficiency", "Unfolding Iteration", "Unfolding Priority"};
  draw_total_uncertainty(h_r06_uncertainty_all_total_up, h_r06_uncertainty_all_total_down, h_r06_uncertainty_all_stat, h_r06_uncertainty_all_up, h_r06_uncertainty_all_down, "r06_all", color, legend, 0.6);
  color.clear(); legend.clear();

  color = {kRed, kYellow+1, kAzure+2, kGreen+2, kMagenta};
  legend = {"Jet Energy Scale", "Jet Energy Resolution", "Jet Trigger Efficiency", "Unfolding Iteration", "Unfolding Priority"};
  draw_total_uncertainty(h_r08_uncertainty_all_total_up, h_r08_uncertainty_all_total_down, h_r08_uncertainty_all_stat, h_r08_uncertainty_all_up, h_r08_uncertainty_all_down, "r08_all", color, legend, 0.8);
  color.clear(); legend.clear();

  // sPHENIX all result
  double g_r02_sphenix_all_x[calibnpt], g_r02_sphenix_all_xerrdown[calibnpt], g_r02_sphenix_all_xerrup[calibnpt], g_r02_sphenix_all_y[calibnpt], g_r02_sphenix_all_yerrdown[calibnpt], g_r02_sphenix_all_yerrup[calibnpt], g_r02_sphenix_all_stat[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    g_r02_sphenix_all_x[ib] = (calibptbins[ib] + calibptbins[ib+1]) / 2.;
    g_r02_sphenix_all_xerrdown[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r02_sphenix_all_xerrup[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r02_sphenix_all_y[ib] = h_r02_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r02_sphenix_all_yerrdown[ib] = -1 * h_r02_uncertainty_all_total_down->GetBinContent(ib + n_underflowbin + 1) * h_r02_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r02_sphenix_all_yerrup[ib] = h_r02_uncertainty_all_total_up->GetBinContent(ib + n_underflowbin + 1) * h_r02_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r02_sphenix_all_stat[ib] = h_r02_uncertainty_all_stat->GetBinContent(ib + n_underflowbin + 1) * h_r02_unfold_all->GetBinContent(ib + n_underflowbin + 1);
  }
  TGraphAsymmErrors* g_r02_sphenix_result = new TGraphAsymmErrors(calibnpt, g_r02_sphenix_all_x, g_r02_sphenix_all_y, g_r02_sphenix_all_xerrdown, g_r02_sphenix_all_xerrup, g_r02_sphenix_all_yerrdown, g_r02_sphenix_all_yerrup);
  g_r02_sphenix_result->SetMarkerStyle(20); g_r02_sphenix_result->SetMarkerSize(1.5); g_r02_sphenix_result->SetMarkerColor(kAzure+2); g_r02_sphenix_result->SetLineWidth(0); g_r02_sphenix_result->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_r02_sphenix_result_point = new TGraphErrors(calibnpt, g_r02_sphenix_all_x, g_r02_sphenix_all_y, g_r02_sphenix_all_xerrdown, g_r02_sphenix_all_stat);
  g_r02_sphenix_result_point->SetMarkerStyle(20); g_r02_sphenix_result_point->SetMarkerSize(1.5); g_r02_sphenix_result_point->SetMarkerColor(kAzure+2); g_r02_sphenix_result_point->SetLineColor(kAzure+2); g_r02_sphenix_result_point->SetLineWidth(2);
  f_out->cd(); g_r02_sphenix_result->Write("g_r02_sphenix_result"); g_r02_sphenix_result_point->Write("g_r02_sphenix_result_point");

  double g_r03_sphenix_all_x[calibnpt], g_r03_sphenix_all_xerrdown[calibnpt], g_r03_sphenix_all_xerrup[calibnpt], g_r03_sphenix_all_y[calibnpt], g_r03_sphenix_all_yerrdown[calibnpt], g_r03_sphenix_all_yerrup[calibnpt], g_r03_sphenix_all_stat[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    g_r03_sphenix_all_x[ib] = (calibptbins[ib] + calibptbins[ib+1]) / 2.;
    g_r03_sphenix_all_xerrdown[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r03_sphenix_all_xerrup[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r03_sphenix_all_y[ib] = h_r03_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r03_sphenix_all_yerrdown[ib] = -1 * h_r03_uncertainty_all_total_down->GetBinContent(ib + n_underflowbin + 1) * h_r03_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r03_sphenix_all_yerrup[ib] = h_r03_uncertainty_all_total_up->GetBinContent(ib + n_underflowbin + 1) * h_r03_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r03_sphenix_all_stat[ib] = h_r03_uncertainty_all_stat->GetBinContent(ib + n_underflowbin + 1) * h_r03_unfold_all->GetBinContent(ib + n_underflowbin + 1);
  }
  TGraphAsymmErrors* g_r03_sphenix_result = new TGraphAsymmErrors(calibnpt, g_r03_sphenix_all_x, g_r03_sphenix_all_y, g_r03_sphenix_all_xerrdown, g_r03_sphenix_all_xerrup, g_r03_sphenix_all_yerrdown, g_r03_sphenix_all_yerrup);
  g_r03_sphenix_result->SetMarkerStyle(20); g_r03_sphenix_result->SetMarkerSize(1.5); g_r03_sphenix_result->SetMarkerColor(kAzure+2); g_r03_sphenix_result->SetLineWidth(0); g_r03_sphenix_result->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_r03_sphenix_result_point = new TGraphErrors(calibnpt, g_r03_sphenix_all_x, g_r03_sphenix_all_y, g_r03_sphenix_all_xerrdown, g_r03_sphenix_all_stat);
  g_r03_sphenix_result_point->SetMarkerStyle(20); g_r03_sphenix_result_point->SetMarkerSize(1.5); g_r03_sphenix_result_point->SetMarkerColor(kAzure+2); g_r03_sphenix_result_point->SetLineColor(kAzure+2); g_r03_sphenix_result_point->SetLineWidth(2);
  f_out->cd(); g_r03_sphenix_result->Write("g_r03_sphenix_result"); g_r03_sphenix_result_point->Write("g_r03_sphenix_result_point");
  f_out_forratio->cd(); g_r03_sphenix_result->Write("g_r03_sphenix_result");

  double g_r04_sphenix_all_x[calibnpt], g_r04_sphenix_all_xerrdown[calibnpt], g_r04_sphenix_all_xerrup[calibnpt], g_r04_sphenix_all_y[calibnpt], g_r04_sphenix_all_yerrdown[calibnpt], g_r04_sphenix_all_yerrup[calibnpt], g_r04_sphenix_all_stat[calibnpt];
  double g_sphenix_yerrup_cancelled[calibnpt], g_sphenix_yerrdown_cancelled[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    g_r04_sphenix_all_x[ib] = (calibptbins[ib] + calibptbins[ib+1]) / 2.;
    g_r04_sphenix_all_xerrdown[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r04_sphenix_all_xerrup[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r04_sphenix_all_y[ib] = h_r04_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r04_sphenix_all_yerrdown[ib] = -1 * h_r04_uncertainty_all_total_down->GetBinContent(ib + n_underflowbin + 1) * h_r04_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r04_sphenix_all_yerrup[ib] = h_r04_uncertainty_all_total_up->GetBinContent(ib + n_underflowbin + 1) * h_r04_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r04_sphenix_all_stat[ib] = h_r04_uncertainty_all_stat->GetBinContent(ib + n_underflowbin + 1) * h_r04_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_sphenix_yerrdown_cancelled[ib] = h_r04_uncertainty_all_total_down_cancelled->GetBinContent(ib + n_underflowbin + 1) * h_r04_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_sphenix_yerrup_cancelled[ib] = h_r04_uncertainty_all_total_up_cancelled->GetBinContent(ib + n_underflowbin + 1) * h_r04_unfold_all->GetBinContent(ib + n_underflowbin + 1);
  }
  TGraphAsymmErrors* g_r04_sphenix_result = new TGraphAsymmErrors(calibnpt, g_r04_sphenix_all_x, g_r04_sphenix_all_y, g_r04_sphenix_all_xerrdown, g_r04_sphenix_all_xerrup, g_r04_sphenix_all_yerrdown, g_r04_sphenix_all_yerrup);
  g_r04_sphenix_result->SetMarkerStyle(20); g_r04_sphenix_result->SetMarkerSize(1.5); g_r04_sphenix_result->SetMarkerColor(kAzure+2); g_r04_sphenix_result->SetLineWidth(0); g_r04_sphenix_result->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_r04_sphenix_result_point = new TGraphErrors(calibnpt, g_r04_sphenix_all_x, g_r04_sphenix_all_y, g_r04_sphenix_all_xerrdown, g_r04_sphenix_all_stat);
  g_r04_sphenix_result_point->SetMarkerStyle(20); g_r04_sphenix_result_point->SetMarkerSize(1.5); g_r04_sphenix_result_point->SetMarkerColor(kAzure+2); g_r04_sphenix_result_point->SetLineColor(kAzure+2); g_r04_sphenix_result_point->SetLineWidth(2);
  TGraphAsymmErrors* g_r04_sphenix_result_cancelled = new TGraphAsymmErrors(calibnpt, g_r04_sphenix_all_x, g_r04_sphenix_all_y, g_r04_sphenix_all_xerrdown, g_r04_sphenix_all_xerrup, g_sphenix_yerrdown_cancelled, g_sphenix_yerrup_cancelled);
  g_r04_sphenix_result_cancelled->SetMarkerStyle(20); g_r04_sphenix_result_cancelled->SetMarkerSize(1.5); g_r04_sphenix_result_cancelled->SetMarkerColor(kAzure+2); g_r04_sphenix_result_cancelled->SetLineWidth(0); g_r04_sphenix_result_cancelled->SetFillColorAlpha(kAzure + 1, 0.5);
  f_out->cd(); g_r04_sphenix_result->Write("g_r04_sphenix_result"); g_r04_sphenix_result_point->Write("g_r04_sphenix_result_point"); g_r04_sphenix_result_cancelled->Write("g_r04_sphenix_result_cancelled");
  
  double g_r05_sphenix_all_x[calibnpt], g_r05_sphenix_all_xerrdown[calibnpt], g_r05_sphenix_all_xerrup[calibnpt], g_r05_sphenix_all_y[calibnpt], g_r05_sphenix_all_yerrdown[calibnpt], g_r05_sphenix_all_yerrup[calibnpt], g_r05_sphenix_all_stat[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    g_r05_sphenix_all_x[ib] = (calibptbins[ib] + calibptbins[ib+1]) / 2.;
    g_r05_sphenix_all_xerrdown[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r05_sphenix_all_xerrup[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r05_sphenix_all_y[ib] = h_r05_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r05_sphenix_all_yerrdown[ib] = -1 * h_r05_uncertainty_all_total_down->GetBinContent(ib + n_underflowbin + 1) * h_r05_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r05_sphenix_all_yerrup[ib] = h_r05_uncertainty_all_total_up->GetBinContent(ib + n_underflowbin + 1) * h_r05_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r05_sphenix_all_stat[ib] = h_r05_uncertainty_all_stat->GetBinContent(ib + n_underflowbin + 1) * h_r05_unfold_all->GetBinContent(ib + n_underflowbin + 1);
  }
  TGraphAsymmErrors* g_r05_sphenix_result = new TGraphAsymmErrors(calibnpt, g_r05_sphenix_all_x, g_r05_sphenix_all_y, g_r05_sphenix_all_xerrdown, g_r05_sphenix_all_xerrup, g_r05_sphenix_all_yerrdown, g_r05_sphenix_all_yerrup);
  g_r05_sphenix_result->SetMarkerStyle(20); g_r05_sphenix_result->SetMarkerSize(1.5); g_r05_sphenix_result->SetMarkerColor(kAzure+2); g_r05_sphenix_result->SetLineWidth(0); g_r05_sphenix_result->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_r05_sphenix_result_point = new TGraphErrors(calibnpt, g_r05_sphenix_all_x, g_r05_sphenix_all_y, g_r05_sphenix_all_xerrdown, g_r05_sphenix_all_stat);
  g_r05_sphenix_result_point->SetMarkerStyle(20); g_r05_sphenix_result_point->SetMarkerSize(1.5); g_r05_sphenix_result_point->SetMarkerColor(kAzure+2); g_r05_sphenix_result_point->SetLineColor(kAzure+2); g_r05_sphenix_result_point->SetLineWidth(2);
  f_out->cd(); g_r05_sphenix_result->Write("g_r05_sphenix_result"); g_r05_sphenix_result_point->Write("g_r05_sphenix_result_point");
  f_out_forratio->cd(); g_r05_sphenix_result->Write("g_r05_sphenix_result");

  double g_r06_sphenix_all_x[calibnpt], g_r06_sphenix_all_xerrdown[calibnpt], g_r06_sphenix_all_xerrup[calibnpt], g_r06_sphenix_all_y[calibnpt], g_r06_sphenix_all_yerrdown[calibnpt], g_r06_sphenix_all_yerrup[calibnpt], g_r06_sphenix_all_stat[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    g_r06_sphenix_all_x[ib] = (calibptbins[ib] + calibptbins[ib+1]) / 2.;
    g_r06_sphenix_all_xerrdown[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r06_sphenix_all_xerrup[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r06_sphenix_all_y[ib] = h_r06_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r06_sphenix_all_yerrdown[ib] = -1 * h_r06_uncertainty_all_total_down->GetBinContent(ib + n_underflowbin + 1) * h_r06_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r06_sphenix_all_yerrup[ib] = h_r06_uncertainty_all_total_up->GetBinContent(ib + n_underflowbin + 1) * h_r06_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r06_sphenix_all_stat[ib] = h_r06_uncertainty_all_stat->GetBinContent(ib + n_underflowbin + 1) * h_r06_unfold_all->GetBinContent(ib + n_underflowbin + 1);
  }
  TGraphAsymmErrors* g_r06_sphenix_result = new TGraphAsymmErrors(calibnpt, g_r06_sphenix_all_x, g_r06_sphenix_all_y, g_r06_sphenix_all_xerrdown, g_r06_sphenix_all_xerrup, g_r06_sphenix_all_yerrdown, g_r06_sphenix_all_yerrup);
  g_r06_sphenix_result->SetMarkerStyle(20); g_r06_sphenix_result->SetMarkerSize(1.5); g_r06_sphenix_result->SetMarkerColor(kAzure+2); g_r06_sphenix_result->SetLineWidth(0); g_r06_sphenix_result->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_r06_sphenix_result_point = new TGraphErrors(calibnpt, g_r06_sphenix_all_x, g_r06_sphenix_all_y, g_r06_sphenix_all_xerrdown, g_r06_sphenix_all_stat);
  g_r06_sphenix_result_point->SetMarkerStyle(20); g_r06_sphenix_result_point->SetMarkerSize(1.5); g_r06_sphenix_result_point->SetMarkerColor(kAzure+2); g_r06_sphenix_result_point->SetLineColor(kAzure+2); g_r06_sphenix_result_point->SetLineWidth(2);
  f_out->cd(); g_r06_sphenix_result->Write("g_r06_sphenix_result"); g_r06_sphenix_result_point->Write("g_r06_sphenix_result_point");

  double g_r08_sphenix_all_x[calibnpt], g_r08_sphenix_all_xerrdown[calibnpt], g_r08_sphenix_all_xerrup[calibnpt], g_r08_sphenix_all_y[calibnpt], g_r08_sphenix_all_yerrdown[calibnpt], g_r08_sphenix_all_yerrup[calibnpt], g_r08_sphenix_all_stat[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    g_r08_sphenix_all_x[ib] = (calibptbins[ib] + calibptbins[ib+1]) / 2.;
    g_r08_sphenix_all_xerrdown[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r08_sphenix_all_xerrup[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r08_sphenix_all_y[ib] = h_r08_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r08_sphenix_all_yerrdown[ib] = -1 * h_r08_uncertainty_all_total_down->GetBinContent(ib + n_underflowbin + 1) * h_r08_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r08_sphenix_all_yerrup[ib] = h_r08_uncertainty_all_total_up->GetBinContent(ib + n_underflowbin + 1) * h_r08_unfold_all->GetBinContent(ib + n_underflowbin + 1);
    g_r08_sphenix_all_stat[ib] = h_r08_uncertainty_all_stat->GetBinContent(ib + n_underflowbin + 1) * h_r08_unfold_all->GetBinContent(ib + n_underflowbin + 1);
  }
  TGraphAsymmErrors* g_r08_sphenix_result = new TGraphAsymmErrors(calibnpt, g_r08_sphenix_all_x, g_r08_sphenix_all_y, g_r08_sphenix_all_xerrdown, g_r08_sphenix_all_xerrup, g_r08_sphenix_all_yerrdown, g_r08_sphenix_all_yerrup);
  g_r08_sphenix_result->SetMarkerStyle(20); g_r08_sphenix_result->SetMarkerSize(1.5); g_r08_sphenix_result->SetMarkerColor(kAzure+2); g_r08_sphenix_result->SetLineWidth(0); g_r08_sphenix_result->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_r08_sphenix_result_point = new TGraphErrors(calibnpt, g_r08_sphenix_all_x, g_r08_sphenix_all_y, g_r08_sphenix_all_xerrdown, g_r08_sphenix_all_stat);
  g_r08_sphenix_result_point->SetMarkerStyle(20); g_r08_sphenix_result_point->SetMarkerSize(1.5); g_r08_sphenix_result_point->SetMarkerColor(kAzure+2); g_r08_sphenix_result_point->SetLineColor(kAzure+2); g_r08_sphenix_result_point->SetLineWidth(2);
  f_out->cd(); g_r08_sphenix_result->Write("g_r08_sphenix_result"); g_r08_sphenix_result_point->Write("g_r08_sphenix_result_point");

  // sPHENIX 60 zvertex result
  double g_r04_sphenix_zvertex60_x[calibnpt], g_r04_sphenix_zvertex60_xerrdown[calibnpt], g_r04_sphenix_zvertex60_xerrup[calibnpt], g_r04_sphenix_zvertex60_y[calibnpt], g_r04_sphenix_zvertex60_yerrdown[calibnpt], g_r04_sphenix_zvertex60_yerrup[calibnpt], g_r04_sphenix_zvertex60_stat[calibnpt];
  double g_sphenix_zvertex60_yerrdown_cancelled[calibnpt], g_sphenix_zvertex60_yerrup_cancelled[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    g_r04_sphenix_zvertex60_x[ib] = (calibptbins[ib] + calibptbins[ib+1]) / 2.;
    g_r04_sphenix_zvertex60_xerrdown[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r04_sphenix_zvertex60_xerrup[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_r04_sphenix_zvertex60_y[ib] = h_r04_unfold_zvertex60->GetBinContent(ib + n_underflowbin + 1);
    g_r04_sphenix_zvertex60_yerrdown[ib] = -1 * h_r04_uncertainty_zvertex60_total_down->GetBinContent(ib + n_underflowbin + 1) * h_r04_unfold_zvertex60->GetBinContent(ib + n_underflowbin + 1);
    g_r04_sphenix_zvertex60_yerrup[ib] = h_r04_uncertainty_zvertex60_total_up->GetBinContent(ib + n_underflowbin + 1) * h_r04_unfold_zvertex60->GetBinContent(ib + n_underflowbin + 1);
    g_r04_sphenix_zvertex60_stat[ib] = h_r04_uncertainty_zvertex60_stat->GetBinContent(ib + n_underflowbin + 1) * h_r04_unfold_zvertex60->GetBinContent(ib + n_underflowbin + 1);
    g_sphenix_zvertex60_yerrdown_cancelled[ib] = -1 * h_r04_uncertainty_zvertex60_total_down_cancelled->GetBinContent(ib + n_underflowbin + 1) * h_r04_unfold_zvertex60->GetBinContent(ib + n_underflowbin + 1);
    g_sphenix_zvertex60_yerrup_cancelled[ib] = h_r04_uncertainty_zvertex60_total_up_cancelled->GetBinContent(ib + n_underflowbin + 1) * h_r04_unfold_zvertex60->GetBinContent(ib + n_underflowbin + 1);
  }
  TGraphAsymmErrors* g_r04_sphenix_zvertex60_result = new TGraphAsymmErrors(calibnpt, g_r04_sphenix_zvertex60_x, g_r04_sphenix_zvertex60_y, g_r04_sphenix_zvertex60_xerrdown, g_r04_sphenix_zvertex60_xerrup, g_r04_sphenix_zvertex60_yerrdown, g_r04_sphenix_zvertex60_yerrup);
  g_r04_sphenix_zvertex60_result->SetMarkerStyle(20); g_r04_sphenix_zvertex60_result->SetMarkerSize(1.5); g_r04_sphenix_zvertex60_result->SetMarkerColor(kAzure+2); g_r04_sphenix_zvertex60_result->SetLineWidth(0); g_r04_sphenix_zvertex60_result->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_r04_sphenix_zvertex60_result_point = new TGraphErrors(calibnpt, g_r04_sphenix_zvertex60_x, g_r04_sphenix_zvertex60_y, g_r04_sphenix_zvertex60_xerrdown, g_r04_sphenix_zvertex60_stat);
  g_r04_sphenix_zvertex60_result_point->SetMarkerStyle(20); g_r04_sphenix_zvertex60_result_point->SetMarkerSize(1.5); g_r04_sphenix_zvertex60_result_point->SetMarkerColor(kAzure+2); g_r04_sphenix_zvertex60_result_point->SetLineColor(kAzure+2); g_r04_sphenix_zvertex60_result_point->SetLineWidth(2); 
  f_out->cd(); g_r04_sphenix_zvertex60_result->Write("g_r04_sphenix_zvertex60_result"); g_r04_sphenix_zvertex60_result_point->Write("g_r04_sphenix_zvertex60_result_point");

  // sPHENIX ALTZ / DEFAULT ratio result
  double g_sphenix_DAratio_x[calibnpt], g_sphenix_DAratio_xerrdown[calibnpt], g_sphenix_DAratio_xerrup[calibnpt], g_sphenix_DAratio_y[calibnpt], g_sphenix_DAratio_yerrdown[calibnpt], g_sphenix_DAratio_yerrup[calibnpt], g_sphenix_DAratio_stat[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    g_sphenix_DAratio_x[ib] = (calibptbins[ib] + calibptbins[ib+1]) / 2.;
    g_sphenix_DAratio_xerrdown[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_sphenix_DAratio_xerrup[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_sphenix_DAratio_y[ib] = g_r04_sphenix_zvertex60_y[ib] / g_r04_sphenix_all_y[ib];
    g_sphenix_DAratio_stat[ib] = g_sphenix_DAratio_y[ib] * sqrt(pow(g_r04_sphenix_zvertex60_stat[ib] / g_r04_sphenix_zvertex60_y[ib], 2) + pow(g_r04_sphenix_all_stat[ib] / g_r04_sphenix_all_y[ib], 2));
    g_sphenix_DAratio_yerrdown[ib] = g_sphenix_DAratio_y[ib] * sqrt(pow(g_sphenix_zvertex60_yerrdown_cancelled[ib] / g_r04_sphenix_zvertex60_y[ib], 2) + pow(g_sphenix_yerrup_cancelled[ib] / g_r04_sphenix_all_y[ib], 2));
    g_sphenix_DAratio_yerrup[ib] = g_sphenix_DAratio_y[ib] * sqrt(pow(g_sphenix_zvertex60_yerrup_cancelled[ib] / g_r04_sphenix_zvertex60_y[ib], 2) + pow(g_sphenix_yerrdown_cancelled[ib] / g_r04_sphenix_all_y[ib], 2));
  }
  TGraphAsymmErrors* g_sphenix_DAratio_result = new TGraphAsymmErrors(calibnpt, g_sphenix_DAratio_x, g_sphenix_DAratio_y, g_sphenix_DAratio_xerrdown, g_sphenix_DAratio_xerrup, g_sphenix_DAratio_yerrdown, g_sphenix_DAratio_yerrup);
  g_sphenix_DAratio_result->SetMarkerStyle(20); g_sphenix_DAratio_result->SetMarkerSize(1.5); g_sphenix_DAratio_result->SetMarkerColor(kAzure+2); g_sphenix_DAratio_result->SetLineWidth(0);  g_sphenix_DAratio_result->SetFillColorAlpha(kAzure+1, 0.5);
  TGraphErrors* g_sphenix_DAratio_result_point = new TGraphErrors(calibnpt, g_sphenix_DAratio_x, g_sphenix_DAratio_y, g_sphenix_DAratio_xerrdown, g_sphenix_DAratio_stat);
  g_sphenix_DAratio_result_point->SetMarkerStyle(20); g_sphenix_DAratio_result_point->SetMarkerSize(1.5); g_sphenix_DAratio_result_point->SetMarkerColor(kAzure+2); g_sphenix_DAratio_result_point->SetLineColor(kAzure+2); g_sphenix_DAratio_result_point->SetLineWidth(2);

  // Pythia result
  TFile *f_sim_specturm = new TFile("output_sim_r04.root", "READ");
  TH1D* h_pythia_spectrum = (TH1D*)f_sim_specturm->Get("h_truthjet_pt_record_all");
  h_pythia_spectrum->SetName("h_pythia_spectrum");
  for (int ib = 1; ib <= h_pythia_spectrum->GetNbinsX(); ++ib) {
    h_pythia_spectrum->SetBinContent(ib, h_pythia_spectrum->GetBinContent(ib)/ (float)(1.4 * h_pythia_spectrum->GetBinWidth(ib)));
    h_pythia_spectrum->SetBinError(ib, h_pythia_spectrum->GetBinError(ib)/ (float)(1.4 * h_pythia_spectrum->GetBinWidth(ib)));
  }
  double g_pythia_x[truthnpt], g_pythia_xerrdown[truthnpt], g_pythia_xerrup[truthnpt], g_pythia_y[truthnpt], g_pythia_yerrdown[truthnpt], g_pythia_yerrup[truthnpt];
  for (int ib = 0; ib < truthnpt; ++ib) {
    g_pythia_x[ib] = (truthptbins[ib] + truthptbins[ib+1]) / 2;
    g_pythia_xerrdown[ib] = (truthptbins[ib+1] - truthptbins[ib]) / 2;
    g_pythia_xerrup[ib] = (truthptbins[ib+1] - truthptbins[ib]) / 2;
    g_pythia_y[ib] = h_pythia_spectrum->GetBinContent(ib+1);
    g_pythia_yerrdown[ib] = h_pythia_spectrum->GetBinError(ib+1);
    g_pythia_yerrup[ib] = h_pythia_spectrum->GetBinError(ib+1);
  }
  TGraphAsymmErrors* g_pythia_result = new TGraphAsymmErrors(truthnpt, g_pythia_x, g_pythia_y, g_pythia_xerrdown, g_pythia_xerrup, g_pythia_yerrdown, g_pythia_yerrup);
  g_pythia_result->SetMarkerStyle(21);
  g_pythia_result->SetMarkerSize(1.5);
  g_pythia_result->SetMarkerColor(kMagenta+2);
  g_pythia_result->SetLineColor(kMagenta+2);

  // Herwig result
  TFile *f_sim_specturm_herwig = new TFile("herwigXsections.root", "READ");
  TH1D* h_herwig_spectrum_04 = (TH1D*)f_sim_specturm_herwig->Get("Herwigxsections");
  double g_herwig_x_04[11], g_herwig_y_04[11], g_herwig_xerrdown_04[11], g_herwig_xerrup_04[11], g_herwig_yerrdown_04[11], g_herwig_yerrup_04[11];
  for (int ib = 0; ib < 11; ++ib) {
    g_herwig_x_04[ib] = h_herwig_spectrum_04->GetBinCenter(ib+1);
    g_herwig_xerrdown_04[ib] = h_herwig_spectrum_04->GetBinWidth(ib+1) / 2;
    g_herwig_xerrup_04[ib] = h_herwig_spectrum_04->GetBinWidth(ib+1) / 2;
    g_herwig_y_04[ib] = h_herwig_spectrum_04->GetBinContent(ib+1);
    g_herwig_yerrdown_04[ib] = h_herwig_spectrum_04->GetBinError(ib+1);
    g_herwig_yerrup_04[ib] = h_herwig_spectrum_04->GetBinError(ib+1);
  }
  TGraphAsymmErrors* g_herwig_result_04 = new TGraphAsymmErrors(11, g_herwig_x_04, g_herwig_y_04, g_herwig_xerrdown_04, g_herwig_xerrup_04, g_herwig_yerrdown_04, g_herwig_yerrup_04);
  g_herwig_result_04->SetMarkerStyle(33);
  g_herwig_result_04->SetMarkerSize(1.5);
  g_herwig_result_04->SetMarkerColor(kGreen+2);
  g_herwig_result_04->SetLineColor(kGreen+2);

  // SCET result
  ifstream myfile_scet;
  TGraphErrors *tjet_scet = new TGraphErrors();
  TGraphErrors *tjet_scet_point = new TGraphErrors();
  myfile_scet.open("sphenixScet.dat");
  double scet_x, scet_xerr, scet_xerr1, scet_xerr2, scet_y, scet_yerr;
  for (int index = 0; index < 10; ++index) {
    myfile_scet >> scet_x >> scet_xerr1 >> scet_xerr2 >> scet_y >> scet_yerr;
    scet_xerr = scet_xerr2 - scet_x;

    tjet_scet->SetPoint(index, scet_x, scet_y);
    tjet_scet->SetPointError(index, scet_xerr, scet_yerr);

    tjet_scet_point->SetPoint(index, scet_x, scet_y);
    tjet_scet_point->SetPointError(index, 0, 0);
  }
  tjet_scet->SetMarkerStyle(20);
  tjet_scet->SetMarkerSize(1.5);
  tjet_scet->SetMarkerColor(kGreen+2);
  tjet_scet->SetLineColor(kGreen+2);
  tjet_scet->SetFillColorAlpha(kGreen+1, 0.5);
  tjet_scet_point->SetMarkerStyle(20);
  tjet_scet_point->SetMarkerSize(1.5);
  tjet_scet_point->SetMarkerColor(kGreen+2);
  tjet_scet_point->SetLineColor(kGreen+2);
  tjet_scet_point->SetLineWidth(2);
  myfile_scet.close();

  // NNLO result
  ifstream myfile_nnlo;
  TGraphErrors *tjet_nnlo = new TGraphErrors();
  TGraphErrors *tjet_nnlo_point = new TGraphErrors();
  myfile_nnlo.open("nnloJet2.dat");
  double nnlo_x, nnlo_xerr1, nnlo_xerr2, nnlo_y, nnlo_yerr, nnlo_drop;
  for (int index = 0; index < 9; ++index) {
    myfile_nnlo >> nnlo_x >> nnlo_xerr1 >> nnlo_xerr2 >> nnlo_y >> nnlo_yerr >> nnlo_drop;
    double nnlo_xerr = nnlo_xerr2 - nnlo_x;

    tjet_nnlo->SetPoint(index, nnlo_x, nnlo_y/1000.); // Divide 1000 to convert from fb to pb
    tjet_nnlo->SetPointError(index, nnlo_xerr, nnlo_yerr/1000.);

    tjet_nnlo_point->SetPoint(index, nnlo_x, nnlo_y/1000.);
    tjet_nnlo_point->SetPointError(index, 0, 0);
  }
  tjet_nnlo->SetMarkerStyle(20);
  tjet_nnlo->SetMarkerSize(1.5);
  tjet_nnlo->SetMarkerColor(kYellow+2);
  tjet_nnlo->SetLineColor(kYellow+2);
  tjet_nnlo->SetFillColorAlpha(kYellow+1, 0.5);
  tjet_nnlo_point->SetMarkerStyle(20);
  tjet_nnlo_point->SetMarkerSize(1.5);
  tjet_nnlo_point->SetMarkerColor(kYellow+2);
  tjet_nnlo_point->SetLineColor(kYellow+2);
  myfile_nnlo.close();

  // NLO result
  ifstream infile1, infile05, infile2;
  infile1.open("sPh-R-sc1.dat"); infile05.open("sPh-R-sc05.dat"); infile2.open("sPh-R-sc2.dat");
  TGraph *t_r02jet = new TGraph(); TGraph *t_r02jet05 = new TGraph(); TGraph *t_r02jet2 = new TGraph();
  TGraph *t_r03jet = new TGraph(); TGraph *t_r03jet05 = new TGraph(); TGraph *t_r03jet2 = new TGraph();
  TGraph *t_r04jet = new TGraph(); TGraph *t_r04jet05 = new TGraph(); TGraph *t_r04jet2 = new TGraph();
  TGraph *t_r05jet = new TGraph(); TGraph *t_r05jet05 = new TGraph(); TGraph *t_r05jet2 = new TGraph();
  TGraph *t_r06jet = new TGraph(); TGraph *t_r06jet05 = new TGraph(); TGraph *t_r06jet2 = new TGraph();
  TGraph *t_r07jet = new TGraph(); TGraph *t_r07jet05 = new TGraph(); TGraph *t_r07jet2 = new TGraph();
  TGraph *t_r08jet = new TGraph(); TGraph *t_r08jet05 = new TGraph(); TGraph *t_r08jet2 = new TGraph();
  double radius, drop4, drop5, voge_pt1, xsection, voge_pt05, xsection05, voge_pt2, xsection2;
  std::vector<double> all_xpoint, r02_ypoint05, r02_ypoint1, r02_ypoint2, r03_ypoint05, r03_ypoint1, r03_ypoint2, r04_ypoint05, r04_ypoint1, r04_ypoint2, r05_ypoint05, r05_ypoint1, r05_ypoint2, r06_ypoint05, r06_ypoint1, r06_ypoint2, r07_ypoint05, r07_ypoint1, r07_ypoint2, r08_ypoint05, r08_ypoint1, r08_ypoint2;
  std::vector<double> r02_ratio2, r02_ratio05, r03_ratio2, r03_ratio05, r04_ratio2, r04_ratio05, r05_ratio2, r05_ratio05, r06_ratio2, r06_ratio05, r07_ratio2, r07_ratio05, r08_ratio2, r08_ratio05;
  int count_r02 = 0, count_r03 = 0, count_r04 = 0, count_r05 = 0, count_r06 = 0, count_r07 = 0, count_r08 = 0;
  for (int line = 0; line < 259; line ++) {
    infile1 >> radius >> drop4 >> drop5 >> voge_pt1 >> xsection;
    infile05 >> radius >> drop4 >> drop5 >> voge_pt05 >> xsection05;
    infile2 >> radius >> drop4 >> drop5 >> voge_pt2 >> xsection2;
    if (voge_pt1 != voge_pt05 || voge_pt1 != voge_pt2) {
      cout << "Error: pt values do not match at line in new vogelsang: " << line << endl;
      continue;
    }

    if (radius == 0.2) {
      t_r02jet->SetPoint(count_r02, voge_pt1, xsection);
      t_r02jet05->SetPoint(count_r02, voge_pt05, xsection05);
      t_r02jet2->SetPoint(count_r02, voge_pt2, xsection2);
      all_xpoint.push_back(voge_pt1);
      r02_ypoint05.push_back(xsection05);
      r02_ypoint1.push_back(xsection);
      r02_ypoint2.push_back(xsection2);
      r02_ratio05.push_back(xsection05/xsection);
      r02_ratio2.push_back(xsection2/xsection);
      count_r02++;
    } else if (radius == 0.3) {
      t_r03jet->SetPoint(count_r03, voge_pt1, xsection);
      t_r03jet05->SetPoint(count_r03, voge_pt05, xsection05);
      t_r03jet2->SetPoint(count_r03, voge_pt2, xsection2);
      r03_ypoint05.push_back(xsection05);
      r03_ypoint1.push_back(xsection);
      r03_ypoint2.push_back(xsection2);
      r03_ratio05.push_back(xsection05/xsection);
      r03_ratio2.push_back(xsection2/xsection);
      count_r03++;
    } else if (radius == 0.4) {
      t_r04jet->SetPoint(count_r04, voge_pt1, xsection);
      t_r04jet05->SetPoint(count_r04, voge_pt05, xsection05);
      t_r04jet2->SetPoint(count_r04, voge_pt2, xsection2);
      r04_ypoint05.push_back(xsection05);
      r04_ypoint1.push_back(xsection);
      r04_ypoint2.push_back(xsection2);
      r04_ratio05.push_back(xsection05/xsection);
      r04_ratio2.push_back(xsection2/xsection);
      count_r04++;
    } else if (radius == 0.5) {
      t_r05jet->SetPoint(count_r05, voge_pt1, xsection);
      t_r05jet05->SetPoint(count_r05, voge_pt05, xsection05);
      t_r05jet2->SetPoint(count_r05, voge_pt2, xsection2);
      r05_ypoint05.push_back(xsection05);
      r05_ypoint1.push_back(xsection);
      r05_ypoint2.push_back(xsection2);
      r05_ratio05.push_back(xsection05/xsection);
      r05_ratio2.push_back(xsection2/xsection);
      count_r05++;
    } else if (radius == 0.6) {
      t_r06jet->SetPoint(count_r06, voge_pt1, xsection);
      t_r06jet05->SetPoint(count_r06, voge_pt05, xsection05);
      t_r06jet2->SetPoint(count_r06, voge_pt2, xsection2);
      r06_ypoint05.push_back(xsection05);
      r06_ypoint1.push_back(xsection);
      r06_ypoint2.push_back(xsection2);
      r06_ratio05.push_back(xsection05/xsection);
      r06_ratio2.push_back(xsection2/xsection);
      count_r06++;
    } else if (radius == 0.7) {
      t_r07jet->SetPoint(count_r07, voge_pt1, xsection);
      t_r07jet05->SetPoint(count_r07, voge_pt05, xsection05);
      t_r07jet2->SetPoint(count_r07, voge_pt2, xsection2);
      r07_ypoint05.push_back(xsection05);
      r07_ypoint1.push_back(xsection);
      r07_ypoint2.push_back(xsection2);
      r07_ratio05.push_back(xsection05/xsection);
      r07_ratio2.push_back(xsection2/xsection);
      count_r07++;
    } else if (radius == 0.8) {
      t_r08jet->SetPoint(count_r08, voge_pt1, xsection);
      t_r08jet05->SetPoint(count_r08, voge_pt05, xsection05);
      t_r08jet2->SetPoint(count_r08, voge_pt2, xsection2);
      r08_ypoint05.push_back(xsection05);
      r08_ypoint1.push_back(xsection);
      r08_ypoint2.push_back(xsection2);
      r08_ratio05.push_back(xsection05/xsection);
      r08_ratio2.push_back(xsection2/xsection);
      count_r08++;
    } else {
      cout << "Error: radius value does not match expected values in new vogelsang: " << radius << " at line " << line << endl;
    }
  }
  if (count_r02 != count_r03 || count_r02 != count_r04 || count_r02 != count_r05 || 
      count_r02 != count_r06 || count_r02 != count_r07 || count_r02 != count_r08) {
    cout << "Error: different number of pT points per radius in Vogelsang file!" << endl;
    return;
  }

  TSpline3 *spline_r02 = new TSpline3("spline_r02", t_r02jet); TSpline3 *spline_r02_05 = new TSpline3("spline_r02_05", t_r02jet05); TSpline3 *spline_r02_2 = new TSpline3("spline_r02_2", t_r02jet2);
  TSpline3 *spline_r03 = new TSpline3("spline_r03", t_r03jet); TSpline3 *spline_r03_05 = new TSpline3("spline_r03_05", t_r03jet05); TSpline3 *spline_r03_2 = new TSpline3("spline_r03_2", t_r03jet2);
  TSpline3 *spline_r04 = new TSpline3("spline_r04", t_r04jet); TSpline3 *spline_r04_05 = new TSpline3("spline_r04_05", t_r04jet05); TSpline3 *spline_r04_2 = new TSpline3("spline_r04_2", t_r04jet2);
  TSpline3 *spline_r05 = new TSpline3("spline_r05", t_r05jet); TSpline3 *spline_r05_05 = new TSpline3("spline_r05_05", t_r05jet05); TSpline3 *spline_r05_2 = new TSpline3("spline_r05_2", t_r05jet2);
  TSpline3 *spline_r06 = new TSpline3("spline_r06", t_r06jet); TSpline3 *spline_r06_05 = new TSpline3("spline_r06_05", t_r06jet05); TSpline3 *spline_r06_2 = new TSpline3("spline_r06_2", t_r06jet2);
  TSpline3 *spline_r07 = new TSpline3("spline_r07", t_r07jet); TSpline3 *spline_r07_05 = new TSpline3("spline_r07_05", t_r07jet05); TSpline3 *spline_r07_2 = new TSpline3("spline_r07_2", t_r07jet2);
  TSpline3 *spline_r08 = new TSpline3("spline_r08", t_r08jet); TSpline3 *spline_r08_05 = new TSpline3("spline_r08_05", t_r08jet05); TSpline3 *spline_r08_2 = new TSpline3("spline_r08_2", t_r08jet2);

  std::vector <double> all_xpoint_shadow, r02_ypoint_shadow, r03_ypoint_shadow, r04_ypoint_shadow, r05_ypoint_shadow, r06_ypoint_shadow, r07_ypoint_shadow, r08_ypoint_shadow;
  std::vector <double> r02_ratio_shadow, r03_ratio_shadow, r04_ratio_shadow, r05_ratio_shadow, r06_ratio_shadow, r07_ratio_shadow, r08_ratio_shadow;
  for (int i = 0; i < all_xpoint.size(); i++) {
    all_xpoint_shadow.push_back(all_xpoint[i]);
    r02_ypoint_shadow.push_back(std::min({r02_ypoint1[i], r02_ypoint05[i], r02_ypoint2[i]}));
    r03_ypoint_shadow.push_back(std::min({r03_ypoint1[i], r03_ypoint05[i], r03_ypoint2[i]}));
    r04_ypoint_shadow.push_back(std::min({r04_ypoint1[i], r04_ypoint05[i], r04_ypoint2[i]}));
    r05_ypoint_shadow.push_back(std::min({r05_ypoint1[i], r05_ypoint05[i], r05_ypoint2[i]}));
    r06_ypoint_shadow.push_back(std::min({r06_ypoint1[i], r06_ypoint05[i], r06_ypoint2[i]}));
    r07_ypoint_shadow.push_back(std::min({r07_ypoint1[i], r07_ypoint05[i], r07_ypoint2[i]}));
    r08_ypoint_shadow.push_back(std::min({r08_ypoint1[i], r08_ypoint05[i], r08_ypoint2[i]}));

    r02_ratio_shadow.push_back(std::min({1.0, r02_ratio05[i], r02_ratio2[i]}));
    r03_ratio_shadow.push_back(std::min({1.0, r03_ratio05[i], r03_ratio2[i]}));
    r04_ratio_shadow.push_back(std::min({1.0, r04_ratio05[i], r04_ratio2[i]}));
    r05_ratio_shadow.push_back(std::min({1.0, r05_ratio05[i], r05_ratio2[i]}));
    r06_ratio_shadow.push_back(std::min({1.0, r06_ratio05[i], r06_ratio2[i]}));
    r07_ratio_shadow.push_back(std::min({1.0, r07_ratio05[i], r07_ratio2[i]}));
    r08_ratio_shadow.push_back(std::min({1.0, r08_ratio05[i], r08_ratio2[i]}));
  }
  for (int i = all_xpoint.size(); i > 0; i--) {
    all_xpoint_shadow.push_back(all_xpoint[i-1]);

    r02_ypoint_shadow.push_back(std::max({r02_ypoint1[i-1], r02_ypoint05[i-1], r02_ypoint2[i-1]}));
    r03_ypoint_shadow.push_back(std::max({r03_ypoint1[i-1], r03_ypoint05[i-1], r03_ypoint2[i-1]}));
    r04_ypoint_shadow.push_back(std::max({r04_ypoint1[i-1], r04_ypoint05[i-1], r04_ypoint2[i-1]}));
    r05_ypoint_shadow.push_back(std::max({r05_ypoint1[i-1], r05_ypoint05[i-1], r05_ypoint2[i-1]}));
    r06_ypoint_shadow.push_back(std::max({r06_ypoint1[i-1], r06_ypoint05[i-1], r06_ypoint2[i-1]}));
    r07_ypoint_shadow.push_back(std::max({r07_ypoint1[i-1], r07_ypoint05[i-1], r07_ypoint2[i-1]}));
    r08_ypoint_shadow.push_back(std::max({r08_ypoint1[i-1], r08_ypoint05[i-1], r08_ypoint2[i-1]}));

    r02_ratio_shadow.push_back(std::max({1.0, r02_ratio05[i-1], r02_ratio2[i-1]}));
    r03_ratio_shadow.push_back(std::max({1.0, r03_ratio05[i-1], r03_ratio2[i-1]}));
    r04_ratio_shadow.push_back(std::max({1.0, r04_ratio05[i-1], r04_ratio2[i-1]}));
    r05_ratio_shadow.push_back(std::max({1.0, r05_ratio05[i-1], r05_ratio2[i-1]}));
    r06_ratio_shadow.push_back(std::max({1.0, r06_ratio05[i-1], r06_ratio2[i-1]}));
    r07_ratio_shadow.push_back(std::max({1.0, r07_ratio05[i-1], r07_ratio2[i-1]}));
    r08_ratio_shadow.push_back(std::max({1.0, r08_ratio05[i-1], r08_ratio2[i-1]}));
  }
  TGraph *t_r02jet_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r02_ypoint_shadow[0]); t_r02jet_shadow->SetLineColor(kGreen); t_r02jet_shadow->SetFillColorAlpha(kGreen+1, 0.5);
  TGraph *t_r03jet_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r03_ypoint_shadow[0]); t_r03jet_shadow->SetLineColor(kGreen); t_r03jet_shadow->SetFillColorAlpha(kGreen+1, 0.5);
  TGraph *t_r04jet_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r04_ypoint_shadow[0]); t_r04jet_shadow->SetLineColor(kGreen); t_r04jet_shadow->SetFillColorAlpha(kGreen+1, 0.5);
  TGraph *t_r05jet_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r05_ypoint_shadow[0]); t_r05jet_shadow->SetLineColor(kGreen); t_r05jet_shadow->SetFillColorAlpha(kGreen+1, 0.5);
  TGraph *t_r06jet_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r06_ypoint_shadow[0]); t_r06jet_shadow->SetLineColor(kGreen); t_r06jet_shadow->SetFillColorAlpha(kGreen+1, 0.5);
  TGraph *t_r07jet_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r07_ypoint_shadow[0]); t_r07jet_shadow->SetLineColor(kGreen); t_r07jet_shadow->SetFillColorAlpha(kGreen+1, 0.5);
  TGraph *t_r08jet_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r08_ypoint_shadow[0]); t_r08jet_shadow->SetLineColor(kGreen); t_r08jet_shadow->SetFillColorAlpha(kGreen+1, 0.5);

  TGraph *t_r02jet_ratio_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r02_ratio_shadow[0]); t_r02jet_ratio_shadow->SetLineWidth(kGreen); t_r02jet_ratio_shadow->SetFillColorAlpha(kGreen+1, 0.5);
  TGraph *t_r03jet_ratio_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r03_ratio_shadow[0]); t_r03jet_ratio_shadow->SetLineWidth(kGreen); t_r03jet_ratio_shadow->SetFillColorAlpha(kGreen+1, 0.5);
  TGraph *t_r04jet_ratio_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r04_ratio_shadow[0]); t_r04jet_ratio_shadow->SetLineWidth(kGreen); t_r04jet_ratio_shadow->SetFillColorAlpha(kGreen+1, 0.5);
  TGraph *t_r05jet_ratio_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r05_ratio_shadow[0]); t_r05jet_ratio_shadow->SetLineWidth(kGreen); t_r05jet_ratio_shadow->SetFillColorAlpha(kGreen+1, 0.5);
  TGraph *t_r06jet_ratio_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r06_ratio_shadow[0]); t_r06jet_ratio_shadow->SetLineWidth(kGreen); t_r06jet_ratio_shadow->SetFillColorAlpha(kGreen+1, 0.5);
  TGraph *t_r07jet_ratio_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r07_ratio_shadow[0]); t_r07jet_ratio_shadow->SetLineWidth(kGreen); t_r07jet_ratio_shadow->SetFillColorAlpha(kGreen+1, 0.5);
  TGraph *t_r08jet_ratio_shadow = new TGraph(2*all_xpoint.size(), &all_xpoint_shadow[0], &r08_ratio_shadow[0]); t_r08jet_ratio_shadow->SetLineWidth(kGreen); t_r08jet_ratio_shadow->SetFillColorAlpha(kGreen+1, 0.5);
  infile1.close(); infile05.close(); infile2.close();

  // New Vogelsang
  ifstream infile_voge, infile_nnpdf, infile_msht;
  infile_voge.open("sPh-R-sc1-band.dat"); infile_nnpdf.open("sPh-R-sc1-nnpdf30.dat"); infile_msht.open("sPh-R-sc1-msht20.dat");
  TGraph *t_r02jet_voge = new TGraph(); TGraph *t_r02jet_voge_low = new TGraph(); TGraph *t_r02jet_voge_high = new TGraph(); std::vector<double> r02_ypoint_voge, r02_ypoint_voge_low, r02_ypoint_voge_high, r02_ypoint_voge_ratiolow, r02_ypoint_voge_ratiohigh;
  TGraph *t_r03jet_voge = new TGraph(); TGraph *t_r03jet_voge_low = new TGraph(); TGraph *t_r03jet_voge_high = new TGraph(); std::vector<double> r03_ypoint_voge, r03_ypoint_voge_low, r03_ypoint_voge_high, r03_ypoint_voge_ratiolow, r03_ypoint_voge_ratiohigh;
  TGraph *t_r04jet_voge = new TGraph(); TGraph *t_r04jet_voge_low = new TGraph(); TGraph *t_r04jet_voge_high = new TGraph(); std::vector<double> r04_ypoint_voge, r04_ypoint_voge_low, r04_ypoint_voge_high, r04_ypoint_voge_ratiolow, r04_ypoint_voge_ratiohigh;
  TGraph *t_r05jet_voge = new TGraph(); TGraph *t_r05jet_voge_low = new TGraph(); TGraph *t_r05jet_voge_high = new TGraph(); std::vector<double> r05_ypoint_voge, r05_ypoint_voge_low, r05_ypoint_voge_high, r05_ypoint_voge_ratiolow, r05_ypoint_voge_ratiohigh;
  TGraph *t_r06jet_voge = new TGraph(); TGraph *t_r06jet_voge_low = new TGraph(); TGraph *t_r06jet_voge_high = new TGraph(); std::vector<double> r06_ypoint_voge, r06_ypoint_voge_low, r06_ypoint_voge_high, r06_ypoint_voge_ratiolow, r06_ypoint_voge_ratiohigh;
  TGraph *t_r07jet_voge = new TGraph(); TGraph *t_r07jet_voge_low = new TGraph(); TGraph *t_r07jet_voge_high = new TGraph(); std::vector<double> r07_ypoint_voge, r07_ypoint_voge_low, r07_ypoint_voge_high, r07_ypoint_voge_ratiolow, r07_ypoint_voge_ratiohigh;
  TGraph *t_r08jet_voge = new TGraph(); TGraph *t_r08jet_voge_low = new TGraph(); TGraph *t_r08jet_voge_high = new TGraph(); std::vector<double> r08_ypoint_voge, r08_ypoint_voge_low, r08_ypoint_voge_high, r08_ypoint_voge_ratiolow, r08_ypoint_voge_ratiohigh;
  TGraph *t_r02jet_nnpdf = new TGraph(); TGraph *t_r02jet_nnpdf_ratio = new TGraph(); std::vector<double> r02_ypoint_nnpdf, r02_ypoint_nnpdf_ratio;
  TGraph *t_r03jet_nnpdf = new TGraph(); TGraph *t_r03jet_nnpdf_ratio = new TGraph(); std::vector<double> r03_ypoint_nnpdf, r03_ypoint_nnpdf_ratio;
  TGraph *t_r04jet_nnpdf = new TGraph(); TGraph *t_r04jet_nnpdf_ratio = new TGraph(); std::vector<double> r04_ypoint_nnpdf, r04_ypoint_nnpdf_ratio;
  TGraph *t_r05jet_nnpdf = new TGraph(); TGraph *t_r05jet_nnpdf_ratio = new TGraph(); std::vector<double> r05_ypoint_nnpdf, r05_ypoint_nnpdf_ratio;
  TGraph *t_r06jet_nnpdf = new TGraph(); TGraph *t_r06jet_nnpdf_ratio = new TGraph(); std::vector<double> r06_ypoint_nnpdf, r06_ypoint_nnpdf_ratio;
  TGraph *t_r07jet_nnpdf = new TGraph(); TGraph *t_r07jet_nnpdf_ratio = new TGraph(); std::vector<double> r07_ypoint_nnpdf, r07_ypoint_nnpdf_ratio;
  TGraph *t_r08jet_nnpdf = new TGraph(); TGraph *t_r08jet_nnpdf_ratio = new TGraph(); std::vector<double> r08_ypoint_nnpdf, r08_ypoint_nnpdf_ratio;
  TGraph *t_r02jet_msht = new TGraph(); TGraph *t_r02jet_msht_ratio = new TGraph(); std::vector<double> r02_ypoint_msht, r02_ypoint_msht_ratio;
  TGraph *t_r03jet_msht = new TGraph(); TGraph *t_r03jet_msht_ratio = new TGraph(); std::vector<double> r03_ypoint_msht, r03_ypoint_msht_ratio;
  TGraph *t_r04jet_msht = new TGraph(); TGraph *t_r04jet_msht_ratio = new TGraph(); std::vector<double> r04_ypoint_msht, r04_ypoint_msht_ratio;
  TGraph *t_r05jet_msht = new TGraph(); TGraph *t_r05jet_msht_ratio = new TGraph(); std::vector<double> r05_ypoint_msht, r05_ypoint_msht_ratio;
  TGraph *t_r06jet_msht = new TGraph(); TGraph *t_r06jet_msht_ratio = new TGraph(); std::vector<double> r06_ypoint_msht, r06_ypoint_msht_ratio;
  TGraph *t_r07jet_msht = new TGraph(); TGraph *t_r07jet_msht_ratio = new TGraph(); std::vector<double> r07_ypoint_msht, r07_ypoint_msht_ratio;
  TGraph *t_r08jet_msht = new TGraph(); TGraph *t_r08jet_msht_ratio = new TGraph(); std::vector<double> r08_ypoint_msht, r08_ypoint_msht_ratio;
  std::vector<double> all_xpoint_new;
  double drop_eta1, drop_eta2;
  double pt_voge, xsection_voge, xsection_low_voge, xsection_high_voge, xsection_nnpdf, xsection_msht;
  int count_new_r02 = 0, count_new_r03 = 0, count_new_r04 = 0, count_new_r05 = 0, count_new_r06 = 0, count_new_r07 = 0, count_new_r08 = 0;
  for (int line = 0; line < 259; line ++) {
    infile_voge >> radius >> drop_eta1 >> drop_eta2 >> pt_voge >> xsection_voge >> xsection_low_voge >> xsection_high_voge;
    infile_nnpdf >> radius >> drop_eta1 >> drop_eta2 >> pt_voge >> xsection_nnpdf;
    infile_msht >> radius >> drop_eta1 >> drop_eta2 >> pt_voge >> xsection_msht;

    if (radius == 0.2) {
      t_r02jet_voge->SetPoint(count_new_r02, pt_voge, xsection_voge); t_r02jet_voge_low->SetPoint(count_new_r02, pt_voge, xsection_low_voge); t_r02jet_voge_high->SetPoint(count_new_r02, pt_voge, xsection_high_voge);
      t_r02jet_nnpdf->SetPoint(count_new_r02, pt_voge, xsection_nnpdf); t_r02jet_nnpdf_ratio->SetPoint(count_new_r02, pt_voge, xsection_nnpdf/xsection_voge);
      t_r02jet_msht->SetPoint(count_new_r02, pt_voge, xsection_msht); t_r02jet_msht_ratio->SetPoint(count_new_r02, pt_voge, xsection_msht/xsection_voge);
      r02_ypoint_voge.push_back(xsection_voge); r02_ypoint_voge_low.push_back(xsection_low_voge); r02_ypoint_voge_high.push_back(xsection_high_voge); r02_ypoint_voge_ratiolow.push_back(xsection_low_voge/xsection_voge); r02_ypoint_voge_ratiohigh.push_back(xsection_high_voge/xsection_voge);
      r02_ypoint_nnpdf.push_back(xsection_nnpdf); r02_ypoint_nnpdf_ratio.push_back(xsection_nnpdf/xsection_voge);
      r02_ypoint_msht.push_back(xsection_msht); r02_ypoint_msht_ratio.push_back(xsection_msht/xsection_voge);
      all_xpoint_new.push_back(pt_voge);
      count_new_r02++;
    } else if (radius == 0.3) {
      t_r03jet_voge->SetPoint(count_new_r03, pt_voge, xsection_voge); t_r03jet_voge_low->SetPoint(count_new_r03, pt_voge, xsection_low_voge); t_r03jet_voge_high->SetPoint(count_new_r03, pt_voge, xsection_high_voge);
      t_r03jet_nnpdf->SetPoint(count_new_r03, pt_voge, xsection_nnpdf); t_r03jet_nnpdf_ratio->SetPoint(count_new_r03, pt_voge, xsection_nnpdf/xsection_voge);
      t_r03jet_msht->SetPoint(count_new_r03, pt_voge, xsection_msht); t_r03jet_msht_ratio->SetPoint(count_new_r03, pt_voge, xsection_msht/xsection_voge);
      r03_ypoint_voge.push_back(xsection_voge); r03_ypoint_voge_low.push_back(xsection_low_voge); r03_ypoint_voge_high.push_back(xsection_high_voge); r03_ypoint_voge_ratiolow.push_back(xsection_low_voge/xsection_voge); r03_ypoint_voge_ratiohigh.push_back(xsection_high_voge/xsection_voge);
      r03_ypoint_nnpdf.push_back(xsection_nnpdf); r03_ypoint_nnpdf_ratio.push_back(xsection_nnpdf/xsection_voge);
      r03_ypoint_msht.push_back(xsection_msht); r03_ypoint_msht_ratio.push_back(xsection_msht/xsection_voge);
      count_new_r03++;
    } else if (radius == 0.4) {
      t_r04jet_voge->SetPoint(count_new_r04, pt_voge, xsection_voge); t_r04jet_voge_low->SetPoint(count_new_r04, pt_voge, xsection_low_voge); t_r04jet_voge_high->SetPoint(count_new_r04, pt_voge, xsection_high_voge);
      t_r04jet_nnpdf->SetPoint(count_new_r04, pt_voge, xsection_nnpdf); t_r04jet_nnpdf_ratio->SetPoint(count_new_r04, pt_voge, xsection_nnpdf/xsection_voge);
      t_r04jet_msht->SetPoint(count_new_r04, pt_voge, xsection_msht); t_r04jet_msht_ratio->SetPoint(count_new_r04, pt_voge, xsection_msht/xsection_voge);
      r04_ypoint_voge.push_back(xsection_voge); r04_ypoint_voge_low.push_back(xsection_low_voge); r04_ypoint_voge_high.push_back(xsection_high_voge); r04_ypoint_voge_ratiolow.push_back(xsection_low_voge/xsection_voge); r04_ypoint_voge_ratiohigh.push_back(xsection_high_voge/xsection_voge);
      r04_ypoint_nnpdf.push_back(xsection_nnpdf); r04_ypoint_nnpdf_ratio.push_back(xsection_nnpdf/xsection_voge);
      r04_ypoint_msht.push_back(xsection_msht); r04_ypoint_msht_ratio.push_back(xsection_msht/xsection_voge);
      count_new_r04++;
    } else if (radius == 0.5) {
      t_r05jet_voge->SetPoint(count_new_r05, pt_voge, xsection_voge); t_r05jet_voge_low->SetPoint(count_new_r05, pt_voge, xsection_low_voge); t_r05jet_voge_high->SetPoint(count_new_r05, pt_voge, xsection_high_voge);
      t_r05jet_nnpdf->SetPoint(count_new_r05, pt_voge, xsection_nnpdf); t_r05jet_nnpdf_ratio->SetPoint(count_new_r05, pt_voge, xsection_nnpdf/xsection_voge);
      t_r05jet_msht->SetPoint(count_new_r05, pt_voge, xsection_msht); t_r05jet_msht_ratio->SetPoint(count_new_r05, pt_voge, xsection_msht/xsection_voge);
      r05_ypoint_voge.push_back(xsection_voge); r05_ypoint_voge_low.push_back(xsection_low_voge); r05_ypoint_voge_high.push_back(xsection_high_voge); r05_ypoint_voge_ratiolow.push_back(xsection_low_voge/xsection_voge); r05_ypoint_voge_ratiohigh.push_back(xsection_high_voge/xsection_voge);
      r05_ypoint_nnpdf.push_back(xsection_nnpdf); r05_ypoint_nnpdf_ratio.push_back(xsection_nnpdf/xsection_voge);
      r05_ypoint_msht.push_back(xsection_msht); r05_ypoint_msht_ratio.push_back(xsection_msht/xsection_voge);
      count_new_r05++;
    } else if (radius == 0.6) {
      t_r06jet_voge->SetPoint(count_new_r06, pt_voge, xsection_voge); t_r06jet_voge_low->SetPoint(count_new_r06, pt_voge, xsection_low_voge); t_r06jet_voge_high->SetPoint(count_new_r06, pt_voge, xsection_high_voge);
      t_r06jet_nnpdf->SetPoint(count_new_r06, pt_voge, xsection_nnpdf); t_r06jet_nnpdf_ratio->SetPoint(count_new_r06, pt_voge, xsection_nnpdf/xsection_voge);
      t_r06jet_msht->SetPoint(count_new_r06, pt_voge, xsection_msht); t_r06jet_msht_ratio->SetPoint(count_new_r06, pt_voge, xsection_msht/xsection_voge);
      r06_ypoint_voge.push_back(xsection_voge); r06_ypoint_voge_low.push_back(xsection_low_voge); r06_ypoint_voge_high.push_back(xsection_high_voge); r06_ypoint_voge_ratiolow.push_back(xsection_low_voge/xsection_voge); r06_ypoint_voge_ratiohigh.push_back(xsection_high_voge/xsection_voge);
      r06_ypoint_nnpdf.push_back(xsection_nnpdf); r06_ypoint_nnpdf_ratio.push_back(xsection_nnpdf/xsection_voge);
      r06_ypoint_msht.push_back(xsection_msht); r06_ypoint_msht_ratio.push_back(xsection_msht/xsection_voge);
      count_new_r06++;
    } else if (radius == 0.7) {
      t_r07jet_voge->SetPoint(count_new_r07, pt_voge, xsection_voge); t_r07jet_voge_low->SetPoint(count_new_r07, pt_voge, xsection_low_voge); t_r07jet_voge_high->SetPoint(count_new_r07, pt_voge, xsection_high_voge);
      t_r07jet_nnpdf->SetPoint(count_new_r07, pt_voge, xsection_nnpdf); t_r07jet_nnpdf_ratio->SetPoint(count_new_r07, pt_voge, xsection_nnpdf/xsection_voge);
      t_r07jet_msht->SetPoint(count_new_r07, pt_voge, xsection_msht); t_r07jet_msht_ratio->SetPoint(count_new_r07, pt_voge, xsection_msht/xsection_voge);
      r07_ypoint_voge.push_back(xsection_voge); r07_ypoint_voge_low.push_back(xsection_low_voge); r07_ypoint_voge_high.push_back(xsection_high_voge); r07_ypoint_voge_ratiolow.push_back(xsection_low_voge/xsection_voge); r07_ypoint_voge_ratiohigh.push_back(xsection_high_voge/xsection_voge);
      r07_ypoint_nnpdf.push_back(xsection_nnpdf); r07_ypoint_nnpdf_ratio.push_back(xsection_nnpdf/xsection_voge);
      r07_ypoint_msht.push_back(xsection_msht); r07_ypoint_msht_ratio.push_back(xsection_msht/xsection_voge);
      count_new_r07++;
    } else if (radius == 0.8) {
      t_r08jet_voge->SetPoint(count_new_r08, pt_voge, xsection_voge); t_r08jet_voge_low->SetPoint(count_new_r08, pt_voge, xsection_low_voge); t_r08jet_voge_high->SetPoint(count_new_r08, pt_voge, xsection_high_voge);
      t_r08jet_nnpdf->SetPoint(count_new_r08, pt_voge, xsection_nnpdf); t_r08jet_nnpdf_ratio->SetPoint(count_new_r08, pt_voge, xsection_nnpdf/xsection_voge);
      t_r08jet_msht->SetPoint(count_new_r08, pt_voge, xsection_msht); t_r08jet_msht_ratio->SetPoint(count_new_r08, pt_voge, xsection_msht/xsection_voge);
      r08_ypoint_voge.push_back(xsection_voge); r08_ypoint_voge_low.push_back(xsection_low_voge); r08_ypoint_voge_high.push_back(xsection_high_voge); r08_ypoint_voge_ratiolow.push_back(xsection_low_voge/xsection_voge); r08_ypoint_voge_ratiohigh.push_back(xsection_high_voge/xsection_voge);
      r08_ypoint_nnpdf.push_back(xsection_nnpdf); r08_ypoint_nnpdf_ratio.push_back(xsection_nnpdf/xsection_voge);
      r08_ypoint_msht.push_back(xsection_msht); r08_ypoint_msht_ratio.push_back(xsection_msht/xsection_voge);
      count_new_r08++;
    } else {
      cout << "Error: radius value does not match expected values in new vogelsang: " << radius << " at line " << line << endl;
    }
  }
  if (count_new_r02 != count_new_r03 || count_new_r02 != count_new_r04 || count_new_r02 != count_new_r05 || 
      count_new_r02 != count_new_r06 || count_new_r02 != count_new_r07 || count_new_r02 != count_new_r08) {
    cout << "Error: different number of pT points per radius in Vogelsang file!" << endl;
    return;
  }

  TSpline3 *spline_r02_voge = new TSpline3("spline_r02_voge", t_r02jet_voge); TSpline3 *spline_r02_voge_high = new TSpline3("spline_r02_voge_high", t_r02jet_voge_high); TSpline3 *spline_r02_voge_low = new TSpline3("spline_r02_voge_low", t_r02jet_voge_low);
  TSpline3 *spline_r03_voge = new TSpline3("spline_r03_voge", t_r03jet_voge); TSpline3 *spline_r03_voge_high = new TSpline3("spline_r03_voge_high", t_r03jet_voge_high); TSpline3 *spline_r03_voge_low = new TSpline3("spline_r03_voge_low", t_r03jet_voge_low);
  TSpline3 *spline_r04_voge = new TSpline3("spline_r04_voge", t_r04jet_voge); TSpline3 *spline_r04_voge_high = new TSpline3("spline_r04_voge_high", t_r04jet_voge_high); TSpline3 *spline_r04_voge_low = new TSpline3("spline_r04_voge_low", t_r04jet_voge_low);
  TSpline3 *spline_r05_voge = new TSpline3("spline_r05_voge", t_r05jet_voge); TSpline3 *spline_r05_voge_high = new TSpline3("spline_r05_voge_high", t_r05jet_voge_high); TSpline3 *spline_r05_voge_low = new TSpline3("spline_r05_voge_low", t_r05jet_voge_low);
  TSpline3 *spline_r06_voge = new TSpline3("spline_r06_voge", t_r06jet_voge); TSpline3 *spline_r06_voge_high = new TSpline3("spline_r06_voge_high", t_r06jet_voge_high); TSpline3 *spline_r06_voge_low = new TSpline3("spline_r06_voge_low", t_r06jet_voge_low);
  TSpline3 *spline_r07_voge = new TSpline3("spline_r07_voge", t_r07jet_voge); TSpline3 *spline_r07_voge_high = new TSpline3("spline_r07_voge_high", t_r07jet_voge_high); TSpline3 *spline_r07_voge_low = new TSpline3("spline_r07_voge_low", t_r07jet_voge_low);
  TSpline3 *spline_r08_voge = new TSpline3("spline_r08_voge", t_r08jet_voge); TSpline3 *spline_r08_voge_high = new TSpline3("spline_r08_voge_high", t_r08jet_voge_high); TSpline3 *spline_r08_voge_low = new TSpline3("spline_r08_voge_low", t_r08jet_voge_low);
  spline_r02_voge->SetLineColor(kRed+1); spline_r02_voge->SetLineWidth(2); spline_r02_voge_high->SetLineColor(kRed+1); spline_r02_voge_high->SetLineWidth(2); spline_r02_voge_low->SetLineColor(kRed+1); spline_r02_voge_low->SetLineWidth(2);
  spline_r03_voge->SetLineColor(kRed+1); spline_r03_voge->SetLineWidth(2); spline_r03_voge_high->SetLineColor(kRed+1); spline_r03_voge_high->SetLineWidth(2); spline_r03_voge_low->SetLineColor(kRed+1); spline_r03_voge_low->SetLineWidth(2);
  spline_r04_voge->SetLineColor(kRed+1); spline_r04_voge->SetLineWidth(2); spline_r04_voge_high->SetLineColor(kRed+1); spline_r04_voge_high->SetLineWidth(2); spline_r04_voge_low->SetLineColor(kRed+1); spline_r04_voge_low->SetLineWidth(2);
  spline_r05_voge->SetLineColor(kRed+1); spline_r05_voge->SetLineWidth(2); spline_r05_voge_high->SetLineColor(kRed+1); spline_r05_voge_high->SetLineWidth(2); spline_r05_voge_low->SetLineColor(kRed+1); spline_r05_voge_low->SetLineWidth(2);
  spline_r06_voge->SetLineColor(kRed+1); spline_r06_voge->SetLineWidth(2); spline_r06_voge_high->SetLineColor(kRed+1); spline_r06_voge_high->SetLineWidth(2); spline_r06_voge_low->SetLineColor(kRed+1); spline_r06_voge_low->SetLineWidth(2);
  spline_r07_voge->SetLineColor(kRed+1); spline_r07_voge->SetLineWidth(2); spline_r07_voge_high->SetLineColor(kRed+1); spline_r07_voge_high->SetLineWidth(2); spline_r07_voge_low->SetLineColor(kRed+1); spline_r07_voge_low->SetLineWidth(2);
  spline_r08_voge->SetLineColor(kRed+1); spline_r08_voge->SetLineWidth(2); spline_r08_voge_high->SetLineColor(kRed+1); spline_r08_voge_high->SetLineWidth(2); spline_r08_voge_low->SetLineColor(kRed+1); spline_r08_voge_low->SetLineWidth(2);

  TSpline3 *spline_r02_nnpdf = new TSpline3("spline_r02_nnpdf", t_r02jet_nnpdf); TSpline3 *spline_r02_nnpdf_ratio = new TSpline3("spline_r02_nnpdf_ratio", t_r02jet_nnpdf_ratio); spline_r02_nnpdf->SetLineColor(kOrange-1); spline_r02_nnpdf->SetLineWidth(2); spline_r02_nnpdf_ratio->SetLineColor(kOrange-1); spline_r02_nnpdf_ratio->SetLineWidth(2);
  TSpline3 *spline_r03_nnpdf = new TSpline3("spline_r03_nnpdf", t_r03jet_nnpdf); TSpline3 *spline_r03_nnpdf_ratio = new TSpline3("spline_r03_nnpdf_ratio", t_r03jet_nnpdf_ratio); spline_r03_nnpdf->SetLineColor(kOrange-1); spline_r03_nnpdf->SetLineWidth(2); spline_r03_nnpdf_ratio->SetLineColor(kOrange-1); spline_r03_nnpdf_ratio->SetLineWidth(2);
  TSpline3 *spline_r04_nnpdf = new TSpline3("spline_r04_nnpdf", t_r04jet_nnpdf); TSpline3 *spline_r04_nnpdf_ratio = new TSpline3("spline_r04_nnpdf_ratio", t_r04jet_nnpdf_ratio); spline_r04_nnpdf->SetLineColor(kOrange-1); spline_r04_nnpdf->SetLineWidth(2); spline_r04_nnpdf_ratio->SetLineColor(kOrange-1); spline_r04_nnpdf_ratio->SetLineWidth(2);
  TSpline3 *spline_r05_nnpdf = new TSpline3("spline_r05_nnpdf", t_r05jet_nnpdf); TSpline3 *spline_r05_nnpdf_ratio = new TSpline3("spline_r05_nnpdf_ratio", t_r05jet_nnpdf_ratio); spline_r05_nnpdf->SetLineColor(kOrange-1); spline_r05_nnpdf->SetLineWidth(2); spline_r05_nnpdf_ratio->SetLineColor(kOrange-1); spline_r05_nnpdf_ratio->SetLineWidth(2);
  TSpline3 *spline_r06_nnpdf = new TSpline3("spline_r06_nnpdf", t_r06jet_nnpdf); TSpline3 *spline_r06_nnpdf_ratio = new TSpline3("spline_r06_nnpdf_ratio", t_r06jet_nnpdf_ratio); spline_r06_nnpdf->SetLineColor(kOrange-1); spline_r06_nnpdf->SetLineWidth(2); spline_r06_nnpdf_ratio->SetLineColor(kOrange-1); spline_r06_nnpdf_ratio->SetLineWidth(2);
  TSpline3 *spline_r07_nnpdf = new TSpline3("spline_r07_nnpdf", t_r07jet_nnpdf); TSpline3 *spline_r07_nnpdf_ratio = new TSpline3("spline_r07_nnpdf_ratio", t_r07jet_nnpdf_ratio); spline_r07_nnpdf->SetLineColor(kOrange-1); spline_r07_nnpdf->SetLineWidth(2); spline_r07_nnpdf_ratio->SetLineColor(kOrange-1); spline_r07_nnpdf_ratio->SetLineWidth(2);
  TSpline3 *spline_r08_nnpdf = new TSpline3("spline_r08_nnpdf", t_r08jet_nnpdf); TSpline3 *spline_r08_nnpdf_ratio = new TSpline3("spline_r08_nnpdf_ratio", t_r08jet_nnpdf_ratio); spline_r08_nnpdf->SetLineColor(kOrange-1); spline_r08_nnpdf->SetLineWidth(2); spline_r08_nnpdf_ratio->SetLineColor(kOrange-1); spline_r08_nnpdf_ratio->SetLineWidth(2);

  TSpline3 *spline_r02_msht = new TSpline3("spline_r02_msht", t_r02jet_msht); TSpline3 *spline_r02_msht_ratio = new TSpline3("spline_r02_msht_ratio", t_r02jet_msht_ratio); spline_r02_msht->SetLineColor(kMagenta+2); spline_r02_msht->SetLineWidth(2); spline_r02_msht_ratio->SetLineColor(kMagenta+2); spline_r02_msht_ratio->SetLineWidth(2);
  TSpline3 *spline_r03_msht = new TSpline3("spline_r03_msht", t_r03jet_msht); TSpline3 *spline_r03_msht_ratio = new TSpline3("spline_r03_msht_ratio", t_r03jet_msht_ratio); spline_r03_msht->SetLineColor(kMagenta+2); spline_r03_msht->SetLineWidth(2); spline_r03_msht_ratio->SetLineColor(kMagenta+2); spline_r03_msht_ratio->SetLineWidth(2);
  TSpline3 *spline_r04_msht = new TSpline3("spline_r04_msht", t_r04jet_msht); TSpline3 *spline_r04_msht_ratio = new TSpline3("spline_r04_msht_ratio", t_r04jet_msht_ratio); spline_r04_msht->SetLineColor(kMagenta+2); spline_r04_msht->SetLineWidth(2); spline_r04_msht_ratio->SetLineColor(kMagenta+2); spline_r04_msht_ratio->SetLineWidth(2);
  TSpline3 *spline_r05_msht = new TSpline3("spline_r05_msht", t_r05jet_msht); TSpline3 *spline_r05_msht_ratio = new TSpline3("spline_r05_msht_ratio", t_r05jet_msht_ratio); spline_r05_msht->SetLineColor(kMagenta+2); spline_r05_msht->SetLineWidth(2); spline_r05_msht_ratio->SetLineColor(kMagenta+2); spline_r05_msht_ratio->SetLineWidth(2);
  TSpline3 *spline_r06_msht = new TSpline3("spline_r06_msht", t_r06jet_msht); TSpline3 *spline_r06_msht_ratio = new TSpline3("spline_r06_msht_ratio", t_r06jet_msht_ratio); spline_r06_msht->SetLineColor(kMagenta+2); spline_r06_msht->SetLineWidth(2); spline_r06_msht_ratio->SetLineColor(kMagenta+2); spline_r06_msht_ratio->SetLineWidth(2);
  TSpline3 *spline_r07_msht = new TSpline3("spline_r07_msht", t_r07jet_msht); TSpline3 *spline_r07_msht_ratio = new TSpline3("spline_r07_msht_ratio", t_r07jet_msht_ratio); spline_r07_msht->SetLineColor(kMagenta+2); spline_r07_msht->SetLineWidth(2); spline_r07_msht_ratio->SetLineColor(kMagenta+2); spline_r07_msht_ratio->SetLineWidth(2);
  TSpline3 *spline_r08_msht = new TSpline3("spline_r08_msht", t_r08jet_msht); TSpline3 *spline_r08_msht_ratio = new TSpline3("spline_r08_msht_ratio", t_r08jet_msht_ratio); spline_r08_msht->SetLineColor(kMagenta+2); spline_r08_msht->SetLineWidth(2); spline_r08_msht_ratio->SetLineColor(kMagenta+2); spline_r08_msht_ratio->SetLineWidth(2);

  std::vector <double> all_xpoint_new_shadow;
  std::vector <double> r02_ypoint_voge_shadow, r03_ypoint_voge_shadow, r04_ypoint_voge_shadow, r05_ypoint_voge_shadow, r06_ypoint_voge_shadow, r07_ypoint_voge_shadow, r08_ypoint_voge_shadow;
  std::vector <double> r02_ypoint_voge_ratioshadow, r03_ypoint_voge_ratioshadow, r04_ypoint_voge_ratioshadow, r05_ypoint_voge_ratioshadow, r06_ypoint_voge_ratioshadow, r07_ypoint_voge_ratioshadow, r08_ypoint_voge_ratioshadow;
  for (int i = 0; i < all_xpoint_new.size(); i++) {
    all_xpoint_new_shadow.push_back(all_xpoint_new[i]);
    r02_ypoint_voge_shadow.push_back(r02_ypoint_voge_low[i]);
    r03_ypoint_voge_shadow.push_back(r03_ypoint_voge_low[i]);
    r04_ypoint_voge_shadow.push_back(r04_ypoint_voge_low[i]);
    r05_ypoint_voge_shadow.push_back(r05_ypoint_voge_low[i]);
    r06_ypoint_voge_shadow.push_back(r06_ypoint_voge_low[i]);
    r07_ypoint_voge_shadow.push_back(r07_ypoint_voge_low[i]);
    r08_ypoint_voge_shadow.push_back(r08_ypoint_voge_low[i]);
    r02_ypoint_voge_ratioshadow.push_back(r02_ypoint_voge_ratiolow[i]);
    r03_ypoint_voge_ratioshadow.push_back(r03_ypoint_voge_ratiolow[i]);
    r04_ypoint_voge_ratioshadow.push_back(r04_ypoint_voge_ratiolow[i]);
    r05_ypoint_voge_ratioshadow.push_back(r05_ypoint_voge_ratiolow[i]);
    r06_ypoint_voge_ratioshadow.push_back(r06_ypoint_voge_ratiolow[i]);
    r07_ypoint_voge_ratioshadow.push_back(r07_ypoint_voge_ratiolow[i]);
    r08_ypoint_voge_ratioshadow.push_back(r08_ypoint_voge_ratiolow[i]);
  }
  for (int i = all_xpoint_new.size(); i > 0; i--) {
    all_xpoint_new_shadow.push_back(all_xpoint_new[i-1]);
    r02_ypoint_voge_shadow.push_back(r02_ypoint_voge_high[i-1]);
    r03_ypoint_voge_shadow.push_back(r03_ypoint_voge_high[i-1]);
    r04_ypoint_voge_shadow.push_back(r04_ypoint_voge_high[i-1]);
    r05_ypoint_voge_shadow.push_back(r05_ypoint_voge_high[i-1]);
    r06_ypoint_voge_shadow.push_back(r06_ypoint_voge_high[i-1]);
    r07_ypoint_voge_shadow.push_back(r07_ypoint_voge_high[i-1]);
    r08_ypoint_voge_shadow.push_back(r08_ypoint_voge_high[i-1]);
    r02_ypoint_voge_ratioshadow.push_back(r02_ypoint_voge_ratiohigh[i-1]);
    r03_ypoint_voge_ratioshadow.push_back(r03_ypoint_voge_ratiohigh[i-1]);
    r04_ypoint_voge_ratioshadow.push_back(r04_ypoint_voge_ratiohigh[i-1]);
    r05_ypoint_voge_ratioshadow.push_back(r05_ypoint_voge_ratiohigh[i-1]);
    r06_ypoint_voge_ratioshadow.push_back(r06_ypoint_voge_ratiohigh[i-1]);
    r07_ypoint_voge_ratioshadow.push_back(r07_ypoint_voge_ratiohigh[i-1]);
    r08_ypoint_voge_ratioshadow.push_back(r08_ypoint_voge_ratiohigh[i-1]);
  }
  TGraph *t_r02jet_voge_shadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r02_ypoint_voge_shadow[0]); t_r02jet_voge_shadow->SetFillColorAlpha(kRed-4, 0.5); t_r02jet_voge_shadow->SetLineColor(kRed-4);
  TGraph *t_r03jet_voge_shadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r03_ypoint_voge_shadow[0]); t_r03jet_voge_shadow->SetFillColorAlpha(kRed-4, 0.5); t_r03jet_voge_shadow->SetLineColor(kRed-4);
  TGraph *t_r04jet_voge_shadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r04_ypoint_voge_shadow[0]); t_r04jet_voge_shadow->SetFillColorAlpha(kRed-4, 0.5); t_r04jet_voge_shadow->SetLineColor(kRed-4);
  TGraph *t_r05jet_voge_shadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r05_ypoint_voge_shadow[0]); t_r05jet_voge_shadow->SetFillColorAlpha(kRed-4, 0.5); t_r05jet_voge_shadow->SetLineColor(kRed-4);
  TGraph *t_r06jet_voge_shadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r06_ypoint_voge_shadow[0]); t_r06jet_voge_shadow->SetFillColorAlpha(kRed-4, 0.5); t_r06jet_voge_shadow->SetLineColor(kRed-4);
  TGraph *t_r07jet_voge_shadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r07_ypoint_voge_shadow[0]); t_r07jet_voge_shadow->SetFillColorAlpha(kRed-4, 0.5); t_r07jet_voge_shadow->SetLineColor(kRed-4);
  TGraph *t_r08jet_voge_shadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r08_ypoint_voge_shadow[0]); t_r08jet_voge_shadow->SetFillColorAlpha(kRed-4, 0.5); t_r08jet_voge_shadow->SetLineColor(kRed-4);
  TGraph *t_r02jet_voge_ratioshadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r02_ypoint_voge_ratioshadow[0]); t_r02jet_voge_ratioshadow->SetFillColorAlpha(kRed-4, 0.5);
  TGraph *t_r03jet_voge_ratioshadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r03_ypoint_voge_ratioshadow[0]); t_r03jet_voge_ratioshadow->SetFillColorAlpha(kRed-4, 0.5);
  TGraph *t_r04jet_voge_ratioshadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r04_ypoint_voge_ratioshadow[0]); t_r04jet_voge_ratioshadow->SetFillColorAlpha(kRed-4, 0.5);
  TGraph *t_r05jet_voge_ratioshadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r05_ypoint_voge_ratioshadow[0]); t_r05jet_voge_ratioshadow->SetFillColorAlpha(kRed-4, 0.5);
  TGraph *t_r06jet_voge_ratioshadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r06_ypoint_voge_ratioshadow[0]); t_r06jet_voge_ratioshadow->SetFillColorAlpha(kRed-4, 0.5);
  TGraph *t_r07jet_voge_ratioshadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r07_ypoint_voge_ratioshadow[0]); t_r07jet_voge_ratioshadow->SetFillColorAlpha(kRed-4, 0.5);
  TGraph *t_r08jet_voge_ratioshadow = new TGraph(2*all_xpoint_new.size(), &all_xpoint_new_shadow[0], &r08_ypoint_voge_ratioshadow[0]); t_r08jet_voge_ratioshadow->SetFillColorAlpha(kRed-4, 0.5);
  infile_voge.close(); infile_nnpdf.close(); infile_msht.close();

  // PHENIX result
  //double phenixpointx[] = {8.45355, 9.45319, 10.8252, 13.0123, 15.7055, 18.739, 22.0908, 26.2576, 31.1771, 37.4739};
  //const static int phenixnpoint = sizeof(phenixpointx) / sizeof(phenixpointx[0]);
  //double phenixpointxerrdown[phenixnpoint] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
  //double phenixpointxerrup[phenixnpoint] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
  //double phenixpointy[phenixnpoint] = {5.86972e+05, 2.99829e+05, 1.3335e+05, 4.07571e+04, 1.24647e+04, 4.28294e+03, 1.39962e+03, 4.11469e+02, 1.01631e+02, 2.49519e+01};
  //double phenixpointysysdown[phenixnpoint] = {2.32513e+05, 1.08448e+05, 4.80741e+04, 1.44942e+04, 4.19428e+03, 1.21361e+03, 3.87209e+02, 9.41715e+01, 1.20866e+01, 1.06009};
  //double phenixpointysysup[phenixnpoint] = {1.90074e+05, 9.64104e+04, 4.43763e+04, 1.42428e+04, 4.28854e+03, 1.30758e+03, 4.68572e+02, 1.15733e+02, 1.5672e+01, 2.62781};
  //double phenixpointystat[phenixnpoint] = {3.5247e+03, 2.16461e+03, 1.15278e+03, 4.6706e+02, 1.86848e+02, 1.04529e+02, 4.38652e+01, 2.17658e+01, 9.3723, 4.46448};
  double phenixpointx[] = {13.01, 15.7, 18.74, 22.08, 26.26, 31.19, 37.47};
  const static int phenixnpoint = sizeof(phenixpointx) / sizeof(phenixpointx[0]);
  double phenixpointxerrdown[phenixnpoint] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
  double phenixpointxerrup[phenixnpoint] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
  double phenixpointy[] = {7.78e+04, 2.46e+04, 8.8e+03, 2.95e+03, 8e+02, 1.9e+02, 4.1e+01};
  double phenixpointysysdown[] = {2.27e+04, 6.9e+03, 3e+03, 1.21e+03, 3.8e+02, 8e+01, 2.1e+01};
  double phenixpointysysup[] = {2.3e+04, 7e+03, 3e+03, 1.21e+03, 3.8e+02, 8e+01, 2.1e+01};
  double phenixpointystat[] = {3.2e+03, 1.3e+03, 7e+02, 3e+02, 1.2e+02, 5e+01, 2.1e+01};
  double phenixpointyerrdown[phenixnpoint];
  double phenixpointyerrup[phenixnpoint];
  for (int i = 0; i < phenixnpoint; ++i) {
    phenixpointyerrdown[i] = TMath::Sqrt(phenixpointystat[i]*phenixpointystat[i] + phenixpointysysdown[i]*phenixpointysysdown[i]);
    phenixpointyerrup[i] = TMath::Sqrt(phenixpointystat[i]*phenixpointystat[i] + phenixpointysysup[i]*phenixpointysysup[i]);
  }
  TGraphAsymmErrors* g_phenix_result = new TGraphAsymmErrors(phenixnpoint, phenixpointx, phenixpointy, phenixpointxerrdown, phenixpointxerrup, phenixpointyerrdown, phenixpointyerrup);
  g_phenix_result->SetMarkerStyle(20);
  g_phenix_result->SetMarkerSize(1.5);
  g_phenix_result->SetMarkerColor(kGray+2);
  g_phenix_result->SetLineColor(kGray+2);
  g_phenix_result->SetFillColorAlpha(kGray+1, 0.5);
  ofstream fout_phenix("phenix_result_point.txt");
  for (int ib = 0; ib < phenixnpoint; ++ib) {
    fout_phenix << phenixpointx[ib] << " " << phenixpointxerrdown[ib] << " " << phenixpointxerrup[ib] << " " << phenixpointy[ib] << " " << phenixpointyerrdown[ib] << " " << phenixpointyerrup[ib] << " " << phenixpointystat[ib] << std::endl;
  }
  TGraphErrors* g_phenix_result_point = new TGraphErrors(phenixnpoint, phenixpointx, phenixpointy, phenixpointxerrdown, phenixpointystat);
  g_phenix_result_point->SetMarkerStyle(20);
  g_phenix_result_point->SetMarkerSize(1.5);
  g_phenix_result_point->SetMarkerColor(kGray+2);
  g_phenix_result_point->SetLineColor(kGray+2);
  g_phenix_result_point->SetLineWidth(2);

  // New PHENIX result
  double phenixnewpointx[] = {8.5, 9.5, 11, 13.25, 16, 19, 22.5, 26.75, 31.5, 38};
  const static int phenixnewnpoint = sizeof(phenixnewpointx) / sizeof(phenixnewpointx[0]);
  double phenixnewpointxerrdown[] = {0.5, 0.5, 1, 1.25, 1.5, 1.5, 2, 2.25, 2.5, 4};
  double phenixnewpointxerrup[] = {0.5, 0.5, 1, 1.25, 1.5, 1.5, 2, 2.25, 2.5, 4};
  double phenixnewpointy[] = {5.870e+5, 2.998e+5, 1.334e+5, 4.08e+4, 1.246e+4, 4.28e+3, 1.40e+3, 4.11e+2, 1.02e+2, 2.50e+1};
  double phenixnewpointysysdown_1[] = {1.901e+5, 9.64e+4, 4.44e+4, 1.42e+4, 4.29e+3, 1.31e+3, 4.7e+2, 1.16e+2, 1.6e+1, 2.6e+0};
  double phenixnewpointysysup_1[] = {2.325e+5, 1.084e+5, 4.81e+4, 1.45e+4, 4.19e+3, 1.21e+3, 3.9e+2, 9.4e+1, 1.2e+1, 1.1e+0};
  double phenixnewpointystat[] = {3.5e+3, 2.2e+3, 1.2e+3, 5e+2, 1.9e+2, 1e+2, 4e+1, 2.2e+1, 9e+0, 4.5e+0};
  double phenixnewpointsys_2[] = {4.11e+4, 2.1e+4, 9.3e+3, 2.9e+3, 8.7e+2, 3e+2, 1e+2, 2.9e+1, 7e+0, 1.7e+0};
  double phenixnewpointyerrdown[phenixnewnpoint];
  double phenixnewpointyerrup[phenixnewnpoint];
  double etacorrection[phenixnewnpoint] = {0.961433, 0.952788, 0.955192, 0.951236, 0.935517, 0.920912, 0.893216, 0.864733, 0.825994, 0.767519};
  double phenixnewpointy_etacorrected[phenixnewnpoint];
  double phenixnewpointyerrdown_etacorrected[phenixnewnpoint];
  double phenixnewpointyerrup_etacorrected[phenixnewnpoint];
  double phenixnewpointystat_etacorrected[phenixnewnpoint];
  for (int i = 0; i < phenixnewnpoint; ++i) {
    phenixnewpointyerrdown[i] = TMath::Sqrt(phenixnewpointystat[i]*phenixnewpointystat[i] + phenixnewpointsys_2[i]*phenixnewpointsys_2[i] + phenixnewpointysysdown_1[i]*phenixnewpointysysdown_1[i]);
    phenixnewpointyerrup[i] = TMath::Sqrt(phenixnewpointystat[i]*phenixnewpointystat[i] + phenixnewpointsys_2[i]*phenixnewpointsys_2[i] + phenixnewpointysysup_1[i]*phenixnewpointysysup_1[i]);
    phenixnewpointy_etacorrected[i] = phenixnewpointy[i] * etacorrection[i];
    phenixnewpointyerrdown_etacorrected[i] = phenixnewpointyerrdown[i] * etacorrection[i];
    phenixnewpointyerrup_etacorrected[i] = phenixnewpointyerrup[i] * etacorrection[i];
    phenixnewpointystat_etacorrected[i] = phenixnewpointystat[i] * etacorrection[i];
  }
  TGraphAsymmErrors* g_phenix_new_result = new TGraphAsymmErrors(phenixnewnpoint, phenixnewpointx, phenixnewpointy, phenixnewpointxerrdown, phenixnewpointxerrup, phenixnewpointyerrdown, phenixnewpointyerrup);
  g_phenix_new_result->SetMarkerStyle(20);
  g_phenix_new_result->SetMarkerSize(1.5);
  g_phenix_new_result->SetMarkerColor(kGray+2);
  g_phenix_new_result->SetLineColor(kGray+2);
  g_phenix_new_result->SetFillColorAlpha(kGray+1, 0.5);
  ofstream fout_phenix_new("phenix_new_result_point.txt");
  for (int ib = 0; ib < phenixnewnpoint; ++ib) {
    fout_phenix_new << phenixnewpointx[ib] << " " << phenixnewpointxerrdown[ib] << " " << phenixnewpointxerrup[ib] << " " << phenixnewpointy[ib] << " " << phenixnewpointyerrdown[ib] << " " << phenixnewpointyerrup[ib] << " " << phenixnewpointystat[ib] << std::endl;
  }
  TGraphErrors* g_phenix_new_result_point = new TGraphErrors(phenixnewnpoint, phenixnewpointx, phenixnewpointy, phenixnewpointxerrdown, phenixnewpointystat);
  g_phenix_new_result_point->SetMarkerStyle(20);
  g_phenix_new_result_point->SetMarkerSize(1.5);
  g_phenix_new_result_point->SetMarkerColor(kGray+2);
  g_phenix_new_result_point->SetLineColor(kGray+2);
  g_phenix_new_result_point->SetLineWidth(2);

  TGraphAsymmErrors* g_phenix_new_result_etacorrected = new TGraphAsymmErrors(phenixnewnpoint, phenixnewpointx, phenixnewpointy_etacorrected, phenixnewpointxerrdown, phenixnewpointxerrup, phenixnewpointyerrdown_etacorrected, phenixnewpointyerrup_etacorrected);
  g_phenix_new_result_etacorrected->SetMarkerStyle(20);
  g_phenix_new_result_etacorrected->SetMarkerSize(1.5);
  g_phenix_new_result_etacorrected->SetMarkerColor(kGray+2);
  g_phenix_new_result_etacorrected->SetLineColor(kGray+2);
  g_phenix_new_result_etacorrected->SetFillColorAlpha(kGray+1, 0.5);
  ofstream fout_phenix_new_etacorrected("phenix_new_result_etacorrected_point.txt");
  for (int ib = 0; ib < phenixnewnpoint; ++ib) {
    fout_phenix_new_etacorrected << phenixnewpointx[ib] << " " << phenixnewpointxerrdown[ib] << " " << phenixnewpointxerrup[ib] << " " << phenixnewpointy_etacorrected[ib] << " " << phenixnewpointyerrdown_etacorrected[ib] << " " << phenixnewpointyerrup_etacorrected[ib] << " " << phenixnewpointystat_etacorrected[ib] << std::endl;
  }
  TGraphErrors* g_phenix_new_result_etacorrected_point = new TGraphErrors(phenixnewnpoint, phenixnewpointx, phenixnewpointy_etacorrected, phenixnewpointxerrdown, phenixnewpointystat_etacorrected);
  g_phenix_new_result_etacorrected_point->SetMarkerStyle(20);
  g_phenix_new_result_etacorrected_point->SetMarkerSize(1.5);
  g_phenix_new_result_etacorrected_point->SetMarkerColor(kGray+2);
  g_phenix_new_result_etacorrected_point->SetLineColor(kGray+2);
  g_phenix_new_result_etacorrected_point->SetLineWidth(2);

  f_out_forratio->cd(); g_phenix_new_result->Write("g_phenix_result"); g_phenix_new_result_etacorrected->Write("g_phenix_result_etacorrected");

  // New STAR result
  double starnewpointx1[] = {7.55, 8.95, 10.6, 12.55, 14.85, 17.55, 20.75, 24.55, 29, 34.3, 40.6, 48};
  const static int starnewnpoint1 = sizeof(starnewpointx1) / sizeof(starnewpointx1[0]);
  double starnewpointxerrdown1[] = {0.65, 0.75, 0.9, 1.05, 1.25, 1.45, 1.75, 2.05, 2.4, 2.9, 3.4, 4.0};
  double starnewpointxerrup1[] = {0.65, 0.75, 0.9, 1.05, 1.25, 1.45, 1.75, 2.05, 2.4, 2.9, 3.4, 4.0};
  double starnewpointy1[] = {5.42e+6, 1.95e+6, 7.03e+5, 2.53e+5, 9.51e+4, 3.14e+4, 1.01e+4, 3.06e+3, 8.31e+2, 2.05e+2, 4.15e+1, 5.59e+0};
  double starnewpointysyspercentage1[] = {0.171, 0.254, 0.125, 0.173, 0.107, 0.088, 0.13, 0.099, 0.123, 0.129, 0.162, 0.264};
  double starnewpointystat1[] = {0.004, 0.005, 0.005, 0.005, 0.006, 0.007, 0.01, 0.014, 0.023, 0.038, 0.073, 0.154};
  double starnewpointyerrdown1[starnewnpoint1];
  double starnewpointyerrup1[starnewnpoint1];
  for (int i = 0; i < starnewnpoint1; ++i) {
    starnewpointyerrdown1[i] = TMath::Sqrt(TMath::Power(starnewpointystat1[i]*starnewpointy1[i], 2) + TMath::Power(starnewpointysyspercentage1[i]*starnewpointy1[i], 2));
    starnewpointyerrup1[i] = TMath::Sqrt(TMath::Power(starnewpointystat1[i]*starnewpointy1[i], 2) + TMath::Power(starnewpointysyspercentage1[i]*starnewpointy1[i], 2));
  }
  TGraphAsymmErrors* g_star_new_result1 = new TGraphAsymmErrors(starnewnpoint1, starnewpointx1, starnewpointy1, starnewpointxerrdown1, starnewpointxerrup1, starnewpointyerrdown1, starnewpointyerrup1);
  g_star_new_result1->SetMarkerStyle(20);
  g_star_new_result1->SetMarkerSize(1.5);
  g_star_new_result1->SetMarkerColor(kGreen+2);
  g_star_new_result1->SetLineColor(kGreen+2);
  g_star_new_result1->SetFillColorAlpha(kGreen+1, 0.5);
  ofstream fout_star_new("star_new_result_point.txt");
  for (int ib = 0; ib < starnewnpoint1; ++ib) {
    fout_star_new << starnewpointx1[ib] << " " << starnewpointxerrdown1[ib] << " " << starnewpointxerrup1[ib] << " " << starnewpointy1[ib] << " " << starnewpointyerrdown1[ib] << " " << starnewpointyerrup1[ib] << " " << starnewpointystat1[ib] << std::endl;
  }
  TGraphErrors* g_star_new_result_point1 = new TGraphErrors(starnewnpoint1, starnewpointx1, starnewpointy1, starnewpointxerrdown1, starnewpointystat1);
  g_star_new_result_point1->SetMarkerStyle(20);
  g_star_new_result_point1->SetMarkerSize(1.5);
  g_star_new_result_point1->SetMarkerColor(kGreen+2);
  g_star_new_result_point1->SetLineColor(kGreen+2);
  g_star_new_result_point1->SetLineWidth(2);

  double starnewpointx2[] = {7.55, 8.95, 10.6, 12.55, 14.85, 17.55, 20.75, 24.55, 29, 34.3, 40.6, 48};
  const static int starnewnpoint2 = sizeof(starnewpointx2) / sizeof(starnewpointx2[0]);
  double starnewpointxerrdown2[] = {0.65, 0.75, 0.9, 1.05, 1.25, 1.45, 1.75, 2.05, 2.4, 2.9, 3.4, 4.0};
  double starnewpointxerrup2[] = {0.65, 0.75, 0.9, 1.05, 1.25, 1.45, 1.75, 2.05, 2.4, 2.9, 3.4, 4.0};
  double starnewpointy2[] = {4.84e+6, 1.67e+6, 6.54e+5, 2.14e+5, 7.99e+4, 2.71e+4, 7.96e+3, 2.26e+3, 5.25e+2, 1.27e+2, 1.82e+1, 1.39e+0};
  double starnewpointysyspercentage2[] = {0.191, 0.203, 0.194, 0.166, 0.099, 0.12, 0.098, 0.128, 0.105, 0.152, 0.179, 0.351};
  double starnewpointystat2[] = {0.004, 0.006, 0.006, 0.007, 0.007, 0.01, 0.013, 0.02, 0.035, 0.057, 0.138, 0.435};
  double starnewpointyerrdown2[starnewnpoint2];
  double starnewpointyerrup2[starnewnpoint2];
  for (int i = 0; i < starnewnpoint2; ++i) {
    starnewpointyerrdown2[i] = TMath::Sqrt(TMath::Power(starnewpointystat2[i]*starnewpointy2[i], 2) + TMath::Power(starnewpointysyspercentage2[i]*starnewpointy2[i], 2));
    starnewpointyerrup2[i] = TMath::Sqrt(TMath::Power(starnewpointystat2[i]*starnewpointy2[i], 2) + TMath::Power(starnewpointysyspercentage2[i]*starnewpointy2[i], 2));
  }
  TGraphAsymmErrors* g_star_new_result2 = new TGraphAsymmErrors(starnewnpoint2, starnewpointx2, starnewpointy2, starnewpointxerrdown2, starnewpointxerrup2, starnewpointyerrdown2, starnewpointyerrup2);
  g_star_new_result2->SetMarkerStyle(20);
  g_star_new_result2->SetMarkerSize(1.5);
  g_star_new_result2->SetMarkerColor(kGreen-3);
  g_star_new_result2->SetLineColor(kGreen-3);
  g_star_new_result2->SetFillColorAlpha(kGreen-9, 0.5);
  for (int ib = 0; ib < starnewnpoint2; ++ib) {
    fout_star_new << starnewpointx2[ib] << " " << starnewpointxerrdown2[ib] << " " << starnewpointxerrup2[ib] << " " << starnewpointy2[ib] << " " << starnewpointyerrdown2[ib] << " " << starnewpointyerrup2[ib] << " " << starnewpointystat2[ib] << std::endl;
  }
  TGraphErrors* g_star_new_result_point2 = new TGraphErrors(starnewnpoint2, starnewpointx2, starnewpointy2, starnewpointxerrdown2, starnewpointystat2);
  g_star_new_result_point2->SetMarkerStyle(20);
  g_star_new_result_point2->SetMarkerSize(1.5);
  g_star_new_result_point2->SetMarkerColor(kGreen-3);
  g_star_new_result_point2->SetLineColor(kGreen-3);
  g_star_new_result_point2->SetLineWidth(2);

  f_out_forratio->cd(); g_star_new_result1->Write("g_star_result_eta_0-0p5"); g_star_new_result2->Write("g_star_result_eta_0p5-0p9");

  // sPHENIX all ratio
  double g_r02_sphenix_ratio_x[calibnpt], g_r02_sphenix_ratio_xerrdown[calibnpt], g_r02_sphenix_ratio_xerrup[calibnpt], g_r02_sphenix_ratio_y[calibnpt], g_r02_sphenix_ratio_yerrdown[calibnpt], g_r02_sphenix_ratio_yerrup[calibnpt], g_r02_sphenix_ratio_yerrforpoint[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    double pt_low = calibptbins[ib];
    double pt_high = calibptbins[ib+1];
    double bin_width = pt_high - pt_low;
    g_r02_sphenix_ratio_x[ib] = (pt_low + pt_high) / 2.;
    g_r02_sphenix_ratio_xerrdown[ib] = bin_width / 2.;
    g_r02_sphenix_ratio_xerrup[ib] = bin_width / 2.;
    double voge_y = IntegrateSpline(spline_r02_voge, pt_low, pt_high) / bin_width;
    double voge_y_high = IntegrateSpline(spline_r02_voge_high, pt_low, pt_high) / bin_width;
    double voge_y_low = IntegrateSpline(spline_r02_voge_low, pt_low, pt_high) / bin_width;
    double voge_y_up = voge_y_high - voge_y;
    double voge_y_down = voge_y_low - voge_y;
    g_r02_sphenix_ratio_y[ib] = g_r02_sphenix_all_y[ib] / voge_y;
    g_r02_sphenix_ratio_yerrdown[ib] = g_r02_sphenix_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r02_sphenix_all_yerrdown[ib]/g_r02_sphenix_all_y[ib], 2) + TMath::Power(voge_y_down/voge_y, 2));
    g_r02_sphenix_ratio_yerrup[ib] = g_r02_sphenix_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r02_sphenix_all_yerrup[ib]/g_r02_sphenix_all_y[ib], 2) + TMath::Power(voge_y_up/voge_y, 2));
    g_r02_sphenix_ratio_yerrforpoint[ib] = 0;
  }
  TGraphAsymmErrors* g_r02_sphenix_ratio_all = new TGraphAsymmErrors(calibnpt, g_r02_sphenix_ratio_x, g_r02_sphenix_ratio_y, g_r02_sphenix_ratio_xerrdown, g_r02_sphenix_ratio_xerrup, g_r02_sphenix_ratio_yerrdown, g_r02_sphenix_ratio_yerrup);
  g_r02_sphenix_ratio_all->SetLineWidth(0);  g_r02_sphenix_ratio_all->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_r02_sphenix_ratio_all_point = new TGraphErrors(calibnpt, g_r02_sphenix_ratio_x, g_r02_sphenix_ratio_y, g_r02_sphenix_ratio_xerrdown, g_r02_sphenix_ratio_yerrforpoint);
  g_r02_sphenix_ratio_all_point->SetMarkerStyle(20); g_r02_sphenix_ratio_all_point->SetMarkerSize(1.5); g_r02_sphenix_ratio_all_point->SetMarkerColor(kAzure+2); g_r02_sphenix_ratio_all_point->SetLineColor(kAzure+2); g_r02_sphenix_ratio_all_point->SetLineWidth(2);

  double g_r03_sphenix_ratio_x[calibnpt], g_r03_sphenix_ratio_xerrdown[calibnpt], g_r03_sphenix_ratio_xerrup[calibnpt], g_r03_sphenix_ratio_y[calibnpt], g_r03_sphenix_ratio_yerrdown[calibnpt], g_r03_sphenix_ratio_yerrup[calibnpt], g_r03_sphenix_ratio_yerrforpoint[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    double pt_low = calibptbins[ib];
    double pt_high = calibptbins[ib+1];
    double bin_width = pt_high - pt_low;
    g_r03_sphenix_ratio_x[ib] = (pt_low + pt_high) / 2.;
    g_r03_sphenix_ratio_xerrdown[ib] = bin_width / 2.;
    g_r03_sphenix_ratio_xerrup[ib] = bin_width / 2.;
    double voge_y = IntegrateSpline(spline_r03_voge, pt_low, pt_high) / bin_width;
    double voge_y_high = IntegrateSpline(spline_r03_voge_high, pt_low, pt_high) / bin_width;
    double voge_y_low = IntegrateSpline(spline_r03_voge_low, pt_low, pt_high) / bin_width;
    double voge_y_up = voge_y_high - voge_y;
    double voge_y_down = voge_y_low - voge_y;
    g_r03_sphenix_ratio_y[ib] = g_r03_sphenix_all_y[ib] / voge_y;
    g_r03_sphenix_ratio_yerrdown[ib] = g_r03_sphenix_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r03_sphenix_all_yerrdown[ib]/g_r03_sphenix_all_y[ib], 2) + TMath::Power(voge_y_down/voge_y, 2));
    g_r03_sphenix_ratio_yerrup[ib] = g_r03_sphenix_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r03_sphenix_all_yerrup[ib]/g_r03_sphenix_all_y[ib], 2) + TMath::Power(voge_y_up/voge_y, 2));
    g_r03_sphenix_ratio_yerrforpoint[ib] = 0;
  }
  TGraphAsymmErrors* g_r03_sphenix_ratio_all = new TGraphAsymmErrors(calibnpt, g_r03_sphenix_ratio_x, g_r03_sphenix_ratio_y, g_r03_sphenix_ratio_xerrdown, g_r03_sphenix_ratio_xerrup, g_r03_sphenix_ratio_yerrdown, g_r03_sphenix_ratio_yerrup);
  g_r03_sphenix_ratio_all->SetLineWidth(0);  g_r03_sphenix_ratio_all->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_r03_sphenix_ratio_all_point = new TGraphErrors(calibnpt, g_r03_sphenix_ratio_x, g_r03_sphenix_ratio_y, g_r03_sphenix_ratio_xerrdown, g_r03_sphenix_ratio_yerrforpoint);
  g_r03_sphenix_ratio_all_point->SetMarkerStyle(20); g_r03_sphenix_ratio_all_point->SetMarkerSize(1.5); g_r03_sphenix_ratio_all_point->SetMarkerColor(kAzure+2); g_r03_sphenix_ratio_all_point->SetLineColor(kAzure+2); g_r03_sphenix_ratio_all_point->SetLineWidth(2);

  double g_r04_sphenix_ratio_x[calibnpt], g_r04_sphenix_ratio_xerrdown[calibnpt], g_r04_sphenix_ratio_xerrup[calibnpt], g_r04_sphenix_ratio_y[calibnpt], g_r04_sphenix_ratio_yerrdown[calibnpt], g_r04_sphenix_ratio_yerrup[calibnpt], g_r04_sphenix_ratio_yerrforpoint[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    double pt_low = calibptbins[ib];
    double pt_high = calibptbins[ib+1];
    double bin_width = pt_high - pt_low;
    g_r04_sphenix_ratio_x[ib] = (pt_low + pt_high) / 2.;
    g_r04_sphenix_ratio_xerrdown[ib] = bin_width / 2.;
    g_r04_sphenix_ratio_xerrup[ib] = bin_width / 2.;
    double voge_y = IntegrateSpline(spline_r04_voge, pt_low, pt_high) / bin_width;
    double voge_y_high = IntegrateSpline(spline_r04_voge_high, pt_low, pt_high) / bin_width;
    double voge_y_low = IntegrateSpline(spline_r04_voge_low, pt_low, pt_high) / bin_width;
    double voge_y_up = voge_y_high - voge_y;
    double voge_y_down = voge_y_low - voge_y;
    g_r04_sphenix_ratio_y[ib] = g_r04_sphenix_all_y[ib] / voge_y;
    g_r04_sphenix_ratio_yerrdown[ib] = g_r04_sphenix_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r04_sphenix_all_yerrdown[ib]/g_r04_sphenix_all_y[ib], 2) + TMath::Power(voge_y_down/voge_y, 2));
    g_r04_sphenix_ratio_yerrup[ib] = g_r04_sphenix_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r04_sphenix_all_yerrup[ib]/g_r04_sphenix_all_y[ib], 2) + TMath::Power(voge_y_up/voge_y, 2));
    g_r04_sphenix_ratio_yerrforpoint[ib] = 0;
  }
  TGraphAsymmErrors* g_r04_sphenix_ratio_all = new TGraphAsymmErrors(calibnpt, g_r04_sphenix_ratio_x, g_r04_sphenix_ratio_y, g_r04_sphenix_ratio_xerrdown, g_r04_sphenix_ratio_xerrup, g_r04_sphenix_ratio_yerrdown, g_r04_sphenix_ratio_yerrup);
  g_r04_sphenix_ratio_all->SetLineWidth(0);  g_r04_sphenix_ratio_all->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_r04_sphenix_ratio_all_point = new TGraphErrors(calibnpt, g_r04_sphenix_ratio_x, g_r04_sphenix_ratio_y, g_r04_sphenix_ratio_xerrdown, g_r04_sphenix_ratio_yerrforpoint);
  g_r04_sphenix_ratio_all_point->SetMarkerStyle(20); g_r04_sphenix_ratio_all_point->SetMarkerSize(1.5); g_r04_sphenix_ratio_all_point->SetMarkerColor(kAzure+2); g_r04_sphenix_ratio_all_point->SetLineColor(kAzure+2); g_r04_sphenix_ratio_all_point->SetLineWidth(2);

  double g_r05_sphenix_ratio_x[calibnpt], g_r05_sphenix_ratio_xerrdown[calibnpt], g_r05_sphenix_ratio_xerrup[calibnpt], g_r05_sphenix_ratio_y[calibnpt], g_r05_sphenix_ratio_yerrdown[calibnpt], g_r05_sphenix_ratio_yerrup[calibnpt], g_r05_sphenix_ratio_yerrforpoint[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    double pt_low = calibptbins[ib];
    double pt_high = calibptbins[ib+1];
    double bin_width = pt_high - pt_low;
    g_r05_sphenix_ratio_x[ib] = (pt_low + pt_high) / 2.;
    g_r05_sphenix_ratio_xerrdown[ib] = bin_width / 2.;
    g_r05_sphenix_ratio_xerrup[ib] = bin_width / 2.;
    double voge_y = IntegrateSpline(spline_r05_voge, pt_low, pt_high) / bin_width;
    double voge_y_high = IntegrateSpline(spline_r05_voge_high, pt_low, pt_high) / bin_width;
    double voge_y_low = IntegrateSpline(spline_r05_voge_low, pt_low, pt_high) / bin_width;
    double voge_y_up = voge_y_high - voge_y;
    double voge_y_down = voge_y_low - voge_y;
    g_r05_sphenix_ratio_y[ib] = g_r05_sphenix_all_y[ib] / voge_y;
    g_r05_sphenix_ratio_yerrdown[ib] = g_r05_sphenix_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r05_sphenix_all_yerrdown[ib]/g_r05_sphenix_all_y[ib], 2) + TMath::Power(voge_y_down/voge_y, 2));
    g_r05_sphenix_ratio_yerrup[ib] = g_r05_sphenix_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r05_sphenix_all_yerrup[ib]/g_r05_sphenix_all_y[ib], 2) + TMath::Power(voge_y_up/voge_y, 2));
    g_r05_sphenix_ratio_yerrforpoint[ib] = 0;
  }
  TGraphAsymmErrors* g_r05_sphenix_ratio_all = new TGraphAsymmErrors(calibnpt, g_r05_sphenix_ratio_x, g_r05_sphenix_ratio_y, g_r05_sphenix_ratio_xerrdown, g_r05_sphenix_ratio_xerrup, g_r05_sphenix_ratio_yerrdown, g_r05_sphenix_ratio_yerrup);
  g_r05_sphenix_ratio_all->SetLineWidth(0);  g_r05_sphenix_ratio_all->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_r05_sphenix_ratio_all_point = new TGraphErrors(calibnpt, g_r05_sphenix_ratio_x, g_r05_sphenix_ratio_y, g_r05_sphenix_ratio_xerrdown, g_r05_sphenix_ratio_yerrforpoint);
  g_r05_sphenix_ratio_all_point->SetMarkerStyle(20); g_r05_sphenix_ratio_all_point->SetMarkerSize(1.5); g_r05_sphenix_ratio_all_point->SetMarkerColor(kAzure+2); g_r05_sphenix_ratio_all_point->SetLineColor(kAzure+2); g_r05_sphenix_ratio_all_point->SetLineWidth(2);

  double g_r06_sphenix_ratio_x[calibnpt], g_r06_sphenix_ratio_xerrdown[calibnpt], g_r06_sphenix_ratio_xerrup[calibnpt], g_r06_sphenix_ratio_y[calibnpt], g_r06_sphenix_ratio_yerrdown[calibnpt], g_r06_sphenix_ratio_yerrup[calibnpt], g_r06_sphenix_ratio_yerrforpoint[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    double pt_low = calibptbins[ib];
    double pt_high = calibptbins[ib+1];
    double bin_width = pt_high - pt_low;
    g_r06_sphenix_ratio_x[ib] = (pt_low + pt_high) / 2.;
    g_r06_sphenix_ratio_xerrdown[ib] = bin_width / 2.;
    g_r06_sphenix_ratio_xerrup[ib] = bin_width / 2.;
    double voge_y = IntegrateSpline(spline_r06_voge, pt_low, pt_high) / bin_width;
    double voge_y_high = IntegrateSpline(spline_r06_voge_high, pt_low, pt_high) / bin_width;
    double voge_y_low = IntegrateSpline(spline_r06_voge_low, pt_low, pt_high) / bin_width;
    double voge_y_up = voge_y_high - voge_y;
    double voge_y_down = voge_y_low - voge_y;
    g_r06_sphenix_ratio_y[ib] = g_r06_sphenix_all_y[ib] / voge_y;
    g_r06_sphenix_ratio_yerrdown[ib] = g_r06_sphenix_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r06_sphenix_all_yerrdown[ib]/g_r06_sphenix_all_y[ib], 2) + TMath::Power(voge_y_down/voge_y, 2));
    g_r06_sphenix_ratio_yerrup[ib] = g_r06_sphenix_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r06_sphenix_all_yerrup[ib]/g_r06_sphenix_all_y[ib], 2) + TMath::Power(voge_y_up/voge_y, 2));
    g_r06_sphenix_ratio_yerrforpoint[ib] = 0;
  }
  TGraphAsymmErrors* g_r06_sphenix_ratio_all = new TGraphAsymmErrors(calibnpt, g_r06_sphenix_ratio_x, g_r06_sphenix_ratio_y, g_r06_sphenix_ratio_xerrdown, g_r06_sphenix_ratio_xerrup, g_r06_sphenix_ratio_yerrdown, g_r06_sphenix_ratio_yerrup);
  g_r06_sphenix_ratio_all->SetLineWidth(0);  g_r06_sphenix_ratio_all->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_r06_sphenix_ratio_all_point = new TGraphErrors(calibnpt, g_r06_sphenix_ratio_x, g_r06_sphenix_ratio_y, g_r06_sphenix_ratio_xerrdown, g_r06_sphenix_ratio_yerrforpoint);
  g_r06_sphenix_ratio_all_point->SetMarkerStyle(20); g_r06_sphenix_ratio_all_point->SetMarkerSize(1.5); g_r06_sphenix_ratio_all_point->SetMarkerColor(kAzure+2); g_r06_sphenix_ratio_all_point->SetLineColor(kAzure+2); g_r06_sphenix_ratio_all_point->SetLineWidth(2);

  double g_r08_sphenix_ratio_x[calibnpt], g_r08_sphenix_ratio_xerrdown[calibnpt], g_r08_sphenix_ratio_xerrup[calibnpt], g_r08_sphenix_ratio_y[calibnpt], g_r08_sphenix_ratio_yerrdown[calibnpt], g_r08_sphenix_ratio_yerrup[calibnpt], g_r08_sphenix_ratio_yerrforpoint[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    double pt_low = calibptbins[ib];
    double pt_high = calibptbins[ib+1];
    double bin_width = pt_high - pt_low;
    g_r08_sphenix_ratio_x[ib] = (pt_low + pt_high) / 2.;
    g_r08_sphenix_ratio_xerrdown[ib] = bin_width / 2.;
    g_r08_sphenix_ratio_xerrup[ib] = bin_width / 2.;
    double voge_y = IntegrateSpline(spline_r08_voge, pt_low, pt_high) / bin_width;
    double voge_y_high = IntegrateSpline(spline_r08_voge_high, pt_low, pt_high) / bin_width;
    double voge_y_low = IntegrateSpline(spline_r08_voge_low, pt_low, pt_high) / bin_width;
    double voge_y_up = voge_y_high - voge_y;
    double voge_y_down = voge_y_low - voge_y;
    g_r08_sphenix_ratio_y[ib] = g_r08_sphenix_all_y[ib] / voge_y;
    g_r08_sphenix_ratio_yerrdown[ib] = g_r08_sphenix_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r08_sphenix_all_yerrdown[ib]/g_r08_sphenix_all_y[ib], 2) + TMath::Power(voge_y_down/voge_y, 2));
    g_r08_sphenix_ratio_yerrup[ib] = g_r08_sphenix_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r08_sphenix_all_yerrup[ib]/g_r08_sphenix_all_y[ib], 2) + TMath::Power(voge_y_up/voge_y, 2));
    g_r08_sphenix_ratio_yerrforpoint[ib] = 0;
  }
  TGraphAsymmErrors* g_r08_sphenix_ratio_all = new TGraphAsymmErrors(calibnpt, g_r08_sphenix_ratio_x, g_r08_sphenix_ratio_y, g_r08_sphenix_ratio_xerrdown, g_r08_sphenix_ratio_xerrup, g_r08_sphenix_ratio_yerrdown, g_r08_sphenix_ratio_yerrup);
  g_r08_sphenix_ratio_all->SetLineWidth(0);  g_r08_sphenix_ratio_all->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_r08_sphenix_ratio_all_point = new TGraphErrors(calibnpt, g_r08_sphenix_ratio_x, g_r08_sphenix_ratio_y, g_r08_sphenix_ratio_xerrdown, g_r08_sphenix_ratio_yerrforpoint);
  g_r08_sphenix_ratio_all_point->SetMarkerStyle(20); g_r08_sphenix_ratio_all_point->SetMarkerSize(1.5); g_r08_sphenix_ratio_all_point->SetMarkerColor(kAzure+2); g_r08_sphenix_ratio_all_point->SetLineColor(kAzure+2); g_r08_sphenix_ratio_all_point->SetLineWidth(2);

  // sPHENIX zvertex60 ratio
  double g_sphenix_zvertex60_ratio_x[calibnpt], g_sphenix_zvertex60_ratio_xerrdown[calibnpt], g_sphenix_zvertex60_ratio_xerrup[calibnpt], g_sphenix_zvertex60_ratio_y[calibnpt], g_sphenix_zvertex60_ratio_yerrdown[calibnpt], g_sphenix_zvertex60_ratio_yerrup[calibnpt], g_sphenix_zvertex60_ratio_yerrforpoint[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ib) {
    double pt_low = calibptbins[ib];
    double pt_high = calibptbins[ib + 1];
    double bin_width = pt_high - pt_low;
    g_sphenix_zvertex60_ratio_x[ib] = (pt_low + pt_high) / 2.;
    g_sphenix_zvertex60_ratio_xerrdown[ib] = bin_width / 2.;
    g_sphenix_zvertex60_ratio_xerrup[ib] = bin_width / 2.;
    double voge_y = IntegrateSpline(spline_r06, pt_low, pt_high) / bin_width;
    double voge_y_high = IntegrateSpline(spline_r06_2, pt_low, pt_high) / bin_width;
    double voge_y_low = IntegrateSpline(spline_r06_05, pt_low, pt_high) / bin_width;
    double voge_y_up = voge_y_high - voge_y;
    double voge_y_down = voge_y_low - voge_y;
    g_sphenix_zvertex60_ratio_y[ib] = g_r04_sphenix_zvertex60_y[ib] / voge_y;
    g_sphenix_zvertex60_ratio_yerrdown[ib] = g_sphenix_zvertex60_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r04_sphenix_zvertex60_yerrdown[ib]/g_r04_sphenix_zvertex60_y[ib], 2) + TMath::Power(voge_y_down/voge_y, 2));
    g_sphenix_zvertex60_ratio_yerrup[ib] = g_sphenix_zvertex60_ratio_y[ib] * TMath::Sqrt(TMath::Power(g_r04_sphenix_zvertex60_yerrup[ib]/g_r04_sphenix_zvertex60_y[ib], 2) + TMath::Power(voge_y_up/voge_y, 2));
    g_sphenix_zvertex60_ratio_yerrforpoint[ib] = 0;
  }
  TGraphAsymmErrors* g_sphenix_ratio_zvertex60 = new TGraphAsymmErrors(calibnpt, g_sphenix_zvertex60_ratio_x, g_sphenix_zvertex60_ratio_y, g_sphenix_zvertex60_ratio_xerrdown, g_sphenix_zvertex60_ratio_xerrup, g_sphenix_zvertex60_ratio_yerrdown, g_sphenix_zvertex60_ratio_yerrup);
  g_sphenix_ratio_zvertex60->SetLineWidth(0);  g_sphenix_ratio_zvertex60->SetFillColorAlpha(kAzure + 1, 0.5);
  TGraphErrors* g_sphenix_ratio_zvertex60_point = new TGraphErrors(calibnpt, g_sphenix_zvertex60_ratio_x, g_sphenix_zvertex60_ratio_y, g_sphenix_zvertex60_ratio_xerrdown, g_sphenix_zvertex60_ratio_yerrforpoint);
  g_sphenix_ratio_zvertex60_point->SetMarkerStyle(20); g_sphenix_ratio_zvertex60_point->SetMarkerSize(1.5); g_sphenix_ratio_zvertex60_point->SetMarkerColor(kAzure+2); g_sphenix_ratio_zvertex60_point->SetLineColor(kAzure+2); g_sphenix_ratio_zvertex60_point->SetLineWidth(2);

  ////////////// Plotting //////////////
  TCanvas *can_all_result = new TCanvas("can_all_result", "", 850, 800);
  gStyle->SetPalette(57);
  can_all_result->SetTopMargin(0.03);
  can_all_result->SetLeftMargin(0.15);
  can_all_result->SetBottomMargin(0.12);
  can_all_result->SetRightMargin(0.05);
  can_all_result->SetLogy();
  TH2F *frame_all_result = new TH2F("frame_all_result", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r02_sphenix_result->GetN(), g_r02_sphenix_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r08_sphenix_result->GetN(), g_r08_sphenix_result->GetY())));
  frame_all_result->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_all_result->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_all_result->GetXaxis()->SetTitleOffset(0.98);
  frame_all_result->GetYaxis()->SetTitleOffset(1.15);
  frame_all_result->GetXaxis()->SetLabelSize(0.045);
  frame_all_result->GetYaxis()->SetLabelSize(0.045);
  frame_all_result->GetXaxis()->CenterTitle();
  frame_all_result->GetYaxis()->CenterTitle();
  frame_all_result->Draw();
  g_r02_sphenix_result->SetFillColorAlpha(kYellow+1, 0.5);
  g_r03_sphenix_result->SetFillColorAlpha(kGreen+1, 0.5);
  g_r04_sphenix_result->SetFillColorAlpha(kAzure+1, 0.5);
  g_r05_sphenix_result->SetFillColorAlpha(kBlue-4, 0.5);
  g_r06_sphenix_result->SetFillColorAlpha(kMagenta+1, 0.5);
  g_r08_sphenix_result->SetFillColorAlpha(kRed+1, 0.5);
  g_r02_sphenix_result->Draw("2 same");
  g_r03_sphenix_result->Draw("2 same");
  g_r04_sphenix_result->Draw("2 same");
  g_r05_sphenix_result->Draw("2 same");
  g_r06_sphenix_result->Draw("2 same");
  g_r08_sphenix_result->Draw("2 same");
  myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200 GeV", 0.04);
  TLegend *leg_all_result = new TLegend(0.15, 0.08, 0.9, 0.33);
  leg_all_result->SetBorderSize(0);
  leg_all_result->SetFillStyle(0);
  leg_all_result->SetTextSize(0.03);
  leg_all_result->SetTextFont(42);
  g_r02_sphenix_result->SetMarkerColor(kYellow+2); g_r02_sphenix_result_point->SetMarkerSize(1.5); g_r02_sphenix_result_point->SetMarkerColor(kYellow+2); g_r02_sphenix_result_point->SetLineColor(kYellow+2);
  g_r03_sphenix_result->SetMarkerColor(kGreen+2); g_r03_sphenix_result_point->SetMarkerSize(1.5); g_r03_sphenix_result_point->SetMarkerColor(kGreen+2); g_r03_sphenix_result_point->SetLineColor(kGreen+2);
  g_r04_sphenix_result->SetMarkerColor(kAzure+2); g_r04_sphenix_result_point->SetMarkerSize(1.5); g_r04_sphenix_result_point->SetMarkerColor(kAzure+2); g_r04_sphenix_result_point->SetLineColor(kAzure+2);
  g_r05_sphenix_result->SetMarkerColor(kBlue); g_r05_sphenix_result_point->SetMarkerSize(1.5); g_r05_sphenix_result_point->SetMarkerColor(kBlue); g_r05_sphenix_result_point->SetLineColor(kBlue);
  g_r06_sphenix_result->SetMarkerColor(kMagenta+2); g_r06_sphenix_result_point->SetMarkerSize(1.5); g_r06_sphenix_result_point->SetMarkerColor(kMagenta+2); g_r06_sphenix_result_point->SetLineColor(kMagenta+2);
  g_r08_sphenix_result->SetMarkerColor(kRed+2); g_r08_sphenix_result_point->SetMarkerSize(1.5); g_r08_sphenix_result_point->SetMarkerColor(kRed+2); g_r08_sphenix_result_point->SetLineColor(kRed+2);
  leg_all_result->AddEntry(g_r02_sphenix_result, "Anti-k_{t} #kern[-0.3]{#it{R}} = 0.2", "pf");
  leg_all_result->AddEntry(g_r03_sphenix_result, "Anti-k_{t} #kern[-0.3]{#it{R}} = 0.3", "pf");
  leg_all_result->AddEntry(g_r04_sphenix_result, "Anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", "pf");
  leg_all_result->AddEntry(g_r05_sphenix_result, "Anti-k_{t} #kern[-0.3]{#it{R}} = 0.5", "pf");
  leg_all_result->AddEntry(g_r06_sphenix_result, "Anti-k_{t} #kern[-0.3]{#it{R}} = 0.6", "pf");
  leg_all_result->AddEntry(g_r08_sphenix_result, "Anti-k_{t} #kern[-0.3]{#it{R}} = 0.8", "pf");
  leg_all_result->Draw();
  g_r02_sphenix_result_point->Draw("P same");
  g_r03_sphenix_result_point->Draw("P same");
  g_r04_sphenix_result_point->Draw("P same");
  g_r05_sphenix_result_point->Draw("P same");
  g_r06_sphenix_result_point->Draw("P same");
  g_r08_sphenix_result_point->Draw("P same");
  can_all_result->SaveAs("figure/spectrum_all_result.png");
  gStyle->SetPalette(55);
  g_r02_sphenix_result->SetFillColorAlpha(kAzure+1, 0.5); g_r02_sphenix_result->SetMarkerColor(kAzure+2);
  g_r03_sphenix_result->SetFillColorAlpha(kAzure+1, 0.5); g_r03_sphenix_result->SetMarkerColor(kAzure+2);
  g_r04_sphenix_result->SetFillColorAlpha(kAzure+1, 0.5); g_r04_sphenix_result->SetMarkerColor(kAzure+2);
  g_r05_sphenix_result->SetFillColorAlpha(kAzure+1, 0.5); g_r05_sphenix_result->SetMarkerColor(kAzure+2);
  g_r06_sphenix_result->SetFillColorAlpha(kAzure+1, 0.5); g_r06_sphenix_result->SetMarkerColor(kAzure+2);
  g_r08_sphenix_result->SetFillColorAlpha(kAzure+1, 0.5); g_r08_sphenix_result->SetMarkerColor(kAzure+2);
  g_r02_sphenix_result_point->SetMarkerColor(kAzure+2); g_r02_sphenix_result_point->SetLineColor(kAzure+2);
  g_r03_sphenix_result_point->SetMarkerColor(kAzure+2); g_r03_sphenix_result_point->SetLineColor(kAzure+2);
  g_r04_sphenix_result_point->SetMarkerColor(kAzure+2); g_r04_sphenix_result_point->SetLineColor(kAzure+2);
  g_r05_sphenix_result_point->SetMarkerColor(kAzure+2); g_r05_sphenix_result_point->SetLineColor(kAzure+2);
  g_r06_sphenix_result_point->SetMarkerColor(kAzure+2); g_r06_sphenix_result_point->SetLineColor(kAzure+2);
  g_r08_sphenix_result_point->SetMarkerColor(kAzure+2); g_r08_sphenix_result_point->SetLineColor(kAzure+2);

  //TCanvas *can_r04_all_result = new TCanvas("can_r04_all_result", "", 850, 800);
  //gStyle->SetPalette(57);
  //can_r04_all_result->SetTopMargin(0.03);
  //can_r04_all_result->SetLeftMargin(0.15);
  //can_r04_all_result->SetBottomMargin(0.12);
  //can_r04_all_result->SetRightMargin(0.05);
  //can_r04_all_result->SetLogy();
  //TH2F *frame_r04_all_result = new TH2F("frame_r04_all_result", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r04_sphenix_result->GetN(), g_r04_sphenix_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r04_sphenix_result->GetN(), g_r04_sphenix_result->GetY())));
  //frame_r04_all_result->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  //frame_r04_all_result->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  //frame_r04_all_result->GetXaxis()->SetTitleOffset(0.98);
  //frame_r04_all_result->GetYaxis()->SetTitleOffset(1.15);
  //frame_r04_all_result->GetXaxis()->SetLabelSize(0.045);
  //frame_r04_all_result->GetYaxis()->SetLabelSize(0.045);
  //frame_r04_all_result->GetXaxis()->CenterTitle();
  //frame_r04_all_result->GetYaxis()->CenterTitle();
  //frame_r04_all_result->Draw();
  //g_r04_sphenix_result->Draw("2 same");
  //myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  //myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200 GeV", 0.04);
  //myText(0.6, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.04);
  //myText(0.6, 0.75, 1, "|#eta^{jet}| < 0.7", 0.04);
  //TLegend *leg_r04_all_result = new TLegend(0.15, 0.08, 0.9, 0.33);
  //leg_r04_all_result->SetBorderSize(0);
  //leg_r04_all_result->SetFillStyle(0);
  //leg_r04_all_result->SetTextSize(0.03);
  //leg_r04_all_result->SetTextFont(42);
  //g_pythia_result->SetMarkerStyle(54);
  //leg_r04_all_result->AddEntry(g_r04_sphenix_result, "sPHENIX Run 2024 data", "pf");
  //leg_r04_all_result->AddEntry(g_pythia_result, "PYTHIA8 truth jet spectrum", "pl");
  //leg_r04_all_result->AddEntry(tjet_nnlo, "NNLO", "pf");
  //leg_r04_all_result->AddEntry(tjet_scet, "SCET", "pf");
  //leg_r04_all_result->AddEntry(tjet_shadow, "NLO pQCD (No Hadronization)", "f");
  //leg_r04_all_result->AddEntry(tjet_shadow, "NLO(parton) + CT18 PDF", "");
  //leg_r04_all_result->Draw();
  //g_r04_sphenix_result->Draw("P same");
  //g_r04_sphenix_result_point->Draw("P same");
  //g_pythia_result->Draw("p same");
  //tjet_shadow->Draw("F same");
  //tjet_nnlo->Draw("2 same");
  //tjet_scet->Draw("2 same");
  //tjet_scet_point->SetMarkerSize(1.5);
  //tjet_scet_point->Draw("p same");
  //tjet_nnlo_point->SetMarkerSize(1.5);
  //tjet_nnlo_point->Draw("p same");
  //can_r04_all_result->SaveAs("figure/spectrum_r04_all_result.png");
  //gStyle->SetPalette(55);

  TCanvas *can_phenix_comparison = new TCanvas("can_phenix_comparison", "", 850, 800);
  gStyle->SetPalette(57);
  can_phenix_comparison->SetTopMargin(0.03);
  can_phenix_comparison->SetLeftMargin(0.15);
  can_phenix_comparison->SetBottomMargin(0.12);
  can_phenix_comparison->SetRightMargin(0.05);
  can_phenix_comparison->SetLogy();
  TH2F *frame_phenix_comparison = new TH2F("frame_phenix_comparison", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r03_sphenix_result->GetN(), g_r03_sphenix_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r03_sphenix_result->GetN(), g_r03_sphenix_result->GetY())));
  frame_phenix_comparison->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_phenix_comparison->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_phenix_comparison->GetXaxis()->SetTitleOffset(0.98);
  frame_phenix_comparison->GetYaxis()->SetTitleOffset(1.15);
  frame_phenix_comparison->GetXaxis()->SetLabelSize(0.045);
  frame_phenix_comparison->GetYaxis()->SetLabelSize(0.045);
  frame_phenix_comparison->GetXaxis()->CenterTitle();
  frame_phenix_comparison->GetYaxis()->CenterTitle();
  frame_phenix_comparison->Draw();
  g_r03_sphenix_result->SetFillColorAlpha(kAzure+1, 0.5); g_r03_sphenix_result->SetMarkerColor(kAzure+2); g_r03_sphenix_result_point->SetMarkerSize(1.5); g_r03_sphenix_result_point->SetMarkerColor(kAzure+2); g_r03_sphenix_result_point->SetLineColor(kAzure+2);
  g_phenix_new_result->SetFillColorAlpha(kGray+1, 0.5); g_phenix_new_result->SetMarkerColor(kGray+3); g_phenix_new_result_point->SetMarkerSize(1.5); g_phenix_new_result_point->SetMarkerColor(kGray+3); g_phenix_new_result_point->SetLineColor(kGray+3);
  g_r03_sphenix_result->Draw("2 same");
  g_phenix_new_result->Draw("2 same");
  myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200 GeV", 0.04);
  myText(0.6, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.3", 0.04);
  TLegend *leg_phenix_comparison = new TLegend(0.15, 0.08, 0.9, 0.33);
  leg_phenix_comparison->SetBorderSize(0);
  leg_phenix_comparison->SetFillStyle(0);
  leg_phenix_comparison->SetTextSize(0.03);
  leg_phenix_comparison->SetTextFont(42);
  leg_phenix_comparison->AddEntry(g_r03_sphenix_result, "sPHENIX Run 2024", "pf");
  leg_phenix_comparison->AddEntry(g_phenix_new_result, "PHENIX Run12 + Run15", "pf");
  leg_phenix_comparison->Draw();
  g_r03_sphenix_result_point->Draw("P same");
  g_phenix_new_result->Draw("P same");
  can_phenix_comparison->SaveAs("figure/spectrum_phenix_comparison.png");
  gStyle->SetPalette(55);

  TCanvas *can_phenix_comparison_etacorrected = new TCanvas("can_phenix_comparison_etacorrected", "", 850, 800);
  gStyle->SetPalette(57);
  can_phenix_comparison_etacorrected->SetTopMargin(0.03);
  can_phenix_comparison_etacorrected->SetLeftMargin(0.15);
  can_phenix_comparison_etacorrected->SetBottomMargin(0.12);
  can_phenix_comparison_etacorrected->SetRightMargin(0.05);
  can_phenix_comparison_etacorrected->SetLogy();
  TH2F *frame_phenix_comparison_etacorrected = new TH2F("frame_phenix_comparison_etacorrected", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r03_sphenix_result->GetN(), g_r03_sphenix_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r03_sphenix_result->GetN(), g_r03_sphenix_result->GetY())));
  frame_phenix_comparison_etacorrected->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_phenix_comparison_etacorrected->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_phenix_comparison_etacorrected->GetXaxis()->SetTitleOffset(0.98);
  frame_phenix_comparison_etacorrected->GetYaxis()->SetTitleOffset(1.15);
  frame_phenix_comparison_etacorrected->GetXaxis()->SetLabelSize(0.045);
  frame_phenix_comparison_etacorrected->GetYaxis()->SetLabelSize(0.045);
  frame_phenix_comparison_etacorrected->GetXaxis()->CenterTitle();
  frame_phenix_comparison_etacorrected->GetYaxis()->CenterTitle();
  frame_phenix_comparison_etacorrected->Draw();
  g_r03_sphenix_result->SetFillColorAlpha(kAzure+1, 0.5); g_r03_sphenix_result->SetMarkerColor(kAzure+2); g_r03_sphenix_result_point->SetMarkerSize(1.5); g_r03_sphenix_result_point->SetMarkerColor(kAzure+2); g_r03_sphenix_result_point->SetLineColor(kAzure+2);
  g_phenix_new_result_etacorrected->SetFillColorAlpha(kGray, 0.5); g_phenix_new_result_etacorrected->SetMarkerColor(kGray+2); g_phenix_new_result_etacorrected_point->SetMarkerSize(1.5); g_phenix_new_result_etacorrected_point->SetMarkerColor(kGray+2); g_phenix_new_result_etacorrected_point->SetLineColor(kGray+2);
  g_r03_sphenix_result->Draw("2 same");
  g_phenix_new_result_etacorrected->Draw("2 same");
  myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200 GeV", 0.04);
  myText(0.6, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.3", 0.04);
  TLegend *leg_phenix_comparison_etacorrected = new TLegend(0.15, 0.08, 0.9, 0.33);
  leg_phenix_comparison_etacorrected->SetBorderSize(0);
  leg_phenix_comparison_etacorrected->SetFillStyle(0);
  leg_phenix_comparison_etacorrected->SetTextSize(0.03);
  leg_phenix_comparison_etacorrected->SetTextFont(42);
  leg_phenix_comparison_etacorrected->AddEntry(g_r03_sphenix_result, "sPHENIX Run 2024", "pf");
  leg_phenix_comparison_etacorrected->AddEntry(g_phenix_new_result_etacorrected, "PHENIX Run12 + Run15", "pf");
  leg_phenix_comparison_etacorrected->Draw();
  g_r03_sphenix_result_point->Draw("P same");
  g_phenix_new_result_etacorrected_point->Draw("P same");
  can_phenix_comparison_etacorrected->SaveAs("figure/spectrum_phenix_comparison_etacorrected.png");
  gStyle->SetPalette(55);

  TCanvas *can_phenix_comparison_both = new TCanvas("can_phenix_comparison_both", "", 850, 800);
  gStyle->SetPalette(57);
  can_phenix_comparison_both->SetTopMargin(0.03);
  can_phenix_comparison_both->SetLeftMargin(0.15);
  can_phenix_comparison_both->SetBottomMargin(0.12);
  can_phenix_comparison_both->SetRightMargin(0.05);
  can_phenix_comparison_both->SetLogy();
  TH2F *frame_phenix_comparison_both = new TH2F("frame_phenix_comparison_both", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r03_sphenix_result->GetN(), g_r03_sphenix_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r03_sphenix_result->GetN(), g_r03_sphenix_result->GetY())));
  frame_phenix_comparison_both->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_phenix_comparison_both->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_phenix_comparison_both->GetXaxis()->SetTitleOffset(0.98);
  frame_phenix_comparison_both->GetYaxis()->SetTitleOffset(1.15);
  frame_phenix_comparison_both->GetXaxis()->SetLabelSize(0.045);
  frame_phenix_comparison_both->GetYaxis()->SetLabelSize(0.045);
  frame_phenix_comparison_both->GetXaxis()->CenterTitle();
  frame_phenix_comparison_both->GetYaxis()->CenterTitle();
  frame_phenix_comparison_both->Draw();
  g_r03_sphenix_result->SetFillColorAlpha(kAzure+1, 0.5);
  g_r03_sphenix_result->SetMarkerColor(kAzure+2); g_r03_sphenix_result_point->SetMarkerSize(1.5); g_r03_sphenix_result_point->SetMarkerColor(kAzure+2); g_r03_sphenix_result_point->SetLineColor(kAzure+2);
  g_phenix_new_result->SetFillColorAlpha(kGray+1, 0.5); g_phenix_new_result->SetMarkerColor(kGray+3); g_phenix_new_result_point->SetMarkerSize(1.5); g_phenix_new_result_point->SetMarkerColor(kGray+3); g_phenix_new_result_point->SetLineColor(kGray+3);
  g_phenix_new_result_etacorrected->SetFillColorAlpha(kGray, 0.5); g_phenix_new_result_etacorrected->SetMarkerColor(kGray+2); g_phenix_new_result_etacorrected_point->SetMarkerSize(1.5); g_phenix_new_result_etacorrected_point->SetMarkerColor(kGray+2); g_phenix_new_result_etacorrected_point->SetLineColor(kGray+2);
  g_r03_sphenix_result->Draw("2 same");
  g_phenix_new_result->Draw("2 same");
  g_phenix_new_result_etacorrected->Draw("2 same");
  myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200GeV", 0.04);
  myText(0.6, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.3", 0.04);
  TLegend *leg_phenix_comparison_both = new TLegend(0.15, 0.08, 0.9, 0.33);
  leg_phenix_comparison_both->SetBorderSize(0);
  leg_phenix_comparison_both->SetFillStyle(0);
  leg_phenix_comparison_both->SetTextSize(0.03);
  leg_phenix_comparison_both->SetTextFont(42);
  leg_phenix_comparison_both->AddEntry(g_r03_sphenix_result, "sPHENIX Run 2024", "pf");
  leg_phenix_comparison_both->AddEntry(g_phenix_new_result, "PHENIX Run12 + Run15", "pf");
  leg_phenix_comparison_both->AddEntry(g_phenix_new_result_etacorrected, "PHENIX Run12 + Run15 (#eta Corrected)", "pf");
  leg_phenix_comparison_both->Draw();
  g_r03_sphenix_result_point->Draw("P same");
  g_phenix_new_result_point->Draw("P same");
  g_phenix_new_result_etacorrected_point->Draw("P same");
  can_phenix_comparison_both->SaveAs("figure/spectrum_phenix_comparison_both.png");
  gStyle->SetPalette(55);

  TCanvas *can_star_comparison = new TCanvas("can_star_comparison", "", 850, 800);
  gStyle->SetPalette(57);
  can_star_comparison->SetTopMargin(0.03);
  can_star_comparison->SetLeftMargin(0.15);
  can_star_comparison->SetBottomMargin(0.12);
  can_star_comparison->SetRightMargin(0.05);
  can_star_comparison->SetLogy();
  TH2F *frame_star_comparison = new TH2F("frame_star_comparison", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r05_sphenix_result->GetN(), g_r05_sphenix_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r05_sphenix_result->GetN(), g_r05_sphenix_result->GetY())));
  frame_star_comparison->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_star_comparison->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_star_comparison->GetXaxis()->SetTitleOffset(0.98);
  frame_star_comparison->GetYaxis()->SetTitleOffset(1.15);
  frame_star_comparison->GetXaxis()->SetLabelSize(0.045);
  frame_star_comparison->GetYaxis()->SetLabelSize(0.045);
  frame_star_comparison->GetXaxis()->CenterTitle();
  frame_star_comparison->GetYaxis()->CenterTitle();
  frame_star_comparison->Draw();
  g_r05_sphenix_result->SetFillColorAlpha(kAzure+1, 0.5);
  g_r05_sphenix_result->Draw("2 same");
  g_star_new_result1->Draw("2 same");
  g_star_new_result2->Draw("2 same");
  myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200 GeV", 0.04);
  myText(0.6, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.5", 0.04);
  TLegend *leg_star_comparison = new TLegend(0.15, 0.08, 0.9, 0.33);
  leg_star_comparison->SetBorderSize(0);
  leg_star_comparison->SetFillStyle(0);
  leg_star_comparison->SetTextSize(0.03);
  leg_star_comparison->SetTextFont(42);
  g_r05_sphenix_result->SetMarkerColor(kAzure+2); g_r05_sphenix_result_point->SetMarkerSize(1.5); g_r05_sphenix_result_point->SetMarkerColor(kAzure+2); g_r05_sphenix_result_point->SetLineColor(kAzure+2);
  g_star_new_result1->SetMarkerColor(kGreen+2); g_star_new_result_point1->SetMarkerSize(1.5); g_star_new_result_point1->SetMarkerColor(kGreen+2); g_star_new_result_point1->SetLineColor(kGreen+2);
  g_star_new_result2->SetMarkerColor(kGreen-3); g_star_new_result_point2->SetMarkerSize(1.5); g_star_new_result_point2->SetMarkerColor(kGreen-3); g_star_new_result_point2->SetLineColor(kGreen-3);
  leg_star_comparison->AddEntry(g_r05_sphenix_result, "sPHENIX Run 2024 (|#eta| < 0.6)", "pf");
  leg_star_comparison->AddEntry(g_star_new_result1, "STAR Run 2012 (|#eta| < 0.5)", "pf");
  leg_star_comparison->AddEntry(g_star_new_result2, "STAR Run 2012 (0.5 < |#eta| < 0.9)", "pf");
  leg_star_comparison->Draw();
  g_r05_sphenix_result_point->Draw("P same");
  g_star_new_result1->Draw("P same");
  g_star_new_result2->Draw("P same");
  can_star_comparison->SaveAs("figure/spectrum_star_comparison.png");
  gStyle->SetPalette(55);

  // Draw ratio plot
  TCanvas *can_r02_ratio_all = new TCanvas("can_r02_ratio_all", "", 850, 1091);
  can_r02_ratio_all->Divide(1, 2);
  gStyle->SetPalette(57);
  TPad *pad_r02_1_all = (TPad*)can_r02_ratio_all->cd(1);
  pad_r02_1_all->SetPad(0, 0.33333, 1, 1);
  pad_r02_1_all->SetTopMargin(0.03);
  pad_r02_1_all->SetLeftMargin(0.15);
  pad_r02_1_all->SetBottomMargin(0.035);
  pad_r02_1_all->SetRightMargin(0.05);
  pad_r02_1_all->SetLogy();
  TH2F *frame_r02_ratio_all = new TH2F("frame_r02_ratio_all", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r02_sphenix_result->GetN(), g_r02_sphenix_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r02_sphenix_result->GetN(), g_r02_sphenix_result->GetY())));
  frame_r02_ratio_all->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_r02_ratio_all->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_r02_ratio_all->GetXaxis()->SetTitleOffset(0.98);
  frame_r02_ratio_all->GetYaxis()->SetTitleOffset(1.15);
  frame_r02_ratio_all->GetXaxis()->SetLabelSize(0.045);
  frame_r02_ratio_all->GetYaxis()->SetLabelSize(0.045);
  frame_r02_ratio_all->GetXaxis()->SetLabelOffset(2);
  frame_r02_ratio_all->GetXaxis()->CenterTitle();
  frame_r02_ratio_all->GetYaxis()->CenterTitle();
  frame_r02_ratio_all->Draw();
  myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200 GeV", 0.04);
  myText(0.6, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.2", 0.04);
  myText(0.6, 0.75, 1, "|#eta^{jet}| < 0.9", 0.04);
  TLegend *leg_r02_ratio_all = new TLegend(0.15, 0.08, 0.9, 0.33);
  leg_r02_ratio_all->SetBorderSize(0);
  leg_r02_ratio_all->SetFillStyle(0);
  leg_r02_ratio_all->SetTextSize(0.03);
  leg_r02_ratio_all->SetTextFont(42);
  leg_r02_ratio_all->AddEntry(g_r02_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg_r02_ratio_all->AddEntry(spline_r02_voge, "NLO(parton) + CT18 PDF", "l");
  leg_r02_ratio_all->AddEntry(t_r02jet_shadow, "NLO(parton) scale uncertainty", "f");
  leg_r02_ratio_all->AddEntry(t_r02jet_voge_shadow, "NLO(parton) PDF 68% CL", "f");
  leg_r02_ratio_all->AddEntry(spline_r02_nnpdf, "NLO(parton) + NNPDF3.0", "l");
  leg_r02_ratio_all->AddEntry(spline_r02_msht, "NLO(parton) + MSHT20", "l");
  leg_r02_ratio_all->Draw();
  g_r02_sphenix_result->Draw("2 same");
  g_r02_sphenix_result_point->Draw("P same");
  t_r02jet_shadow->Draw("F same");
  t_r02jet_voge_shadow->Draw("F same");
  spline_r02_voge->Draw("L same");
  spline_r02_nnpdf->Draw("L same");
  spline_r02_msht->Draw("L same");
  TPad *pad_r02_2_all = (TPad*)can_r02_ratio_all->cd(2);
  pad_r02_2_all->SetPad(0, 0, 1, 0.333333);
  pad_r02_2_all->SetTopMargin(0.02);
  pad_r02_2_all->SetLeftMargin(0.15);
  pad_r02_2_all->SetBottomMargin(0.25);
  pad_r02_2_all->SetRightMargin(0.05);
  TH2F *frame_r02_ratio_all2 = new TH2F("frame_r02_ratio_all2", "", 10, calibptbins[0], calibptbins[calibnpt], 2, 0., 1.6);
  frame_r02_ratio_all2->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_r02_ratio_all2->GetYaxis()->SetTitle("Ratio");
  frame_r02_ratio_all2->GetXaxis()->SetTitleOffset(0.98);
  frame_r02_ratio_all2->GetYaxis()->SetTitleOffset(0.5);
  frame_r02_ratio_all2->GetXaxis()->SetLabelSize(0.045);
  frame_r02_ratio_all2->GetYaxis()->SetLabelSize(0.045);
  frame_r02_ratio_all2->GetXaxis()->CenterTitle();
  frame_r02_ratio_all2->GetYaxis()->CenterTitle();
  frame_r02_ratio_all2->GetYaxis()->SetTitleOffset(frame_r02_ratio_all->GetYaxis()->GetTitleOffset()*1/2.);
  frame_r02_ratio_all2->GetYaxis()->SetLabelOffset(frame_r02_ratio_all->GetYaxis()->GetLabelOffset()*1/2.);
  frame_r02_ratio_all2->GetXaxis()->SetLabelSize(frame_r02_ratio_all->GetXaxis()->GetLabelSize()*2);
  frame_r02_ratio_all2->GetYaxis()->SetLabelSize(frame_r02_ratio_all->GetYaxis()->GetLabelSize()*2);
  frame_r02_ratio_all2->GetXaxis()->SetTitleSize(frame_r02_ratio_all->GetXaxis()->GetTitleSize()*2);
  frame_r02_ratio_all2->GetYaxis()->SetTitleSize(frame_r02_ratio_all->GetYaxis()->GetTitleSize()*2);
  frame_r02_ratio_all2->Draw();
  g_r02_sphenix_ratio_all->Draw("2 same");
  g_r02_sphenix_ratio_all->Draw("P same");
  g_r02_sphenix_ratio_all_point->Draw("P same");
  t_r02jet_ratio_shadow->Draw("F same");
  t_r02jet_voge_ratioshadow->Draw("F same");
  spline_r02_nnpdf_ratio->Draw("L same");
  spline_r02_msht_ratio->Draw("L same");
  TLine *line_r02_ratio_all = new TLine(calibptbins[0], 1, calibptbins[calibnpt], 1);
  line_r02_ratio_all->SetLineStyle(3);
  line_r02_ratio_all->Draw("same");
  can_r02_ratio_all->SaveAs("figure/spectrum_r02_vogeratio_all.png");

  TCanvas *can_r03_ratio_all = new TCanvas("can_r03_ratio_all", "", 850, 1091);
  can_r03_ratio_all->Divide(1, 2);
  gStyle->SetPalette(57);
  TPad *pad_r03_1_all = (TPad*)can_r03_ratio_all->cd(1);
  pad_r03_1_all->SetPad(0, 0.33333, 1, 1);
  pad_r03_1_all->SetTopMargin(0.03);
  pad_r03_1_all->SetLeftMargin(0.15);
  pad_r03_1_all->SetBottomMargin(0.035);
  pad_r03_1_all->SetRightMargin(0.05);
  pad_r03_1_all->SetLogy();
  TH2F *frame_r03_ratio_all = new TH2F("frame_r03_ratio_all", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r03_sphenix_result->GetN(), g_r03_sphenix_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r03_sphenix_result->GetN(), g_r03_sphenix_result->GetY())));
  frame_r03_ratio_all->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_r03_ratio_all->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_r03_ratio_all->GetXaxis()->SetTitleOffset(0.98);
  frame_r03_ratio_all->GetYaxis()->SetTitleOffset(1.15);
  frame_r03_ratio_all->GetXaxis()->SetLabelSize(0.045);
  frame_r03_ratio_all->GetYaxis()->SetLabelSize(0.045);
  frame_r03_ratio_all->GetXaxis()->SetLabelOffset(2);
  frame_r03_ratio_all->GetXaxis()->CenterTitle();
  frame_r03_ratio_all->GetYaxis()->CenterTitle();
  frame_r03_ratio_all->Draw();
  myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200 GeV", 0.04);
  myText(0.6, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.3", 0.04);
  myText(0.6, 0.75, 1, "|#eta^{jet}| < 0.8", 0.04);
  TLegend *leg_r03_ratio_all = new TLegend(0.15, 0.08, 0.9, 0.33);
  leg_r03_ratio_all->SetBorderSize(0);
  leg_r03_ratio_all->SetFillStyle(0);
  leg_r03_ratio_all->SetTextSize(0.03);
  leg_r03_ratio_all->SetTextFont(42);
  leg_r03_ratio_all->AddEntry(g_r03_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg_r03_ratio_all->AddEntry(spline_r03_voge, "NLO(parton) + CT18 PDF", "l");
  leg_r03_ratio_all->AddEntry(t_r03jet_shadow, "NLO(parton) scale uncertainty", "f");
  leg_r03_ratio_all->AddEntry(t_r03jet_voge_shadow, "NLO(parton) PDF 68% CL", "f");
  leg_r03_ratio_all->AddEntry(spline_r03_nnpdf, "NLO(parton) + NNPDF3.0", "l");
  leg_r03_ratio_all->AddEntry(spline_r03_msht, "NLO(parton) + MSHT20", "l");
  leg_r03_ratio_all->Draw();
  g_r03_sphenix_result->Draw("2 same");
  g_r03_sphenix_result_point->Draw("P same");
  t_r03jet_shadow->Draw("F same");
  t_r03jet_voge_shadow->Draw("F same");
  spline_r03_voge->Draw("L same");
  spline_r03_nnpdf->Draw("L same");
  spline_r03_msht->Draw("L same");
  TPad *pad_r03_2_all = (TPad*)can_r03_ratio_all->cd(2);
  pad_r03_2_all->SetPad(0, 0, 1, 0.333333);
  pad_r03_2_all->SetTopMargin(0.02);
  pad_r03_2_all->SetLeftMargin(0.15);
  pad_r03_2_all->SetBottomMargin(0.25);
  pad_r03_2_all->SetRightMargin(0.05);
  TH2F *frame_r03_ratio_all2 = new TH2F("frame_r03_ratio_all2", "", 10, calibptbins[0], calibptbins[calibnpt], 2, 0., 1.6);
  frame_r03_ratio_all2->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_r03_ratio_all2->GetYaxis()->SetTitle("Ratio");
  frame_r03_ratio_all2->GetXaxis()->SetTitleOffset(0.98);
  frame_r03_ratio_all2->GetYaxis()->SetTitleOffset(0.5);
  frame_r03_ratio_all2->GetXaxis()->SetLabelSize(0.045);
  frame_r03_ratio_all2->GetYaxis()->SetLabelSize(0.045);
  frame_r03_ratio_all2->GetXaxis()->CenterTitle();
  frame_r03_ratio_all2->GetYaxis()->CenterTitle();
  frame_r03_ratio_all2->GetYaxis()->SetTitleOffset(frame_r03_ratio_all->GetYaxis()->GetTitleOffset()*1/2.);
  frame_r03_ratio_all2->GetYaxis()->SetLabelOffset(frame_r03_ratio_all->GetYaxis()->GetLabelOffset()*1/2.);
  frame_r03_ratio_all2->GetXaxis()->SetLabelSize(frame_r03_ratio_all->GetXaxis()->GetLabelSize()*2);
  frame_r03_ratio_all2->GetYaxis()->SetLabelSize(frame_r03_ratio_all->GetYaxis()->GetLabelSize()*2);
  frame_r03_ratio_all2->GetXaxis()->SetTitleSize(frame_r03_ratio_all->GetXaxis()->GetTitleSize()*2);
  frame_r03_ratio_all2->GetYaxis()->SetTitleSize(frame_r03_ratio_all->GetYaxis()->GetTitleSize()*2);
  frame_r03_ratio_all2->Draw();
  g_r03_sphenix_ratio_all->Draw("2 same");
  g_r03_sphenix_ratio_all->Draw("P same");
  g_r03_sphenix_ratio_all_point->Draw("P same");
  t_r03jet_ratio_shadow->Draw("F same");
  t_r03jet_voge_ratioshadow->Draw("F same");
  spline_r03_nnpdf_ratio->Draw("L same");
  spline_r03_msht_ratio->Draw("L same");
  TLine *line_r03_ratio_all = new TLine(calibptbins[0], 1, calibptbins[calibnpt], 1);
  line_r03_ratio_all->SetLineStyle(3);
  line_r03_ratio_all->Draw("same");
  can_r03_ratio_all->SaveAs("figure/spectrum_r03_vogeratio_all.png");

  TCanvas *can_r04_ratio_all = new TCanvas("can_r04_ratio_all", "", 850, 1091);
  can_r04_ratio_all->Divide(1, 2);
  gStyle->SetPalette(57);
  TPad *pad_r04_1_all = (TPad*)can_r04_ratio_all->cd(1);
  pad_r04_1_all->SetPad(0, 0.33333, 1, 1);
  pad_r04_1_all->SetTopMargin(0.03);
  pad_r04_1_all->SetLeftMargin(0.15);
  pad_r04_1_all->SetBottomMargin(0.035);
  pad_r04_1_all->SetRightMargin(0.05);
  pad_r04_1_all->SetLogy();
  TH2F *frame_r04_ratio_all = new TH2F("frame_r04_ratio_all", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r04_sphenix_result->GetN(), g_r04_sphenix_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r04_sphenix_result->GetN(), g_r04_sphenix_result->GetY())));
  frame_r04_ratio_all->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_r04_ratio_all->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_r04_ratio_all->GetXaxis()->SetTitleOffset(0.98);
  frame_r04_ratio_all->GetYaxis()->SetTitleOffset(1.15);
  frame_r04_ratio_all->GetXaxis()->SetLabelSize(0.045);
  frame_r04_ratio_all->GetYaxis()->SetLabelSize(0.045);
  frame_r04_ratio_all->GetXaxis()->SetLabelOffset(2);
  frame_r04_ratio_all->GetXaxis()->CenterTitle();
  frame_r04_ratio_all->GetYaxis()->CenterTitle();
  frame_r04_ratio_all->Draw();
  myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200 GeV", 0.04);
  myText(0.6, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.04);
  myText(0.6, 0.75, 1, "|#eta^{jet}| < 0.7", 0.04);
  TLegend *leg_r04_ratio_all = new TLegend(0.15, 0.08, 0.9, 0.33);
  leg_r04_ratio_all->SetBorderSize(0);
  leg_r04_ratio_all->SetFillStyle(0);
  leg_r04_ratio_all->SetTextSize(0.03);
  leg_r04_ratio_all->SetTextFont(42);
  leg_r04_ratio_all->AddEntry(g_r04_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg_r04_ratio_all->AddEntry(spline_r04_voge, "NLO(parton) + CT18 PDF", "l");
  leg_r04_ratio_all->AddEntry(t_r04jet_shadow, "NLO(parton) scale uncertainty", "f");
  leg_r04_ratio_all->AddEntry(t_r04jet_voge_shadow, "NLO(parton) PDF 68% CL", "f");
  leg_r04_ratio_all->AddEntry(spline_r04_nnpdf, "NLO(parton) + NNPDF3.0", "l");
  leg_r04_ratio_all->AddEntry(spline_r04_msht, "NLO(parton) + MSHT20", "l");
  leg_r04_ratio_all->Draw();
  g_r04_sphenix_result->Draw("2 same");
  g_r04_sphenix_result_point->Draw("P same");
  t_r04jet_shadow->Draw("F same");
  t_r04jet_voge_shadow->Draw("F same");
  spline_r04_voge->Draw("L same");
  spline_r04_nnpdf->Draw("L same");
  spline_r04_msht->Draw("L same");
  TPad *pad_r04_2_all = (TPad*)can_r04_ratio_all->cd(2);
  pad_r04_2_all->SetPad(0, 0, 1, 0.333333);
  pad_r04_2_all->SetTopMargin(0.02);
  pad_r04_2_all->SetLeftMargin(0.15);
  pad_r04_2_all->SetBottomMargin(0.25);
  pad_r04_2_all->SetRightMargin(0.05);
  TH2F *frame_r04_ratio_all2 = new TH2F("frame_r04_ratio_all2", "", 10, calibptbins[0], calibptbins[calibnpt], 2, 0., 1.6);
  frame_r04_ratio_all2->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_r04_ratio_all2->GetYaxis()->SetTitle("Ratio");
  frame_r04_ratio_all2->GetXaxis()->SetTitleOffset(0.98);
  frame_r04_ratio_all2->GetYaxis()->SetTitleOffset(0.5);
  frame_r04_ratio_all2->GetXaxis()->SetLabelSize(0.045);
  frame_r04_ratio_all2->GetYaxis()->SetLabelSize(0.045);
  frame_r04_ratio_all2->GetXaxis()->CenterTitle();
  frame_r04_ratio_all2->GetYaxis()->CenterTitle();
  frame_r04_ratio_all2->GetYaxis()->SetTitleOffset(frame_r04_ratio_all->GetYaxis()->GetTitleOffset()*1/2.);
  frame_r04_ratio_all2->GetYaxis()->SetLabelOffset(frame_r04_ratio_all->GetYaxis()->GetLabelOffset()*1/2.);
  frame_r04_ratio_all2->GetXaxis()->SetLabelSize(frame_r04_ratio_all->GetXaxis()->GetLabelSize()*2);
  frame_r04_ratio_all2->GetYaxis()->SetLabelSize(frame_r04_ratio_all->GetYaxis()->GetLabelSize()*2);
  frame_r04_ratio_all2->GetXaxis()->SetTitleSize(frame_r04_ratio_all->GetXaxis()->GetTitleSize()*2);
  frame_r04_ratio_all2->GetYaxis()->SetTitleSize(frame_r04_ratio_all->GetYaxis()->GetTitleSize()*2);
  frame_r04_ratio_all2->Draw();
  g_r04_sphenix_ratio_all->Draw("2 same");
  g_r04_sphenix_ratio_all->Draw("P same");
  g_r04_sphenix_ratio_all_point->Draw("P same");
  t_r04jet_ratio_shadow->Draw("F same");
  t_r04jet_voge_ratioshadow->Draw("F same");
  spline_r04_nnpdf_ratio->Draw("L same");
  spline_r04_msht_ratio->Draw("L same");
  TLine *line_r04_ratio_all = new TLine(calibptbins[0], 1, calibptbins[calibnpt], 1);
  line_r04_ratio_all->SetLineStyle(3);
  line_r04_ratio_all->Draw("same");
  can_r04_ratio_all->SaveAs("figure/spectrum_r04_vogeratio_all.png");

  TCanvas *can_r05_ratio_all = new TCanvas("can_r05_ratio_all", "", 850, 1091);
  can_r05_ratio_all->Divide(1, 2);
  gStyle->SetPalette(57);
  TPad *pad_r05_1_all = (TPad*)can_r05_ratio_all->cd(1);
  pad_r05_1_all->SetPad(0, 0.33333, 1, 1);
  pad_r05_1_all->SetTopMargin(0.03);
  pad_r05_1_all->SetLeftMargin(0.15);
  pad_r05_1_all->SetBottomMargin(0.035);
  pad_r05_1_all->SetRightMargin(0.05);
  pad_r05_1_all->SetLogy();
  TH2F *frame_r05_ratio_all = new TH2F("frame_r05_ratio_all", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r05_sphenix_result->GetN(), g_r05_sphenix_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r05_sphenix_result->GetN(), g_r05_sphenix_result->GetY())));
  frame_r05_ratio_all->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_r05_ratio_all->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_r05_ratio_all->GetXaxis()->SetTitleOffset(0.98);
  frame_r05_ratio_all->GetYaxis()->SetTitleOffset(1.15);
  frame_r05_ratio_all->GetXaxis()->SetLabelSize(0.045);
  frame_r05_ratio_all->GetYaxis()->SetLabelSize(0.045);
  frame_r05_ratio_all->GetXaxis()->SetLabelOffset(2);
  frame_r05_ratio_all->GetXaxis()->CenterTitle();
  frame_r05_ratio_all->GetYaxis()->CenterTitle();
  frame_r05_ratio_all->Draw();
  myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200 GeV", 0.04);
  myText(0.6, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.5", 0.04);
  myText(0.6, 0.75, 1, "|#eta^{jet}| < 0.6", 0.04);
  TLegend *leg_r05_ratio_all = new TLegend(0.15, 0.08, 0.9, 0.33);
  leg_r05_ratio_all->SetBorderSize(0);
  leg_r05_ratio_all->SetFillStyle(0);
  leg_r05_ratio_all->SetTextSize(0.03);
  leg_r05_ratio_all->SetTextFont(42);
  leg_r05_ratio_all->AddEntry(g_r05_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg_r05_ratio_all->AddEntry(spline_r05_voge, "NLO(parton) + CT18 PDF", "l");
  leg_r05_ratio_all->AddEntry(t_r05jet_shadow, "NLO(parton) scale uncertainty", "f");
  leg_r05_ratio_all->AddEntry(t_r05jet_voge_shadow, "NLO(parton) PDF 68% CL", "f");
  leg_r05_ratio_all->AddEntry(spline_r05_nnpdf, "NLO(parton) + NNPDF3.0", "l");
  leg_r05_ratio_all->AddEntry(spline_r05_msht, "NLO(parton) + MSHT20", "l");
  leg_r05_ratio_all->Draw();
  g_r05_sphenix_result->Draw("2 same");
  g_r05_sphenix_result_point->Draw("P same");
  t_r05jet_shadow->Draw("F same");
  t_r05jet_voge_shadow->Draw("F same");
  spline_r05_voge->Draw("L same");
  spline_r05_nnpdf->Draw("L same");
  spline_r05_msht->Draw("L same");
  TPad *pad_r05_2_all = (TPad*)can_r05_ratio_all->cd(2);
  pad_r05_2_all->SetPad(0, 0, 1, 0.333333);
  pad_r05_2_all->SetTopMargin(0.02);
  pad_r05_2_all->SetLeftMargin(0.15);
  pad_r05_2_all->SetBottomMargin(0.25);
  pad_r05_2_all->SetRightMargin(0.05);
  TH2F *frame_r05_ratio_all2 = new TH2F("frame_r05_ratio_all2", "", 10, calibptbins[0], calibptbins[calibnpt], 2, 0., 1.6);
  frame_r05_ratio_all2->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_r05_ratio_all2->GetYaxis()->SetTitle("Ratio");
  frame_r05_ratio_all2->GetXaxis()->SetTitleOffset(0.98);
  frame_r05_ratio_all2->GetYaxis()->SetTitleOffset(0.5);
  frame_r05_ratio_all2->GetXaxis()->SetLabelSize(0.045);
  frame_r05_ratio_all2->GetYaxis()->SetLabelSize(0.045);
  frame_r05_ratio_all2->GetXaxis()->CenterTitle();
  frame_r05_ratio_all2->GetYaxis()->CenterTitle();
  frame_r05_ratio_all2->GetYaxis()->SetTitleOffset(frame_r05_ratio_all->GetYaxis()->GetTitleOffset()*1/2.);
  frame_r05_ratio_all2->GetYaxis()->SetLabelOffset(frame_r05_ratio_all->GetYaxis()->GetLabelOffset()*1/2.);
  frame_r05_ratio_all2->GetXaxis()->SetLabelSize(frame_r05_ratio_all->GetXaxis()->GetLabelSize()*2);
  frame_r05_ratio_all2->GetYaxis()->SetLabelSize(frame_r05_ratio_all->GetYaxis()->GetLabelSize()*2);
  frame_r05_ratio_all2->GetXaxis()->SetTitleSize(frame_r05_ratio_all->GetXaxis()->GetTitleSize()*2);
  frame_r05_ratio_all2->GetYaxis()->SetTitleSize(frame_r05_ratio_all->GetYaxis()->GetTitleSize()*2);
  frame_r05_ratio_all2->Draw();
  g_r05_sphenix_ratio_all->Draw("2 same");
  g_r05_sphenix_ratio_all->Draw("P same");
  g_r05_sphenix_ratio_all_point->Draw("P same");
  t_r05jet_ratio_shadow->Draw("F same");
  t_r05jet_voge_ratioshadow->Draw("F same");
  spline_r05_nnpdf_ratio->Draw("L same");
  spline_r05_msht_ratio->Draw("L same");
  TLine *line_r05_ratio_all = new TLine(calibptbins[0], 1, calibptbins[calibnpt], 1);
  line_r05_ratio_all->SetLineStyle(3);
  line_r05_ratio_all->Draw("same");
  can_r05_ratio_all->SaveAs("figure/spectrum_r05_vogeratio_all.png");

  TCanvas *can_r06_ratio_all = new TCanvas("can_r06_ratio_all", "", 850, 1091);
  can_r06_ratio_all->Divide(1, 2);
  gStyle->SetPalette(57);
  TPad *pad_r06_1_all = (TPad*)can_r06_ratio_all->cd(1);
  pad_r06_1_all->SetPad(0, 0.33333, 1, 1);
  pad_r06_1_all->SetTopMargin(0.03);
  pad_r06_1_all->SetLeftMargin(0.15);
  pad_r06_1_all->SetBottomMargin(0.035);
  pad_r06_1_all->SetRightMargin(0.05);
  pad_r06_1_all->SetLogy();
  TH2F *frame_r06_ratio_all = new TH2F("frame_r06_ratio_all", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r06_sphenix_result->GetN(), g_r06_sphenix_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r06_sphenix_result->GetN(), g_r06_sphenix_result->GetY())));
  frame_r06_ratio_all->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_r06_ratio_all->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_r06_ratio_all->GetXaxis()->SetTitleOffset(0.98);
  frame_r06_ratio_all->GetYaxis()->SetTitleOffset(1.15);
  frame_r06_ratio_all->GetXaxis()->SetLabelSize(0.045);
  frame_r06_ratio_all->GetYaxis()->SetLabelSize(0.045);
  frame_r06_ratio_all->GetXaxis()->SetLabelOffset(2);
  frame_r06_ratio_all->GetXaxis()->CenterTitle();
  frame_r06_ratio_all->GetYaxis()->CenterTitle();
  frame_r06_ratio_all->Draw();
  myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200 GeV", 0.04);
  myText(0.6, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.6", 0.04);
  myText(0.6, 0.75, 1, "|#eta^{jet}| < 0.5", 0.04);
  TLegend *leg_r06_ratio_all = new TLegend(0.15, 0.08, 0.9, 0.33);
  leg_r06_ratio_all->SetBorderSize(0);
  leg_r06_ratio_all->SetFillStyle(0);
  leg_r06_ratio_all->SetTextSize(0.03);
  leg_r06_ratio_all->SetTextFont(42);
  leg_r06_ratio_all->AddEntry(g_r06_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg_r06_ratio_all->AddEntry(spline_r06_voge, "NLO(parton) + CT18 PDF", "l");
  leg_r06_ratio_all->AddEntry(t_r06jet_shadow, "NLO(parton) scale uncertainty", "f");
  leg_r06_ratio_all->AddEntry(t_r06jet_voge_shadow, "NLO(parton) PDF 68% CL", "f");
  leg_r06_ratio_all->AddEntry(spline_r06_nnpdf, "NLO(parton) + NNPDF3.0", "l");
  leg_r06_ratio_all->AddEntry(spline_r06_msht, "NLO(parton) + MSHT20", "l");
  leg_r06_ratio_all->Draw();
  g_r06_sphenix_result->Draw("2 same");
  g_r06_sphenix_result_point->Draw("P same");
  t_r06jet_shadow->Draw("F same");
  t_r06jet_voge_shadow->Draw("F same");
  spline_r06_voge->Draw("L same");
  spline_r06_nnpdf->Draw("L same");
  spline_r06_msht->Draw("L same");
  TPad *pad_r06_2_all = (TPad*)can_r06_ratio_all->cd(2);
  pad_r06_2_all->SetPad(0, 0, 1, 0.333333);
  pad_r06_2_all->SetTopMargin(0.02);
  pad_r06_2_all->SetLeftMargin(0.15);
  pad_r06_2_all->SetBottomMargin(0.25);
  pad_r06_2_all->SetRightMargin(0.05);
  TH2F *frame_r06_ratio_all2 = new TH2F("frame_r06_ratio_all2", "", 10, calibptbins[0], calibptbins[calibnpt], 2, 0., 1.6);
  frame_r06_ratio_all2->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_r06_ratio_all2->GetYaxis()->SetTitle("Ratio");
  frame_r06_ratio_all2->GetXaxis()->SetTitleOffset(0.98);
  frame_r06_ratio_all2->GetYaxis()->SetTitleOffset(0.5);
  frame_r06_ratio_all2->GetXaxis()->SetLabelSize(0.045);
  frame_r06_ratio_all2->GetYaxis()->SetLabelSize(0.045);
  frame_r06_ratio_all2->GetXaxis()->CenterTitle();
  frame_r06_ratio_all2->GetYaxis()->CenterTitle();
  frame_r06_ratio_all2->GetYaxis()->SetTitleOffset(frame_r06_ratio_all->GetYaxis()->GetTitleOffset()*1/2.);
  frame_r06_ratio_all2->GetYaxis()->SetLabelOffset(frame_r06_ratio_all->GetYaxis()->GetLabelOffset()*1/2.);
  frame_r06_ratio_all2->GetXaxis()->SetLabelSize(frame_r06_ratio_all->GetXaxis()->GetLabelSize()*2);
  frame_r06_ratio_all2->GetYaxis()->SetLabelSize(frame_r06_ratio_all->GetYaxis()->GetLabelSize()*2);
  frame_r06_ratio_all2->GetXaxis()->SetTitleSize(frame_r06_ratio_all->GetXaxis()->GetTitleSize()*2);
  frame_r06_ratio_all2->GetYaxis()->SetTitleSize(frame_r06_ratio_all->GetYaxis()->GetTitleSize()*2);
  frame_r06_ratio_all2->Draw();
  g_r06_sphenix_ratio_all->Draw("2 same");
  g_r06_sphenix_ratio_all->Draw("P same");
  g_r06_sphenix_ratio_all_point->Draw("P same");
  t_r06jet_ratio_shadow->Draw("F same");
  t_r06jet_voge_ratioshadow->Draw("F same");
  spline_r06_nnpdf_ratio->Draw("L same");
  spline_r06_msht_ratio->Draw("L same");
  TLine *line_r06_ratio_all = new TLine(calibptbins[0], 1, calibptbins[calibnpt], 1);
  line_r06_ratio_all->SetLineStyle(3);
  line_r06_ratio_all->Draw("same");
  can_r06_ratio_all->SaveAs("figure/spectrum_r06_vogeratio_all.png");

  TCanvas *can_r08_ratio_all = new TCanvas("can_r08_ratio_all", "", 850, 1091);
  can_r08_ratio_all->Divide(1, 2);
  gStyle->SetPalette(57);
  TPad *pad_r08_1_all = (TPad*)can_r08_ratio_all->cd(1);
  pad_r08_1_all->SetPad(0, 0.33333, 1, 1);
  pad_r08_1_all->SetTopMargin(0.03);
  pad_r08_1_all->SetLeftMargin(0.15);
  pad_r08_1_all->SetBottomMargin(0.035);
  pad_r08_1_all->SetRightMargin(0.05);
  pad_r08_1_all->SetLogy();
  TH2F *frame_r08_ratio_all = new TH2F("frame_r08_ratio_all", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r08_sphenix_result->GetN(), g_r08_sphenix_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r08_sphenix_result->GetN(), g_r08_sphenix_result->GetY())));
  frame_r08_ratio_all->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_r08_ratio_all->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_r08_ratio_all->GetXaxis()->SetTitleOffset(0.98);
  frame_r08_ratio_all->GetYaxis()->SetTitleOffset(1.15);
  frame_r08_ratio_all->GetXaxis()->SetLabelSize(0.045);
  frame_r08_ratio_all->GetYaxis()->SetLabelSize(0.045);
  frame_r08_ratio_all->GetXaxis()->SetLabelOffset(2);
  frame_r08_ratio_all->GetXaxis()->CenterTitle();
  frame_r08_ratio_all->GetYaxis()->CenterTitle();
  frame_r08_ratio_all->Draw();
  myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200 GeV", 0.04);
  myText(0.6, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.8", 0.04);
  myText(0.6, 0.75, 1, "|#eta^{jet}| < 0.3", 0.04);
  TLegend *leg_r08_ratio_all = new TLegend(0.15, 0.08, 0.9, 0.33);
  leg_r08_ratio_all->SetBorderSize(0);
  leg_r08_ratio_all->SetFillStyle(0);
  leg_r08_ratio_all->SetTextSize(0.03);
  leg_r08_ratio_all->SetTextFont(42);
  leg_r08_ratio_all->AddEntry(g_r08_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg_r08_ratio_all->AddEntry(spline_r08_voge, "NLO(parton) + CT18 PDF", "l");
  leg_r08_ratio_all->AddEntry(t_r08jet_shadow, "NLO(parton) scale uncertainty", "f");
  leg_r08_ratio_all->AddEntry(t_r08jet_voge_shadow, "NLO(parton) PDF 68% CL", "f");
  leg_r08_ratio_all->AddEntry(spline_r08_nnpdf, "NLO(parton) + NNPDF3.0", "l");
  leg_r08_ratio_all->AddEntry(spline_r08_msht, "NLO(parton) + MSHT20", "l");
  leg_r08_ratio_all->Draw();
  g_r08_sphenix_result->Draw("2 same");
  g_r08_sphenix_result_point->Draw("P same");
  t_r08jet_shadow->Draw("F same");
  t_r08jet_voge_shadow->Draw("F same");
  spline_r08_voge->Draw("L same");
  spline_r08_nnpdf->Draw("L same");
  spline_r08_msht->Draw("L same");
  TPad *pad_r08_2_all = (TPad*)can_r08_ratio_all->cd(2);
  pad_r08_2_all->SetPad(0, 0, 1, 0.333333);
  pad_r08_2_all->SetTopMargin(0.02);
  pad_r08_2_all->SetLeftMargin(0.15);
  pad_r08_2_all->SetBottomMargin(0.25);
  pad_r08_2_all->SetRightMargin(0.05);
  TH2F *frame_r08_ratio_all2 = new TH2F("frame_r08_ratio_all2", "", 10, calibptbins[0], calibptbins[calibnpt], 2, 0., 1.6);
  frame_r08_ratio_all2->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_r08_ratio_all2->GetYaxis()->SetTitle("Ratio");
  frame_r08_ratio_all2->GetXaxis()->SetTitleOffset(0.98);
  frame_r08_ratio_all2->GetYaxis()->SetTitleOffset(0.5);
  frame_r08_ratio_all2->GetXaxis()->SetLabelSize(0.045);
  frame_r08_ratio_all2->GetYaxis()->SetLabelSize(0.045);
  frame_r08_ratio_all2->GetXaxis()->CenterTitle();
  frame_r08_ratio_all2->GetYaxis()->CenterTitle();
  frame_r08_ratio_all2->GetYaxis()->SetTitleOffset(frame_r08_ratio_all->GetYaxis()->GetTitleOffset()*1/2.);
  frame_r08_ratio_all2->GetYaxis()->SetLabelOffset(frame_r08_ratio_all->GetYaxis()->GetLabelOffset()*1/2.);
  frame_r08_ratio_all2->GetXaxis()->SetLabelSize(frame_r08_ratio_all->GetXaxis()->GetLabelSize()*2);
  frame_r08_ratio_all2->GetYaxis()->SetLabelSize(frame_r08_ratio_all->GetYaxis()->GetLabelSize()*2);
  frame_r08_ratio_all2->GetXaxis()->SetTitleSize(frame_r08_ratio_all->GetXaxis()->GetTitleSize()*2);
  frame_r08_ratio_all2->GetYaxis()->SetTitleSize(frame_r08_ratio_all->GetYaxis()->GetTitleSize()*2);
  frame_r08_ratio_all2->Draw();
  g_r08_sphenix_ratio_all->Draw("2 same");
  g_r08_sphenix_ratio_all->Draw("P same");
  g_r08_sphenix_ratio_all_point->Draw("P same");
  t_r08jet_ratio_shadow->Draw("F same");
  t_r08jet_voge_ratioshadow->Draw("F same");
  spline_r08_nnpdf_ratio->Draw("L same");
  spline_r08_msht_ratio->Draw("L same");
  TLine *line_r08_ratio_all = new TLine(calibptbins[0], 1, calibptbins[calibnpt], 1);
  line_r08_ratio_all->SetLineStyle(3);
  line_r08_ratio_all->Draw("same");
  can_r08_ratio_all->SaveAs("figure/spectrum_r08_vogeratio_all.png");

  //TCanvas *can_ratio_zvertex60 = new TCanvas("can_ratio_zvertex60", "", 850, 1091);
  //can_ratio_zvertex60->Divide(1, 2);
  //gStyle->SetPalette(57);
  //TPad *pad_1_zvertex60 = (TPad*)can_ratio_zvertex60->cd(1);
  //pad_1_zvertex60->SetPad(0, 0.33333, 1, 1);
  //pad_1_zvertex60->SetTopMargin(0.03);
  //pad_1_zvertex60->SetLeftMargin(0.15);
  //pad_1_zvertex60->SetBottomMargin(0.035);
  //pad_1_zvertex60->SetRightMargin(0.05);
  //pad_1_zvertex60->SetLogy();
  //TH2F *frame_ratio_zvertex60 = new TH2F("frame_ratio_zvertex60", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r04_sphenix_zvertex60_result->GetN(), g_r04_sphenix_zvertex60_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r04_sphenix_zvertex60_result->GetN(), g_r04_sphenix_zvertex60_result->GetY())));
  //frame_ratio_zvertex60->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  //frame_ratio_zvertex60->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  //frame_ratio_zvertex60->GetXaxis()->SetTitleOffset(0.98);
  //frame_ratio_zvertex60->GetYaxis()->SetTitleOffset(1.15);
  //frame_ratio_zvertex60->GetXaxis()->SetLabelSize(0.045);
  //frame_ratio_zvertex60->GetYaxis()->SetLabelSize(0.045);
  //frame_ratio_zvertex60->GetXaxis()->SetLabelOffset(2);
  //frame_ratio_zvertex60->GetXaxis()->CenterTitle();
  //frame_ratio_zvertex60->GetYaxis()->CenterTitle();
  //frame_ratio_zvertex60->Draw();
  //myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  //myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200 GeV", 0.04);
  //myText(0.6, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.04);
  //myText(0.6, 0.75, 1, "|#eta^{jet}| < 0.7", 0.04);
  //TLegend *leg_ratio_zvertex60 = new TLegend(0.15, 0.08, 0.9, 0.33);
  //leg_ratio_zvertex60->SetBorderSize(0);
  //leg_ratio_zvertex60->SetFillStyle(0);
  //leg_ratio_zvertex60->SetTextSize(0.03);
  //leg_ratio_zvertex60->SetTextFont(42);
  //leg_ratio_zvertex60->AddEntry(g_r04_sphenix_zvertex60_result, "sPHENIX Run 2024 data", "pf");
  //leg_ratio_zvertex60->AddEntry(tjet_shadow, "NLO pQCD (No Hadronization)", "f");
  //leg_ratio_zvertex60->AddEntry(tjet_shadow, "NLO(parton) + CT18 PDF", "");
  //leg_ratio_zvertex60->Draw();
  //g_r04_sphenix_zvertex60_result->Draw("2 same");
  //g_r04_sphenix_zvertex60_result->Draw("P same");
  //g_r04_sphenix_zvertex60_result_point->Draw("P same");
  //t_r04jet_shadow->Draw("F same");
  //TPad *pad_2_zvertex60 = (TPad*)can_ratio_zvertex60->cd(2);
  //pad_2_zvertex60->SetPad(0, 0, 1, 0.333333);
  //pad_2_zvertex60->SetTopMargin(0.02);
  //pad_2_zvertex60->SetLeftMargin(0.15);
  //pad_2_zvertex60->SetBottomMargin(0.25);
  //pad_2_zvertex60->SetRightMargin(0.05);
  //TH2F *frame_ratio_zvertex602 = new TH2F("frame_ratio_zvertex602", "", 10, calibptbins[0], calibptbins[calibnpt], 2, 0., 1.1);
  //frame_ratio_zvertex602->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  //frame_ratio_zvertex602->GetYaxis()->SetTitle("Ratio");
  //frame_ratio_zvertex602->GetXaxis()->SetTitleOffset(0.98);
  //frame_ratio_zvertex602->GetYaxis()->SetTitleOffset(0.5);
  //frame_ratio_zvertex602->GetXaxis()->SetLabelSize(0.045);
  //frame_ratio_zvertex602->GetYaxis()->SetLabelSize(0.045);
  //frame_ratio_zvertex602->GetXaxis()->CenterTitle();
  //frame_ratio_zvertex602->GetYaxis()->CenterTitle();
  //frame_ratio_zvertex602->GetYaxis()->SetTitleOffset(frame_ratio_zvertex60->GetYaxis()->GetTitleOffset()*1/2.);
  //frame_ratio_zvertex602->GetYaxis()->SetLabelOffset(frame_ratio_zvertex60->GetYaxis()->GetLabelOffset()*1/2.);
  //frame_ratio_zvertex602->GetXaxis()->SetLabelSize(frame_ratio_zvertex60->GetXaxis()->GetLabelSize()*2);
  //frame_ratio_zvertex602->GetYaxis()->SetLabelSize(frame_ratio_zvertex60->GetYaxis()->GetLabelSize()*2);
  //frame_ratio_zvertex602->GetXaxis()->SetTitleSize(frame_ratio_zvertex60->GetXaxis()->GetTitleSize()*2);
  //frame_ratio_zvertex602->GetYaxis()->SetTitleSize(frame_ratio_zvertex60->GetYaxis()->GetTitleSize()*2);
  //frame_ratio_zvertex602->Draw();
  //g_sphenix_ratio_zvertex60->Draw("2 same");
  //g_sphenix_ratio_zvertex60->Draw("P same");
  //g_sphenix_ratio_zvertex60_point->Draw("P same");
  //TLine *line_ratio_zvertex60 = new TLine(calibptbins[0], 1, calibptbins[calibnpt], 1);
  //line_ratio_zvertex60->SetLineStyle(3);
  //line_ratio_zvertex60->Draw("same");
  //can_ratio_zvertex60->SaveAs("figure/spectrum_r04_ratio_zvertex60.png");

  TCanvas *can_DAratio_all = new TCanvas("can_DAratio_all", "", 850, 1091);
  can_DAratio_all->Divide(1, 2);
  gStyle->SetPalette(57);
  TPad *pad_DA1_all = (TPad*)can_DAratio_all->cd(1);
  pad_DA1_all->SetPad(0, 0.33333, 1, 1);
  pad_DA1_all->SetTopMargin(0.03);
  pad_DA1_all->SetLeftMargin(0.15);
  pad_DA1_all->SetBottomMargin(0.035);
  pad_DA1_all->SetRightMargin(0.05);
  pad_DA1_all->SetLogy();
  TH2F *frame_DAratio_all = new TH2F("frame_DAratio_all", "", 10, calibptbins[0], calibptbins[calibnpt], 10, low_edge_yaxis(TMath::MinElement(g_r04_sphenix_result->GetN(), g_r04_sphenix_result->GetY())), high_edge_yaxis(TMath::MaxElement(g_r04_sphenix_result->GetN(), g_r04_sphenix_result->GetY())));
  frame_DAratio_all->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_DAratio_all->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_DAratio_all->GetXaxis()->SetTitleOffset(0.98);
  frame_DAratio_all->GetYaxis()->SetTitleOffset(1.15);
  frame_DAratio_all->GetXaxis()->SetLabelSize(0.045);
  frame_DAratio_all->GetYaxis()->SetLabelSize(0.045);
  frame_DAratio_all->GetXaxis()->SetLabelOffset(2);
  frame_DAratio_all->GetXaxis()->CenterTitle();
  frame_DAratio_all->GetYaxis()->CenterTitle();
  frame_DAratio_all->Draw();
  myText(0.6, 0.9, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  myText(0.6, 0.85, 1, "p+p #sqrt{s} = 200 GeV", 0.04);
  myText(0.6, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.04);
  myText(0.6, 0.75, 1, "|#eta^{jet}| < 0.7", 0.04);
  TLegend *leg_DAratio_all = new TLegend(0.15, 0.08, 0.9, 0.33);
  leg_DAratio_all->SetBorderSize(0);
  leg_DAratio_all->SetFillStyle(0);
  leg_DAratio_all->SetTextSize(0.03);
  leg_DAratio_all->SetTextFont(42);
  leg_DAratio_all->AddEntry(g_r04_sphenix_result, "sPHENIX Run 2024 data (DEFAULT)", "pf");
  leg_DAratio_all->AddEntry(g_r04_sphenix_zvertex60_result, "sPHENIX Run 2024 data (ALTZ)", "pf");
  leg_DAratio_all->Draw();
  g_r04_sphenix_zvertex60_result->SetMarkerColor(kRed+1);
  g_r04_sphenix_zvertex60_result->SetLineColor(kRed+1);
  g_r04_sphenix_zvertex60_result->SetFillColorAlpha(kRed-4, 0.5);
  g_r04_sphenix_zvertex60_result_point->SetMarkerColor(kRed+1);
  g_r04_sphenix_zvertex60_result_point->SetLineColor(kRed+1);
  g_r04_sphenix_zvertex60_result->Draw("2 same");
  g_r04_sphenix_zvertex60_result->Draw("P same");
  g_r04_sphenix_zvertex60_result_point->Draw("P same");
  g_r04_sphenix_result->Draw("2 same");
  g_r04_sphenix_result->Draw("P same");
  g_r04_sphenix_result_point->Draw("P same");
  TPad *pad_DA2_all = (TPad*)can_DAratio_all->cd(2);
  pad_DA2_all->SetPad(0, 0, 1, 0.333333);
  pad_DA2_all->SetTopMargin(0.02);
  pad_DA2_all->SetLeftMargin(0.15);
  pad_DA2_all->SetBottomMargin(0.25);
  pad_DA2_all->SetRightMargin(0.05);
  TH2F *frame_DAratio_all2 = new TH2F("frame_DAratio_all2", "", 10, calibptbins[0], calibptbins[calibnpt], 2, 0.0, 2.0);
  frame_DAratio_all2->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_DAratio_all2->GetYaxis()->SetTitle("Ratio");
  frame_DAratio_all2->GetXaxis()->SetTitleOffset(0.98);
  frame_DAratio_all2->GetYaxis()->SetTitleOffset(0.5);
  frame_DAratio_all2->GetXaxis()->SetLabelSize(0.045);
  frame_DAratio_all2->GetYaxis()->SetLabelSize(0.045);
  frame_DAratio_all2->GetXaxis()->CenterTitle();
  frame_DAratio_all2->GetYaxis()->CenterTitle();
  frame_DAratio_all2->GetYaxis()->SetTitleOffset(frame_DAratio_all->GetYaxis()->GetTitleOffset()*1/2.);
  frame_DAratio_all2->GetYaxis()->SetLabelOffset(frame_DAratio_all->GetYaxis()->GetLabelOffset()*1/2.);
  frame_DAratio_all2->GetXaxis()->SetLabelSize(frame_DAratio_all->GetXaxis()->GetLabelSize()*2);
  frame_DAratio_all2->GetYaxis()->SetLabelSize(frame_DAratio_all->GetYaxis()->GetLabelSize()*2);
  frame_DAratio_all2->GetXaxis()->SetTitleSize(frame_DAratio_all->GetXaxis()->GetTitleSize()*2);
  frame_DAratio_all2->GetYaxis()->SetTitleSize(frame_DAratio_all->GetYaxis()->GetTitleSize()*2);
  frame_DAratio_all2->Draw();
  g_sphenix_DAratio_result->Draw("2 same");
  g_sphenix_DAratio_result_point->Draw("P same");
  TLine *line_DAratio_all = new TLine(calibptbins[0], 1, calibptbins[calibnpt], 1);
  line_DAratio_all->SetLineStyle(3);
  line_DAratio_all->Draw("same");
  can_DAratio_all->SaveAs("figure/default_altz_spectrum_ratio.png");
}
