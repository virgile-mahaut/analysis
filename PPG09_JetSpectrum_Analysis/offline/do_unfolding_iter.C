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
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "unfold_Def.h"
#include "/sphenix/user/hanpuj/CaloDataAna24_skimmed/src/draw_template.C"              
#include "BootstrapGenerator/BootstrapGenerator.h"
#include "BootstrapGenerator/TH1DBootstrap.h"
#include "BootstrapGenerator/TH2DBootstrap.h"

R__LOAD_LIBRARY(libBootstrapGenerator.so)

void do_unfolding_iteration_check(TH2D* h_respmatrix, TH1D* h_spectrum, string figure_name, int& optimal_iter) {
  TH1D* h_meas = (TH1D*)h_respmatrix->ProjectionX(Form("h_meas_%s", figure_name.c_str()));
  TH1D* h_truth = (TH1D*)h_respmatrix->ProjectionY(Form("h_truth_%s", figure_name.c_str()));
  const int niter = 20;
  TH1D* h_unfolded[niter];
  for (int i = 0; i < niter; ++i) {
    RooUnfoldResponse* response = new RooUnfoldResponse(h_meas, h_truth, h_respmatrix, "response", "");
    RooUnfoldBayes unfold(response, h_spectrum, i + 1);
    h_unfolded[i] = (TH1D*)unfold.Hunfold(RooUnfold::kErrors);
    h_unfolded[i]->SetName(Form("h_unfolded_iter_%d", i + 1));
    delete response;
  }

  TGraph* g_var_iter = new TGraph(niter-1);
  TGraph* g_var_error = new TGraph(niter-1);
  TGraph* g_var_total = new TGraph(niter-1);
  for (int i = 1; i < niter; ++i) {
    double var_iter = 0;
    double var_error = 0;
    double var_total = 0;
    for (int j = 1+n_underflowbin; j <= h_unfolded[i]->GetNbinsX()-n_overflowbin; ++j) {
      double diff = h_unfolded[i]->GetBinContent(j) - h_unfolded[i-1]->GetBinContent(j);
      double cont = (h_unfolded[i]->GetBinContent(j) + h_unfolded[i-1]->GetBinContent(j)) / 2.0;
      var_iter += diff * diff / (double)(cont * cont);
      var_error += h_unfolded[i]->GetBinError(j) * h_unfolded[i]->GetBinError(j) / (double)(h_unfolded[i]->GetBinContent(j) * h_unfolded[i]->GetBinContent(j));
    }
    var_iter = std::sqrt(var_iter / (double)(h_unfolded[i]->GetNbinsX()-n_overflowbin-n_underflowbin));
    var_error = std::sqrt(var_error / (double)(h_unfolded[i]->GetNbinsX()-n_overflowbin-n_underflowbin));
    var_total = std::sqrt(var_iter * var_iter + var_error * var_error);

    g_var_iter->SetPoint(i-1, i+1, var_iter);
    g_var_error->SetPoint(i-1, i+1, var_error);
    g_var_total->SetPoint(i-1, i+1, var_total);
  }

  // Find optimal iteration (where total variation is minimized)
  int nPoints = g_var_total->GetN();
  double* yValues = g_var_total->GetY();
  double* xValues = g_var_total->GetX();
  int minIndex = TMath::LocMin(nPoints, yValues);
  optimal_iter = std::round(xValues[minIndex]);

  TCanvas* can = new TCanvas(figure_name.c_str(), "Unfolding Iteration Check", 800, 600);
  g_var_iter->SetMarkerStyle(20);
  g_var_iter->SetMarkerColor(kBlue);
  g_var_iter->SetLineColor(kBlue);
  g_var_error->SetMarkerStyle(20);
  g_var_error->SetMarkerColor(kRed);
  g_var_error->SetLineColor(kRed);
  g_var_total->SetMarkerStyle(20);
  g_var_total->SetMarkerColor(kBlack);
  g_var_total->SetLineColor(kBlack);
  g_var_total->SetMinimum(0);
  g_var_total->GetXaxis()->SetTitle("Number of Iterations");
  g_var_total->GetYaxis()->SetTitle("RMS of Relative Variation");
  g_var_total->SetTitle("Unfolding Iteration Optimization");
  g_var_total->Draw("AP");
  g_var_iter->Draw("P same");
  g_var_error->Draw("P same");
  TLegend* leg = new TLegend(0.55, 0.65, 0.85, 0.85);
  leg->AddEntry(g_var_iter, "Iteration variation", "p");
  leg->AddEntry(g_var_error, "Statistical uncertainty", "p");
  leg->AddEntry(g_var_total, "Total (quadrature sum)", "p");
  leg->Draw();
  g_var_total->Draw("AP");
  g_var_iter->Draw("P same");
  g_var_error->Draw("P same");
  g_var_iter->GetXaxis()->SetTitle("Iteration");
  can->SaveAs((figure_name + ".png").c_str());

  delete h_meas;
  delete h_truth;
  for (int i = 0; i < niter; ++i) {
    delete h_unfolded[i];
  }
  delete g_var_iter;
  delete g_var_error;
  delete g_var_total;
  delete leg;
  delete can; 
}

void get_unfolded_spectrum(TH2D* h_respmatrix, TH1D* h_spectrum, TH1D*& h_unfolded, int niter, std::string hist_name, TH1D* h_eff) {
  TH1D* h_meas = (TH1D*)h_respmatrix->ProjectionX(Form("h_meas_%s", hist_name.c_str()));
  TH1D* h_truth = (TH1D*)h_respmatrix->ProjectionY(Form("h_truth_%s", hist_name.c_str()));
  RooUnfoldResponse* response = new RooUnfoldResponse(h_meas, h_truth, h_respmatrix, "response", "");
  RooUnfoldBayes unfold(response, h_spectrum, niter);
  h_unfolded = (TH1D*)unfold.Hunfold(RooUnfold::kErrors);
  h_unfolded->SetName(hist_name.c_str());
  h_unfolded->Divide(h_eff);
  delete h_meas;
  delete h_truth;
  delete response;
}

void get_unfolded_spectrum_stat(TH2DBootstrap* h_respmatrix, TH1DBootstrap* h_spectrum, TH1DBootstrap*& h_unfolded, int niter, std::string hist_name, TH1DBootstrap* h_eff) {
  TH1D* h_meas_nom,* h_truth_nom;
  TH1D** h_meas_rep = new TH1D*[nrep];
  TH1D** h_truth_rep = new TH1D*[nrep];

  TH1D* h_unfolded_nom;
  TH1D** h_unfolded_rep = new TH1D*[nrep];
  
  h_meas_nom = (TH1D*)((TH2D*)h_respmatrix->GetNominal())->ProjectionX("h_meas_nom");
  h_truth_nom = (TH1D*)((TH2D*)h_respmatrix->GetNominal())->ProjectionY("h_truth_nom");
  for(int i=0; i<nrep; ++i) {
    h_meas_rep[i] = (TH1D*)((TH2D*)h_respmatrix->GetReplica(i))->ProjectionX(("h_meas_rep_"+to_string(i)).c_str());
    h_truth_rep[i] = (TH1D*)((TH2D*)h_respmatrix->GetReplica(i))->ProjectionY(("h_truth_rep_"+to_string(i)).c_str());
  }

  RooUnfoldResponse* response = new RooUnfoldResponse(h_meas_nom, h_truth_nom, (TH2D*)h_respmatrix->GetNominal(), "response", "");
  RooUnfoldBayes unfold(response, (TH1D*)h_spectrum->GetNominal(), niter);
  h_unfolded_nom = (TH1D*)unfold.Hunfold(RooUnfold::kErrors);
  h_unfolded_nom->SetName((hist_name+"_nom").c_str());
  h_unfolded_nom->Divide((TH1D*)h_eff->GetNominal());
  delete response;
  for(int i=0; i<nrep; ++i) {
    RooUnfoldResponse* response = new RooUnfoldResponse(h_meas_rep[i], h_truth_rep[i], (TH2D*)h_respmatrix->GetReplica(i),"response","");
    RooUnfoldBayes unfold(response, (TH1D*)h_spectrum->GetReplica(i), niter);
    h_unfolded_rep[i] = (TH1D*)unfold.Hunfold(RooUnfold::kErrors);
    h_unfolded_rep[i]->SetName((hist_name+"_rep_"+to_string(i)).c_str());
    h_unfolded_rep[i]->Divide((TH1D*)h_eff->GetReplica(i));
    delete response;
  }

  h_unfolded = new TH1DBootstrap(hist_name.c_str(), hist_name.c_str(), h_unfolded_nom, h_unfolded_rep, nrep);
  h_unfolded->SetErrBootstrapRMS();

  for(int i=0; i<nrep; ++i) {
    delete h_truth_rep[i];
    delete h_meas_rep[i];
    delete h_unfolded_rep[i];
  }
  delete h_meas_nom;
  delete h_truth_nom;
  delete h_unfolded_nom;
  delete[] h_truth_rep;
  delete[] h_meas_rep;
  delete[] h_unfolded_rep;
}

////////////////////////////////////////// Main Function //////////////////////////////////////////
void do_unfolding_iter(int radius_index = 4) {
  //********** General Set up **********//
  float jet_radius = 0.1 * radius_index;
 
  //********** Files **********//
  TFile *f_out = new TFile(Form("output_unfolded_r0%d.root", radius_index), "RECREATE");
  TFile *f_data = new TFile(Form("output_data_puritycorr_r0%d.root", radius_index), "READ");
  TFile *f_in_rm = new TFile(Form("output_sim_r0%d.root", radius_index), "READ");
  TFile *f_efficiency = new TFile(Form("output_purityefficiency_r0%d.root", radius_index), "READ");

  ofstream log_file(Form("unfolding_iteration_check_r0%d.txt", radius_index));

  // Get data spectrum
  TH1D* h_calibjet_pt_puritycorr_all = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_all");
  TH1DBootstrap* h_calibjet_pt_puritycorr_all_stat = (TH1DBootstrap*)f_data->Get("h_calibjet_pt_puritycorr_all_stat");
  TH1D* h_calibjet_pt_puritycorr_all_unreweighted = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_all_unreweighted");
  TH1D* h_calibjet_pt_puritycorr_all_jesup = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_all_jesup");
  TH1D* h_calibjet_pt_puritycorr_all_jesdown = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_all_jesdown");
  TH1D* h_calibjet_pt_puritycorr_all_jerup = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_all_jerup");
  TH1D* h_calibjet_pt_puritycorr_all_jerdown = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_all_jerdown");
  TH1D* h_calibjet_pt_puritycorr_all_jetup = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_all_jetup");
  TH1D* h_calibjet_pt_puritycorr_all_jetdown = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_all_jetdown");

  TH1D* h_calibjet_pt_puritycorr_zvertex60 = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_zvertex60");
  TH1DBootstrap* h_calibjet_pt_puritycorr_zvertex60_stat = (TH1DBootstrap*)f_data->Get("h_calibjet_pt_puritycorr_zvertex60_stat");
  TH1D* h_calibjet_pt_puritycorr_zvertex60_unreweighted = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_zvertex60_unreweighted");
  TH1D* h_calibjet_pt_puritycorr_zvertex60_jesup = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_zvertex60_jesup");
  TH1D* h_calibjet_pt_puritycorr_zvertex60_jesdown = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_zvertex60_jesdown");
  TH1D* h_calibjet_pt_puritycorr_zvertex60_jerup = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_zvertex60_jerup");
  TH1D* h_calibjet_pt_puritycorr_zvertex60_jerdown = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_zvertex60_jerdown");
  TH1D* h_calibjet_pt_puritycorr_zvertex60_jetup = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_zvertex60_jetup");
  TH1D* h_calibjet_pt_puritycorr_zvertex60_jetdown = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr_zvertex60_jetdown");

  // Response matrix
  TH2D* h_respmatrix_all = (TH2D*)f_in_rm->Get("h_respmatrix_all");
  TH2DBootstrap* h_respmatrix_all_stat = (TH2DBootstrap*)f_in_rm->Get("h_respmatrix_all_stat");
  TH2D* h_respmatrix_all_unreweighted = (TH2D*)f_in_rm->Get("h_respmatrix_all_unreweighted");
  TH2D* h_respmatrix_all_jesup = (TH2D*)f_in_rm->Get("h_respmatrix_all_jesup");
  TH2D* h_respmatrix_all_jesdown = (TH2D*)f_in_rm->Get("h_respmatrix_all_jesdown");
  TH2D* h_respmatrix_all_jerup = (TH2D*)f_in_rm->Get("h_respmatrix_all_jerup");
  TH2D* h_respmatrix_all_jerdown = (TH2D*)f_in_rm->Get("h_respmatrix_all_jerdown");
  TH2D* h_respmatrix_all_jetup = (TH2D*)f_in_rm->Get("h_respmatrix_all_jetup");
  TH2D* h_respmatrix_all_jetdown = (TH2D*)f_in_rm->Get("h_respmatrix_all_jetdown");

  TH2D* h_respmatrix_zvertex60 = (TH2D*)f_in_rm->Get("h_respmatrix_zvertex60");
  TH2DBootstrap* h_respmatrix_zvertex60_stat = (TH2DBootstrap*)f_in_rm->Get("h_respmatrix_zvertex60_stat");
  TH2D* h_respmatrix_zvertex60_unreweighted = (TH2D*)f_in_rm->Get("h_respmatrix_zvertex60_unreweighted");
  TH2D* h_respmatrix_zvertex60_jesup = (TH2D*)f_in_rm->Get("h_respmatrix_zvertex60_jesup");
  TH2D* h_respmatrix_zvertex60_jesdown = (TH2D*)f_in_rm->Get("h_respmatrix_zvertex60_jesdown");
  TH2D* h_respmatrix_zvertex60_jerup = (TH2D*)f_in_rm->Get("h_respmatrix_zvertex60_jerup");
  TH2D* h_respmatrix_zvertex60_jerdown = (TH2D*)f_in_rm->Get("h_respmatrix_zvertex60_jerdown");
  TH2D* h_respmatrix_zvertex60_jetup = (TH2D*)f_in_rm->Get("h_respmatrix_zvertex60_jetup");
  TH2D* h_respmatrix_zvertex60_jetdown = (TH2D*)f_in_rm->Get("h_respmatrix_zvertex60_jetdown");

  // Get efficiency
  TH1D* h_efficiency_all = (TH1D*)f_efficiency->Get("h_efficiency_all");
  TH1DBootstrap* h_efficiency_all_stat = (TH1DBootstrap*)f_efficiency->Get("h_efficiency_all_stat");
  TH1D* h_efficiency_all_unreweighted = (TH1D*)f_efficiency->Get("h_efficiency_all_unreweighted");
  TH1D* h_efficiency_all_jesup = (TH1D*)f_efficiency->Get("h_efficiency_all_jesup");
  TH1D* h_efficiency_all_jesdown = (TH1D*)f_efficiency->Get("h_efficiency_all_jesdown");
  TH1D* h_efficiency_all_jerup = (TH1D*)f_efficiency->Get("h_efficiency_all_jerup");
  TH1D* h_efficiency_all_jerdown = (TH1D*)f_efficiency->Get("h_efficiency_all_jerdown");
  TH1D* h_efficiency_all_jetup = (TH1D*)f_efficiency->Get("h_efficiency_all_jetup");
  TH1D* h_efficiency_all_jetdown = (TH1D*)f_efficiency->Get("h_efficiency_all_jetdown");

  TH1D* h_efficiency_zvertex60 = (TH1D*)f_efficiency->Get("h_efficiency_zvertex60");
  TH1DBootstrap* h_efficiency_zvertex60_stat = (TH1DBootstrap*)f_efficiency->Get("h_efficiency_zvertex60_stat");
  TH1D* h_efficiency_zvertex60_unreweighted = (TH1D*)f_efficiency->Get("h_efficiency_zvertex60_unreweighted");
  TH1D* h_efficiency_zvertex60_jesup = (TH1D*)f_efficiency->Get("h_efficiency_zvertex60_jesup");
  TH1D* h_efficiency_zvertex60_jesdown = (TH1D*)f_efficiency->Get("h_efficiency_zvertex60_jesdown");
  TH1D* h_efficiency_zvertex60_jerup = (TH1D*)f_efficiency->Get("h_efficiency_zvertex60_jerup");
  TH1D* h_efficiency_zvertex60_jerdown = (TH1D*)f_efficiency->Get("h_efficiency_zvertex60_jerdown");
  TH1D* h_efficiency_zvertex60_jetup = (TH1D*)f_efficiency->Get("h_efficiency_zvertex60_jetup");
  TH1D* h_efficiency_zvertex60_jetdown = (TH1D*)f_efficiency->Get("h_efficiency_zvertex60_jetdown");

  //********** Response Matrix **********//
  int optimal_iter = 1;

  do_unfolding_iteration_check(h_respmatrix_all, h_calibjet_pt_puritycorr_all, Form("figure_unfolding/unfolding_iteration_r0%d_all", radius_index), optimal_iter);
  TH1D* h_unfold_all; get_unfolded_spectrum(h_respmatrix_all, h_calibjet_pt_puritycorr_all, h_unfold_all, optimal_iter, "h_unfold_all", h_efficiency_all);
  TH1DBootstrap* h_unfold_all_stat; get_unfolded_spectrum_stat(h_respmatrix_all_stat, h_calibjet_pt_puritycorr_all_stat, h_unfold_all_stat, optimal_iter, "h_unfold_all_stat", h_efficiency_all_stat);
  log_file << "Optimal iteration for r0" << radius_index << "_all: " << optimal_iter << std::endl;
  TH1D* h_unfold_all_unfoldunc; get_unfolded_spectrum(h_respmatrix_all, h_calibjet_pt_puritycorr_all, h_unfold_all_unfoldunc, optimal_iter+1, "h_unfold_all_unfoldunc", h_efficiency_all);
  log_file << "Optimal iteration for r0" << radius_index << "_all_unfoldunc: " << optimal_iter+1 << std::endl;
  TH1D* h_unfold_all_jesup; get_unfolded_spectrum(h_respmatrix_all_jesup, h_calibjet_pt_puritycorr_all_jesup, h_unfold_all_jesup, optimal_iter, "h_unfold_all_jesup", h_efficiency_all_jesup);
  TH1D* h_unfold_all_jesdown; get_unfolded_spectrum(h_respmatrix_all_jesdown, h_calibjet_pt_puritycorr_all_jesdown, h_unfold_all_jesdown, optimal_iter, "h_unfold_all_jesdown", h_efficiency_all_jesdown);
  TH1D* h_unfold_all_jerup; get_unfolded_spectrum(h_respmatrix_all_jerup, h_calibjet_pt_puritycorr_all_jerup, h_unfold_all_jerup, optimal_iter, "h_unfold_all_jerup", h_efficiency_all_jerup);
  TH1D* h_unfold_all_jerdown; get_unfolded_spectrum(h_respmatrix_all_jerdown, h_calibjet_pt_puritycorr_all_jerdown, h_unfold_all_jerdown, optimal_iter, "h_unfold_all_jerdown", h_efficiency_all_jerdown);
  TH1D* h_unfold_all_jetup; get_unfolded_spectrum(h_respmatrix_all_jetup, h_calibjet_pt_puritycorr_all_jetup, h_unfold_all_jetup, optimal_iter, "h_unfold_all_jetup", h_efficiency_all_jetup);
  TH1D* h_unfold_all_jetdown; get_unfolded_spectrum(h_respmatrix_all_jetdown, h_calibjet_pt_puritycorr_all_jetdown, h_unfold_all_jetdown, optimal_iter, "h_unfold_all_jetdown", h_efficiency_all_jetdown);
  do_unfolding_iteration_check(h_respmatrix_all_unreweighted, h_calibjet_pt_puritycorr_all_unreweighted, Form("figure_unfolding/unfolding_iteration_r0%d_all_unreweighted", radius_index), optimal_iter);
  TH1D* h_unfold_all_unreweighted; get_unfolded_spectrum(h_respmatrix_all_unreweighted, h_calibjet_pt_puritycorr_all_unreweighted, h_unfold_all_unreweighted, optimal_iter, "h_unfold_all_unreweighted", h_efficiency_all_unreweighted);
  log_file << "Optimal iteration for r0" << radius_index << "_all_unreweighted: " << optimal_iter << std::endl;

  do_unfolding_iteration_check(h_respmatrix_zvertex60, h_calibjet_pt_puritycorr_zvertex60, Form("figure_unfolding/unfolding_iteration_r0%d_zvertex60", radius_index), optimal_iter);
  TH1D* h_unfold_zvertex60; get_unfolded_spectrum(h_respmatrix_zvertex60, h_calibjet_pt_puritycorr_zvertex60, h_unfold_zvertex60, optimal_iter, "h_unfold_zvertex60", h_efficiency_zvertex60);
  TH1DBootstrap* h_unfold_zvertex60_stat; get_unfolded_spectrum_stat(h_respmatrix_zvertex60_stat, h_calibjet_pt_puritycorr_zvertex60_stat, h_unfold_zvertex60_stat, optimal_iter, "h_unfold_zvertex60_stat", h_efficiency_zvertex60_stat);
  log_file << "Optimal iteration for r0" << radius_index << "_zvertex60: " << optimal_iter << std::endl;
  TH1D* h_unfold_zvertex60_unfoldunc; get_unfolded_spectrum(h_respmatrix_zvertex60, h_calibjet_pt_puritycorr_zvertex60, h_unfold_zvertex60_unfoldunc, optimal_iter+1, "h_unfold_zvertex60_unfoldunc", h_efficiency_zvertex60);
  log_file << "Optimal iteration for r0" << radius_index << "_zvertex60_unfoldunc: " << optimal_iter+1 << std::endl;
  TH1D* h_unfold_zvertex60_jesup; get_unfolded_spectrum(h_respmatrix_zvertex60_jesup, h_calibjet_pt_puritycorr_zvertex60_jesup, h_unfold_zvertex60_jesup, optimal_iter, "h_unfold_zvertex60_jesup", h_efficiency_zvertex60_jesup);
  TH1D* h_unfold_zvertex60_jesdown; get_unfolded_spectrum(h_respmatrix_zvertex60_jesdown, h_calibjet_pt_puritycorr_zvertex60_jesdown, h_unfold_zvertex60_jesdown, optimal_iter, "h_unfold_zvertex60_jesdown", h_efficiency_zvertex60_jesdown);
  TH1D* h_unfold_zvertex60_jerup; get_unfolded_spectrum(h_respmatrix_zvertex60_jerup, h_calibjet_pt_puritycorr_zvertex60_jerup, h_unfold_zvertex60_jerup, optimal_iter, "h_unfold_zvertex60_jerup", h_efficiency_zvertex60_jerup);
  TH1D* h_unfold_zvertex60_jerdown; get_unfolded_spectrum(h_respmatrix_zvertex60_jerdown, h_calibjet_pt_puritycorr_zvertex60_jerdown, h_unfold_zvertex60_jerdown, optimal_iter, "h_unfold_zvertex60_jerdown", h_efficiency_zvertex60_jerdown);
  TH1D* h_unfold_zvertex60_jetup; get_unfolded_spectrum(h_respmatrix_zvertex60_jetup, h_calibjet_pt_puritycorr_zvertex60_jetup, h_unfold_zvertex60_jetup, optimal_iter, "h_unfold_zvertex60_jetup", h_efficiency_zvertex60_jetup);
  TH1D* h_unfold_zvertex60_jetdown; get_unfolded_spectrum(h_respmatrix_zvertex60_jetdown, h_calibjet_pt_puritycorr_zvertex60_jetdown, h_unfold_zvertex60_jetdown, optimal_iter, "h_unfold_zvertex60_jetdown", h_efficiency_zvertex60_jetdown);
  do_unfolding_iteration_check(h_respmatrix_zvertex60_unreweighted, h_calibjet_pt_puritycorr_zvertex60_unreweighted, Form("figure_unfolding/unfolding_iteration_r0%d_zvertex60_unreweighted", radius_index), optimal_iter);
  TH1D* h_unfold_zvertex60_unreweighted; get_unfolded_spectrum(h_respmatrix_zvertex60_unreweighted, h_calibjet_pt_puritycorr_zvertex60_unreweighted, h_unfold_zvertex60_unreweighted, optimal_iter, "h_unfold_zvertex60_unreweighted", h_efficiency_zvertex60_unreweighted);
  log_file << "Optimal iteration for r0" << radius_index << "_zvertex60_unreweighted: " << optimal_iter << std::endl;

  //********** Writing **********//
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();
  h_unfold_all->Write();
  h_unfold_all_stat->Write();
  h_unfold_all_unfoldunc->Write();
  h_unfold_all_unreweighted->Write();
  h_unfold_all_jesup->Write();
  h_unfold_all_jesdown->Write();
  h_unfold_all_jerup->Write();
  h_unfold_all_jerdown->Write();
  h_unfold_all_jetup->Write();
  h_unfold_all_jetdown->Write();

  h_unfold_zvertex60->Write();
  h_unfold_zvertex60_stat->Write();
  h_unfold_zvertex60_unfoldunc->Write();
  h_unfold_zvertex60_unreweighted->Write();
  h_unfold_zvertex60_jesup->Write();
  h_unfold_zvertex60_jesdown->Write();
  h_unfold_zvertex60_jerup->Write();
  h_unfold_zvertex60_jerdown->Write();
  h_unfold_zvertex60_jetup->Write();
  h_unfold_zvertex60_jetdown->Write();

  f_out->Close();
  std::cout << "All done!" << std::endl;
}
