#include "AnNeutralMeson_nano.h"
#include <fun4all/Fun4AllReturnCodes.h>

// Spin DB
#include <uspin/SpinDBContent.h>
#include <uspin/SpinDBContentv1.h>
#include <uspin/SpinDBOutput.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip> // for setprecision
#include <random>
#include <sstream>

AnNeutralMeson_nano::AnNeutralMeson_nano(const std::string &name,
                                         const std::string &inputlistname,
                                         const std::string &inputfiletemplate,
                                         const std::string &outputfiletemplate)
  : SubsysReco(name)
  , inlistname(inputlistname)
  , infiletemplate(inputfiletemplate)
  , treename("tree_diphoton_compact")
  , outfiletemplate(outputfiletemplate)
{
}

AnNeutralMeson_nano::~AnNeutralMeson_nano()
{
}

int AnNeutralMeson_nano::Init(PHCompositeNode *)
{
  std::ifstream inlistfile;
  inlistfile.open(inlistname.c_str());

  std::string runLine;
  while (std::getline(inlistfile, runLine))
  {
    int runnumber = std::stoi(runLine);
    runList.push_back(runnumber);
  }
  nRuns = runList.size();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnNeutralMeson_nano::process_event(PHCompositeNode *)
{
  std::stringstream infilename;
  std::stringstream outfilename;
  std::cout << "nRuns = " << nRuns << std::endl;
  for (int irun = 0; irun < nRuns; irun++)
  {
    std::cout << "run " << runList[irun] << std::endl;
    infilename.str("");
    infilename << infiletemplate << runList[irun] << ".root";
    outfilename.str("");
    outfilename << outfiletemplate << runList[irun] << ".root";

    // Book histograms
    BookHistos(outfilename.str());

    // Open TTree
    TFile *infile = TFile::Open(infilename.str().c_str(), "READ");

    if (infile == nullptr)
    {
      std::cerr << "Error. Could not open file " << infilename.str()
                << std::endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
      continue;
    }

    TTree *nanoDST = nullptr;
    infile->GetObject(treename.c_str(), nanoDST);
    if (nanoDST == nullptr)
    {
      std::cerr << "Error: tree " << treename << " was not found in the file "
                << infilename.str() << std::endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
      continue;
    }

    // Optimization for reading
    nanoDST->SetBranchStatus("*", 0);
    nanoDST->SetBranchStatus("diphoton_vertex_z", 1);
    nanoDST->SetBranchStatus("diphoton_bunchnumber", 1);
    nanoDST->SetBranchStatus("diphoton_mass", 1);
    nanoDST->SetBranchStatus("diphoton_eta", 1);
    nanoDST->SetBranchStatus("diphoton_phi", 1);
    nanoDST->SetBranchStatus("diphoton_pt", 1);
    nanoDST->SetBranchStatus("diphoton_xf", 1);

    nanoDST->SetBranchAddress(
        "diphoton_vertex_z", &diphoton_vertex_z);
    nanoDST->SetBranchAddress(
        "diphoton_bunchnumber", &diphoton_bunchnumber);
    nanoDST->SetBranchAddress(
        "diphoton_mass", &diphoton_mass);
    nanoDST->SetBranchAddress(
        "diphoton_eta", &diphoton_eta);
    nanoDST->SetBranchAddress(
        "diphoton_phi", &diphoton_phi);
    nanoDST->SetBranchAddress(
        "diphoton_pt", &diphoton_pt);
    nanoDST->SetBranchAddress(
        "diphoton_xf", &diphoton_xf);

    // Spin pattern from SpinDB
    // Get spin pattern from SpinDB

    // Use the default QA level
    SpinDBOutput spin_out(
        "phnxrc");
    SpinDBContent *spin_cont = new SpinDBContentv1();
    spin_out.StoreDBContent(runList[irun], runList[irun]);
    spin_out.GetDBContentStore(spin_cont, runList[irun]);

    crossingshift = spin_cont->GetCrossingShift();

    for (int i = 0; i < nBunches; i++)
    {
      beamspinpat[0][i] = -1 * spin_cont->GetSpinPatternYellow(i);
      beamspinpat[1][i] = -1 * spin_cont->GetSpinPatternBlue(i);
      // spin pattern at sPHENIX is -1*(CDEV pattern)
      // The spin pattern corresponds to the expected spin at IP12, before the
      // siberian snake
    }

    if (seednb != 0) shuffle_spin_pattern(irun);

    // Loop over all recorded events;
    Long64_t nentries = nanoDST->GetEntries();
    std::cout << "nentries = " << nentries << std::endl;
    for (Long64_t ientry = 0; ientry < nentries; ientry++)
    {
      nanoDST->GetEntry(ientry);

      if (diphoton_pt < pTCutMin || diphoton_pt > pTCutMax) continue;

      if (require_low_vtx_cut) {
        if (std::abs(diphoton_vertex_z) > vertex_min) continue;
      }
      if (require_high_vtx_cut) {
        if (std::abs(diphoton_vertex_z) <= vertex_max) continue;
      }

      if (require_phenix_cut) {
        if (std::abs(diphoton_eta) > 0.35) continue;
      }
      
      // Select the right pt bin:
      int iPt = FindBinBinary(diphoton_pt, pTBins, nPtBins);
      if (!((iPt < nPtBins) && (iPt >= 0)))
        continue;

      // Select the right eta bin:
      int ietaBin = FindBinBinary(diphoton_eta, etaBins, nEtaBins);
      if (!((ietaBin < nEtaBins) && (ietaBin >= 0)))
        continue;

      // Select the right xf bin:
      int ixfBin = FindBinBinary(diphoton_xf, xfBins, nXfBins);
      if (!((ixfBin < nXfBins) && (ixfBin >= 0)))
        continue;

      int izvtxBin = FindBinBinary(diphoton_vertex_z, zvtxBins, nZvtxBins);

      // Select particle and region index;
      int ival = FindBinBinary(diphoton_mass, band_limits,
                               nParticles * (nRegions + 1) * 2);
      
      // Ignore any diphoton which is not within the side band or the peak band
      int iP = -1;            // particle index
      int iR = -1;  // region index
      if (!(ival < 0 || ival % 2 == 1))
      {
        iP = ival / 2 / (nRegions + 1);            // particle index
        iR = (ival / 2 % (nRegions + 1) + 1) % 2;  // region index
      }

      // Store QA
      if (iR == 0 && iP != -1) {
        h_pair_meson_zvtx[iP]->Fill(diphoton_vertex_z);
        h_pair_meson_pt_eta[iP][iPt]->Fill(diphoton_eta);
        h_pair_meson_pt_xf[iP][iPt]->Fill(diphoton_xf);
        h_pair_meson_eta_pt[iP][ietaBin]->Fill(diphoton_pt);
        h_pair_meson_eta_xf[iP][ietaBin]->Fill(diphoton_xf);
        h_pair_meson_xf_pt[iP][ixfBin]->Fill(diphoton_pt);
        h_pair_meson_xf_eta[iP][ixfBin]->Fill(diphoton_eta);

        h_pair_meson_2D_xf_pt[iP]->Fill(diphoton_xf, diphoton_pt);
        h_pair_meson_2D_xf_eta[iP]->Fill(diphoton_xf, diphoton_eta);
        h_pair_meson_2D_xf_zvtx[iP]->Fill(diphoton_xf, diphoton_vertex_z);
        h_pair_meson_2D_pt_eta[iP]->Fill(diphoton_pt, diphoton_eta);
        h_pair_meson_2D_pt_zvtx[iP]->Fill(diphoton_pt, diphoton_vertex_z);
        h_pair_meson_2D_eta_zvtx[iP]->Fill(diphoton_eta, diphoton_vertex_z);
      }
      
      // Compute yield
      float phi_beam = 0;
      for (int iB = 0; iB < nBeams; iB++)
      {
        // Forward vs backward asymmetry
        int iDir = (diphoton_eta * beamDirection[iB] > 0 ? 0 : 1);
        int ietaBinRelative =
            (beamDirection[iB] > 0 ? ietaBin : nEtaBins - 1 - ietaBin);
        int ixfBinRelative =
            (beamDirection[iB] > 0 ? ixfBin : nXfBins - 1 - ixfBin);

        if (require_high_xf_cut) {
          if (ixfBinRelative < 6) continue; // i.e. last most forward pT bins 6 and 7 -> xF > 0.035
        }

        if (require_low_xf_cut) {
          if (ixfBinRelative >= 6) continue; // i.e. last most forward pT bins 6 and 7 -> xF > 0.035
        }

        // Store symmetrized invariant mass distributions
        // Only eta and xF distributions are symmetrized!
        h_pair_mass->Fill(diphoton_mass);
        if (iPt >= 0 && iPt < nPtBins && beamDirection[iB] > 0)
          h_pair_mass_pt[iPt]->Fill(diphoton_mass);
        if (izvtxBin >= 0 && izvtxBin < nZvtxBins && beamDirection[iB] > 0)
          h_pair_mass_zvtx[izvtxBin]->Fill(diphoton_mass);
        if (ietaBinRelative >= 0 && ietaBinRelative < nEtaBins)
          h_pair_mass_eta[ietaBinRelative]->Fill(diphoton_mass);
        if (ixfBinRelative >= 0 && ixfBinRelative < nXfBins)
          h_pair_mass_xf[ixfBinRelative]->Fill(diphoton_mass);

        if (iP == -1 || iR == -1) continue;

        // Store pT, eta and xF values in order to compute the average value of
        // each bin
        if (iP != -1 && iR == 0) // Really the signal band
        {
          if (iPt >= 0 && iPt < nPtBins)
          {
            h_average_pt[iP]->Fill(iPt, diphoton_pt);
            h_norm_pt[iP]->Fill(iPt, 1);
          }
          if (ietaBinRelative >= 0 && ietaBinRelative < nEtaBins)
          {
            h_average_eta[iP]->Fill(ietaBinRelative, (beamDirection[iB] > 0 ? diphoton_eta : -diphoton_eta));
            h_norm_eta[iP]->Fill(ietaBinRelative, 1);
          }
          if (ixfBinRelative >= 0 && ixfBinRelative < nXfBins)
          {
            h_average_xf[iP]->Fill(ixfBinRelative, (beamDirection[iB] > 0 ? diphoton_xf : -diphoton_xf));
            h_norm_xf[iP]->Fill(ixfBinRelative, 1);
          }
        }

        // phi global to phi yellow/blue
        phi_beam = WrapAngle(diphoton_phi + phi_shift[iB]); // i.e. phi_beam between - M_PI and +M_PI

        // select spin value;
        int iBunch = (crossingshift + diphoton_bunchnumber) % nBunches;
        int iS = -1;
        if (beamspinpat[iB][iBunch] == 1)
          iS = 0;
        else if (beamspinpat[iB][iBunch] == -1)
          iS = 1;
        if (iS == -1)
          continue;

        // Determine phi index
        if ((phi_beam < phiMin) || (phi_beam > phiMax))
          continue;

        // Increment kinematic-dependent yields
        h_yield_pt[iB][iP][iR][iPt][iDir][iS]->Fill(phi_beam);
        h_yield_eta[iB][iP][iR][ietaBinRelative][iS]->Fill(phi_beam);
        h_yield_xf[iB][iP][iR][ixfBinRelative][iS]->Fill(phi_beam);

        if (store_bunch_yields) {
          // Different pT thresholds for the pi0 and eta trigger
          if ((trigger_mbd && iP == 0 && iPt <= 1) ||
              (trigger_mbd && iP == 1 && iPt <= 2) ||
              (trigger_photon && iP == 0 && iPt >= 2) ||
              (trigger_photon && iP == 1 && iPt >= 3))
          {
            int iPhi = FindBinDirect(phi_beam, -M_PI, M_PI, nPhiBins);
            array_yield_pt[iB][iP][iR][iPt][iDir][iS][iPhi][iBunch]++;
            array_yield_eta[iB][iP][iR][ietaBinRelative][iS][iPhi][iBunch]++;
            array_yield_xf[iB][iP][iR][ixfBinRelative][iS][iPhi][iBunch]++;
          }
        }
      }
    }
    infile->Close();

    // Store histogram yields
    outfile->cd();
    outfile->Write();
    outfile->Close();
    delete outfile;
    outfile = nullptr; // To prevent double deletion with the destructor

    if (store_bunch_yields)
    {
      outfilename.str("");
      outfilename << outbunchtemplate << runList[irun] << ".bin";
      SaveBunchYields(outfilename.str());
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnNeutralMeson_nano::End(PHCompositeNode *)
{
  return 0;
}

void AnNeutralMeson_nano::BookHistos(const std::string &outputfilename)
{
  // Book histograms
  outfile = new TFile(outputfilename.c_str(), "RECREATE");
  outfile->cd();

  // diphoton invariant mass
  h_pair_mass = new TH1F(
      "h_pair_mass",
      ";M_{#gamma} [GeV];counts", 500, 0, 1);

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

  for (int izvtxBin = 0; izvtxBin < nZvtxBins; izvtxBin++)
  {
    std::stringstream h_pair_mass_name;
    h_pair_mass_name << std::fixed << std::setprecision(0) << "h_pair_mass_zvtx_"
                     << izvtxBin;
    std::stringstream h_pair_mass_title;
    h_pair_mass_title << std::fixed << std::setprecision(2) << "diphoton mass ["
                      << zvtxBins[izvtxBin] << " < x_{F} < " << zvtxBins[izvtxBin + 1]
                      << "];M_{#gamma#gamma};counts";
    h_pair_mass_zvtx[izvtxBin] =
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
  h_average_xf[0] = new TH1F(
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
  h_norm_xf[0] = new TH1F(
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
  h_average_xf[1] = new TH1F(
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
  h_norm_xf[1] = new TH1F(
      "h_norm_xf_eta",
      ";ixf; pT norm", nXfBins, 0, nXfBins);

  for (int iP = 0; iP < 2; iP++) {
    std::stringstream h_name;
    h_name << "h_pair_" << particle[iP] << "_zvtx";
    h_pair_meson_zvtx[iP] = new TH1F(
      h_name.str().c_str(),
      ";z_{vtx} [cm]; Counts / [2 cm]",
      200, -200, 200);
  }

  for (int iP = 0; iP < 2; iP++) {
    {
      std::stringstream h_name;
      h_name << "h_pair_" << particle[iP] << "_2D_xf_pt";
      h_pair_meson_2D_xf_pt[iP] = new TH2F(
        h_name.str().c_str(),
        ";x_{F}; p_{T} [GeV]",
        200, -0.2, 0.2,
        200, 0.0, 20.0);
    }
    {
      std::stringstream h_name;
      h_name << "h_pair_" << particle[iP] << "_2D_xf_eta";
      h_pair_meson_2D_xf_eta[iP] = new TH2F(
        h_name.str().c_str(),
        ";x_{F}; #eta",
        200, -0.2, 0.2,
        200, -2.0, 2.0);
    }
    {
      std::stringstream h_name;
      h_name << "h_pair_" << particle[iP] << "_2D_xf_zvtx";
      h_pair_meson_2D_xf_zvtx[iP] = new TH2F(
        h_name.str().c_str(),
        ";x_{F}; z_{vtx} [cm]",
        200, -0.2, 0.2,
        200, -200, 200);
    }
    {
      std::stringstream h_name;
      h_name << "h_pair_" << particle[iP] << "_2D_pt_eta";
      h_pair_meson_2D_pt_eta[iP] = new TH2F(
        h_name.str().c_str(),
        ";p_{T} [GeV]; #eta",
        200, 0.0, 20.0,
        200, -2.0, 2.0);
    }
    {
      std::stringstream h_name;
      h_name << "h_pair_" << particle[iP] << "_2D_pt_zvtx";
      h_pair_meson_2D_pt_zvtx[iP] = new TH2F(
        h_name.str().c_str(),
        ";p_{T} [GeV]; z_{vtx} [cm]",
        200, 0.0, 20.0,
        200, -200, 200);
    }
    {
      std::stringstream h_name;
      h_name << "h_pair_" << particle[iP] << "_2D_eta_zvtx";
      h_pair_meson_2D_eta_zvtx[iP] = new TH2F(
        h_name.str().c_str(),
        ";#eta; z_{vtx} [cm]",
        200, -2.0, 2.0,
        200, -200, 200);
    }
    
  }
  
  // Kinematic relation pT vs eta vs xF
  // diphoton eta pT-bin by pT-bin
  for (int iP = 0; iP < 2; iP++) {
    for (int iPt = 1; iPt <= 9; iPt++) {
      std::stringstream h_name;
      h_name << "h_pair_" << particle[iP] << "_eta_pt_" << iPt;
      h_pair_meson_pt_eta[iP][iPt-1] = new TH1F(
        h_name.str().c_str(),
        ";#eta; Counts / [20 mrad]",
        200, -2.0, 2.0);
    }
  }

  // Kinematic relation pT vs eta vs xF
  // diphoton xF pT-bin by pT-bin
  for (int iP = 0; iP < 2; iP++) {
    for (int iPt = 1; iPt <= 9; iPt++) {
      std::stringstream h_name;
      h_name << "h_pair_" << particle[iP] << "_xf_pt_" << iPt;
      h_pair_meson_pt_xf[iP][iPt-1] = new TH1F(
        h_name.str().c_str(),
        ";x_{F}; Counts / [2 mrad]",
        200, -0.2, 0.2);
    }
  }

  // Kinematic relation pT vs eta vs xF
  // diphoton pT eta-bin by eta-bin
  for (int iP = 0; iP < 2; iP++) {
    for (int iEta = 1; iEta <= 8; iEta++) {
      std::stringstream h_name;
      h_name << "h_pair_" << particle[iP] << "_pt_eta_" << iEta;
      h_pair_meson_eta_pt[iP][iEta-1] = new TH1F(
        h_name.str().c_str(),
        ";p_{T} [GeV]; Counts / [100 MeV]",
        200, 0.0, 20);
    }
  }

  // Kinematic relation pT vs eta vs xF
  // diphoton xF eta-bin by eta-bin
  for (int iP = 0; iP < 2; iP++) {
    for (int iEta = 1; iEta <= 8; iEta++) {
      std::stringstream h_name;
      h_name << "h_pair_" << particle[iP] << "_xf_eta_" << iEta;
      h_pair_meson_eta_xf[iP][iEta-1] = new TH1F(
        h_name.str().c_str(),
        ";x_{F}; Counts / [2 mrad]",
        200, -0.2, 0.2);
    }
  }

  // Kinematic relation pT vs eta vs xF
  // diphoton pT xf-bin by xf-bin
  for (int iP = 0; iP < 2; iP++) {
    for (int iXf = 1; iXf <= 8; iXf++) {
      std::stringstream h_name;
      h_name << "h_pair_" << particle[iP] << "_pt_xf_" << iXf;
      h_pair_meson_xf_pt[iP][iXf-1] = new TH1F(
        h_name.str().c_str(),
        ";p_{T} [GeV]; Counts / [100 MeV]",
        200, 0.0, 20);
    }
  }

  // Kinematic relation pT vs eta vs xF
  // diphoton eta xf-bin by xf-bin
  for (int iP = 0; iP < 2; iP++) {
    for (int iXf = 1; iXf <= 8; iXf++) {
      std::stringstream h_name;
      h_name << "h_pair_" << particle[iP] << "_eta_xf_" << iXf;
      h_pair_meson_xf_eta[iP][iXf-1] = new TH1F(
        h_name.str().c_str(),
        ";#eta; Counts / [20 mrad]",
        200, -2.0, 2.0);
    }
  }
  
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
            for (int iDir = 0; iDir < nDirections; iDir++)
            {
              std::stringstream h_yield_name;
              h_yield_name << "h_yield_" << beams[iB] << "_" << particle[iP]
                           << "_" << regions[iR] << "_"
                           << "pT_" << iPt << "_"
                           << directions[iDir]
                           << "_" << spins[iS];
              std::stringstream h_yield_title;
              h_yield_title << std::fixed << std::setprecision(1)
                            << "P_{T} #in [" << pTBins[iPt] << ", "
                            << pTBins[iPt + 1] << "] " << directions[iDir]
                            << " " << (iP == 0 ? "(#pi^{0} " : "(#eta ")
                            << regions[iR] << " range);"
                            << "#phi_{" << (iB == 0 ? "yellow" : "blue")
                            << "};counts";
              h_yield_pt[iB][iP][iR][iPt][iDir][iS] =
                new TH1F(h_yield_name.str().c_str(),
                         h_yield_title.str().c_str(), 12, -M_PI, M_PI);
            }
          }
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
}

void AnNeutralMeson_nano::shuffle_spin_pattern(const int irun)
{
  // Shuffle the spin pattern if the seed number is non zero
  // Only consider the first 111/(nBunches = 120) bunches, the last 9 are
  // slots are empty.
  const int nFilledBunches = 111;
  int shuffled_indices[nBunches] = {0};
  std::iota(shuffled_indices, shuffled_indices + nFilledBunches, 0);

  std::mt19937 rng(seednb + 100000 * irun);
  // Pseudo-random generator: Mersenne twister with a period of 2^19937 - 1
  std::shuffle(shuffled_indices, shuffled_indices + nFilledBunches, rng);

  // With a probability 1/2, flip all the bunch spins
  std::uniform_int_distribution<std::mt19937::result_type> dist2(0, 1);
  int factor = 2 * (int) dist2(rng) - 1;
  int beamspinpat_tmp[2][nBunches] = {0};
  for (int ibunch = 0; ibunch < nBunches; ibunch++)
  {
    int ishuffle = shuffled_indices[ibunch];
    beamspinpat_tmp[0][ibunch] = factor * beamspinpat[0][ishuffle];
    beamspinpat_tmp[1][ibunch] = factor * beamspinpat[1][ishuffle];
  }
  for (int ibunch = 0; ibunch < nBunches; ibunch++)
  {
    beamspinpat[0][ibunch] = beamspinpat_tmp[0][ibunch];
    beamspinpat[1][ibunch] = beamspinpat_tmp[1][ibunch];
  }
}

float AnNeutralMeson_nano::WrapAngle(const float phi)
{
  float phi_out = phi;
  while (phi_out < -M_PI)
  {
    phi_out += 2 * M_PI;
  }
  while (phi_out > M_PI)
  {
    phi_out -= 2 * M_PI;
  }
  return phi_out;
}

void AnNeutralMeson_nano::SaveBunchYields(const std::string &outputfilename)
{
  // Binary output
  std::ofstream outbinary(outputfilename, std::ios::binary);
  outbinary.write(reinterpret_cast<const char *>(array_yield_pt),
                  sizeof(array_yield_pt));
  outbinary.write(reinterpret_cast<const char *>(array_yield_eta),
                  sizeof(array_yield_eta));
  outbinary.write(reinterpret_cast<const char *>(array_yield_xf),
                  sizeof(array_yield_xf));
  outbinary.close();

  // Reset arrays to zero
  std::memset(array_yield_pt, 0, sizeof(array_yield_pt));
  std::memset(array_yield_eta, 0, sizeof(array_yield_eta));
  std::memset(array_yield_xf, 0, sizeof(array_yield_xf));
}
