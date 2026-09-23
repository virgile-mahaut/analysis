#include "MC_calo_pythia_pi0.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <globalvertex/GlobalVertex.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawCluster.h>
#include <CLHEP/Vector/ThreeVector.h>

#include "Math/VectorUtil.h"

MC_calo_pythia_pi0::MC_calo_pythia_pi0(const std::string &name,
                                       const std::string &outputfilename,
                                       const int seednb)
  : SubsysReco(name)
  , outfilename(outputfilename)
  , seednumber(seednb)
{
}

MC_calo_pythia_pi0::~MC_calo_pythia_pi0()
{
  delete outfile;
  outfile = nullptr;
}

int MC_calo_pythia_pi0::Init(PHCompositeNode * /*topNode*/)
{
  // Efficiency trigger
  if (require_efficiency_matching)
  {
    trigger_efficiency_file = TFile::Open(trigger_efficiency_path.c_str());
    trigger_turnon_curve = (TF1*) trigger_efficiency_file->Get("fit_4_GeV");
  }

  rnd = new TRandom3(seednumber);

  if (do_ultra_smearing)
  {
    // Very large second term
    f_energy_ultra = new TF1("f_energy_ultra", "sqrt([0]/x*[0]/x)", 0.1, 20);
    f_energy_ultra->SetParameter(0, 0.61);
  }
  
  // From Blair pi0EtaByEta
  if (do_smearing)
  {
    const std::string kSmearParFile =
      "/sphenix/u/virgilemahaut/fromBlair/function_compare_wide_05102026.root";
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
  
  // book histograms
  book_histograms();

  _eventcounter = -1;

  return Fun4AllReturnCodes::EVENT_OK;
}

int MC_calo_pythia_pi0::InitRun(PHCompositeNode *topNode)
{
  truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthinfo)
  {
    std::cerr << "Error: No truth info" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _truth_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_Truth_r04");
  if (!_truth_jets && _mc_index >= 0 && _mc_index < nMCSamples)
  {
    std::cerr << "Error: No AntiKt_Truth_r04 jets found." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _reco_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_unsubtracted_r04");
  if (!_reco_jets)
  {
    std::cout << "Info: No AntiKt_unsubtracted_r04 jets found." << std::endl;
    std::cout << "No stitching cut applied on max reco jet pT" << std::endl;
  }

  vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    std::cerr << "Error: No vertex map found." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  clusterContainer =
    findNode::getClass<RawClusterContainer>(topNode,
                                            "CLUSTERINFO_CEMC");
  if (!clusterContainer)
  {
    std::cout << PHWHERE << "AnNeutralMeson - Fatal Error - "
      "CLUSTER_CEMC node is missing. "
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MC_calo_pythia_pi0::process_event(PHCompositeNode * /*topNode*/)
{
  _eventcounter++;

  if ((_eventcounter % 1000) == 0) std::cout << "MC_calo_pythia event " << _eventcounter << std::endl;

  if (Verbosity() >= VERBOSITY_SOME)
  {
    if ((_eventcounter % 1) == 0) std::cout << "MC_calo_pythia event (debug)" << _eventcounter << std::endl;
  }

  /////////////////////////////////////////////////
  //// Truth info

  // In nopileup MC data, we expect one vertex
  PHG4TruthInfoContainer::VtxRange vtxrange = truthinfo->GetPrimaryVtxRange();
  int nb_primary_vertices = 0;
  for (PHG4TruthInfoContainer::ConstVtxIterator iter = vtxrange.first; iter != vtxrange.second; ++iter)
  {
    nb_primary_vertices ++;
    PHG4VtxPoint* vtx = iter->second;
    vertex_z_true = vtx->get_z();
    // Assume (x, y)_vtx = (0, 0) as usual
  }
  if (nb_primary_vertices >= 2) {
    std::cout << "Unexpected number of primary vertices: " << nb_primary_vertices << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Stitch MC sample
  if (!stitch_mc_sample_truth_pt()) {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (do_vertex_reweighting) {
    m_event_weight = std::exp(
          -0.5 * (vertex_z_true * vertex_z_true) *
            (1.0 / (m_sigma_vertex * m_sigma_vertex)
           - 1.0 / (m_sigma_reference * m_sigma_reference))); 
  }
  
  meson_by_track_id.clear();
  photon_by_track_id.clear();
  photons_by_parent.clear();
  nondecay_photon_by_track_id.clear();
  nonphoton_by_track_id.clear();
  
  auto make_truth_particle = [](const PHG4Particle* truth)
  {
    TruthParticleInfo out;

    out.pid = truth->get_pid();
    out.track_id = truth->get_track_id();
    out.parent_id = truth->get_parent_id();

    out.p4 = ROOT::Math::PxPyPzEVector(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());

    return out;
  };
    
  auto collect_meson = [&](const PHG4Particle* truth)
  {
    if (!truth) return;
    
    const int pid = truth->get_pid();
    if (pid != 111 && pid != 221) return; // Either pi0 or eta meson

    // Remove very low energy mesons (less than 100 MeV)
    //if (truth->get_e() < 1) return;

    // Remove photons well outside of EMCal acceptance
    
    auto meson = make_truth_particle(truth);
    
    meson_by_track_id[meson.track_id] = meson;
    //decaying_mesons.push_back(meson);
  };
  
  auto collect_decay_photon = [&](const PHG4Particle* truth)
  {
    if (!truth) return;
    if (truth->get_pid() != 22) {
      //TruthParticleInfo nonphoton = make_truth_particle(truth);
      //nonphoton_by_track_id[nonphoton.track_id] = nonphoton;
      return; // Only photons
    }

    PHG4Particle* parent = truthinfo->GetParticle(truth->get_parent_id());
    if (!parent) {
      //TruthParticleInfo nondecay_photon = make_truth_particle(truth);
      //nondecay_photon_by_track_id[nondecay_photon.track_id] = nondecay_photon;
      return; // only decays
    }

    const int parent_pid = parent->get_pid();

    if (parent_pid != 111 && parent_pid != 221) {
      //TruthParticleInfo nondecay_photon = make_truth_particle(truth);
      //nondecay_photon_by_track_id[nondecay_photon.track_id] = nondecay_photon;
      return; // Either pi0 or eta decay
    }

    // Remove very low energy photons (less than 1000 MeV)
    //if (truth->get_e() < 1) return;
    
    TruthParticleInfo decay_photon = make_truth_particle(truth);

    //truth_decay_photons.push_back(decay_photon);
    photon_by_track_id[decay_photon.track_id] = decay_photon;
    photons_by_parent[parent->get_track_id()].push_back(decay_photon);
  };  

  // Mesons in primary range
  PHG4TruthInfoContainer::Range primary_range = truthinfo->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = primary_range.first; iter != primary_range.second; ++iter)
  {
    const PHG4Particle *truth = iter->second;
    collect_meson(truth);
    collect_decay_photon(truth);
  }
  
  // Secondary range
  PHG4TruthInfoContainer::Range second_range = truthinfo->GetSecondaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator siter = second_range.first; siter != second_range.second; ++siter)
  {
    const PHG4Particle *truth = siter->second;
    collect_meson(truth);
    collect_decay_photon(truth);
  }

  // Fill truth distributions
  fill_truth_histograms();

  // Extract reconstructed vertex from DST_GLOBAL
  if (vertexmap && !vertexmap->empty())
  {
    GlobalVertex *vtx = vertexmap->begin()->second;
    if (vtx)
    {
      vertex_z_reco = vtx->get_z();
    }
    else
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if (Verbosity() > VERBOSITY_SOME)
  {
    std::cout << "vertex = (" << vertex_z_true << ", " << vertex_z_reco << ")" << std::endl;
    std::cout << "(m_sigma_vertex, m_sigma_reference) = (" << m_sigma_vertex << ", " << m_sigma_reference << ")" << std::endl;
    std::cout << "m_event_weight = " << m_event_weight << std::endl;
  }

  // Remove 5 sigma vertex difference (reco_true gives sigma = 2cm)
  // This domain is not stat. significant but can be enhanced with event reweighting
  if (std::abs(vertex_z_reco - vertex_z_true) > 10) return Fun4AllReturnCodes::ABORTEVENT;
  
  h_true_zvtx->Fill(vertex_z_true, m_event_weight);
  h_reco_zvtx->Fill(vertex_z_reco, m_event_weight);
  h_reco_true_zvtx->Fill(vertex_z_reco - vertex_z_true, m_event_weight);
  
  RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
  RawClusterContainer::ConstIterator clusterIter;
  RawClusterContainer::ConstIterator clusterIter2;

  good_photons.clear();
  num_photons = 0;
  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
  {
    //cluster_index++;
    RawCluster* recoCluster = clusterIter->second;

    CLHEP::Hep3Vector vertex(0, 0, vertex_z_reco);
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);
    CLHEP::Hep3Vector pos_vec_cluster = recoCluster->get_position() - vertex;

    float clus_E = E_vec_cluster.mag();
    float clus_chisq = recoCluster->get_chi2();

    // As in data analysis, cut all photons below a certain energy threshold
    // Nominally, we have a cut of 1 GeV
    if (clus_E < clus_E_cut) continue;
    if (clus_chisq > clus_chisq_cut) continue;

    // Then match remaining photons to the true pi0
    float clus_eta = E_vec_cluster.pseudoRapidity();
    float clus_phi = E_vec_cluster.phi();
    float clus_pt = E_vec_cluster.perp();

    // Fill reco vector
    ROOT::Math::PtEtaPhiMVector photon_cluster(clus_pt, clus_eta, clus_phi, 0);
    customCluster reco_cluster;
    reco_cluster.p4 = photon_cluster;
    good_photons.push_back(reco_cluster);
    num_photons++;
  }

  if (require_efficiency_matching)
  {
    // Determine how likely the simulated event would have fired the photon trigger
    double max_cluster_energy = get_max_energy(good_photons);
    double trigger_efficiency = trigger_turnon_curve->Eval(max_cluster_energy);
    m_event_weight *= trigger_efficiency;
  }
  
  // Match truth to reco decay photons (when possible)
  match_truth_index(good_photons);

  // Apply extra smearing to account for resolution differences between MC and data (after truth matching!)
  if (do_smearing)
  {
    for (int iclus = 0; iclus < num_photons; iclus++)
    {
      float clus_eta = good_photons[iclus].p4.Eta();
      float clus_phi = good_photons[iclus].p4.Phi();
      float clus_E = good_photons[iclus].p4.E();
      float eta_nominal = clus_eta;
      float phi_nominal = clus_phi;
      float E_nominal = clus_E;
      smear_photon(clus_eta, clus_phi, clus_E, _smear_index); // cEp05

      if (do_ultra_smearing) {
        smear_ultra(clus_E);
      }

      if (Verbosity() > VERBOSITY_MORE) {
        std::cout << "reco cluster before extra smearing: ("
                  << good_photons[iclus].p4.Pt() << ", "
                  << good_photons[iclus].p4.Eta() << ", "
                  << good_photons[iclus].p4.Phi() << ")" << std::endl;
      }

      good_photons[iclus].p4.SetEta(clus_eta);
      good_photons[iclus].p4.SetPhi(clus_phi);
      good_photons[iclus].p4.SetPt(clus_E / cosh(clus_eta));

      if (Verbosity() > VERBOSITY_MORE) {
        std::cout << "reco cluster after extra smearing: ("
                  << good_photons[iclus].p4.Pt() << ", "
                  << good_photons[iclus].p4.Eta() << ", "
                  << good_photons[iclus].p4.Phi() << ")" << std::endl;
      }

      // Just to check that the reco clusters are really smeared.
      h_smear_E->Fill(clus_E / E_nominal, m_event_weight);
      h_smear_eta->Fill(clus_eta - eta_nominal, m_event_weight);
      h_smear_phi->Fill(clus_phi - phi_nominal, m_event_weight);
    
      h_smear_E_dE->Fill(E_nominal, clus_E / E_nominal, m_event_weight);
      h_smear_E_deta->Fill(E_nominal, clus_eta - eta_nominal, m_event_weight);
      h_smear_E_dphi->Fill(E_nominal, clus_phi - phi_nominal, m_event_weight);
    }
  }

  // Shift clusters according to absolute energy scale uncertainty
  for (int iclus = 0; iclus < num_photons; iclus++)
  {
    float clus_eta = good_photons[iclus].p4.Eta();
    float clus_E = good_photons[iclus].p4.E();
    float clus_E_scale = clus_E + clus_E * m_scale_diff / 100.0;
    float clus_pt_scale = clus_E_scale / cosh(clus_eta);
    good_photons[iclus].p4.SetPt(clus_pt_scale);
  }
  
  // Select the photon pairs
  float diphoton_pt_max = 0;
  std::vector<customCluster> photon1_accepted;
  std::vector<customCluster> photon2_accepted;
  std::vector<ROOT::Math::PtEtaPhiMVector> diphoton_accepted;
  for (int iclus1 = 0; iclus1 < num_photons - 1; iclus1++)
  {
    customCluster photon1 = good_photons[iclus1];
    for (int iclus2 = iclus1 + 1; iclus2 < num_photons; iclus2++)
    {
      customCluster photon2 = good_photons[iclus2];
      
      ROOT::Math::PtEtaPhiMVector diphoton = photon1.p4 + photon2.p4;

      // General diphoton cut
      if (diphoton_cut(photon1.p4, photon2.p4, diphoton)) continue;

      // Distinct diphoton selection between MBD and Photon Trigger

      // For a photon-triggered event
      // Keep candidates at high-pT which satisfy the trigger matching criterion

      // Trigger efficiency matching cut ("Proxy method")
      // Not for the Min-Bias selection
      if (require_efficiency_matching)
      {
        efficiency_match = trigger_efficiency_matching(photon1.p4, photon2.p4, diphoton);
        if (!(efficiency_match))
        {
          continue;
        }
      }

      if (Verbosity() > VERBOSITY_MORE) {
        std::cout << "diphoton (" << diphoton.Pt() << "," << diphoton.Eta() << "," << diphoton.Phi() << ")" << std::endl;
      }
      
      if (diphoton.Pt() > diphoton_pt_max) {
        diphoton_pt_max = diphoton.Pt();
      }

      photon1_accepted.push_back(photon1);
      photon2_accepted.push_back(photon2);
      diphoton_accepted.push_back(diphoton);

      // Start filling histograms with any reconstructed photon pair
      //fill_histograms(photon1, photon2, diphoton);
    }
  }

  if (!stitch_mc_sample_diphoton(diphoton_pt_max)) {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  for (size_t i = 0; i < photon1_accepted.size(); i++) {
    fill_histograms(photon1_accepted[i],
                    photon2_accepted[i],
                    diphoton_accepted[i]);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MC_calo_pythia_pi0::End(PHCompositeNode *)
{
  if (outfile)
  {
    outfile->cd();
    outfile->Write();
    outfile->Close();
    delete outfile;
    outfile = nullptr;
  }
  return 0;
}


void MC_calo_pythia_pi0::book_histograms()
{
  TH1::SetDefaultSumw2(true); // event-wise weight
  TH2::SetDefaultSumw2(true);
  
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  outfile->cd();

  // Vertices
  h_true_zvtx = new TH1F(
    "h_true_zvtx",
    ";z_{vtx} [cm]; Counts / [2 cm]",
    200, -200, 200);
  h_reco_zvtx = new TH1F(
    "h_reco_zvtx",
    ";z_{vtx} [cm]; Counts / [2 cm]",
    200, -200, 200);                       
  h_reco_true_zvtx = new TH1F(
    "h_reco_true_zvtx",
    ";#Deltaz_{vtx} [cm]; Counts / [1 mm]",
    200, -10, 10);                              

  // Smearing 
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

  // Invariant mass decomposition
  h_pair_mass_total = new TH1F(
    "h_pair_mass_total",
    ";M_{#gamma#gamma} [GeV]; Counts / [2 MeV]",
    500, 0, 1);
  for (int iType = 0; iType < nMassTypes; iType++) {
    h_pair_mass[iType] = new TH1F(
      ("h_pair_mass_" + type_suffix[iType]).c_str(),
      ";M_{#gamma#gamma} [GeV]; Counts / [2 MeV]",
      500, 0, 1);
  
    for (int izvtxBin = 0; izvtxBin < nZvtxBins; izvtxBin++)
    {
      std::stringstream h_pair_mass_name;
      h_pair_mass_name << std::fixed << std::setprecision(0)
                       << "h_pair_mass_" << type_suffix[iType] << "_zvtx_"
                       << izvtxBin;
      
      std::stringstream h_pair_mass_title;
      h_pair_mass_title << std::fixed << std::setprecision(2) << "diphoton mass ["
                        << zvtxBins[izvtxBin] << " < z_{vtx} [cm] < " << zvtxBins[izvtxBin + 1]
                        << "];M_{#gamma#gamma}; Counts / [2 MeV]";

      h_pair_mass_zvtx[iType][izvtxBin] = new TH1F(
        h_pair_mass_name.str().c_str(),
        h_pair_mass_title.str().c_str(),
        500, 0, 1);
    }

    for (int iEBin = 0; iEBin < nEBins; iEBin++)
    {
      std::stringstream h_pair_mass_name;
      h_pair_mass_name << std::fixed << std::setprecision(0)
                       << "h_pair_mass_" << type_suffix[iType] << "_E_"
                       << iEBin;

      std::stringstream h_pair_mass_title;
      h_pair_mass_title << std::fixed << std::setprecision(2) << "diphoton mass ["
                        << EBins[iEBin] << " < E [GeV] < " << EBins[iEBin + 1]
                        << "];M_{#gamma#gamma}; Counts / [2 MeV]";

      h_pair_mass_E[iType][iEBin] = new TH1F(
        h_pair_mass_name.str().c_str(),
        h_pair_mass_title.str().c_str(),
        500, 0, 1);
    }

    for (int iptBin = 0; iptBin < nPtBins; iptBin++)
    {
      std::stringstream h_pair_mass_name;
      h_pair_mass_name << std::fixed << std::setprecision(0)
                       << "h_pair_mass_" << type_suffix[iType] << "_pt_"
                       << iptBin;
      std::stringstream h_pair_mass_title;
      h_pair_mass_title << std::fixed << std::setprecision(2) << "diphoton mass ["
                        << pTBins[iptBin] << " < p_{T} [GeV] < " << pTBins[iptBin + 1]
                        << "];M_{#gamma#gamma}; Counts / [2 MeV]";

      h_pair_mass_pt[iType][iptBin] = new TH1F(
        h_pair_mass_name.str().c_str(),
        h_pair_mass_title.str().c_str(),
        500, 0, 1);
    }

    for (int ixfBin = 0; ixfBin < nXfBins; ixfBin++)
    {
      std::stringstream h_pair_mass_name;
      h_pair_mass_name << std::fixed << std::setprecision(0)
                       << "h_pair_mass_" << type_suffix[iType] << "_xf_"
                       << ixfBin;
      std::stringstream h_pair_mass_title;
      h_pair_mass_title << std::fixed << std::setprecision(2) << "diphoton mass ["
                        << xfBins[ixfBin] << " < x_{F} < " << xfBins[ixfBin + 1]
                        << "];M_{#gamma#gamma}; Counts / [2 MeV]";

      h_pair_mass_xf[iType][ixfBin] = new TH1F(
        h_pair_mass_name.str().c_str(),
        h_pair_mass_title.str().c_str(),
        500, 0, 1);
    }

    for (int ietaBin = 0; ietaBin < nEtaBins; ietaBin++)
    {
      std::stringstream h_pair_mass_name;
      h_pair_mass_name << std::fixed << std::setprecision(0)
                       << "h_pair_mass_" << type_suffix[iType] << "_eta_"
                       << ietaBin;
      std::stringstream h_pair_mass_title;
      h_pair_mass_title << std::fixed << std::setprecision(2) << "diphoton mass ["
                        << etaBins[ietaBin] << " < x_{F} < " << etaBins[ietaBin + 1]
                        << "];M_{#gamma#gamma}; Counts / [2 MeV]";

      h_pair_mass_eta[iType][ietaBin] = new TH1F(
        h_pair_mass_name.str().c_str(),
        h_pair_mass_title.str().c_str(),
        500, 0, 1);
    }

    // 1D reco distributions
    h_reco_E[iType] = new TH1F(
      ("h_reco_" + type_suffix[iType] + "_E").c_str(),
      ";E [GeV]; Counts / [100 MeV]",
      200, 0, 20);
    h_reco_pt[iType] = new TH1F(
      ("h_reco_" + type_suffix[iType] + "_pt").c_str(),
      ";p_{T} [GeV]; Counts / [100 MeV]",
      200, 0, 20);
    h_reco_eta[iType] = new TH1F(
      ("h_reco_" + type_suffix[iType] + "_eta").c_str(),
      ";#eta; Counts / [0.02]",
      200, -2.0, 2.0);
    h_reco_xf[iType] = new TH1F(
       ("h_reco_" + type_suffix[iType] + "_xf").c_str(),
      ";x_{F}; Counts / [0.002]",
      200, -0.2, 0.2);
    h_reco_phi[iType] = new TH1F(
      ("h_reco_" + type_suffix[iType] + "_phi").c_str(),
      ";#phi_{G} [rad]; Counts / [#pi/128]",
      256, -M_PI, M_PI);
  }

  // Truth-Matched Distributions
  for (int iType = 0; iType < nMassTypes - 1; iType++) {
    
    // 1D generated distributions
    h_true_E[iType] = new TH1F(
      ("h_true_" + type_suffix[iType] + "_E").c_str(),
      ";E [GeV]; Counts / [100 MeV]",
      200, 0, 20);
    h_true_pt[iType] = new TH1F(
      ("h_true_" + type_suffix[iType] + "_pt").c_str(),
      ";p_{T} [GeV]; Counts / [100 MeV]",
      200, 0, 20);
    h_true_eta[iType] = new TH1F(
      ("h_true_" + type_suffix[iType] + "_eta").c_str(),
      ";#eta; Counts / [0.02]",
      200, -2.0, 2.0);
    h_true_xf[iType] = new TH1F(
       ("h_true_" + type_suffix[iType] + "_xf").c_str(),
      ";x_{F}; Counts / [0.002]",
      200, -0.2, 0.2);
    h_true_phi[iType] = new TH1F(
      ("h_true_" + type_suffix[iType] + "_phi").c_str(),
      ";#phi_{G} [rad]; Counts / [#pi/128]",
      256, -M_PI, M_PI);
    
    // Differences Reco/True
    h_reco_true_E[iType] = new TH1F(
      ("h_reco_true_" + type_suffix[iType] + "_E").c_str(),
      ";E_{Reco}/E_{True}; Counts / [0.02]",
      100, 0.0, 2.0);
    h_reco_true_pt[iType] = new TH1F(
      ("h_reco_true_" + type_suffix[iType] + "_pt").c_str(),
      ";p_{T,Reco}/p_{T,True}; Counts / [0.02]",
      100, 0.0, 2.0);
    h_reco_true_eta[iType] = new TH1F(
      ("h_reco_true_" + type_suffix[iType] + "_eta").c_str(),
      ";#Delta#eta; Counts / [0.0005]",
      100, -0.025, 0.025);
    h_reco_true_xf[iType] = new TH1F(
      ("h_reco_true_" + type_suffix[iType] + "_xf_").c_str(),
      ";#Deltax_{F}; Counts / [0.00005]",
      100, -0.0025, 0.0025);
    h_reco_true_phi[iType] = new TH1F(
      ("h_reco_true_" + type_suffix[iType] + "_phi").c_str(),
      ";#Delta#phi_{G} [rad]; Counts / [0.0005]",
      100, -0.025, 0.025);

    // Response matrices for unfolding
    h_response_pt[iType] = new TH2F(
      ("h_response_" + type_suffix[iType] + "_pt").c_str(),
      ";p_{T,Reco} [GeV]; p_{T,True} [GeV]",
      nPtBins, pTBins,
      nPtBins, pTBins);
    h_response_eta[iType] = new TH2F(
      ("h_response_" + type_suffix[iType] + "_eta").c_str(),
      ";#eta_{Reco}; #eta_{True}",
      nEtaBins, etaBins,
      nEtaBins, etaBins);
    h_response_xf[iType] = new TH2F(
      ("h_response_" + type_suffix[iType] + "_xf_").c_str(),
      ";x_{F,Reco}; x_{F,True}",
      nXfBins, xfBins,
      nXfBins, xfBins);
    h_response_phi[iType] = new TH2F(
      ("h_response_" + type_suffix[iType] + "_phi").c_str(),
      ";#phi_{Reco} [rad]; #phi_{True} [rad]",
      nPhiBins, phiBins,
      nPhiBins, phiBins);

    h_response_fine_pt[iType] = new TH2F(
      ("h_response_fine_" + type_suffix[iType] + "_pt").c_str(),
      ";p_{T,Reco} [GeV]; p_{T,True} [GeV]",
      200, 0, 20,
      200, 0, 20);
    h_response_fine_eta[iType] = new TH2F(
      ("h_response_fine_" + type_suffix[iType] + "_eta").c_str(),
      ";#eta_{Reco}; #eta_{True}",
      200, -2.0, 2.0,
      200, -2.0, 2.0);
    h_response_fine_xf[iType] = new TH2F(
      ("h_response_fine_" + type_suffix[iType] + "_xf_").c_str(),
      ";x_{F,Reco}; x_{F,True}",
      200, -0.2, 0.2,
      200, -0.2, 0.2);
    h_response_fine_phi[iType] = new TH2F(
      ("h_response_fine_" + type_suffix[iType] + "_phi").c_str(),
      ";#phi_{Reco} [rad]; #phi_{True} [rad]",
      200, -M_PI, M_PI,
      200, -M_PI, M_PI);
    
    h_response_pt_phi[iType] = new TH2F(
      ("h_response_" + type_suffix[iType] + "_pt_phi").c_str(),
      ";(p_{T}, #phi)_{Reco} Index; (p_{T}, #phi)_{True} Index",
      nPtBins * nPhiBins, -0.5, static_cast<float>(nPtBins * nPhiBins) - 0.5,
      nPtBins * nPhiBins, -0.5, static_cast<float>(nPtBins * nPhiBins) - 0.5);
    h_response_eta_phi[iType] = new TH2F(
      ("h_response_" + type_suffix[iType] + "_eta_phi").c_str(),
      ";(#eta, #phi)_{Reco} Index; (#eta, #phi)_{True} Index",
      nEtaBins * nPhiBins, -0.5, static_cast<float>(nEtaBins * nPhiBins) - 0.5,
      nEtaBins * nPhiBins, -0.5, static_cast<float>(nEtaBins * nPhiBins) - 0.5);
    h_response_xf_phi[iType] = new TH2F(
      ("h_response_" + type_suffix[iType] + "_xf_phi").c_str(),
      ";(x_{F}, #phi)_{Reco} Index; (x_{F}, #phi)_{True} Index",
      nXfBins * nPhiBins, -0.5, static_cast<float>(nXfBins * nPhiBins) - 0.5,
      nXfBins * nPhiBins, -0.5, static_cast<float>(nXfBins * nPhiBins) - 0.5);

    // Matrices for meson reconstruction efficiency per true bin
    h_true_pt_n[iType] = new TH1F(
      ("h_true_" + type_suffix[iType] + "_pt_n").c_str(),
      ";p_{T,True} [GeV];Count",
      nPtBins, pTBins);
    h_true_eta_n[iType] = new TH1F(
      ("h_true_" + type_suffix[iType] + "_eta_n").c_str(),
      ";#eta_{True};Count",
      nEtaBins, etaBins);
    h_true_xf_n[iType] = new TH1F(
      ("h_true_" + type_suffix[iType] + "_xf_n").c_str(),
      ";x_{F,True};Count",
      nXfBins, xfBins);
    h_true_pt_phi_n[iType] = new TH1F(
      ("h_true_" + type_suffix[iType] + "_pt_phi_n").c_str(),
      ";(p_{T}, #phi)_{True} Index;Count",
      nPtBins * nPhiBins, -0.5, static_cast<float>(nPtBins * nPhiBins) - 0.5);
    h_true_eta_phi_n[iType] = new TH1F(
      ("h_true_" + type_suffix[iType] + "_eta_phi_n").c_str(),
      ";(#eta, #phi)_{True} Index;Count",
      nEtaBins * nPhiBins, -0.5, static_cast<float>(nEtaBins * nPhiBins) - 0.5);
    h_true_xf_phi_n[iType] = new TH1F(
      ("h_true_" + type_suffix[iType] + "_xf_phi_n").c_str(),
      ";(x_{F}, #phi)_{True} Index;Count",
      nXfBins * nPhiBins, -0.5, static_cast<float>(nXfBins * nPhiBins) - 0.5);

    h_reco_pt_n_matched[iType] = new TH1F(
      ("h_reco_" + type_suffix[iType] + "_pt_n_matched").c_str(),
      ";p_{T,Reco} [GeV];Count",
      nPtBins, pTBins);
    h_reco_eta_n_matched[iType] = new TH1F(
      ("h_reco_" + type_suffix[iType] + "_eta_n_matched").c_str(),
      ";#eta_{Reco};Count",
      nEtaBins, etaBins);
    h_reco_xf_n_matched[iType] = new TH1F(
      ("h_reco_" + type_suffix[iType] + "_xf_n_matched").c_str(),
      ";x_{F,Reco};Count",
      nXfBins, xfBins);
    h_reco_pt_phi_n_matched[iType] = new TH1F(
      ("h_reco_" + type_suffix[iType] + "_pt_phi_n_matched").c_str(),
      ";(p_{T}, #phi)_{Reco} Index;Count",
      nPtBins * nPhiBins, -0.5, static_cast<float>(nPtBins * nPhiBins) - 0.5);
    h_reco_eta_phi_n_matched[iType] = new TH1F(
      ("h_reco_" + type_suffix[iType] + "_eta_phi_n_matched").c_str(),
      ";(#eta, #phi)_{Reco} Index;Count",
      nEtaBins * nPhiBins, -0.5, static_cast<float>(nEtaBins * nPhiBins) - 0.5);
    h_reco_xf_phi_n_matched[iType] = new TH1F(
      ("h_reco_" + type_suffix[iType] + "_xf_phi_n_matched").c_str(),
      ";(x_{F}, #phi)_{Reco} Index;Count",
      nXfBins * nPhiBins, -0.5, static_cast<float>(nXfBins * nPhiBins) - 0.5);  

    // Individual Reco/True histograms in coarse binning
    h_true_coarse_pt[iType] = new TH1F(
      ("h_true_coarse_" + type_suffix[iType] + "_pt").c_str(),
      ";p_{T,Reco}; Count",
      nPtBins, pTBins);
    h_true_coarse_eta[iType] = new TH1F(
      ("h_true_coarse_" + type_suffix[iType] + "_eta").c_str(),
      ";#eta_{Reco}; Count",
      nEtaBins, etaBins);
    h_true_coarse_xf[iType] = new TH1F(
      ("h_true_coarse_" + type_suffix[iType] + "_xf_").c_str(),
      ";x_{FReco}; Count",
      nXfBins, xfBins);
    h_true_coarse_phi[iType] = new TH1F(
      ("h_true_coarse_" + type_suffix[iType] + "_phi").c_str(),
      ";#phi_{Reco}; Count",
      nPhiBins, phiBins);
    h_true_coarse_pt_phi[iType] = new TH1F(
      ("h_true_coarse_" + type_suffix[iType] + "_pt_phi").c_str(),
      ";(p_{T}, #phi)_{Reco} Index;Count",
      nPtBins * nPhiBins, -0.5, static_cast<float>(nPtBins * nPhiBins) - 0.5);
    h_true_coarse_eta_phi[iType] = new TH1F(
      ("h_true_coarse_" + type_suffix[iType] + "_eta_phi").c_str(),
      ";(#eta, #phi)_{Reco} Index;Count",
      nEtaBins * nPhiBins, -0.5, static_cast<float>(nEtaBins * nPhiBins) - 0.5);
    h_true_coarse_xf_phi[iType] = new TH1F(
      ("h_true_coarse_" + type_suffix[iType] + "_xf_phi").c_str(),
      ";(x_{F}, #phi)_{Reco} Index;Count",
      nXfBins * nPhiBins, -0.5, static_cast<float>(nXfBins * nPhiBins) - 0.5);

    h_reco_coarse_pt[iType] = new TH1F(
      ("h_reco_coarse_" + type_suffix[iType] + "_pt").c_str(),
      ";p_{T,Reco}; Count",
      nPtBins, pTBins);
    h_reco_coarse_eta[iType] = new TH1F(
      ("h_reco_coarse_" + type_suffix[iType] + "_eta").c_str(),
      ";#eta_{Reco}; Count",
      nEtaBins, etaBins);
    h_reco_coarse_xf[iType] = new TH1F(
      ("h_reco_coarse_" + type_suffix[iType] + "_xf_").c_str(),
      ";x_{FReco}; Count",
      nXfBins, xfBins);
    h_reco_coarse_phi[iType] = new TH1F(
      ("h_reco_coarse_" + type_suffix[iType] + "_phi").c_str(),
      ";#phi_{Reco}; Count",
      nPhiBins, phiBins);
    h_reco_coarse_pt_phi[iType] = new TH1F(
      ("h_reco_coarse_" + type_suffix[iType] + "_pt_phi").c_str(),
      ";(p_{T}, #phi)_{Reco} Index;Count",
      nPtBins * nPhiBins, -0.5, static_cast<float>(nPtBins * nPhiBins) - 0.5);
    h_reco_coarse_eta_phi[iType] = new TH1F(
      ("h_reco_coarse_" + type_suffix[iType] + "_eta_phi").c_str(),
      ";(#eta, #phi)_{Reco} Index;Count",
      nEtaBins * nPhiBins, -0.5, static_cast<float>(nEtaBins * nPhiBins) - 0.5);
    h_reco_coarse_xf_phi[iType] = new TH1F(
      ("h_reco_coarse_" + type_suffix[iType] + "_xf_phi").c_str(),
      ";(x_{F}, #phi)_{Reco} Index;Count",
      nXfBins * nPhiBins, -0.5, static_cast<float>(nXfBins * nPhiBins) - 0.5);
  }
}

void MC_calo_pythia_pi0::fill_truth_histograms()
{
  outfile->cd();
  
  for (const auto& [truth_meson_id, truth_meson_particle] : meson_by_track_id)
  {
    const ROOT::Math::PxPyPzEVector& truth_meson = truth_meson_particle.p4;

    const float truth_meson_pt = truth_meson.Pt();
    const float truth_meson_E = truth_meson.E();
    const float truth_meson_eta = truth_meson.Eta();
    float truth_meson_phi = truth_meson.Phi();
    WrapAngle(truth_meson_phi); // between -pi and +pi
    const float truth_meson_xf = 2 * truth_meson.Pz() / Vs;

    // Select the pt bin:
    int truth_iPt = FindBinBinary(truth_meson_pt, pTBins, nPtBins + 1);
    
    // Select the phi bin:
    int truth_iphiBin = FindBinDirect(truth_meson_phi, -M_PI, M_PI, nPhiBins);

    // Select the eta bin (symmetrized)
    int truth_ietaBin = FindBinBinary(truth_meson_eta, etaBins, nEtaBins + 1);
    //if (truth_ietaBin >= nEtaBins / 2) truth_ietaBin = nEtaBins - 1 - truth_ietaBin;

    // Select the xf bin (symmetrized)
    int truth_ixfBin = FindBinBinary(truth_meson_xf, xfBins, nXfBins + 1);
    //if (truth_ixfBin >= nXfBins / 2) truth_ixfBin = nXfBins - 1 - truth_ixfBin;

    // Warning FindBinBinary starts from 0 and overflow -1 !!!
    bool good_pt_truth_index = (truth_iPt >= 0 && truth_iPt < nPtBins);
    bool good_eta_truth_index = (truth_ietaBin >= 0 && truth_ietaBin < nEtaBins);
    bool good_xf_truth_index = (truth_ixfBin >= 0 && truth_ixfBin < nXfBins);
    bool good_phi_truth_index = (truth_iphiBin >= 0 && truth_iphiBin < nPhiBins);

    int truth_pt_phi_bin = -1;
    int truth_eta_phi_bin = -1;
    int truth_xf_phi_bin = -1;
    if (good_phi_truth_index)
    {
      if (good_pt_truth_index) truth_pt_phi_bin = global_index(truth_iPt, truth_iphiBin);
      if (good_eta_truth_index) truth_eta_phi_bin = global_index(truth_ietaBin, truth_iphiBin);
      if (good_xf_truth_index) truth_xf_phi_bin = global_index(truth_ixfBin, truth_iphiBin);
    }

    if (truth_meson_particle.pid == 111)
    {
      h_true_pt[0]->Fill(truth_meson_pt, m_event_weight);
      h_true_E[0]->Fill(truth_meson_E, m_event_weight);
      h_true_eta[0]->Fill(truth_meson_eta, m_event_weight);
      h_true_xf[0]->Fill(truth_meson_xf, m_event_weight);
      h_true_phi[0]->Fill(truth_meson_phi, m_event_weight);
      
      // N_true for reconstruction efficiency
      h_true_pt_n[0]->Fill(truth_meson_pt, m_event_weight);
      h_true_eta_n[0]->Fill(truth_meson_eta, m_event_weight);
      h_true_xf_n[0]->Fill(truth_meson_xf, m_event_weight);
      if (truth_pt_phi_bin != -1) h_true_pt_phi_n[0]->Fill(truth_pt_phi_bin, m_event_weight);
      if (truth_eta_phi_bin != -1) h_true_eta_phi_n[0]->Fill(truth_eta_phi_bin, m_event_weight);
      if (truth_xf_phi_bin != -1) h_true_xf_phi_n[0]->Fill(truth_xf_phi_bin, m_event_weight);
    }
    else if (truth_meson_particle.pid == 221)
    {
      h_true_pt[1]->Fill(truth_meson_pt, m_event_weight);
      h_true_E[1]->Fill(truth_meson_E, m_event_weight);
      h_true_eta[1]->Fill(truth_meson_eta, m_event_weight);
      h_true_xf[1]->Fill(truth_meson_xf, m_event_weight);
      h_true_phi[1]->Fill(truth_meson_phi, m_event_weight);
      
      // N_true for reconstruction efficiency
      h_true_pt_n[1]->Fill(truth_meson_pt, m_event_weight);
      h_true_eta_n[1]->Fill(truth_meson_eta, m_event_weight);
      h_true_xf_n[1]->Fill(truth_meson_xf, m_event_weight);
      if (truth_pt_phi_bin != -1) h_true_pt_phi_n[1]->Fill(truth_pt_phi_bin, m_event_weight);
      if (truth_eta_phi_bin != -1) h_true_eta_phi_n[1]->Fill(truth_eta_phi_bin, m_event_weight);
      if (truth_xf_phi_bin != -1) h_true_xf_phi_n[1]->Fill(truth_xf_phi_bin, m_event_weight);
    }
  }
}

void MC_calo_pythia_pi0::fill_histograms(const customCluster& photon1,
                                         const customCluster& photon2,
                                         const ROOT::Math::PtEtaPhiMVector& diphoton)
{
  outfile->cd();
  
  // If a pi0 or an eta meson decays into a photon pair
  // it is extremely likely a direct decay (see PDG)
  const float diphoton_E = diphoton.E() / 2; // Average cluster energy
  const float diphoton_pt = diphoton.Pt();
  const float diphoton_eta = diphoton.Eta();
  float diphoton_phi = diphoton.Phi();
  WrapAngle(diphoton_phi); // between -pi and +pi
  const float diphoton_mass = diphoton.mag();
  const float diphoton_xf = 2 * diphoton.Pz() / Vs;

  h_pair_mass_total->Fill(diphoton_mass, m_event_weight);

  // Select the energy bin:
  int iE = FindBinBinary(diphoton_E, pTBins, nPtBins + 1);

  // Select the pt bin:
  int iPt = FindBinBinary(diphoton_pt, pTBins, nPtBins + 1);

  // Select the zvtx bin:
  int izvtx = FindBinBinary(std::abs(vertex_z_reco), zvtxBins, nZvtxBins + 1);

  // Select the phi bin:
  int iphiBin = FindBinDirect(diphoton_phi, -M_PI, M_PI, nPhiBins);

  // Select the eta bin (symmetrized)
  int ietaBin = FindBinBinary(diphoton_eta, etaBins, nEtaBins + 1);
  //if (ietaBin > nEtaBins / 2) ietaBin = nEtaBins - 1 - ietaBin;

  // Select the xf bin (symmetrized)
  int ixfBin = FindBinBinary(diphoton_xf, xfBins, nXfBins + 1);
  //if (ixfBin > nXfBins / 2) ixfBin = nXfBins - 1 - ixfBin;

  // Warning FindBinBinary starts from 0 and overflow -1 !!!
  bool good_pt_index = (iPt >= 0 && iPt < nPtBins);
  bool good_eta_index = (ietaBin >= 0 && ietaBin < nEtaBins);
  bool good_xf_index = (ixfBin >= 0 && ixfBin < nXfBins);
  bool good_phi_index = (iphiBin >= 0 && iphiBin < nPhiBins);

  int reco_pt_phi_bin = -1;
  int reco_eta_phi_bin = -1;
  int reco_xf_phi_bin = -1;
  if (good_phi_index)
  {
    if (good_pt_index) reco_pt_phi_bin = global_index(iPt, iphiBin);
    if (good_eta_index) reco_eta_phi_bin = global_index(ietaBin, iphiBin);
    if (good_xf_index) reco_xf_phi_bin = global_index(ixfBin, iphiBin);
  }

  if (photon1.match_truth_track_id == -1 || photon2.match_truth_track_id == -1 ||
      photon1.match_truth_track_id == photon2.match_truth_track_id)
  {
    // Fill unmatched histograms

    // Reco distributions
    h_reco_pt[MassType::Unmatched]->Fill(diphoton_pt, m_event_weight);
    h_reco_eta[MassType::Unmatched]->Fill(diphoton_eta, m_event_weight);
    h_reco_xf[MassType::Unmatched]->Fill(diphoton_xf, m_event_weight);
    h_reco_phi[MassType::Unmatched]->Fill(diphoton_phi, m_event_weight);

    // Reco Invariant mass
    h_pair_mass[MassType::Unmatched]->Fill(diphoton_mass, m_event_weight);
    if (izvtx >= 0 && izvtx < nZvtxBins)
      for (int izvtxtmp = izvtx; izvtxtmp < nZvtxBins; izvtxtmp++)
        h_pair_mass_zvtx[MassType::Unmatched][izvtxtmp]->Fill(diphoton_mass, m_event_weight);
    if (iPt >= 0 && iPt < nPtBins)
      h_pair_mass_pt[MassType::Unmatched][iPt]->Fill(diphoton_mass, m_event_weight);
    if (iE >= 0 && iE < nEBins)
      h_pair_mass_E[MassType::Unmatched][iE]->Fill(diphoton_mass, m_event_weight);
    if (ietaBin >= 0 && ietaBin < nEtaBins)
      h_pair_mass_eta[MassType::Unmatched][ietaBin]->Fill(diphoton_mass, m_event_weight);
    if (ixfBin >= 0 && ixfBin < nXfBins)
      h_pair_mass_xf[MassType::Unmatched][ixfBin]->Fill(diphoton_mass, m_event_weight);
    
    return;
  }

  // Truth diphoton
  const TruthParticleInfo& truth_photon1 = photon_by_track_id[photon1.match_truth_track_id];
  const TruthParticleInfo& truth_photon2 = photon_by_track_id[photon2.match_truth_track_id];
  const ROOT::Math::PxPyPzEVector truth_diphoton = truth_photon1.p4 + truth_photon2.p4;

  const float truth_diphoton_pt = truth_diphoton.Pt();
  const float truth_diphoton_E = truth_diphoton.E();
  const float truth_diphoton_eta = truth_diphoton.Eta();
  float truth_diphoton_phi = truth_diphoton.Phi();
  WrapAngle(truth_diphoton_phi);// between -pi and +pi
  const float truth_diphoton_xf = 2 * truth_diphoton.Pz() / Vs;

  // Select the energy bin:
  //int truth_iE = FindBinBinary(truth_diphoton_E, pTBins, nPtBins + 1);

  // Select the pt bin:
  int truth_iPt = FindBinBinary(truth_diphoton_pt, pTBins, nPtBins + 1);

  // Select the zvtx bin:
  //int truth_izvtx = FindBinBinary(std::abs(vertex_z_true), zvtxBins, nZvtxBins + 1);

  // Select the phi bin:
  int truth_iphiBin = FindBinDirect(truth_diphoton_phi, -M_PI, M_PI, nPhiBins);

  // Select the eta bin (symmetrized)
  int truth_ietaBin = FindBinBinary(truth_diphoton_eta, etaBins, nEtaBins + 1);
  //if (truth_ietaBin > nEtaBins / 2) truth_ietaBin = nEtaBins - 1 - truth_ietaBin;

  // Select the xf bin (symmetrized)
  int truth_ixfBin = FindBinBinary(truth_diphoton_xf, xfBins, nXfBins + 1);
  //if (truth_ixfBin > nXfBins / 2) truth_ixfBin = nXfBins - 1 - truth_ixfBin;

  // Warning FindBinBinary starts from 0 and overflow -1 !!!
  bool good_pt_truth_index = (truth_iPt >= 0 && truth_iPt < nPtBins);
  bool good_eta_truth_index = (truth_ietaBin >= 0 && truth_ietaBin < nEtaBins);
  bool good_xf_truth_index = (truth_ixfBin >= 0 && truth_ixfBin < nXfBins);
  bool good_phi_truth_index = (truth_iphiBin >= 0 && truth_iphiBin < nPhiBins);

  int truth_pt_phi_bin = -1;
  int truth_eta_phi_bin = -1;
  int truth_xf_phi_bin = -1;
  if (good_phi_truth_index)
  {
    if (good_pt_truth_index) truth_pt_phi_bin = global_index(truth_iPt, truth_iphiBin);
    if (good_eta_truth_index) truth_eta_phi_bin = global_index(truth_ietaBin, truth_iphiBin);
    if (good_xf_truth_index) truth_xf_phi_bin = global_index(truth_ixfBin, truth_iphiBin);
  }
  
  // First case photon-pair matches a true pi0
  if (photon1.match_parent_pid == 111 &&
      photon1.match_parent_track_id == photon2.match_parent_track_id)
  {
    // Fill pi0 histograms

    if (Verbosity() > VERBOSITY_MORE) {
      std::cout << "matched pi0 signal" << std::endl;
      std::cout << "cluster 1: ("
                  << photon1.match_truth_track_id << ", "
                  << photon1.match_truth_pid << ", "
                  << photon1.match_truth_type << ", "
                  << photon1.match_parent_track_id << ", "
                  << photon1.match_parent_pid << ", "
                  << photon1.match_truth_dR << ", "
                  << photon1.match_truth_Eratio << ")" << std::endl;
      std::cout << "cluster 2: ("
                  << photon2.match_truth_track_id << ", "
                  << photon2.match_truth_pid << ", "
                  << photon2.match_truth_type << ", "
                  << photon2.match_parent_track_id << ", "
                  << photon2.match_parent_pid << ", "
                  << photon2.match_truth_dR << ", "
                  << photon2.match_truth_Eratio << ")" << std::endl;
    }

    // Reco distributions
    h_reco_E[MassType::Pi0]->Fill(diphoton_E, m_event_weight);
    h_reco_pt[MassType::Pi0]->Fill(diphoton_pt, m_event_weight);
    h_reco_eta[MassType::Pi0]->Fill(diphoton_eta, m_event_weight);
    h_reco_xf[MassType::Pi0]->Fill(diphoton_xf, m_event_weight);
    h_reco_phi[MassType::Pi0]->Fill(diphoton_phi, m_event_weight);

    // Differences versus MC
    h_reco_true_E[MassType::Pi0]->Fill(diphoton_E / truth_diphoton_E, m_event_weight);
    h_reco_true_pt[MassType::Pi0]->Fill(diphoton_pt / truth_diphoton_pt, m_event_weight);
    h_reco_true_eta[MassType::Pi0]->Fill(diphoton_eta - truth_diphoton_eta, m_event_weight);
    h_reco_true_xf[MassType::Pi0]->Fill(diphoton_xf - truth_diphoton_xf, m_event_weight);
    h_reco_true_phi[MassType::Pi0]->Fill(diphoton_phi - truth_diphoton_phi, m_event_weight);

    // Reco Invariant mass
    h_pair_mass[MassType::Pi0]->Fill(diphoton_mass, m_event_weight);
    if (izvtx >= 0 && izvtx < nZvtxBins)
      for (int izvtxtmp = izvtx; izvtxtmp < nZvtxBins; izvtxtmp++)
        h_pair_mass_zvtx[MassType::Pi0][izvtxtmp]->Fill(diphoton_mass, m_event_weight);
    if (iPt >= 0 && iPt < nPtBins)
      h_pair_mass_pt[MassType::Pi0][iPt]->Fill(diphoton_mass, m_event_weight);
    if (iE >= 0 && iE < nEBins)
      h_pair_mass_E[MassType::Pi0][iE]->Fill(diphoton_mass, m_event_weight);
    if (ietaBin >= 0 && ietaBin < nEtaBins)
      h_pair_mass_eta[MassType::Pi0][ietaBin]->Fill(diphoton_mass, m_event_weight);
    if (ixfBin >= 0 && ixfBin < nXfBins)
      h_pair_mass_xf[MassType::Pi0][ixfBin]->Fill(diphoton_mass, m_event_weight);

    // Efficiency response
    h_reco_pt_n_matched[MassType::Pi0]->Fill(truth_diphoton_pt, m_event_weight);
    h_reco_eta_n_matched[MassType::Pi0]->Fill(truth_diphoton_eta, m_event_weight);
    h_reco_xf_n_matched[MassType::Pi0]->Fill(truth_diphoton_xf, m_event_weight);
    if (truth_pt_phi_bin != -1) h_reco_pt_phi_n_matched[MassType::Pi0]->Fill(truth_pt_phi_bin, m_event_weight);
    if (truth_eta_phi_bin != -1) h_reco_eta_phi_n_matched[MassType::Pi0]->Fill(truth_eta_phi_bin, m_event_weight);
    if (truth_xf_phi_bin != -1) h_reco_xf_phi_n_matched[MassType::Pi0]->Fill(truth_xf_phi_bin, m_event_weight);

    // Response Matrix
    h_response_pt[MassType::Pi0]->Fill(diphoton_pt, truth_diphoton_pt, m_event_weight);
    h_response_eta[MassType::Pi0]->Fill(diphoton_eta, truth_diphoton_eta, m_event_weight);
    h_response_xf[MassType::Pi0]->Fill(diphoton_xf, truth_diphoton_xf, m_event_weight);
    h_response_phi[MassType::Pi0]->Fill(diphoton_phi, truth_diphoton_phi, m_event_weight);
    h_response_fine_pt[MassType::Pi0]->Fill(diphoton_pt, truth_diphoton_pt, m_event_weight);
    h_response_fine_eta[MassType::Pi0]->Fill(diphoton_eta, truth_diphoton_eta, m_event_weight);
    h_response_fine_xf[MassType::Pi0]->Fill(diphoton_xf, truth_diphoton_xf, m_event_weight);
    h_response_fine_phi[MassType::Pi0]->Fill(diphoton_phi, truth_diphoton_phi, m_event_weight);

    
    if (truth_pt_phi_bin != -1 && reco_pt_phi_bin != -1) h_response_pt_phi[MassType::Pi0]->Fill(reco_pt_phi_bin, truth_pt_phi_bin, m_event_weight);
    if (truth_eta_phi_bin != -1 && reco_eta_phi_bin != -1) h_response_eta_phi[MassType::Pi0]->Fill(reco_eta_phi_bin, truth_eta_phi_bin, m_event_weight);
    if (truth_xf_phi_bin != -1 && reco_xf_phi_bin != -1) h_response_xf_phi[MassType::Pi0]->Fill(reco_xf_phi_bin, truth_xf_phi_bin, m_event_weight);

    // Individual reco/true histograms
    h_reco_coarse_pt[MassType::Pi0]->Fill(diphoton_pt, m_event_weight);
    h_true_coarse_pt[MassType::Pi0]->Fill(truth_diphoton_pt, m_event_weight);
    h_reco_coarse_eta[MassType::Pi0]->Fill(diphoton_eta, m_event_weight);
    h_true_coarse_eta[MassType::Pi0]->Fill(truth_diphoton_eta, m_event_weight);
    h_reco_coarse_xf[MassType::Pi0]->Fill(diphoton_xf, m_event_weight);
    h_true_coarse_xf[MassType::Pi0]->Fill(truth_diphoton_xf, m_event_weight);
    h_reco_coarse_phi[MassType::Pi0]->Fill(diphoton_phi, m_event_weight);
    h_true_coarse_phi[MassType::Pi0]->Fill(truth_diphoton_phi, m_event_weight);
    h_reco_coarse_pt_phi[MassType::Pi0]->Fill(reco_pt_phi_bin, m_event_weight);
    h_true_coarse_pt_phi[MassType::Pi0]->Fill(truth_pt_phi_bin, m_event_weight);
    h_reco_coarse_eta_phi[MassType::Pi0]->Fill(reco_eta_phi_bin, m_event_weight);
    h_true_coarse_eta_phi[MassType::Pi0]->Fill(truth_eta_phi_bin, m_event_weight);
    h_reco_coarse_xf_phi[MassType::Pi0]->Fill(reco_xf_phi_bin, m_event_weight);
    h_true_coarse_xf_phi[MassType::Pi0]->Fill(truth_xf_phi_bin, m_event_weight);
  }
  // Second case photon-pair matches a true eta
  else if (photon1.match_parent_pid == 221 &&
           photon1.match_parent_track_id == photon2.match_parent_track_id)
  {
    // Fill eta histograms

    // Reco distributions
    h_reco_pt[MassType::Eta]->Fill(diphoton_pt, m_event_weight);
    h_reco_E[MassType::Eta]->Fill(diphoton_E, m_event_weight);
    h_reco_eta[MassType::Eta]->Fill(diphoton_eta, m_event_weight);
    h_reco_xf[MassType::Eta]->Fill(diphoton_xf, m_event_weight);
    h_reco_phi[MassType::Eta]->Fill(diphoton_phi, m_event_weight);

    // Differences versus MC
    h_reco_true_pt[MassType::Eta]->Fill(diphoton_pt / truth_diphoton_pt, m_event_weight);
    h_reco_true_E[MassType::Eta]->Fill(diphoton_E / truth_diphoton_E, m_event_weight);
    h_reco_true_eta[MassType::Eta]->Fill(diphoton_eta - truth_diphoton_eta, m_event_weight);
    h_reco_true_xf[MassType::Eta]->Fill(diphoton_xf - truth_diphoton_xf, m_event_weight);
    h_reco_true_phi[MassType::Eta]->Fill(diphoton_phi - truth_diphoton_phi, m_event_weight);

    // Reco Invariant mass
    h_pair_mass[MassType::Eta]->Fill(diphoton_mass, m_event_weight);
    if (izvtx >= 0 && izvtx < nZvtxBins)
      for (int izvtxtmp = izvtx; izvtxtmp < nZvtxBins; izvtxtmp++)
        h_pair_mass_zvtx[MassType::Eta][izvtxtmp]->Fill(diphoton_mass, m_event_weight);
    if (iPt >= 0 && iPt < nPtBins)
      h_pair_mass_pt[MassType::Eta][iPt]->Fill(diphoton_mass, m_event_weight);
    if (iE >= 0 && iE < nEBins)
      h_pair_mass_E[MassType::Eta][iE]->Fill(diphoton_mass, m_event_weight);
    if (ietaBin >= 0 && ietaBin < nEtaBins)
      h_pair_mass_eta[MassType::Eta][ietaBin]->Fill(diphoton_mass, m_event_weight);
    if (ixfBin >= 0 && ixfBin < nXfBins)
      h_pair_mass_xf[MassType::Eta][ixfBin]->Fill(diphoton_mass, m_event_weight);

    // Efficiency response
    h_reco_pt_n_matched[MassType::Eta]->Fill(truth_diphoton_pt, m_event_weight);
    h_reco_eta_n_matched[MassType::Eta]->Fill(truth_diphoton_eta, m_event_weight);
    h_reco_xf_n_matched[MassType::Eta]->Fill(truth_diphoton_xf, m_event_weight);
    if (truth_pt_phi_bin != -1) h_reco_pt_phi_n_matched[MassType::Eta]->Fill(truth_pt_phi_bin, m_event_weight);
    if (truth_eta_phi_bin != -1) h_reco_eta_phi_n_matched[MassType::Eta]->Fill(truth_eta_phi_bin, m_event_weight);
    if (truth_xf_phi_bin != -1) h_reco_xf_phi_n_matched[MassType::Eta]->Fill(truth_xf_phi_bin, m_event_weight);

    // Response Matrix
    h_response_pt[MassType::Eta]->Fill(diphoton_pt, truth_diphoton_pt, m_event_weight);
    h_response_eta[MassType::Eta]->Fill(diphoton_eta, truth_diphoton_eta, m_event_weight);
    h_response_xf[MassType::Eta]->Fill(diphoton_xf, truth_diphoton_xf, m_event_weight);
    h_response_phi[MassType::Eta]->Fill(diphoton_phi, truth_diphoton_phi, m_event_weight);

    h_response_fine_pt[MassType::Eta]->Fill(diphoton_pt, truth_diphoton_pt, m_event_weight);
    h_response_fine_eta[MassType::Eta]->Fill(diphoton_eta, truth_diphoton_eta, m_event_weight);
    h_response_fine_xf[MassType::Eta]->Fill(diphoton_xf, truth_diphoton_xf, m_event_weight);
    h_response_fine_phi[MassType::Eta]->Fill(diphoton_phi, truth_diphoton_phi, m_event_weight);
    
    if (truth_pt_phi_bin != -1 && reco_pt_phi_bin != -1) h_response_pt_phi[MassType::Eta]->Fill(reco_pt_phi_bin, truth_pt_phi_bin, m_event_weight);
    if (truth_eta_phi_bin != -1 && reco_eta_phi_bin != -1) h_response_eta_phi[MassType::Eta]->Fill(reco_eta_phi_bin, truth_eta_phi_bin, m_event_weight);
    if (truth_xf_phi_bin != -1 && reco_xf_phi_bin != -1) h_response_xf_phi[MassType::Eta]->Fill(reco_xf_phi_bin, truth_xf_phi_bin, m_event_weight);

    // Individual reco/true histograms
    h_reco_coarse_pt[MassType::Eta]->Fill(diphoton_pt, m_event_weight);
    h_true_coarse_pt[MassType::Eta]->Fill(truth_diphoton_pt, m_event_weight);
    h_reco_coarse_eta[MassType::Eta]->Fill(diphoton_eta, m_event_weight);
    h_true_coarse_eta[MassType::Eta]->Fill(truth_diphoton_eta, m_event_weight);
    h_reco_coarse_xf[MassType::Eta]->Fill(diphoton_xf, m_event_weight);
    h_true_coarse_xf[MassType::Eta]->Fill(truth_diphoton_xf, m_event_weight);
    h_reco_coarse_phi[MassType::Eta]->Fill(diphoton_phi, m_event_weight);
    h_true_coarse_phi[MassType::Eta]->Fill(truth_diphoton_phi, m_event_weight);
    h_reco_coarse_pt_phi[MassType::Eta]->Fill(reco_pt_phi_bin, m_event_weight);
    h_true_coarse_pt_phi[MassType::Eta]->Fill(truth_pt_phi_bin, m_event_weight);
    h_reco_coarse_eta_phi[MassType::Eta]->Fill(reco_eta_phi_bin, m_event_weight);
    h_true_coarse_eta_phi[MassType::Eta]->Fill(truth_eta_phi_bin, m_event_weight);
    h_reco_coarse_xf_phi[MassType::Eta]->Fill(reco_xf_phi_bin, m_event_weight);
    h_true_coarse_xf_phi[MassType::Eta]->Fill(truth_xf_phi_bin, m_event_weight);
  }
  else if (photon1.match_truth_track_id != -1 && photon2.match_truth_track_id != -1 &&
           photon1.match_truth_track_id != photon2.match_truth_track_id)
  {
    // Fill combinatorial histogram

    // Reco distributions
    h_reco_pt[MassType::Combinatorial]->Fill(diphoton_pt, m_event_weight);
    h_reco_E[MassType::Combinatorial]->Fill(diphoton_E, m_event_weight);
    h_reco_eta[MassType::Combinatorial]->Fill(diphoton_eta, m_event_weight);
    h_reco_xf[MassType::Combinatorial]->Fill(diphoton_xf, m_event_weight);
    h_reco_phi[MassType::Combinatorial]->Fill(diphoton_phi, m_event_weight);

    // Differences versus MC
    h_reco_true_pt[MassType::Combinatorial]->Fill(diphoton_pt / truth_diphoton_pt, m_event_weight);
    h_reco_true_E[MassType::Combinatorial]->Fill(diphoton_E / truth_diphoton_E, m_event_weight);
    h_reco_true_eta[MassType::Combinatorial]->Fill(diphoton_eta - truth_diphoton_eta, m_event_weight);
    h_reco_true_xf[MassType::Combinatorial]->Fill(diphoton_xf - truth_diphoton_xf, m_event_weight);
    h_reco_true_phi[MassType::Combinatorial]->Fill(diphoton_phi - truth_diphoton_phi, m_event_weight);

    // Reco Invariant mass
    h_pair_mass[MassType::Combinatorial]->Fill(diphoton_mass, m_event_weight);
    if (izvtx >= 0 && izvtx < nZvtxBins)
      for (int izvtxtmp = izvtx; izvtxtmp < nZvtxBins; izvtxtmp++)
        h_pair_mass_zvtx[MassType::Combinatorial][izvtxtmp]->Fill(diphoton_mass, m_event_weight);
    if (iPt >= 0 && iPt < nPtBins)
      h_pair_mass_pt[MassType::Combinatorial][iPt]->Fill(diphoton_mass, m_event_weight);
    if (iE >= 0 && iE < nEBins)
      h_pair_mass_E[MassType::Combinatorial][iE]->Fill(diphoton_mass, m_event_weight);
    if (ietaBin >= 0 && ietaBin < nEtaBins)
      h_pair_mass_eta[MassType::Combinatorial][ietaBin]->Fill(diphoton_mass, m_event_weight);
    if (ixfBin >= 0 && ixfBin < nXfBins)
      h_pair_mass_xf[MassType::Combinatorial][ixfBin]->Fill(diphoton_mass, m_event_weight);

    // Efficiency response
    h_reco_pt_n_matched[MassType::Combinatorial]->Fill(truth_diphoton_pt, m_event_weight);
    h_reco_eta_n_matched[MassType::Combinatorial]->Fill(truth_diphoton_eta, m_event_weight);
    h_reco_xf_n_matched[MassType::Combinatorial]->Fill(truth_diphoton_xf, m_event_weight);
    if (truth_pt_phi_bin != -1) h_reco_pt_phi_n_matched[MassType::Combinatorial]->Fill(truth_pt_phi_bin, m_event_weight);
    if (truth_eta_phi_bin != -1) h_reco_eta_phi_n_matched[MassType::Combinatorial]->Fill(truth_eta_phi_bin, m_event_weight);
    if (truth_xf_phi_bin != -1) h_reco_xf_phi_n_matched[MassType::Combinatorial]->Fill(truth_xf_phi_bin, m_event_weight);

    // Response Matrix
    h_response_pt[MassType::Combinatorial]->Fill(diphoton_pt, truth_diphoton_pt, m_event_weight);
    h_response_eta[MassType::Combinatorial]->Fill(diphoton_eta, truth_diphoton_eta, m_event_weight);
    h_response_xf[MassType::Combinatorial]->Fill(diphoton_xf, truth_diphoton_xf, m_event_weight);
    h_response_phi[MassType::Combinatorial]->Fill(diphoton_phi, truth_diphoton_phi, m_event_weight);

    h_response_fine_pt[MassType::Combinatorial]->Fill(diphoton_pt, truth_diphoton_pt, m_event_weight);
    h_response_fine_eta[MassType::Combinatorial]->Fill(diphoton_eta, truth_diphoton_eta, m_event_weight);
    h_response_fine_xf[MassType::Combinatorial]->Fill(diphoton_xf, truth_diphoton_xf, m_event_weight);
    h_response_fine_phi[MassType::Combinatorial]->Fill(diphoton_phi, truth_diphoton_phi, m_event_weight);
    
    if (truth_pt_phi_bin != -1 && reco_pt_phi_bin != -1) h_response_pt_phi[MassType::Combinatorial]->Fill(reco_pt_phi_bin, truth_pt_phi_bin, m_event_weight);
    if (truth_eta_phi_bin != -1 && reco_eta_phi_bin != -1) h_response_eta_phi[MassType::Combinatorial]->Fill(reco_eta_phi_bin, truth_eta_phi_bin, m_event_weight);
    if (truth_xf_phi_bin != -1 && reco_xf_phi_bin != -1) h_response_xf_phi[MassType::Combinatorial]->Fill(reco_xf_phi_bin, truth_xf_phi_bin, m_event_weight);

    // Individual reco/true histograms
    h_reco_coarse_pt[MassType::Combinatorial]->Fill(diphoton_pt, m_event_weight);
    h_true_coarse_pt[MassType::Combinatorial]->Fill(truth_diphoton_pt, m_event_weight);
    h_reco_coarse_eta[MassType::Combinatorial]->Fill(diphoton_eta, m_event_weight);
    h_true_coarse_eta[MassType::Combinatorial]->Fill(truth_diphoton_eta, m_event_weight);
    h_reco_coarse_xf[MassType::Combinatorial]->Fill(diphoton_xf, m_event_weight);
    h_true_coarse_xf[MassType::Combinatorial]->Fill(truth_diphoton_xf, m_event_weight);
    h_reco_coarse_phi[MassType::Combinatorial]->Fill(diphoton_phi, m_event_weight);
    h_true_coarse_phi[MassType::Combinatorial]->Fill(truth_diphoton_phi, m_event_weight);
    h_reco_coarse_pt_phi[MassType::Combinatorial]->Fill(reco_pt_phi_bin, m_event_weight);
    h_true_coarse_pt_phi[MassType::Combinatorial]->Fill(truth_pt_phi_bin, m_event_weight);
    h_reco_coarse_eta_phi[MassType::Combinatorial]->Fill(reco_eta_phi_bin, m_event_weight);
    h_true_coarse_eta_phi[MassType::Combinatorial]->Fill(truth_eta_phi_bin, m_event_weight);
    h_reco_coarse_xf_phi[MassType::Combinatorial]->Fill(reco_xf_phi_bin, m_event_weight);
    h_true_coarse_xf_phi[MassType::Combinatorial]->Fill(truth_xf_phi_bin, m_event_weight);
  }
  else
  {
    std::cerr << "fill_histograms: This case should not happen. Please take a look back into the code." << std::endl;
  }
}

bool MC_calo_pythia_pi0::trigger_efficiency_matching(const ROOT::Math::PtEtaPhiMVector& photon1,
                                                     const ROOT::Math::PtEtaPhiMVector& photon2,
                                                     const ROOT::Math::PtEtaPhiMVector& diphoton)
{
  float energy_threshold = 0;

  //if (trigger_mbd_photon_3) energy_threshold = energy_threshold_3; // Photon 3 GeV trigger
  //else if (trigger_mbd_photon_4) energy_threshold = energy_threshold_4; // Photon 4 GeV trigger
  energy_threshold = energy_threshold_4; // No trigger in MC simulation
  
  // else 
  // {
  //   std::cerr << "Error: either photon 3 GeV or photon 4 GeV scaled trigger bit should be fired at this point.\n";
  //   exit(1);
  // }
  
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

float MC_calo_pythia_pi0::get_max_energy(std::vector<customCluster>& reco_clusters)
{
  float max_energy = 0;
  for (const auto& reco_clus : reco_clusters)
  {
    float energy = reco_clus.p4.E();
    if (energy > max_energy) max_energy = energy;
  }
  return max_energy;
}

void MC_calo_pythia_pi0::match_truth_index(std::vector<customCluster>& reco_clusters)
{
  // Minimum cut for matching
  const double delR_cut = 0.07; // More than 2 towers!
  const double res_min_cut = 0.3;
  const double res_max_cut = 2.0;

  // First, clusters are unmatched

  if (Verbosity() >= VERBOSITY_MORE) std::cout << "reco clusters" << std::endl;
  for (auto& reco_clus : reco_clusters)
  {
    reco_clus.match_truth_track_id = -1; // What truth photon matches best
    reco_clus.match_parent_track_id = -1; // If truth match, what parent meson
    reco_clus.match_parent_pid = 0; // pi0, eta or else
    
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << "(" << reco_clus.p4.E() << ", " << reco_clus.p4.Pt() << ", " << reco_clus.p4.Eta() << ", " << reco_clus.p4.Phi() << ")" << std::endl;
    }
  }

  if (Verbosity() >= VERBOSITY_MORE)
  {
    std::cout << "truth mesons" << std::endl;
    for (const auto& [track_id, meson] : meson_by_track_id)
    {
      std::cout << "(" << track_id << ", " << meson.pid << ", " << meson.p4.E() << ", " << meson.p4.Pt() << ", " << meson.p4.Eta() << ", " << meson.p4.Phi() << ")" << std::endl;
    }
    std::cout << "truth decay photons" << std::endl;
    for (const auto& [track_id, photon] : photon_by_track_id)
    {
      std::cout << "(" << track_id << ", " << photon.parent_id << ", " << photon.p4.E() << ", " << photon.p4.Pt() << ", " << photon.p4.Eta() << ", " << photon.p4.Phi() << ")" << " i.e.";
      std::cout << "(" << photon.p4.Px() << ", " << photon.p4.Py() << ", " << photon.p4.Pz() << ", " << photon.p4.E() << ")" << std::endl;
    }
  }

  auto try_match_truth = [&](const TruthParticleInfo& truth,
                             TruthMatchType type,
                             int parent_track_id,
                             int parent_pid,
                             int& best_track_id,
                             int& best_pid,
                             int& best_parent_id,
                             int& best_parent_pid,
                             TruthMatchType& best_type,
                             double& best_score,
                             double& best_dR,
                             double& best_Eratio,
                             const customCluster& reco_clus)
  {
    const double delR = ROOT::Math::VectorUtil::DeltaR(reco_clus.p4, truth.p4);
    if (delR > delR_cut) return;

    const double Eratio = reco_clus.p4.E() / truth.p4.E();
    if (Eratio < res_min_cut || Eratio > res_max_cut) return;

    const double score = delR / delR_cut;

    if (score < best_score)
    {
      best_score = score;
      best_track_id = truth.track_id;
      best_pid = truth.pid;
      best_parent_id = parent_track_id;
      best_parent_pid = parent_pid;
      best_type = type;
      best_dR = delR;
      best_Eratio = Eratio;
    }
  };

  // Find the truth photon that best matches reco cluster (if any)
  for (auto& reco_clus : reco_clusters)
  {
    int best_track_id = -1;
    int best_pid = 0;
    int best_parent_id = -1;
    int best_parent_pid = 0;
    TruthMatchType best_type = TruthMatchType::None;

    double best_score = 1000;
    double best_dR = -1;
    double best_Eratio = -1;

    // 1. Decay photons from pi0 / eta
    for (const auto& [parent_id, truth_photons] : photons_by_parent)
    {
      auto meson_it = meson_by_track_id.find(parent_id);
      if (meson_it == meson_by_track_id.end()) continue;

      const int parent_pid = meson_it->second.pid;

      for (const auto& truth_photon : truth_photons)
      {
        try_match_truth(truth_photon,
                        TruthMatchType::DecayPhoton,
                        parent_id,
                        parent_pid,
                        best_track_id,
                        best_pid,
                        best_parent_id,
                        best_parent_pid,
                        best_type,
                        best_score,
                        best_dR,
                        best_Eratio,
                        reco_clus);
      }
    }

    // // 2. Nondecay photons
    // for (const auto& [track_id, truth_photon] : nondecay_photon_by_track_id)
    // {
    //   try_match_truth(truth_photon,
    //                   TruthMatchType::NonDecayPhoton,
    //                   truth_photon.parent_id,
    //                   0,
    //                   best_track_id,
    //                   best_pid,
    //                   best_parent_id,
    //                   best_parent_pid,
    //                   best_type,
    //                   best_score,
    //                   best_dR,
    //                   best_Eratio,
    //                   reco_clus);
    // }

    // // 3. Nonphotons
    // for (const auto& [track_id, truth_particle] : nonphoton_by_track_id)
    // {
    //   try_match_truth(truth_particle,
    //                   TruthMatchType::NonPhoton,
    //                   truth_particle.parent_id,
    //                   0,
    //                   best_track_id,
    //                   best_pid,
    //                   best_parent_id,
    //                   best_parent_pid,
    //                   best_type,
    //                   best_score,
    //                   best_dR,
    //                   best_Eratio,
    //                   reco_clus);
    // }

    if (best_track_id != -1)
    {
      reco_clus.match_truth_track_id = best_track_id;
      reco_clus.match_truth_pid = best_pid;
      reco_clus.match_parent_track_id = best_parent_id;
      reco_clus.match_parent_pid = best_parent_pid;
      reco_clus.match_truth_type = best_type;
      reco_clus.match_truth_dR = best_dR;
      reco_clus.match_truth_Eratio = best_Eratio;

      if (Verbosity() >= VERBOSITY_MORE)
      {
        std::cout << "matched cluster: ("
                  << reco_clus.match_truth_track_id << ", "
                  << reco_clus.match_truth_pid << ", "
                  << reco_clus.match_truth_type << ", "
                  << reco_clus.match_parent_track_id << ", "
                  << reco_clus.match_parent_pid << ", "
                  << reco_clus.match_truth_dR << ", "
                  << reco_clus.match_truth_Eratio << ")" << std::endl;
      }
      
    }
  }
}

void MC_calo_pythia_pi0::smear_ultra(float& ecore)
{
  if (!f_energy_ultra) {
    std::cerr << "Error: f_energy_ultra not defined" << std::endl;
    return;
  }

  const double E_res = f_energy_ultra->Eval(ecore);
  float smeared_energy = std::max(0.0, rnd->Gaus(ecore, ecore * E_res));
  ecore = smeared_energy;
}

void MC_calo_pythia_pi0::smear_photon(float& eta,
                                      float& phi,
                                      float& ecore,
                                      int smear_idx)
{
  // Defaults if parameter file missing or Eval invalid
  float E_res = 0.08f;
  float P_res = 0.0001f;

  
  const double se = f_energy_smear[smear_idx]->Eval(ecore);
  if (std::isfinite(se) && se > 0.0)
  {
    E_res = se;
  }

  const double sp = f_position_smear[smear_idx]->Eval(ecore);
  if (std::isfinite(sp) && sp > 0.0)
  {
    P_res = sp;
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

void MC_calo_pythia_pi0::WrapAngle(float& phi)
{
  while (phi > static_cast<float>(M_PI)) {
    phi -= 2.0 * M_PI;
  }
  while (phi < - static_cast<float>(M_PI)) {
    phi += 2.0 * M_PI;
  }
}

float MC_calo_pythia_pi0::WrapAngleDifference(const float& phi1, const float& phi2)
{
  float difference = std::fmod(phi1 - phi2, 2 * M_PI);
  if (difference < 0)
  {
    difference += 2 * M_PI;
  }
  return std::min(difference, (float)(2 * M_PI - difference));
}

bool MC_calo_pythia_pi0::diphoton_cut(ROOT::Math::PtEtaPhiMVector p1,
                                      ROOT::Math::PtEtaPhiMVector p2,
                                      ROOT::Math::PtEtaPhiMVector ppair)
{
  float alpha = std::abs(p1.E() - p2.E()) / (p1.E() + p2.E());
  float pt = ppair.Pt();
  return (alpha > alphaCut || pt < pTCutMin || pt > pTCutMax);
}

bool MC_calo_pythia_pi0::stitch_mc_sample_truth_pt()
{
  bool truth_cut_valid = false;
  bool reco_cut_valid = false;
  
  if (_mc_index < 0 || _mc_index > nMCSamples) {
    if (Verbosity() > VERBOSITY_SOME)
    {
      std::cout << "_mc_index = " << _mc_index << " out of bounds. No stitching cut." << std::endl;
    }
    return true;
  }
  
  max_truth_jet_pt = 0; 
  for (auto truthjet : *_truth_jets)
  {
    const float truth_jet_pt = truthjet->get_pt();
    if (truth_jet_pt > max_truth_jet_pt) max_truth_jet_pt = truth_jet_pt;
  }
  
  if (_truth_pt_min[_mc_index] <= max_truth_jet_pt && max_truth_jet_pt <= _truth_pt_max[_mc_index])
  {
    if (Verbosity() > VERBOSITY_SOME)
    {
      std::cout << "max truth jet pT = " << max_truth_jet_pt << " in [" << _truth_pt_min[_mc_index] << ", " << _truth_pt_max[_mc_index] << "]" << std::endl;
      std::cout << "Truth-accepted event." << std::endl;
    }
    truth_cut_valid = true; 
  }
  else if (Verbosity() > VERBOSITY_SOME)
  {
    std::cout << "max truth jet pT = " << max_truth_jet_pt << " outside of [" << _truth_pt_min[_mc_index] << ", " << _truth_pt_max[_mc_index] << "]" << std::endl;
    std::cout << "Event does not pass." << std::endl;
  }

  if (!_reco_jets) {
    return truth_cut_valid;
  }

  float max_reco_jet_pt = 0; 
  for (auto recojet : *_reco_jets)
  {
    const float reco_jet_pt = recojet->get_pt();
    if (reco_jet_pt > max_reco_jet_pt) max_reco_jet_pt = reco_jet_pt;
  }
  
  if (max_reco_jet_pt <= _reco_pt_max[_mc_index])
  {
    if (Verbosity() > VERBOSITY_SOME)
    {
      std::cout << "max reco jet pT = " << max_reco_jet_pt << " lower than " << _reco_pt_max[_mc_index] << std::endl;
      std::cout << "Reco-accepted event." << std::endl;
    }
    reco_cut_valid = true; 
  }
  else if (Verbosity() > VERBOSITY_SOME)
  {
    std::cout << "max reco jet pT = " << max_reco_jet_pt << " greater than " << _reco_pt_max[_mc_index] << std::endl;
    std::cout << "Event does not pass." << std::endl;
  }
  
  return (truth_cut_valid && reco_cut_valid); 

  return false;
}

bool MC_calo_pythia_pi0::stitch_mc_sample_diphoton(float diphoton_pt)
{
  //if (_diphoton_pt_min[_mc_index] < diphoton_pt && diphoton_pt < _diphoton_pt_max[_mc_index]) {
  if (diphoton_pt < max_truth_jet_pt) {
    if (Verbosity() > VERBOSITY_MORE) {
      std::cout << "max diphoton pT " << diphoton_pt << " below max truth jet pT " << max_truth_jet_pt << std::endl;
      std::cout << "Diphoton-accepted event." << std::endl;
    }
    return true;
  }

  if (Verbosity() > VERBOSITY_MORE) {
    std::cout << "max diphoton pT " << diphoton_pt << " above max truth jet pT " << max_truth_jet_pt << std::endl;
    std::cout << "Diphoton-rejected event." << std::endl;
  }
  return false;
}
