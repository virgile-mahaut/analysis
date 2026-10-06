#include "UPCMeson.h"

/// Tracking includes
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
//#include <trackbase/TrackVertexCrossingAssoc.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

/// Truth evaluation includes
#include <g4eval/SvtxEvalStack.h>

/// HEPMC truth includes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

/// Fun4All includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <ffaobjects/EventHeader.h>
#include <ffarawobjects/Gl1Packet.h>

/// ROOT includes
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TTree.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <TDatabasePDG.h>


/// C++ includes
#include <cassert>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
#include <array>

/**
 * Constructor of module
 */
UPCMeson::UPCMeson(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , m_outfilename(filename)
  , m_analyzeTracks(true)
  , m_analyzeTruth(true)
{
  /// Initialize variables and trees so we don't accidentally access
  /// memory that was never allocated
  initializeVariables();
}

/**
 * Destructor of module
 */
UPCMeson::~UPCMeson()
{
  std::cout << PHWHERE << "In ~UPCMeson()" << std::endl;
  /*
  if ( m_hepmctree!=nullptr )  delete m_hepmctree;
  if ( m_tracktree!=nullptr )  delete m_tracktree;
  if ( m_truthtree!=nullptr )  delete m_truthtree;
  if ( m_pairtree!=nullptr )   delete m_pairtree;
  if ( m_globaltree!=nullptr ) delete m_globaltree;
  */
}

/**
 * Initialize the module and prepare looping over events
 */
int UPCMeson::Init(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning Init in UPCMeson" << std::endl;
  }

  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  initializeTrees();

  h_phi[0] = new TH1F("h_phi", "#phi [rad]", 60, -M_PI, M_PI);
  h2_eta_phi[0] = new TH2F("h2_phi_eta", ";#eta;#phi [rad]", 24, -5.0, 5.0, 60, -M_PI, M_PI);
  h_mass[0] = new TH1F("h_mass", "mass [GeV]", 1200, 0, 6);
  h_pt[0] = new TH1F("h_pt", "p_{T}", 200, 0, 2);
  h_y[0] = new TH1F("h_y", "y", 24, -1.2, 1.2);
  h_eta[0] = new TH1F("h_eta", "#eta", 24, -5.0, 5.0);

  // like-sign pairs
  h_phi[1] = new TH1F("h_phi_ls", "#phi [rad]", 60, -M_PI, M_PI);
  h2_eta_phi[1] = new TH2F("h2_phi_eta_ls", ";#eta;#phi [rad]", 24, -5.0, 5.0, 60, -M_PI, M_PI);
  h_mass[1] = new TH1F("h_mass_ls", "mass [GeV]", 1200, 0, 6);
  h_pt[1] = new TH1F("h_pt_ls", "p_{T}", 200, 0, 2);
  h_y[1] = new TH1F("h_y_ls", "y", 24, -1.2, 1.2);
  h_eta[1] = new TH1F("h_eta_ls", "#eta", 24, -5.0, 5.0);
 
  h_trig = new TH1F("h_trig", "trig", 65, -0.5, 64.5);
  h_ntracks = new TH1F("h_ntracks", "num tracks", 2000, 0, 2000);
  h2_ntrksvsb = new TH2F("h2_ntrksvsb", "num tracks vs b", 220, 0, 22, 2001, -0.5, 2000.5);
  h2_ntrksvsb->SetXTitle("b [fm]");
  h2_ntrksvsb->SetYTitle("N_{TRKS}");

  h_cross_evt = new TH1F("h_cross_evt", "cross;cross", 700, -CROSS_OFFSET-0.5, 700-CROSS_OFFSET-0.5);
  h_cross = new TH1F("h_cross", "cross;cross", 700, -CROSS_OFFSET-0.5, 700-CROSS_OFFSET-0.5);
  h_bunch = new TH1F("h_bunch", "bunch", 120, -0.5, 119.5);

  h_b_mb = new TH1F("h_b_mb", "b, MB events", 200, 0, 20);
  h_npart_mb = new TH1F("h_npart_mb", "npart, MB events", 401, -0.5, 400.5);
  h_ncoll_mb = new TH1F("h_ncoll_mb", "ncoll, MB events", 1301, -0.5, 1300.5);
  h_b = new TH1F("h_b", "b", 200, 0, 20);
  h_npart = new TH1F("h_npart", "npart", 401, -0.5, 400.5);
  h_ncoll = new TH1F("h_ncoll", "ncoll", 1301, -0.5, 1300.5);

  return 0;
}

/**
 * Main workhorse function where each event is looped over and
 * data from each event is collected from the node tree for analysis
 */
int UPCMeson::GetNodes(PHCompositeNode *topNode)
{
  /// EventHeader node
  _evthdr = findNode::getClass<EventHeader>(topNode, "EventHeader");

  /// GL1 node
  _gl1raw = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  if (!_gl1raw)
  {
    static int ctr = 0;
    if ( ctr<4 )
    {
      std::cout << PHWHERE << "GL1Packet node is missing, no trigger info" << std::endl;
      ctr++;
    }
  }

  /// SVTX tracks node
  _trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_trackmap)
  {
    std::cout << PHWHERE << "SvtxTrackMap node is missing, can't collect tracks" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /*
  _track_vertex_crossing_map = findNode::getClass<TrackVertexCrossingAssoc>(topNode,"TrackVertexCrossingAssocMap");
  if(!m_track_vertex_crossing_map)
  {
    std::cout << PHWHERE << "No TrackVertexCrossingAssocMap on node tree, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  */

  /// G4 truth particle node
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truthinfo)
  {
    static int ctr = 0;
    if ( ctr<4 )
    {
      std::cout << PHWHERE << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles" << std::endl;
      ctr++;
    }
  }

  /// HEPMC info
  _genevent_map = findNode::getClass<PHHepMCGenEventMap>(topNode,"PHHepMCGenEventMap");
  if (!_genevent_map)
  {
    static int ctr = 0;
    if ( ctr<4 )
    {
      std::cout << PHWHERE << "PHHepMCGenEventMap node is missing, can't collect HEPMC info" << std::endl;
      ctr++;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * Main workhorse function where each event is looped over and
 * data from each event is collected from the node tree for analysis
 */
int UPCMeson::process_event(PHCompositeNode *topNode)
{
  if ( m_evt%1000 == 1 )
  {
    std::cout << "m_evt " << m_evt << std::endl;
  }

  Reset(topNode);

  if (Verbosity() > 5)
  {
    std::cout << "Beginning process_event in UPCMeson, " << m_evt << std::endl;
  }

  h_trig->Fill( 0 );  // processed event counter

  if ( _gl1raw )
  {
    m_bunch = _gl1raw->getBunchNumber();
    h_bunch->Fill( m_bunch );

    m_strig = _gl1raw->getScaledVector();
    //std::cout << "strig " << std::hex << m_strig << std::dec << std::endl;

    uint64_t trigbit = 0x1UL;
    for (int ibit=1; ibit<=64; ibit++)
    {
      if ( (m_strig&trigbit) != 0 )
      {
        h_trig->Fill( ibit );
      }

      trigbit = trigbit<<1;
    }
  }

  /// Get all the data nodes
  int status = GetNodes(topNode);
  if ( status != Fun4AllReturnCodes::EVENT_OK )
  {
    return status;
  }

  /// Get the run and eventnumber
  if ( _evthdr )
  {
    m_run = _evthdr->get_RunNumber();
    m_evt = _evthdr->get_EvtSequence();
  }

  if ( m_ntrks!=0 || m_ntrk_sphenix!= 0 )
  {
    std::cout << PHWHERE << " ERROR, evt m_ntrks m_ntrk_sphenix = " << m_evt << "\t"
      << m_ntrks << "\t" << m_ntrk_sphenix << std::endl;
  }

  /// Get global info

  if ( _genevent_map )
  {
    if (Verbosity()>5)
    {
      std::cout << PHWHERE << "processing global info from sim" << std::endl;
    }
    PHHepMCGenEvent *genevent = (_genevent_map->begin())->second; 
    HepMC::GenEvent *hepmc_event{nullptr};
    HepMC::HeavyIon *hepmc_hi{nullptr};
    if (genevent)
    {
      hepmc_event = genevent->getEvent();
      hepmc_hi = hepmc_event->heavy_ion();
      if ( hepmc_hi )
      {
        m_npart_targ =  hepmc_hi->Npart_targ();
        m_npart_proj =  hepmc_hi->Npart_proj();
        m_npart = m_npart_targ + m_npart_proj;
        m_ncoll =  hepmc_hi->Ncoll();
        m_ncoll_hard =  hepmc_hi->Ncoll_hard();
        //std::cout << "ncoll " << m_ncoll << "\t" << m_ncoll_hard << std::endl;
        m_bimpact =  hepmc_hi->impact_parameter();
        //std::cout << "b ntracks " << m_bimpact << "\t" << m_ntrks << std::endl;

        //m_globaltree->Fill();

        h_b_mb->Fill( m_bimpact );
        h_npart_mb->Fill( m_npart );
        h_ncoll_mb->Fill( m_ncoll );
      }
    }
  }

  /// Get the tracks
  if (m_analyzeTracks)
  {
    status = getTracks(topNode);
    if ( status != Fun4AllReturnCodes::EVENT_OK )
    {
      return status;
    }
  }

  /// Get the truth track information
  if (m_analyzeTruth && _truthinfo )
  {
    getPHG4Truth();
  }

  // Fill Global info for events that pass cuts
  if ( _genevent_map )
  {
    h2_ntrksvsb->Fill( m_bimpact, m_ntrks );

    if ( m_ntrks<=3 )
    {
      h_b->Fill( m_bimpact );
      h_npart->Fill( m_npart );
      h_ncoll->Fill( m_ncoll );
    }

    m_globaltree->Fill();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * End the module and finish any data collection. Clean up any remaining
 * loose ends
 */
int UPCMeson::End(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 5)
  {
    std::cout << "Ending UPCMeson analysis package" << std::endl;
  }

  /// Write and close the outfile
  m_outfile->Write();
  m_outfile->Close();

  if (Verbosity() > 1)
  {
    std::cout << "Finished UPCMeson analysis package" << std::endl;
  }

  return 0;
}

/**
 * This method gets all of the HEPMC truth particles from the node tree
 * and stores them in a ROOT TTree. The HEPMC truth particles are what,
 * for example, directly comes out of PYTHIA and thus gives you all of
 * the associated parton information
 */
void UPCMeson::getHEPMCTruth()
{
  /// Could have some print statements for debugging with verbosity
  if (Verbosity() > 1)
  {
    std::cout << "Getting HEPMC truth particles " << std::endl;
  }

  /// You can iterate over the number of events in a hepmc event
  /// for pile up events where you have multiple hard scatterings per bunch crossing
  for (PHHepMCGenEventMap::ConstIter eventIter = _genevent_map->begin(); eventIter != _genevent_map->end(); ++eventIter)
  {
    /// Get the event
    PHHepMCGenEvent *hepmcevent = eventIter->second;

    if (hepmcevent)
    {
      /// Get the event characteristics, inherited from HepMC classes
      HepMC::GenEvent *truthevent = hepmcevent->getEvent();
      if (!truthevent)
      {
        std::cout << PHWHERE << "no evt pointer under phhepmcgeneventmap found " << std::endl;
        return;
      }

      /// Get the parton info
      HepMC::PdfInfo *pdfinfo = truthevent->pdf_info();

      /// Get the parton info as determined from HEPMC
      m_partid1 = pdfinfo->id1();
      m_partid2 = pdfinfo->id2();
      m_x1 = pdfinfo->x1();
      m_x2 = pdfinfo->x2();

      /// Are there multiple partonic intercations in a p+p event
      m_mpi = truthevent->mpi();

      /// Get the PYTHIA signal process id identifying the 2-to-2 hard process
      m_process_id = truthevent->signal_process_id();

      if (Verbosity() > 2)
      {
        std::cout << " Iterating over an event" << std::endl;
      }
      /// Loop over all the truth particles and get their information
      for (HepMC::GenEvent::particle_const_iterator iter = truthevent->particles_begin(); iter != truthevent->particles_end(); ++iter)
      {
        /// Get each pythia particle characteristics
        m_truthenergy = (*iter)->momentum().e();
        m_truthpid = (*iter)->pdg_id();

        m_trutheta = (*iter)->momentum().pseudoRapidity();
        m_truthphi = (*iter)->momentum().phi();
        m_truthpx = (*iter)->momentum().px();
        m_truthpy = (*iter)->momentum().py();
        m_truthpz = (*iter)->momentum().pz();
        m_truthpt = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);

        /// Fill the truth tree
        m_hepmctree->Fill();
        m_numparticlesinevent++;
      }
    }
  }
}

/**
 * This function collects the truth PHG4 stable particles that
 * are produced from PYTHIA, or some other event generator. These
 * are the stable particles, e.g. there are not any (for example)
 * partons here.
 */
void UPCMeson::getPHG4Truth()
{
  /// Get the primary particle range
  PHG4TruthInfoContainer::Range range = _truthinfo->GetPrimaryParticleRange();

  ROOT::Math::XYZTVector v1;

  /// Loop over the G4 truth (stable) particles
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    /// Get this truth particle
    const PHG4Particle *truth = iter->second;

    /// Get this particles momentum, etc.
    m_truthpx = truth->get_px();
    m_truthpy = truth->get_py();
    m_truthpz = truth->get_pz();
    m_truthp = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy + m_truthpz * m_truthpz);
    m_truthenergy = truth->get_e();

    v1.SetPxPyPzE( m_truthpx, m_truthpy, m_truthpz, m_truthenergy );
    m_truthpt = v1.Pt();
    m_truthphi = v1.Phi();
    m_trutheta = v1.Eta();

    /// Check for nans
    if (!std::isfinite(m_trutheta))
    {
      m_trutheta = -99;
    }
    m_truthpid = truth->get_pid();

    // Get the singleton instance of the PDG database
    TDatabasePDG* pdgDB = TDatabasePDG::Instance();

    // Get particle information by PDG code
    TParticlePDG* particle = pdgDB->GetParticle(m_truthpid);

    if (particle)
    {
      m_truthcharge = particle->Charge();
      if ( m_truthcharge != 0 )
      {
        m_ntrk_mc++;  // found charged track in sphenix accept
      }
      if ( (fabs(m_trutheta) < 1.1) && (m_truthpt>0.4) && (m_truthcharge!=0) )
      {
        m_ntrk_sphenix++;
        m_truthtree->Fill();
      }
    }
    else
    {
      std::cout << "Particle not found!" << std::endl;
    }

    //std::cout << "pid ch\t" << m_truthpid << "\t" << m_truthcharge << std::endl;

  }
}

/**
 * This method gets the tracks as reconstructed from the tracker.
 */
int UPCMeson::getTracks(PHCompositeNode *topNode)
{
  // make a cut on low ntracks
  m_ntrks = _trackmap->size();
  h_ntracks->Fill( m_ntrks );
  if (Verbosity() > 1)
  {
    std::cout << "ntracks " << m_ntrks << std::endl;
  }

  /*
  if ( m_ntrks > 3 || m_ntrks < 2 )
  //if ( m_ntrks != 2 )
  {
    return Fun4AllReturnCodes::DISCARDEVENT;
  }
  */

  /// EvalStack for truth track matching
  if ( !m_svtxEvalStack && _truthinfo )
  {
    std::cout << "getting svtx eval stack" << std::endl;
    m_svtxEvalStack = new SvtxEvalStack(topNode);
    m_svtxEvalStack->set_verbosity(Verbosity());
  }

  SvtxTrackEval *trackeval = nullptr;
  if ( m_svtxEvalStack )
  {
    m_svtxEvalStack->next_event(topNode);

    // Get the track evaluator (only available if eval is run during sim production)
    trackeval = m_svtxEvalStack->get_track_eval();
  }

  if (Verbosity() > 1)
  {
    std::cout << "Get the SVTX tracks " << m_ntrks << std::endl;
  }

  // get ntracks for each crossing (add CROSS_OFFSET = 150 to cross since it can be negative)
  std::array<int,700> ntracks_in_cross{0};
  std::array< std::vector<SvtxTrack*>, 700 > tracks_in_cross;
 
 
  for (auto &iter : *_trackmap)
  {
    SvtxTrack *track = iter.second;

    // Get ntracks per crossing
    m_cross = track->get_crossing();
    h_cross->Fill( m_cross );
    if ( static_cast<size_t>(m_cross+CROSS_OFFSET) > ntracks_in_cross.size() )
    {
      std::cout << "ERROR, cross too large " << m_cross << std::endl;
      continue;
    }
    ++ntracks_in_cross[m_cross+CROSS_OFFSET];
    tracks_in_cross[m_cross+CROSS_OFFSET].push_back( track );

    if ( Verbosity()>0 )
    {
      std::cout << "cross " << m_cross << "\ttrkid\t" << track->get_id() << std::endl;
    }

    /// Get the reconstructed track info
    m_tr_px = track->get_px();
    m_tr_py = track->get_py();
    m_tr_pz = track->get_pz();
    m_tr_p = sqrt(m_tr_px * m_tr_px + m_tr_py * m_tr_py + m_tr_pz * m_tr_pz);

    m_tr_pt = sqrt(m_tr_px * m_tr_px + m_tr_py * m_tr_py);

    // Make some cuts on the track to clean up sample
    if (m_tr_pt < 0.5)
    {
      continue;
    }
    m_tr_phi = track->get_phi();
    m_tr_eta = track->get_eta();

    m_charge = track->get_charge();
    m_chisq = track->get_chisq();
    m_ndf = track->get_ndf();
    m_dca = track->get_dca();
    m_tr_x = track->get_x();
    m_tr_y = track->get_y();
    m_tr_z = track->get_z();

    /// Get truth track info that matches this reconstructed track
    PHG4Particle *truthtrack = nullptr;
    if ( trackeval != nullptr )
    {
      truthtrack = trackeval->max_truth_particle_by_nclusters(track);
    }
    if ( truthtrack != nullptr )
    {
      m_truth_is_primary = _truthinfo->is_primary(truthtrack);

      m_truthtrackpx = truthtrack->get_px();
      m_truthtrackpy = truthtrack->get_py();
      m_truthtrackpz = truthtrack->get_pz();
      m_truthtrackp = std::sqrt(m_truthtrackpx * m_truthtrackpx + m_truthtrackpy * m_truthtrackpy + m_truthtrackpz * m_truthtrackpz);

      m_truthtracke = truthtrack->get_e();

      m_truthtrackpt = sqrt(m_truthtrackpx * m_truthtrackpx + m_truthtrackpy * m_truthtrackpy);
      m_truthtrackphi = atan(m_truthtrackpy / m_truthtrackpx);
      m_truthtracketa = atanh(m_truthtrackpz / m_truthtrackp);
      m_truthtrackpid = truthtrack->get_pid();
    }
    else
    {
      // Seems that we often miss the truth track?
      //std::cout << "Missing truth track" << std::endl;
      m_truth_is_primary = -9999;

      m_truthtrackpx = 0.;
      m_truthtrackpy = 0.;
      m_truthtrackpz = 0.;
      m_truthtrackp = 0.;

      m_truthtracke = 0.;

      m_truthtrackpt = 0.;
      m_truthtrackphi = 0.;
      m_truthtracketa = 0.;
      m_truthtrackpid = 0;
    }

    //m_tracktree->Fill();
  }

  ROOT::Math::XYZTVector v1, v2;

  // make pairs
  /*
  for (auto iter1 = _trackmap->begin(); iter1 != _trackmap->end(); iter1++)
  {
    for (auto iter2 = iter1; iter2 != _trackmap->end(); iter2++)
    {
      if ( iter2 == iter1 ) continue;

      //SvtxTrack *track2 = iter2.second;
      //std::cout << "XXX " << iter1->first << "\t" << iter2->first << std::endl;
      SvtxTrack *track1 = iter1->second;
      SvtxTrack *track2 = iter2->second;

      // same sign or opposite
      m_pq1 = track1->get_charge();
      m_pq2 = track2->get_charge();
      //std::cout << "charge " << m_pq1 << "\t" << m_pq2 << std::endl;
      int type = 0;
      if ( m_pq1*m_pq2 > 0 )
      {
        type = 1;
      }

      double px1 = track1->get_px();
      double py1 = track1->get_py();
      double pz1 = track1->get_pz();
      double e1 = sqrt( _mguess*_mguess + px1*px1 + py1*py1 + pz1*pz1 );
      v1.SetPxPyPzE( px1, py1, pz1, e1 );

      double px2 = track2->get_px();
      double py2 = track2->get_py();
      double pz2 = track2->get_pz();
      double e2 = sqrt( _mguess*_mguess + px2*px2 + py2*py2 + pz2*pz2 );
      v2.SetPxPyPzE( px2, py2, pz2, e2 );

      //TLorentzVector sum = v1 + v2;
      ROOT::Math::XYZTVector sum = v1 + v2;
      m_pm = sum.M();
      m_ppt = sum.Pt();
      m_pphi = sum.Phi();
      m_py = sum.Rapidity();
      m_peta = sum.Eta();
      //m_pdphi = ROOT::Math::VectorUtil::DeltaPhi(v1,v2);
      m_pdphi = ROOT::Math::VectorUtil::DeltaPhi(v1,v2);
      m_ppt1 = v1.Pt();
      m_ppz1 = v1.Pz();
      m_pphi1 = v1.Phi();
      m_peta1 = v1.Eta();
      m_ppt2 = v2.Pt();
      m_ppz2 = v2.Pz();
      m_pphi2 = v2.Phi();
      m_peta2 = v2.Eta();

      h_mass[type]->Fill( m_pm );
      h_pt[type]->Fill( m_ppt );
      h_y[type]->Fill( m_py );
      h_eta[type]->Fill( m_peta );
      h2_eta_phi[type]->Fill( m_peta, m_pphi );
      h_phi[type]->Fill( m_pphi );

      m_pairtree->Fill();
    }
  }
  */

  for (size_t icross=0; icross<ntracks_in_cross.size(); icross++)
  {
    m_ntrks_cross = ntracks_in_cross[icross];
    h_cross_evt->Fill( icross - CROSS_OFFSET );

    //std::cout << "icross " << icross << "\t" << ntrks << std::endl;
    if ( (m_ntrks_cross==2) || (m_ntrks_cross==3) )
    {
      m_cross = icross - CROSS_OFFSET;

      std::vector<SvtxTrack*> &t = tracks_in_cross[icross];
      //std::cout << t.size() << std::endl;
      SvtxTrack *track1 = t[0];
      SvtxTrack *track2 = t[1];

      // same sign or opposite
      m_pq1 = track1->get_charge();
      m_pq2 = track2->get_charge();
      //std::cout << "charge " << m_pq1 << "\t" << m_pq2 << std::endl;
      int type = 0;   // unlike sign = 0, like sign = 1
      if ( m_pq1*m_pq2 > 0 )
      {
        type = 1; // like sign
      }

      double px1 = track1->get_px();
      double py1 = track1->get_py();
      double pz1 = track1->get_pz();
      double e1 = sqrt( _mguess*_mguess + px1*px1 + py1*py1 + pz1*pz1 );
      v1.SetPxPyPzE( px1, py1, pz1, e1 );

      double px2 = track2->get_px();
      double py2 = track2->get_py();
      double pz2 = track2->get_pz();
      double e2 = sqrt( _mguess*_mguess + px2*px2 + py2*py2 + pz2*pz2 );
      v2.SetPxPyPzE( px2, py2, pz2, e2 );

      //TLorentzVector sum = v1 + v2;
      ROOT::Math::XYZTVector sum = v1 + v2;
      m_pm = sum.M();
      m_ppt = sum.Pt();
      m_pphi = sum.Phi();
      m_py = sum.Rapidity();
      m_peta = sum.Eta();
      //m_pdphi = ROOT::Math::VectorUtil::DeltaPhi(v1,v2);
      m_pdphi = ROOT::Math::VectorUtil::DeltaPhi(v1,v2);
      m_ppt1 = v1.Pt();
      m_ppz1 = v1.Pz();
      m_pphi1 = v1.Phi();
      m_peta1 = v1.Eta();
      m_ppt2 = v2.Pt();
      m_ppz2 = v2.Pz();
      m_pphi2 = v2.Phi();
      m_peta2 = v2.Eta();

      if ( m_cross != 0 )
      {
        h_mass[type]->Fill( m_pm );
        h_pt[type]->Fill( m_ppt );
        h_y[type]->Fill( m_py );
        h_eta[type]->Fill( m_peta );
        h2_eta_phi[type]->Fill( m_peta, m_pphi );
        h_phi[type]->Fill( m_pphi );
      }

      m_pairtree->Fill();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}



/**
 * This function puts all of the tree branch assignments in one place so as to not
 * clutter up the UPCMeson::Init function.
 */
void UPCMeson::initializeTrees()
{
  m_tracktree = new TTree("tracktree", "A tree with svtx tracks");
  m_tracktree->Branch("px", &m_tr_px, "m_tr_px/D");
  m_tracktree->Branch("py", &m_tr_py, "m_tr_py/D");
  m_tracktree->Branch("pz", &m_tr_pz, "m_tr_pz/D");
  m_tracktree->Branch("p", &m_tr_p, "m_tr_p/D");
  m_tracktree->Branch("pt", &m_tr_pt, "m_tr_pt/D");
  m_tracktree->Branch("phi", &m_tr_phi, "m_tr_phi/D");
  m_tracktree->Branch("eta", &m_tr_eta, "m_tr_eta/D");
  m_tracktree->Branch("q", &m_charge, "m_charge/I");
  m_tracktree->Branch("chisq", &m_chisq, "m_chisq/D");
  m_tracktree->Branch("ndf", &m_ndf, "m_ndf/I");
  m_tracktree->Branch("dca", &m_dca, "m_dca/D");
  m_tracktree->Branch("x", &m_tr_x, "m_tr_x/D");
  m_tracktree->Branch("y", &m_tr_y, "m_tr_y/D");
  m_tracktree->Branch("z", &m_tr_z, "m_tr_z/D");
  m_tracktree->Branch("cross", &m_cross, "m_cross/S");
  m_tracktree->Branch("truth_is_primary", &m_truth_is_primary, "m_truth_is_primary/I");
  m_tracktree->Branch("trupx", &m_truthtrackpx, "m_truthtrackpx/D");
  m_tracktree->Branch("trupy", &m_truthtrackpy, "m_truthtrackpy/D");
  m_tracktree->Branch("trupz", &m_truthtrackpz, "m_truthtrackpz/D");
  m_tracktree->Branch("trup", &m_truthtrackp, "m_truthtrackp/D");
  m_tracktree->Branch("true", &m_truthtracke, "m_truthtracke/D");
  m_tracktree->Branch("trupt", &m_truthtrackpt, "m_truthtrackpt/D");
  m_tracktree->Branch("truphi", &m_truthtrackphi, "m_truthtrackphi/D");
  m_tracktree->Branch("trueta", &m_truthtracketa, "m_truthtracketa/D");
  m_tracktree->Branch("trupid", &m_truthtrackpid, "m_truthtrackpid/I");

  m_globaltree = new TTree("globaltree", "Global Info");
  m_globaltree->Branch("run", &m_run, "run/I");
  m_globaltree->Branch("evt", &m_evt, "evt/I");
  m_globaltree->Branch("npart", &m_npart, "npart/I");
  m_globaltree->Branch("ncoll", &m_ncoll, "ncoll/I");
  m_globaltree->Branch("b", &m_bimpact, "b/F");
  m_globaltree->Branch("totntrks", &m_ntrks, "tottrks/I");
  m_globaltree->Branch("sphntrks", &m_ntrk_sphenix, "sphntrks/I");
  m_globaltree->Branch("hijntrks", &m_ntrk_mc, "mcntrks/I");

  m_truthtree = new TTree("truthg4tree", "A tree with truth g4 particles");
  m_truthtree->Branch("evt", &m_evt, "evt/I");
  m_truthtree->Branch("te", &m_truthenergy, "m_truthe/D");
  m_truthtree->Branch("tp", &m_truthp, "m_truthp/D");
  m_truthtree->Branch("tpx", &m_truthpx, "m_truthpx/D");
  m_truthtree->Branch("tpy", &m_truthpy, "m_truthpy/D");
  m_truthtree->Branch("tpz", &m_truthpz, "m_truthpz/D");
  m_truthtree->Branch("tpt", &m_truthpt, "m_truthpt/D");
  m_truthtree->Branch("tphi", &m_truthphi, "m_truthphi/D");
  m_truthtree->Branch("teta", &m_trutheta, "m_trutheta/D");
  m_truthtree->Branch("tpid", &m_truthpid, "m_truthpid/I");
  m_truthtree->Branch("tq", &m_truthcharge, "m_truthcharge/I");

  m_pairtree = new TTree("pairs", "pairs");
  m_pairtree->Branch("run", &m_run, "run/I");
  m_pairtree->Branch("evt", &m_evt, "evt/I");
  m_pairtree->Branch("cross", &m_cross, "cross/S");
  m_pairtree->Branch("bunch", &m_bunch, "bunch/S");
  m_pairtree->Branch("strig", &m_strig, "strig/l");
  m_pairtree->Branch("ntrks", &m_ntrks_cross, "ntrks/S"); // ntrks in crossing
  m_pairtree->Branch("m", &m_pm, "m/F");
  m_pairtree->Branch("pt", &m_ppt, "pt/F");
  m_pairtree->Branch("phi", &m_pphi, "phi/F");
  m_pairtree->Branch("y", &m_py, "y/F");
  m_pairtree->Branch("eta", &m_peta, "eta/F");
  m_pairtree->Branch("dphi", &m_pdphi, "dphi/F");
  m_pairtree->Branch("pt1", &m_ppt1, "pt1/F");
  m_pairtree->Branch("pz1", &m_ppz1, "pz1/F");
  m_pairtree->Branch("phi1", &m_pphi1, "phi1/F");
  m_pairtree->Branch("eta1", &m_peta1, "eta1/F");
  m_pairtree->Branch("pt2", &m_ppt2, "pt2/F");
  m_pairtree->Branch("pz2", &m_ppz2, "pz2/F");
  m_pairtree->Branch("phi2", &m_pphi2, "phi2/F");
  m_pairtree->Branch("eta2", &m_peta2, "eta2/F");
  m_pairtree->Branch("q1", &m_pq1, "q1/S");
  m_pairtree->Branch("q2", &m_pq2, "q2/S");

}

/**
 * This function initializes all of the member variables in this class so that there
 * are no variables that might not be set before e.g. writing them to the output
 * trees.
 */
void UPCMeson::initializeVariables()
{
  m_partid1 = -99;
  m_partid2 = -99;
  m_x1 = -99;
  m_x2 = -99;
  m_mpi = -99;
  m_process_id = -99;
  m_truthenergy = -99;
  m_trutheta = -99;
  m_truthphi = -99;
  m_truthp = -99;
  m_truthpx = -99;
  m_truthpy = -99;
  m_truthpz = -99;
  m_truthpt = -99;
  m_numparticlesinevent = -99;
  m_truthpid = -99;

  m_tr_px = -99;
  m_tr_py = -99;
  m_tr_pz = -99;
  m_tr_p = -99;
  m_tr_pt = -99;
  m_tr_phi = -99;
  m_tr_eta = -99;
  m_charge = -99;
  m_chisq = -99;
  m_ndf = -99;
  m_dca = -99;
  m_tr_x = -99;
  m_tr_y = -99;
  m_tr_z = -99;
  m_truth_is_primary = -99;
  m_truthtrackpx = -99;
  m_truthtrackpy = -99;
  m_truthtrackpz = -99;
  m_truthtrackp = -99;
  m_truthtracke = -99;
  m_truthtrackpt = -99;
  m_truthtrackphi = -99;
  m_truthtracketa = -99;
  m_truthtrackpid = -99;

  m_ntrks = 0;
  m_ntrk_sphenix = 0;
  m_ntrk_mc = 0;

}

int UPCMeson::Reset(PHCompositeNode * /*topNode*/)
{
  //std::cout << "In Reset()" << std::endl;
  initializeVariables();
  return 0;
}

