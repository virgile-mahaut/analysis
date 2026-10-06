#include "lwtrackntuplizer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackbase_historic/TrackSeed.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

#include <cmath>
#include <climits>
#include <iostream>
#include <limits>

namespace
{
template <class Container>
void Clean(Container& c)
{
  Container().swap(c);
}

float nan_value()
{
  return std::numeric_limits<float>::quiet_NaN();
}
}  // namespace

lwtrackntuplizer::lwtrackntuplizer(const std::string& name)
  : SubsysReco(name)
{
}

lwtrackntuplizer::~lwtrackntuplizer() = default;

int lwtrackntuplizer::Init(PHCompositeNode* /*topNode*/)
{
  m_outFile = new TFile(m_outFileName.c_str(), "RECREATE");
  if (!m_outFile || m_outFile->IsZombie())
  {
    std::cout << PHWHERE << " failed to create output file " << m_outFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_outFile->SetCompressionLevel(7);
  SetupTree();
  Cleanup();

  return Fun4AllReturnCodes::EVENT_OK;
}

int lwtrackntuplizer::InitRun(PHCompositeNode* topNode)
{
  const int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int lwtrackntuplizer::process_event(PHCompositeNode* topNode)
{
  const int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  Cleanup();

  if (m_trackMap)
  {
    for (const auto& entry : *m_trackMap)
    {
      const auto* track = entry.second;
      if (!track)
      {
        continue;
      }
      FillTrack(track);
    }
  }

  m_nTracks = m_trackID.size();
  if (m_outTree)
  {
    m_outTree->Fill();
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}

int lwtrackntuplizer::ResetEvent(PHCompositeNode* /*topNode*/)
{
  Cleanup();
  return Fun4AllReturnCodes::EVENT_OK;
}

int lwtrackntuplizer::EndRun(const int /*runnumber*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int lwtrackntuplizer::End(PHCompositeNode* /*topNode*/)
{
  if (m_outFile)
  {
    m_outFile->cd();
    if (m_outTree)
    {
      m_outTree->Write("", TObject::kOverwrite);
    }
    m_outFile->Close();
    delete m_outFile;
    m_outFile = nullptr;
    m_outTree = nullptr;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int lwtrackntuplizer::Reset(PHCompositeNode* /*topNode*/)
{
  Cleanup();
  return Fun4AllReturnCodes::EVENT_OK;
}

void lwtrackntuplizer::Print(const std::string& what) const
{
  std::cout << "lwtrackntuplizer::Print - " << what
            << ", output=" << m_outFileName
            << ", tree=" << m_treeName
            << ", trackmap=" << m_trackMapName
            << std::endl;
}

int lwtrackntuplizer::GetNodes(PHCompositeNode* topNode)
{
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!m_trackMap)
  {
    std::cout << PHWHERE << " missing required node " << m_trackMapName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  if (!m_clusterMap)
  {
    std::cout << PHWHERE << " missing required node " << m_clusterContainerName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, m_geometryNodeName);
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << " missing required node " << m_geometryNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tpcGeomContainer = findNode::getClass<PHG4TpcGeomContainer>(topNode, m_tpcGeomNodeName);
  if (!m_tpcGeomContainer)
  {
    std::cout << PHWHERE << " missing required node " << m_tpcGeomNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName);

  return Fun4AllReturnCodes::EVENT_OK;
}

void lwtrackntuplizer::SetupTree()
{
  if (!m_outFile)
  {
    return;
  }

  m_outFile->cd();
  m_outTree = new TTree(m_treeName.c_str(), m_treeName.c_str());

  m_outTree->Branch("event", &m_event, "event/i");
  m_outTree->Branch("nTracks", &m_nTracks, "nTracks/i");

  m_outTree->Branch("trackID", &m_trackID);
  m_outTree->Branch("crossing", &m_crossing);
  m_outTree->Branch("px", &m_px);
  m_outTree->Branch("py", &m_py);
  m_outTree->Branch("pz", &m_pz);
  m_outTree->Branch("pt", &m_pt);
  m_outTree->Branch("eta", &m_eta);
  m_outTree->Branch("phi", &m_phi);
  m_outTree->Branch("deltapt", &m_deltapt);
  m_outTree->Branch("deltaeta", &m_deltaeta);
  m_outTree->Branch("deltaphi", &m_deltaphi);
  m_outTree->Branch("charge", &m_charge);
  m_outTree->Branch("quality", &m_quality);
  m_outTree->Branch("chisq", &m_chisq);
  m_outTree->Branch("ndf", &m_ndf);
  m_outTree->Branch("nhits", &m_nhits);
  m_outTree->Branch("nmaps", &m_nmaps);
  m_outTree->Branch("nintt", &m_nintt);
  m_outTree->Branch("ntpc", &m_ntpc);
  m_outTree->Branch("nmms", &m_nmms);
  m_outTree->Branch("ntpc1", &m_ntpc1);
  m_outTree->Branch("ntpc11", &m_ntpc11);
  m_outTree->Branch("ntpc2", &m_ntpc2);
  m_outTree->Branch("ntpc3", &m_ntpc3);
  m_outTree->Branch("dedx", &m_dedx);
  m_outTree->Branch("pidedx", &m_pidedx);
  m_outTree->Branch("kdedx", &m_kdedx);
  m_outTree->Branch("prdedx", &m_prdedx);
  m_outTree->Branch("vertexID", &m_vertexID);
  m_outTree->Branch("vx", &m_vx);
  m_outTree->Branch("vy", &m_vy);
  m_outTree->Branch("vz", &m_vz);
  m_outTree->Branch("dca2d", &m_dca2d);
  m_outTree->Branch("dca2dsigma", &m_dca2dsigma);
  m_outTree->Branch("dca3dxy", &m_dca3dxy);
  m_outTree->Branch("dca3dxysigma", &m_dca3dxysigma);
  m_outTree->Branch("dca3dz", &m_dca3dz);
  m_outTree->Branch("dca3dzsigma", &m_dca3dzsigma);
  m_outTree->Branch("pcax", &m_pcax);
  m_outTree->Branch("pcay", &m_pcay);
  m_outTree->Branch("pcaz", &m_pcaz);
  m_outTree->Branch("hlxpt", &m_hlxpt);
  m_outTree->Branch("hlxeta", &m_hlxeta);
  m_outTree->Branch("hlxphi", &m_hlxphi);
  m_outTree->Branch("hlxX0", &m_hlxX0);
  m_outTree->Branch("hlxY0", &m_hlxY0);
  m_outTree->Branch("hlxZ0", &m_hlxZ0);
  m_outTree->Branch("hlxcharge", &m_hlxcharge);
}

void lwtrackntuplizer::FillTrack(const SvtxTrack* track)
{
  const float nan = nan_value();

  const auto* tpcseed = track->get_tpc_seed();
  const auto* silseed = track->get_silicon_seed();

  int nhits = 0;
  int nmaps = 0;
  int nintt = 0;
  int ntpc = 0;
  int nmms = 0;
  int ntpc1 = 0;
  int ntpc11 = 0;
  int ntpc2 = 0;
  int ntpc3 = 0;

  const auto count_seed = [&](const TrackSeed* seed)
  {
    if (!seed)
    {
      return;
    }

    nhits += seed->size_cluster_keys();
    for (auto iter = seed->begin_cluster_keys(); iter != seed->end_cluster_keys(); ++iter)
    {
      const auto clusterKey = *iter;
      const unsigned int layer = TrkrDefs::getLayer(clusterKey);
      switch (TrkrDefs::getTrkrId(clusterKey))
      {
      case TrkrDefs::TrkrId::mvtxId:
        ++nmaps;
        break;
      case TrkrDefs::TrkrId::inttId:
        ++nintt;
        break;
      case TrkrDefs::TrkrId::tpcId:
        ++ntpc;
        if ((layer - 7U) < 8U)
        {
          ++ntpc11;
        }
        if ((layer - 7U) < 16U)
        {
          ++ntpc1;
        }
        else if ((layer - 7U) < 32U)
        {
          ++ntpc2;
        }
        else if ((layer - 7U) < 48U)
        {
          ++ntpc3;
        }
        break;
      case TrkrDefs::TrkrId::micromegasId:
        ++nmms;
        break;
      default:
        break;
      }
    }
  };

  count_seed(tpcseed);
  count_seed(silseed);

  float dedx = nan;
  if (tpcseed)
  {
    auto* inner1 = m_tpcGeomContainer->GetLayerCellGeom(7);
    auto* inner2 = m_tpcGeomContainer->GetLayerCellGeom(8);
    auto* middle = m_tpcGeomContainer->GetLayerCellGeom(27);
    auto* outer = m_tpcGeomContainer->GetLayerCellGeom(50);

    if (inner1 && inner2 && middle && outer)
    {
      float layerThicknesses[4] = {
        static_cast<float>(inner1->get_thickness()),
        static_cast<float>(inner2->get_thickness()),
        static_cast<float>(middle->get_thickness()),
        static_cast<float>(outer->get_thickness())};
      dedx = TrackAnalysisUtils::calc_dedx(
          const_cast<TrackSeed*>(tpcseed), m_clusterMap, m_tGeometry, layerThicknesses);
    }
  }

  const float px = track->get_px();
  const float py = track->get_py();
  const float pz = track->get_pz();
  const TVector3 momentum(px, py, pz);
  const float pt = momentum.Pt();
  const float eta = momentum.Eta();
  const float phi = momentum.Phi();

  const float cvxx = track->get_error(3, 3);
  const float cvxy = track->get_error(3, 4);
  const float cvxz = track->get_error(3, 5);
  const float cvyy = track->get_error(4, 4);
  const float cvyz = track->get_error(4, 5);
  const float cvzz = track->get_error(5, 5);

  const double pt2 = static_cast<double>(px) * px + static_cast<double>(py) * py;
  const double p2 = pt2 + static_cast<double>(pz) * pz;

  float deltapt = nan;
  if (pt2 > 0.)
  {
    const double arg = (static_cast<double>(cvxx) * px * px + 2. * static_cast<double>(cvxy) * px * py + static_cast<double>(cvyy) * py * py) / pt2;
    if (std::isfinite(arg) && arg >= 0.)
    {
      deltapt = std::sqrt(arg);
    }
  }

  float deltaeta = nan;
  if (pt2 > 0. && p2 > 0.)
  {
    const double numerator =
        static_cast<double>(cvzz) * pt2 * pt2 +
        static_cast<double>(pz) *
            (-2. * (static_cast<double>(cvxz) * px + static_cast<double>(cvyz) * py) * pt2 +
             static_cast<double>(cvxx) * px * px * pz +
             static_cast<double>(cvyy) * py * py * pz +
             2. * static_cast<double>(cvxy) * px * py * pz);
    const double denominator = pt2 * pt2 * p2;
    const double arg = numerator / denominator;
    if (std::isfinite(arg) && arg >= 0.)
    {
      deltaeta = std::sqrt(arg);
    }
  }

  float deltaphi = nan;
  if (pt2 > 0.)
  {
    const double denominator = pt2 * pt2;
    const double arg =
        (static_cast<double>(cvyy) * px * px - 2. * static_cast<double>(cvxy) * px * py + static_cast<double>(cvxx) * py * py) /
        denominator;
    if (std::isfinite(arg) && arg >= 0.)
    {
      deltaphi = std::sqrt(arg);
    }
  }

  const int vertexID = track->get_vertex_id();
  float vx = nan;
  float vy = nan;
  float vz = nan;
  if (vertexID >= 0 && m_vertexMap)
  {
    auto vertexIter = m_vertexMap->find(vertexID);
    if (vertexIter != m_vertexMap->end() && vertexIter->second)
    {
      vx = vertexIter->second->get_x();
      vy = vertexIter->second->get_y();
      vz = vertexIter->second->get_z();
    }
  }

  float hlxpt = nan;
  float hlxeta = nan;
  float hlxphi = nan;
  float hlxX0 = nan;
  float hlxY0 = nan;
  float hlxZ0 = nan;
  int hlxcharge = 0;
  if (tpcseed)
  {
    // The TPC seed stores the helix-fit parameters that TrkrNtuplizer labels as hlx* fields
    hlxpt = tpcseed->get_pt();
    hlxeta = tpcseed->get_eta();
    hlxphi = tpcseed->get_phi();
    hlxX0 = tpcseed->get_X0();
    hlxY0 = tpcseed->get_Y0();
    hlxZ0 = tpcseed->get_Z0();
    if (std::isfinite(tpcseed->get_qOverR()) && tpcseed->get_qOverR() != 0.)
    {
      hlxcharge = (tpcseed->get_qOverR() > 0.) ? 1 : -1;
    }
  }

  m_trackID.push_back(track->get_id());
  m_crossing.push_back(track->get_crossing() == SHRT_MAX ? nan : static_cast<float>(track->get_crossing()));
  m_px.push_back(px);
  m_py.push_back(py);
  m_pz.push_back(pz);
  m_pt.push_back(pt);
  m_eta.push_back(eta);
  m_phi.push_back(phi);
  m_deltapt.push_back(deltapt);
  m_deltaeta.push_back(deltaeta);
  m_deltaphi.push_back(deltaphi);
  m_charge.push_back(track->get_charge());
  m_quality.push_back(track->get_quality());
  m_chisq.push_back(track->get_chisq());
  m_ndf.push_back(track->get_ndf());
  m_nhits.push_back(nhits);
  m_nmaps.push_back(nmaps);
  m_nintt.push_back(nintt);
  m_ntpc.push_back(ntpc);
  m_nmms.push_back(nmms);
  m_ntpc1.push_back(ntpc1);
  m_ntpc11.push_back(ntpc11);
  m_ntpc2.push_back(ntpc2);
  m_ntpc3.push_back(ntpc3);
  m_dedx.push_back(dedx);
  m_pidedx.push_back(nan);
  m_kdedx.push_back(nan);
  m_prdedx.push_back(nan);
  m_vertexID.push_back(vertexID);
  m_vx.push_back(vx);
  m_vy.push_back(vy);
  m_vz.push_back(vz);
  m_dca2d.push_back(track->get_dca2d());
  m_dca2dsigma.push_back(track->get_dca2d_error());
  m_dca3dxy.push_back(track->get_dca3d_xy());
  m_dca3dxysigma.push_back(track->get_dca3d_xy_error());
  m_dca3dz.push_back(track->get_dca3d_z());
  m_dca3dzsigma.push_back(track->get_dca3d_z_error());
  m_pcax.push_back(track->get_x());
  m_pcay.push_back(track->get_y());
  m_pcaz.push_back(track->get_z());
  m_hlxpt.push_back(hlxpt);
  m_hlxeta.push_back(hlxeta);
  m_hlxphi.push_back(hlxphi);
  m_hlxX0.push_back(hlxX0);
  m_hlxY0.push_back(hlxY0);
  m_hlxZ0.push_back(hlxZ0);
  m_hlxcharge.push_back(hlxcharge);
}

void lwtrackntuplizer::Cleanup()
{
  m_nTracks = 0;

  Clean(m_trackID);
  Clean(m_crossing);
  Clean(m_px);
  Clean(m_py);
  Clean(m_pz);
  Clean(m_pt);
  Clean(m_eta);
  Clean(m_phi);
  Clean(m_deltapt);
  Clean(m_deltaeta);
  Clean(m_deltaphi);
  Clean(m_charge);
  Clean(m_quality);
  Clean(m_chisq);
  Clean(m_ndf);
  Clean(m_nhits);
  Clean(m_nmaps);
  Clean(m_nintt);
  Clean(m_ntpc);
  Clean(m_nmms);
  Clean(m_ntpc1);
  Clean(m_ntpc11);
  Clean(m_ntpc2);
  Clean(m_ntpc3);
  Clean(m_dedx);
  Clean(m_pidedx);
  Clean(m_kdedx);
  Clean(m_prdedx);
  Clean(m_vertexID);
  Clean(m_vx);
  Clean(m_vy);
  Clean(m_vz);
  Clean(m_dca2d);
  Clean(m_dca2dsigma);
  Clean(m_dca3dxy);
  Clean(m_dca3dxysigma);
  Clean(m_dca3dz);
  Clean(m_dca3dzsigma);
  Clean(m_pcax);
  Clean(m_pcay);
  Clean(m_pcaz);
  Clean(m_hlxpt);
  Clean(m_hlxeta);
  Clean(m_hlxphi);
  Clean(m_hlxX0);
  Clean(m_hlxY0);
  Clean(m_hlxZ0);
  Clean(m_hlxcharge);
}
