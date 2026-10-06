// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef LWTRACKNTUPLIZER_H
#define LWTRACKNTUPLIZER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class PHG4TpcGeomContainer;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class TFile;
class TTree;
class TrkrClusterContainer;

class lwtrackntuplizer : public SubsysReco
{
 public:
  explicit lwtrackntuplizer(const std::string& name = "lwtrackntuplizer");
  ~lwtrackntuplizer() override;

  int Init(PHCompositeNode* topNode) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int ResetEvent(PHCompositeNode* topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode* topNode) override;
  int Reset(PHCompositeNode* topNode) override;

  void Print(const std::string& what = "ALL") const override;

  void setOutputName(const std::string& name) { m_outFileName = name; }
  void setTreeName(const std::string& name) { m_treeName = name; }
  void setClusterContainerName(const std::string& name) { m_clusterContainerName = name; }
  void setTrackMapName(const std::string& name) { m_trackMapName = name; }
  void setVertexMapName(const std::string& name) { m_vertexMapName = name; }

 private:
  int GetNodes(PHCompositeNode* topNode);
  void SetupTree();
  void FillTrack(const SvtxTrack* track);
  void Cleanup();

  TFile* m_outFile = nullptr;
  TTree* m_outTree = nullptr;

  std::string m_outFileName = "lwtrackntuple.root";
  std::string m_treeName = "TrackTree";

  std::string m_trackMapName = "SvtxTrackMap";
  std::string m_vertexMapName = "SvtxVertexMap";
  std::string m_clusterContainerName = "TRKR_CLUSTER_SEED";
  std::string m_geometryNodeName = "ActsGeometry";
  std::string m_tpcGeomNodeName = "TPCGEOMCONTAINER";

  SvtxTrackMap* m_trackMap = nullptr;
  SvtxVertexMap* m_vertexMap = nullptr;
  TrkrClusterContainer* m_clusterMap = nullptr;
  ActsGeometry* m_tGeometry = nullptr;
  PHG4TpcGeomContainer* m_tpcGeomContainer = nullptr;

  unsigned int m_event = 0;
  unsigned int m_nTracks = 0;

  std::vector<int> m_trackID;
  std::vector<float> m_crossing;
  std::vector<float> m_px;
  std::vector<float> m_py;
  std::vector<float> m_pz;
  std::vector<float> m_pt;
  std::vector<float> m_eta;
  std::vector<float> m_phi;
  std::vector<float> m_deltapt;
  std::vector<float> m_deltaeta;
  std::vector<float> m_deltaphi;
  std::vector<int> m_charge;
  std::vector<float> m_quality;
  std::vector<float> m_chisq;
  std::vector<int> m_ndf;
  std::vector<int> m_nhits;
  std::vector<int> m_nmaps;
  std::vector<int> m_nintt;
  std::vector<int> m_ntpc;
  std::vector<int> m_nmms;
  std::vector<int> m_ntpc1;
  std::vector<int> m_ntpc11;
  std::vector<int> m_ntpc2;
  std::vector<int> m_ntpc3;
  std::vector<float> m_dedx;
  std::vector<float> m_pidedx;
  std::vector<float> m_kdedx;
  std::vector<float> m_prdedx;
  std::vector<int> m_vertexID;
  std::vector<float> m_vx;
  std::vector<float> m_vy;
  std::vector<float> m_vz;
  std::vector<float> m_dca2d;
  std::vector<float> m_dca2dsigma;
  std::vector<float> m_dca3dxy;
  std::vector<float> m_dca3dxysigma;
  std::vector<float> m_dca3dz;
  std::vector<float> m_dca3dzsigma;
  std::vector<float> m_pcax;
  std::vector<float> m_pcay;
  std::vector<float> m_pcaz;
  std::vector<float> m_hlxpt;
  std::vector<float> m_hlxeta;
  std::vector<float> m_hlxphi;
  std::vector<float> m_hlxX0;
  std::vector<float> m_hlxY0;
  std::vector<float> m_hlxZ0;
  std::vector<int> m_hlxcharge;
};

#endif  // LWTRACKNTUPLIZER_H
