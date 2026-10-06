#ifndef VERTEX_H
#define VERTEX_H

#include <Riostream.h>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdlib.h>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TRandom3.h>
#include <TText.h>
#include <TTree.h>
#include <TVector3.h>

using namespace std;

// structure to hold vertex information (MBD and silicon)
struct vertex
{
    bool trigger12;
    float beamspot_x;
    float beamspot_y;
    float mbd_z;
    float silicon_z_crossing0;
    vector<float> silicon_z_nonzeroCrossing;

    // constructor
    vertex() : trigger12(false), beamspot_x(-9999.0), beamspot_y(-9999.0), mbd_z(-9999.0), silicon_z_crossing0(-9999.0) {}

    void reset()
    {
        trigger12 = false;
        beamspot_x = -9999.0;
        beamspot_y = -9999.0;
        mbd_z = -9999.0;
        silicon_z_crossing0 = -9999.0;
        silicon_z_nonzeroCrossing.clear();
    }

    void identify()
    {
        std::cout << "----- Vertex Information -----" << std::endl;
        std::cout << "trigger12: " << trigger12 << std::endl;
        std::cout << "beamspot_x: " << beamspot_x << std::endl;
        std::cout << "beamspot_y: " << beamspot_y << std::endl;
        std::cout << "mbd_z: " << mbd_z << std::endl;
        std::cout << "silicon_z_crossing0: " << silicon_z_crossing0 << std::endl;
        std::cout << "silicon_z_nonzeroCrossing: ";
        for (const auto &z : silicon_z_nonzeroCrossing)
        {
            std::cout << z << " ";
        }
        std::cout << std::endl;
        std::cout << "------------------------------" << std::endl;
    }

};


struct VtxData
{
    bool isdata, is_min_bias, isGoodVtx;
    vector<int> firedTriggers;
    bool firedTrig12_vtxle10cm;
    int event, NClusLayer1, NTruthVtx;
    uint64_t INTT_BCO, gl1bco, bcotr;
    float TruthPV_x, TruthPV_y, TruthPV_z;
    float PV_x, PV_y, PV_z;
    float sigmaGaus_PVz, sigmaGaus_err_PVz; // from the Gaussian fit
    float maxGroup_ratio, maxGroup_width;   // from the max group method
    float PV_z_rand; // For acceptance correction
    float Centrality_bimp, Centrality_impactparam, Centrality_mbd;
    float mbd_south_charge_sum, mbd_north_charge_sum, mbd_charge_sum, mbd_charge_asymm, mbd_z_vtx;
    float SiliconVertexZ_crossing0;
};

void SetVtxMinitree(TTree *outTree, VtxData &vtxdata)
{
    outTree->Branch("counter", &vtxdata.event);
    outTree->Branch("INTT_BCO", &vtxdata.INTT_BCO);
    outTree->Branch("gl1bco", &vtxdata.gl1bco);
    outTree->Branch("bcotr", &vtxdata.bcotr);
    outTree->Branch("is_min_bias", &vtxdata.is_min_bias);
    outTree->Branch("isGoodVtx", &vtxdata.isGoodVtx);
    outTree->Branch("firedTriggers", &vtxdata.firedTriggers);
    outTree->Branch("firedTrig12_vtxle10cm", &vtxdata.firedTrig12_vtxle10cm);
    outTree->Branch("NClusLayer1", &vtxdata.NClusLayer1);
    if (!vtxdata.isdata)
    {
        outTree->Branch("NTruthVtx", &vtxdata.NTruthVtx);
        outTree->Branch("TruthPV_x", &vtxdata.TruthPV_x);
        outTree->Branch("TruthPV_y", &vtxdata.TruthPV_y);
        outTree->Branch("TruthPV_z", &vtxdata.TruthPV_z);
        outTree->Branch("Centrality_bimp", &vtxdata.Centrality_bimp);
        outTree->Branch("Centrality_impactparam", &vtxdata.Centrality_impactparam);
    }
    outTree->Branch("MBD_centrality", &vtxdata.Centrality_mbd);
    outTree->Branch("mbd_south_charge_sum", &vtxdata.mbd_south_charge_sum);
    outTree->Branch("mbd_north_charge_sum", &vtxdata.mbd_north_charge_sum);
    outTree->Branch("mbd_charge_sum", &vtxdata.mbd_charge_sum);
    outTree->Branch("mbd_charge_asymm", &vtxdata.mbd_charge_asymm);
    outTree->Branch("mbd_z_vtx", &vtxdata.mbd_z_vtx);
    outTree->Branch("PV_x", &vtxdata.PV_x);
    outTree->Branch("PV_y", &vtxdata.PV_y);
    outTree->Branch("PV_z", &vtxdata.PV_z);
    outTree->Branch("sigmaGaus_PVz", &vtxdata.sigmaGaus_PVz);
    outTree->Branch("sigmaGaus_err_PVz", &vtxdata.sigmaGaus_err_PVz);
    outTree->Branch("maxGroup_ratio", &vtxdata.maxGroup_ratio);
    outTree->Branch("maxGroup_width", &vtxdata.maxGroup_width);
    outTree->Branch("SiliconVertexZ_crossing0", &vtxdata.SiliconVertexZ_crossing0);
}

std::map<int, vertex> EvtVtx_map_event(const char *vtxfname)
{
    std::map<int, vertex> EvtVtx_map;

    vertex vtx;

    TFile *f = new TFile(vtxfname, "READ");
    TTree *t = (TTree *)f->Get("vertex");
    int event;
    float PV_x, PV_y, PV_z;
    t->SetBranchAddress("counter", &event);
    t->SetBranchAddress("beamspotX", &PV_x);
    t->SetBranchAddress("beamspotY", &PV_y);
    t->SetBranchAddress("MbdVertexZ", &PV_z);
    for (int ev = 0; ev < t->GetEntriesFast(); ev++)
    {
        t->GetEntry(ev);
        vtx.beamspot_x = PV_x;
        vtx.beamspot_y = PV_y;
        vtx.mbd_z = PV_z; //! for time being
        vtx.silicon_z_crossing0 = PV_z; //! for time being

        EvtVtx_map[event] = vtx;
    }

    return EvtVtx_map;
}

std::map<uint64_t, vertex> EvtVtx_map_inttbco(const char *vtxfname)
{
    std::map<uint64_t, vertex> EvtVtx_map;

    vertex vtx;

    TFile *f = new TFile(vtxfname, "READ");
    TTree *t = (TTree *)f->Get("vertex");
    uint64_t gl1bco, bcotr;
    bool trigger12;
    float beamspotX, beamspotY, MbdVertexZ, SiliconVertexZ_crossing0;
    t->SetBranchAddress("gl1bco", &gl1bco);
    t->SetBranchAddress("bcotr", &bcotr);
    t->SetBranchAddress("trigger12", &trigger12);
    t->SetBranchAddress("beamspotX", &beamspotX);
    t->SetBranchAddress("beamspotY", &beamspotY);
    t->SetBranchAddress("MbdVertexZ", &MbdVertexZ);
    t->SetBranchAddress("SiliconVertexZ_crossing0", &SiliconVertexZ_crossing0);
    for (int ev = 0; ev < t->GetEntriesFast(); ev++)
    {
        t->GetEntry(ev);
        
        vtx.beamspot_x = beamspotX;
        vtx.beamspot_y = beamspotY;
        vtx.mbd_z = MbdVertexZ;
        vtx.silicon_z_crossing0 = SiliconVertexZ_crossing0;
        vtx.trigger12 = trigger12;

        EvtVtx_map[gl1bco] = vtx;
        vtx.reset();
    }

    return EvtVtx_map;
}

std::vector<float> getBeamspot(const char *fname, const char *tname)
{
    TString mergename = Form("%s", fname);
    // TString cmd = Form("hadd -f -j 20 %s %s/minitree_00*.root", mergename.Data(), fpath);
    // system(cmd.Data());

    float beamspotx, beamspoty;

    TFile *f = new TFile(mergename, "READ");
    TTree *t = (TTree *)f->Get(tname);
    float bs_x, bs_y;
    t->SetBranchAddress("beamspotX", &bs_x);
    t->SetBranchAddress("beamspotY", &bs_y);
    for (int i = 0; i < 1; i++) // the same for all entries
    {
        t->GetEntry(i);
        beamspotx = bs_x;
        beamspoty = bs_y;
    }
    f->Close();

    return {beamspotx, beamspoty};
}

TH1F *VtxZ_ReweiHist(const char *filename = "/sphenix/user/hjheng/sPHENIXRepo/analysis/dNdEta_Run2023/analysis_INTT/plot/RecoPV_ana/VtxZ_reweight_HIJING_ana419_20240910.root", const char *histname = "VtxZ_reweight_VtxZm10to10")
{
    TFile *f = new TFile(filename, "READ");
    TH1F *h = (TH1F *)f->Get(histname);
    h->SetDirectory(0);
    f->Close();
    
    return h;
}

#endif