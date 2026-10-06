#include "analysisHelper.h"

// ─────────────────────────────────────────────────────────────────────────────
//  fillData.C
//
//  Fills measured (reco-only) histograms from real data.
//  No truth information, no response matrix, no pT-hat cut.
//
//  Output — dijet pT:
//    hJetPt_meas       — measured reco dijet pT (flat index, unfolding input)
//
//  Output — wEEC 2D per reco dijet pT bin fr:
//    hWEEC2D_data_fr   — (ΔΦ, reco pair weight) histogram
//                        pair weight = pT_i * pT_j / <pT>²
//                        using vertex-shifted eta for reco towers
// ─────────────────────────────────────────────────────────────────────────────
void fillData(int runnumber = 47722, int seg = 0, const char* outDir = ".", float towCut = 0.25)
{
    TH1D* hJetPt_meas = new TH1D("hJetPt_meas", "", nRecoFlat(), 0, nRecoFlat());
    hJetPt_meas->Sumw2();

    // Vtx distribution from data — used by makeVtxWeights.C to form the
    // data/MC ratio that is stored in analysisHelper.h.
    TH1D* hVtxData = new TH1D("hVtxData",
        "Data vtx_z;v_{tz} (cm);events",
        nVtxBins, vtxZMin, vtxZMax);
    // No Sumw2: data events are unweighted counts.

    // 3D data histograms: one per ΔΦ bin, in reco-index space
    std::vector<TH1D*> hWEEC3D_data(nDphi, nullptr);
    for (int k = 0; k < nDphi; ++k) {
        hWEEC3D_data[k] = new TH1D(
            std::format("hWEEC3D_data_{}", k).c_str(),
            std::format("Data wEEC 3D k{};flat(iL,iS,iPw);pairs", k).c_str(),
            nRecoFlat3D(), 0, nRecoFlat3D());
        hWEEC3D_data[k]->Sumw2();
    }

    // ── open input ────────────────────────────────────────────────────────
    std::string dataFile = std::format(
        "/sphenix/tg/tg01/jets/bkimelman/VandyDSTs_Apr4_2026/"
        "VandyDSTs_run2pp_ana521_2025p007_v001-{:08}-{}_{}.root",
        runnumber, seg, seg + 24);
    TFile* fIn = new TFile(dataFile.c_str(), "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "Cannot open " << dataFile << "\n";
        if (fIn) delete fIn; return;
    }
    TTree* T = (TTree*)fIn->Get("T");
    if (!T) { fIn->Close(); delete fIn; return; }

    EventInfo*          eventInfo = nullptr;
    std::vector<Tower>* towers    = nullptr;
    T->SetBranchAddress("EventInfo", &eventInfo);
    T->SetBranchAddress("TowerInfo", &towers);

    long nEvents = T->GetEntries();
    for (long ev = 0; ev < nEvents; ++ev)
    {
        std::cout << "working on event " << ev << std::endl;
        T->GetEntry(ev);

        double vtx_z = eventInfo->get_z_vtx();
        if (std::abs(vtx_z) > 10.0) continue;

        // Fill vtx distribution for all events passing the vtx cut,
        // regardless of whether they pass the dijet selection.
        hVtxData->Fill(vtx_z);

        if (!eventInfo->is_dijet_event(2)) continue;
        double rLeadPT = eventInfo->get_lead_pT(2);
        double rSublPT = eventInfo->get_sublead_pT(2);
        if (rLeadPT < recoLeadMin || rSublPT < recoSublMin) continue;
        

        int iLreco = FindBin(rLeadPT, recoLeadPtBins);
        int iSreco = FindBin(rSublPT, recoSublPtBins);
        if (iLreco < 0 || iSreco < 0) continue;

        int    fr      = RecoFlatIndex(iLreco, iSreco);
        double rFlat   = RecoFlatBinCenter(fr);
        hJetPt_meas->Fill(rFlat);   // no MC weight for data

        double rPTmean2 = 0.25 * (rLeadPT + rSublPT) * (rLeadPT + rSublPT);

        double recoMap[24][64] = {};
        shiftEtaEdges(vtx_z);

        for (auto& tow : *towers) {
            if (tow.get_calo() < 1 || tow.get_calo() > 3) continue;
            fastjet::PseudoJet tj(tow.px(), tow.py(), tow.pz(), tow.e());
            int etaBin = getEtaBin(tj.pseudorapidity(), vtx_z);
            int phiBin = getPhiBin(tj.phi());
            if (etaBin >= 0 && phiBin >= 0)
                recoMap[etaBin][phiBin] += tow.e();
        }

        for (int i = 0; i < 24*64; ++i) {
            int ei = i / 64, pi = i % 64;
            double phi_i = getPhiCenter(pi);
            double rpT_i = recoMap[ei][pi] / cosh(getEtaCenter(ei, vtx_z));
            if (rpT_i <= towCut || rpT_i >= 80) continue;

            for (int j = i+1; j < 24*64; ++j) {
                int ej = j / 64, pj = j % 64;
                double phi_j = getPhiCenter(pj);
                double rpT_j = recoMap[ej][pj] / cosh(getEtaCenter(ej, vtx_z));
                if (rpT_j <= towCut || rpT_j >= 80) continue;

                double dphi   = DeltaPhi(phi_i, phi_j);
                double rPairW = rpT_i * rpT_j / rPTmean2;
                if (rPairW < pairWeightMin || rPairW > pairWeightMax) continue;
                int iDphi = FindBin(dphi, dPhiBins);
                int iPw   = FindBin(rPairW, pairWeightBins);
                if (iDphi < 0 || iPw < 0) continue;
                hWEEC3D_data[iDphi]->Fill(
                    RecoFlat3DBinCenter(RecoFlat3DIndex(iLreco, iSreco, iPw)));
            }
        }
    }

    fIn->Close(); delete fIn;

    std::string outName = std::format(
        "{}/Data/data_measured_{:08}_seg{:06d}_to_{:06d}.root",
        outDir, runnumber, seg, seg + 24);

    std::cout << "writing file " << outName << std::endl;

    TFile* fOut = new TFile(outName.c_str(), "RECREATE");
    hJetPt_meas->Write();
    for (int k = 0; k < nDphi; ++k)
        hWEEC3D_data[k]->Write();
    hVtxData->Write();
    fOut->Close();
    std::cout << "Written: " << outName << "\n";
}
