#include "analysisHelper.h"
#include "RooUnfoldResponse.h"

// ─────────────────────────────────────────────────────────────────────────────
//  fillResponse.C
//
//  Fills the dijet pT response matrix and all wEEC histograms needed for
//  both the 2D dijet pT unfolding and the new 3D wEEC unfolding.
//
//  Truth towers are constructed from TruthParticles using the UNSHIFTED tower
//  grid (same eta/phi indices as reco) so that matched pairs share indices.
//  Individual particles must have E > 0.5 GeV to contribute to a truth tower.
//  Tower thresholds (0.25–80 GeV) are applied after summing.
//
//  Mode::kFull  — all events fill both response and measured histograms
//  Mode::kHalf  — random ~50% split: training half fills response, measured
//                 half fills hJetPt_meas and wEEC measured histograms
//
//  ── 2D dijet pT unfolding outputs ──────────────────────────────────────────
//  response_jet_pT          — RooUnfoldResponse for dijet pT unfolding
//  hJetPt_truth_matched     — matched truth per-event
//  hJetPt_truth             — all truth (matched + missed) per-event
//  hJetPt_reco_matched      — matched reco per-event
//  hJetPt_reco              — all reco (matched + fake) per-event
//  hJetPt_meas              — measured reco per-event
//  hJetPt_truth_meas        — truth from measured half (kHalf diagnostic)
//
//  ── 3D wEEC unfolding outputs, per ΔΦ bin k ────────────────────────────────
//  response_wEEC3D_{k}      — RooUnfoldResponse in 3D flat space
//                             axes: (lead pT bin, subl pT bin, pair weight bin)
//                             filled per tower pair in matched dijet events;
//                             dijet-level fakes/misses also fill this response
//  hWEEC3D_truth_{k}        — truth-tower pairs for all truth-valid events
//                             (closure target: what unfolding should recover)
//  hWEEC3D_meas_{k}         — reco-tower pairs for measured events
//  hWEEC3D_data_{k}         — reco-tower pairs from real data (fillData.C)
//
//  ── retained 2D wEEC histograms (for old unfolding path) ───────────────────
//  hWEEC_fake_{fr}          — flat 2D, fake dijet reco pairs
//  hWEEC_miss_{ft}          — flat 2D, missed dijet truth pairs
//  hWEEC_truth_{ft}         — flat 2D, truth pairs per truth dijet pT bin
//  hWEEC_meas_{fr}          — flat 2D, measured reco pairs per reco dijet pT bin
//  response_wEEC_{ft}       — old per-ft flat-2D response (kept for reference)
// ─────────────────────────────────────────────────────────────────────────────
void fillResponse(int jetSample = 20, int seg = 0, const char* outDir = ".",
                  Mode mode = Mode::kFull, float towCut = 0.25)
{
    const bool halfClosure = (mode == Mode::kHalf);
    const bool vtxOnly     = (mode == Mode::kVtx);

    std::string outName = std::format(
        "{}/Jet{}/response_Jet{}_seg{:06d}_to_{:06d}-{}.root",
        outDir, jetSample, jetSample, seg, seg + 24, ModeLabel(mode));
        
    std::cout << "making file " << outName << std::endl;
    std::cout.flush();

    TFile* fOut = new TFile(outName.c_str(), "RECREATE");

    std::cout << "made file " << outName << std::endl;
    std::cout.flush();

    // ── vtx histogram (all modes) ─────────────────────────────────────────
    // hVtxMC is filled in every mode so the vtx distribution can always be
    // inspected.  In kVtx mode this is the only histogram filled — all
    // tower-pair and response-matrix work is skipped entirely.
    TH1D* hVtxMC = new TH1D("hVtxMC",
        "MC vtx_z (xsec weighted);v_{z} (cm);weighted events",
        nVtxBins, vtxZMin, vtxZMax);
    hVtxMC->Sumw2();

    if (vtxOnly)
    {
        // ── kVtx fast path ────────────────────────────────────────────────
        // Open the MC file and loop over events, filling hVtxMC with the
        // cross-section weight.  No tower maps, no jet matching, no pair loops.
        TFile* simFile = new TFile(SimFilePath(jetSample, seg).c_str(), "READ");
        if (!simFile || simFile->IsZombie()) {
            std::cerr << "Cannot open sim file\n";
            fOut->Close(); return;
        }

        TTree* T = (TTree*)simFile->Get("T");
        if (!T) { simFile->Close(); fOut->Close(); return; }

        EventInfo*                  eventInfo  = nullptr;
        std::vector<JetInfo>* truthJets = nullptr;
        T->SetBranchAddress("EventInfo", &eventInfo);
        T->SetBranchAddress("TruthJetInfo_r04", &truthJets);

        long nEvents = T->GetEntries();
        for (long ev = 0; ev < nEvents; ++ev) {
            std::cout << "kVtx: working on event " << ev << std::endl;
            T->GetEntry(ev);

            double vtx_z = eventInfo->get_z_vtx();
            if (std::abs(vtx_z) > 10.0) continue;

            // pThat slice cut — same guard as the full loop
            double maxTruthPT = -1;
            for (auto& j : *truthJets)
                if (j.pt() > maxTruthPT) maxTruthPT = j.pt();
            if (maxTruthPT < get_jet_pTLow(jetSample) ||
                maxTruthPT > get_jet_pTHigh(jetSample)) continue;

            double evtWeight = eventInfo->get_cross_section();
            hVtxMC->Fill(vtx_z, evtWeight);
        }

        simFile->Close(); delete simFile;

        fOut->cd();
        hVtxMC->Write();
        fOut->Close();
        std::cout << "Written (vtxOnly): " << outName << "\n";
        return;
    }

    // ── dijet pT histograms ───────────────────────────────────────────────
    TH1D* hJetPt_truth_matched = new TH1D("hJetPt_truth_matched", "", nTrueFlat(), 0, nTrueFlat());
    TH1D* hJetPt_truth         = new TH1D("hJetPt_truth",         "", nTrueFlat(), 0, nTrueFlat());
    TH1D* hJetPt_reco_matched  = new TH1D("hJetPt_reco_matched",  "", nRecoFlat(), 0, nRecoFlat());
    TH1D* hJetPt_reco          = new TH1D("hJetPt_reco",          "", nRecoFlat(), 0, nRecoFlat());
    TH1D* hJetPt_meas          = new TH1D("hJetPt_meas",          "", nRecoFlat(), 0, nRecoFlat());
    TH1D* hJetPt_truth_meas    = new TH1D("hJetPt_truth_meas",    "", nTrueFlat(), 0, nTrueFlat());
    hJetPt_truth_matched->Sumw2(); hJetPt_truth       ->Sumw2();
    hJetPt_reco_matched ->Sumw2(); hJetPt_reco        ->Sumw2();
    hJetPt_meas         ->Sumw2(); hJetPt_truth_meas  ->Sumw2();

    std::cout << "made jet hists" << std::endl;
    std::cout.flush();

    RooUnfoldResponse* respJetPt = nullptr;
    {
        //TH1D hJPReco ("hJetPt_reco_tmpl",  "", nRecoFlat, 0, nRecoFlat);
        //TH1D hJPTruth("hJetPt_truth_tmpl", "", nTrueFlat, 0, nTrueFlat);
        respJetPt = new RooUnfoldResponse(nRecoFlat(), 0, nRecoFlat(), nTrueFlat(), 0, nTrueFlat());
    }

    std::cout << "made jet response matrix" << std::endl;
    std::cout.flush();

    // ── 3D wEEC response matrices: one per ΔΦ bin k ───────────────────────
    // Each maps flat 3D (lead pT bin, subl pT bin, pair weight bin) on the
    // reco side to the same space on the truth side.
    //   iFlat3D = iL * (nSubl*nPairWeight) + iS * nPairWeight + iPw
    // Dijet-level fakes/misses are included directly in these response
    // matrices (all pairs from fake/miss dijet events call Fake/Miss).
    // Tower-level fakes/misses within matched events are also included.
    std::vector<RooUnfoldResponse*> respWEEC3D(nDphi, nullptr);
    for (int k = 0; k < nDphi; ++k)
        {
            //TH1D hReco ("hWEEC_reco_tmpl",  "", nRecoFlat3D, 0, nRecoFlat3D);
            //TH1D hTruth("hWEEC_truth_tmpl", "", nTrueFlat3D, 0, nTrueFlat3D);
        respWEEC3D[k] = new RooUnfoldResponse(nRecoFlat3D(), 0, nRecoFlat3D(), nTrueFlat3D(), 0, nTrueFlat3D());
        }

    std::cout << "made 3D response matrices" << std::endl;
    std::cout.flush();

    // ── 3D wEEC truth histograms (per ΔΦ bin k, training half) ───────────
    // Filled for ALL truth-valid events (matched + missed) per truth-tower pair.
    // Closure target: after unfolding hWEEC3D_meas, the projected result
    // should match the projection of this histogram.
    std::vector<TH1D*> hWEEC3D_truth(nDphi, nullptr);
    std::vector<TH1D*> hWEEC3D_misses(nDphi, nullptr);
    for (int k = 0; k < nDphi; ++k) {
        hWEEC3D_truth[k] = new TH1D(
            std::format("hWEEC3D_truth_{}", k).c_str(),
            std::format("Truth wEEC 3D k{};flat(iL,iS,iPw);pairs", k).c_str(),
            nTrueFlat3D(), 0, nTrueFlat3D());
        hWEEC3D_truth[k]->Sumw2();

        hWEEC3D_misses[k] = new TH1D(
            std::format("hWEEC3D_misses_{}", k).c_str(),
            std::format("Misses wEEC 3D k{};flat(iL,iS,iPw);pairs", k).c_str(),
            nTrueFlat3D(), 0, nTrueFlat3D());
        hWEEC3D_misses[k]->Sumw2();
    }

    std::cout << "made 3D truth hists" << std::endl;
    std::cout.flush();

    // ── Pair weight mean numerator/denominator (per ΔΦ bin k) ────────────
    // hWEEC3D_pwNum_{k}: sum of (evtWeight * pairWeight) per truth flat3D bin
    // hWEEC3D_pwDen_{k}: sum of evtWeight per truth flat3D bin
    // Dividing num/den at projection time gives the luminosity-weighted mean
    // pair weight within each bin, which is more accurate than the bin-center
    // approximation when the pair weight distribution is not flat within a bin.
    // Both histograms are filled at every site that fills hWEEC3D_truth.
    std::vector<TH1D*> hWEEC3D_pwNum(nDphi, nullptr);
    std::vector<TH1D*> hWEEC3D_pwDen(nDphi, nullptr);
    for (int k = 0; k < nDphi; ++k) {
        hWEEC3D_pwNum[k] = new TH1D(
            std::format("hWEEC3D_pwNum_{}", k).c_str(),
            std::format("PW num k{};flat(iL,iS,iPw);sum(w*evtW)", k).c_str(),
            nTrueFlat3D(), 0, nTrueFlat3D());
        hWEEC3D_pwDen[k] = new TH1D(
            std::format("hWEEC3D_pwDen_{}", k).c_str(),
            std::format("PW den k{};flat(iL,iS,iPw);sum(evtW)", k).c_str(),
            nTrueFlat3D(), 0, nTrueFlat3D());
        hWEEC3D_pwNum[k]->Sumw2();
        hWEEC3D_pwDen[k]->Sumw2();
    }

    // ── 3D wEEC measured histograms (per ΔΦ bin k, measured half) ─────────
    std::vector<TH1D*> hWEEC3D_meas(nDphi, nullptr);
    std::vector<TH1D*> hWEEC3D_fakes(nDphi, nullptr);
    std::vector<TH1D*> hWEEC3D_meas_truth(nDphi, nullptr);
    for (int k = 0; k < nDphi; ++k) {
        hWEEC3D_meas[k] = new TH1D(
            std::format("hWEEC3D_meas_{}", k).c_str(),
            std::format("Meas wEEC 3D k{};flat(iL,iS,iPw);pairs", k).c_str(),
            nRecoFlat3D(), 0, nRecoFlat3D());
        hWEEC3D_meas[k]->Sumw2();

        hWEEC3D_fakes[k] = new TH1D(
            std::format("hWEEC3D_fakes_{}", k).c_str(),
            std::format("Fakes wEEC 3D k{};flat(iL,iS,iPw);pairs", k).c_str(),
            nRecoFlat3D(), 0, nRecoFlat3D());
        hWEEC3D_fakes[k]->Sumw2();

        hWEEC3D_meas_truth[k] = new TH1D(
            std::format("hWEEC3D_meas_truth_{}", k).c_str(),
            std::format("Truth Meas Half wEEC 3D k{};flat(iL,iS,iPw);pairs", k).c_str(),
            nTrueFlat3D(), 0, nTrueFlat3D());
        hWEEC3D_meas_truth[k]->Sumw2();
    }

    std::cout << "made 3D meas hists" << std::endl;
    std::cout.flush();

    // ── 3D wEEC counts histograms (per ΔΦ bin k, response half) ──────────
    // 2D histogram matching the response matrix exactly: reco flat3D on X,
    // truth flat3D on Y.  Filled with unity weight (not the cross-section
    // weight) for every tower pair that enters a response matrix Fill call.
    // Fake and Miss calls are not counted here — they have no (reco, truth)
    // cell to trim.  Used by wEEC_doUnfolding.C for cell-by-cell trimming.
    std::vector<TH2D*> hWEEC3D_counts(nDphi, nullptr);
    for (int k = 0; k < nDphi; ++k) {
        hWEEC3D_counts[k] = new TH2D(
            std::format("hWEEC3D_counts_{}", k).c_str(),
            std::format("Counts wEEC 3D k{};reco flat3D;truth flat3D", k).c_str(),
            nRecoFlat3D(), 0, nRecoFlat3D(),
            nTrueFlat3D(), 0, nTrueFlat3D());
        // No Sumw2: raw event counts, not weighted sums.
    }

    std::cout << "made 3D counts hists" << std::endl;
    std::cout.flush();

    // ── Direct wEEC cross-check histograms ───────────────────────────────
    // Filled as direct ΔΦ histograms for cross-checks only — not unfolding
    // inputs.  Both filled for every trueInRange event (matched + missed),
    // mirroring hWEEC3D_truth coverage.
    //
    // hWEEC_particle:   particle-level wEEC — pairs over final-state particles
    //                   (calo==4) before tower binning.  No tower thresholds.
    //                   The physical observable the measurement corrects to.
    //
    // hWEEC_truthTower: truth-tower wEEC — pairs over truthMap with the same
    //                   0.25–80 GeV thresholds as the unfolding truth.
    //                   Quantifies tower discretization effects vs particles.
    TH1D* hWEEC_particle   = new TH1D("hWEEC_particle",
        "Particle-level wEEC;#Delta#phi;wEEC", nDphi, dPhiBins.data());
    TH1D* hWEEC_truthTower = new TH1D("hWEEC_truthTower",
        "Truth-tower wEEC;#Delta#phi;wEEC", nDphi, dPhiBins.data());
    hWEEC_particle  ->Sumw2();
    hWEEC_truthTower->Sumw2();

    TFile* simFile = new TFile(SimFilePath(jetSample, seg).c_str(), "READ");
    if (!simFile || simFile->IsZombie()) {
        std::cerr << "Cannot open " << SimFilePath(jetSample, seg) << "\n";
        if (simFile) delete simFile; return;
    }
    TTree* T = (TTree*)simFile->Get("T");
    if (!T) { simFile->Close(); delete simFile; return; }

    TRandom3 rng(jetSample * 100000 + seg);

    // ── branches ──────────────────────────────────────────────────────────
    EventInfo*            eventInfo = nullptr;
    std::vector<JetInfo>* recoJets  = nullptr;
    std::vector<JetInfo>* truthJets = nullptr;
    std::vector<Tower>*   towers    = nullptr;
    std::vector<Tower>*   particles = nullptr;

    T->SetBranchAddress("EventInfo",        &eventInfo);
    T->SetBranchAddress("JetInfo_r04",      &recoJets);
    T->SetBranchAddress("TruthJetInfo_r04", &truthJets);
    T->SetBranchAddress("TowerInfo",        &towers);
    T->SetBranchAddress("TruthParticles",   &particles);

    // ── event loop ────────────────────────────────────────────────────────
    long nEvents = T->GetEntries();
    long nGood = 0, nGoodReco = 0, nGoodTruth = 0, nMatched = 0;
    for (long ev = 0; ev < nEvents; ++ev)
    {
        T->GetEntry(ev);

        const bool isEven   = (rng.Rndm() < 0.5);
        const bool fillResp = !halfClosure || isEven;
        const bool fillMeas = !halfClosure || !isEven;

        double vtx_z = eventInfo->get_z_vtx();
     
        std::cout << std::format("Working on event {:5d}: vtx_z={:.2f}\n", ev, vtx_z);
        std::cout.flush();
        
        if (std::abs(vtx_z) > 10.0) continue;

        double evtWeight = eventInfo->get_cross_section();

        // Apply vtx_z reweighting only when building the response for real
        // data (kData mode).  In closure modes (kFull, kHalf) the MC vtx
        // distribution is self-consistent and reweighting would introduce
        // a spurious bias.
        if (mode == Mode::kData)
            evtWeight *= GetVtxWeight(vtx_z);

        double maxTruthPT = -1;
        for (auto& j : *truthJets)
            if (j.pt() > maxTruthPT) maxTruthPT = j.pt();
        if (maxTruthPT < get_jet_pTLow(jetSample) ||
            maxTruthPT > get_jet_pTHigh(jetSample)) continue;

        bool rValid = eventInfo->is_dijet_event(2)
                   && eventInfo->get_lead_pT(2)    >= recoLeadMin
                   && eventInfo->get_sublead_pT(2) >= recoSublMin;
        bool tValid = eventInfo->is_dijetTruth_event(2)
                   && eventInfo->get_leadTruth_pT(2)    >= trueLeadMin
                   && eventInfo->get_subleadTruth_pT(2) >= trueSublMin;

        if(!rValid && !tValid) continue;

        // Fill vtx histogram with the raw (pre-vtx-reweight) cross-section
        // so the stored MC distribution always reflects the true MC shape,
        // usable as the denominator in makeVtxWeights.C regardless of mode.
        hVtxMC->Fill(vtx_z, eventInfo->get_cross_section());

        double rLeadPT = rValid ? eventInfo->get_lead_pT(2)         : -1;
        double rSublPT = rValid ? eventInfo->get_sublead_pT(2)      : -1;
        double tLeadPT = tValid ? eventInfo->get_leadTruth_pT(2)    : -1;
        double tSublPT = tValid ? eventInfo->get_subleadTruth_pT(2) : -1;

        int iLreco = rValid ? FindBin(rLeadPT, recoLeadPtBins) : -1;
        int iSreco = rValid ? FindBin(rSublPT, recoSublPtBins) : -1;
        int iLtrue = tValid ? FindBin(tLeadPT, trueLeadPtBins) : -1;
        int iStrue = tValid ? FindBin(tSublPT, trueSublPtBins) : -1;

        bool recoInRange = rValid && iLreco >= 0 && iSreco >= 0;
        bool trueInRange = tValid && iLtrue >= 0 && iStrue >= 0;

        bool geoMatch = false;
        if (rValid && tValid)
            geoMatch = GetGeoMatch(*recoJets, *truthJets,
                                   rLeadPT, rSublPT, tLeadPT, tSublPT);

        // ── build reco tower map (vertex-shifted eta) ─────────────────────
        double recoMap[24][64] = {};
        if (recoInRange) {
            shiftEtaEdges(vtx_z);
            for (auto& tow : *towers) {
                if (tow.get_calo() < 1 || tow.get_calo() > 3) continue;
                fastjet::PseudoJet tj(tow.px(), tow.py(), tow.pz(), tow.e());
                int etaBin = getEtaBin(tj.pseudorapidity(), vtx_z);
                int phiBin = getPhiBin(tj.phi());
                if (etaBin >= 0 && phiBin >= 0)
                    recoMap[etaBin][phiBin] += tow.e();
            }
        }

        // ── build truth tower map (unshifted eta — shares indices with reco) ──
        // Particles without calo 4 are not used since they are not truth particles
        // Particles with E > 0.5 GeV are assigned to the unshifted grid.
        // Tower thresholds applied after summing, same as reco.
        double truthMap[24][64] = {};
        if (trueInRange) {
            for (auto& part : *particles) {
                if(part.get_calo() != 4) continue;
                //if (part.e() <= 0.5) continue;
                fastjet::PseudoJet pj(part.px(), part.py(), part.pz(), part.e());
                int etaBin = getEtaBin(pj.pseudorapidity(), 0.0);  // unshifted
                int phiBin = getPhiBin(pj.phi());
                if (etaBin >= 0 && phiBin >= 0)
                    truthMap[etaBin][phiBin] += part.e();
            }
        }

        // ── fill dijet pT histograms ──────────────────────────────────────
        if (fillResp) {
            if (recoInRange && trueInRange && geoMatch) {
                double rFlat = RecoFlatBinCenter(RecoFlatIndex(iLreco, iSreco));
                double tFlat = TrueFlatBinCenter(TrueFlatIndex(iLtrue, iStrue));
                respJetPt->Fill(rFlat, tFlat, evtWeight);
                hJetPt_reco_matched ->Fill(rFlat, evtWeight);
                hJetPt_reco         ->Fill(rFlat, evtWeight);
                hJetPt_truth_matched->Fill(tFlat, evtWeight);
                hJetPt_truth        ->Fill(tFlat, evtWeight);
            } else {
                if (recoInRange) {
                    double rFlat = RecoFlatBinCenter(RecoFlatIndex(iLreco, iSreco));
                    respJetPt->Fake(rFlat, evtWeight);
                    hJetPt_reco->Fill(rFlat, evtWeight);
                }
                if (trueInRange) {
                    double tFlat = TrueFlatBinCenter(TrueFlatIndex(iLtrue, iStrue));
                    respJetPt->Miss(tFlat, evtWeight);
                    hJetPt_truth->Fill(tFlat, evtWeight);
                }
            }
        }
        if (fillMeas && recoInRange) {
            hJetPt_meas->Fill(RecoFlatBinCenter(RecoFlatIndex(iLreco, iSreco)), evtWeight);
            if (halfClosure && trueInRange)
                hJetPt_truth_meas->Fill(TrueFlatBinCenter(TrueFlatIndex(iLtrue, iStrue)), evtWeight);
        }

        // ══════════════════════════════════════════════════════════════════════
        //  3D wEEC tower-pair loop
        //
        //  The event is classified into exactly one of: matched, fake, missed.
        //  Each class has one tower-pair loop that computes all quantities once
        //  and fills every relevant histogram in a single pass.
        // ══════════════════════════════════════════════════════════════════════

        const bool isMatched = recoInRange && trueInRange && geoMatch;
        const bool isFake    = recoInRange && !isMatched;
        const bool isMiss    = trueInRange && !isMatched;

        ++nGood;
        if (recoInRange) ++nGoodReco;
        if (trueInRange) ++nGoodTruth;
        if (isMatched)   ++nMatched;

        std::cout << std::format(
            "ev {:5d}: reco={} truth={} match={} | "
            "good={} goodReco={} goodTruth={} matched={}\n",
            ev, (int)recoInRange, (int)trueInRange, (int)isMatched,
            nGood, nGoodReco, nGoodTruth, nMatched);
        std::cout.flush();

        // ── direct wEEC cross-checks (all trueInRange events) ────────────
        // Both histograms are filled here once per event, independently of
        // the matched/fake/miss classification, so their coverage matches
        // hWEEC3D_truth exactly (matched + missed).
        if (fillResp && trueInRange)
        {
            const double tPTmean2 = 0.25*(tLeadPT+tSublPT)*(tLeadPT+tSublPT);

            // ── particle-level wEEC ───────────────────────────────────────
            // Iterate over all calo==4 particles directly — no tower thresholds.
            // Build a lightweight list of (pT, phi) for this event first to
            // avoid recomputing pseudorapidity inside the double loop.
            struct PInfo { double pT, phi; };
            std::vector<PInfo> plist;
            plist.reserve(particles->size());
            for (auto& part : *particles) {
                if (part.get_calo() != 4) continue;
                fastjet::PseudoJet pj(part.px(), part.py(), part.pz(), part.e());
                double eta = pj.pseudorapidity();
                // Restrict to calorimeter acceptance — same η range as towers
                if (std::abs(eta) > 1.1) continue;
                double phi = pj.phi();
                double pT = part.e() / std::cosh(eta);
                if (pT <= 0) continue;
                plist.push_back({pT, phi});
            }
            for (int i = 0; i < (int)plist.size(); ++i)
            for (int j = i+1; j < (int)plist.size(); ++j) {
                double dphi   = DeltaPhi(plist[i].phi, plist[j].phi);
                double pairW  = plist[i].pT * plist[j].pT / tPTmean2;
                if (pairW < pairWeightMin || pairW > pairWeightMax) continue;
                hWEEC_particle->Fill(dphi, evtWeight * pairW);
            }

            // ── truth-tower wEEC ──────────────────────────────────────────
            // Same tower thresholds and eta/phi lookup as the unfolding truth.
            for (int i = 0; i < 24*64; ++i) {
                const int ei = i/64, pi = i%64;
                double phi_i = getPhiCenter(pi);
                const double tpT_i = truthMap[ei][pi] / cosh(getEtaCenter(ei, 0.0));
                if (tpT_i <= towCut || tpT_i >= 80) continue;

                for (int j = i+1; j < 24*64; ++j) {
                    const int ej = j/64, pj = j%64;
                    double phi_j = getPhiCenter(pj);
                    const double tpT_j = truthMap[ej][pj] / cosh(getEtaCenter(ej, 0.0));
                    if (tpT_j <= towCut || tpT_j >= 80) continue;

                    double dphi  = DeltaPhi(phi_i, phi_j);
                    double pairW = tpT_i * tpT_j / tPTmean2;
                    if (pairW < pairWeightMin || pairW > pairWeightMax) continue;
                    hWEEC_truthTower->Fill(dphi, evtWeight * pairW);
                }
            }
        }

        // ── matched dijet ─────────────────────────────────────────────────
        // Handles all four tower-pair cases in one loop:
        //   hasReco && hasTruth  → response Fill + truth/meas hists
        //   hasReco only         → Fake + meas hist
        //   hasTruth only        → Miss + truth hist
        //   neither              → skipped
        // fillMeas fills (hWEEC3D_meas) are included here so no separate
        // reco-pair loop is needed for matched events.
        if (isMatched)
        {
            const double rPTmean2 = 0.25*(rLeadPT+rSublPT)*(rLeadPT+rSublPT);
            const double tPTmean2 = 0.25*(tLeadPT+tSublPT)*(tLeadPT+tSublPT);

            for (int i = 0; i < 24*64; ++i) {
                const int ei = i/64, pi = i%64;
                const double rE_i = recoMap[ei][pi];
                const double tE_i = truthMap[ei][pi];

                double phi_i = getPhiCenter(pi);
                const double rpT_i = rE_i/cosh(getEtaCenter(ei, vtx_z));
                const double tpT_i = tE_i/cosh(getEtaCenter(ei, 0.0));
                
                const bool rOk_i = (rpT_i > towCut && rpT_i < 80);
                const bool tOk_i = (tpT_i > towCut && tpT_i < 80);
                if (!rOk_i && !tOk_i) continue;

                for (int j = i+1; j < 24*64; ++j) {
                    const int ej = j/64, pj = j%64;
                    const double rE_j = recoMap[ej][pj];
                    const double tE_j = truthMap[ej][pj];

                    double phi_j = getPhiCenter(pj);
                    const double rpT_j = rE_j/cosh(getEtaCenter(ej, vtx_z));
                    const double tpT_j = tE_j/cosh(getEtaCenter(ej, 0.0));

                    const bool rOk_j = (rpT_j > towCut && rpT_j < 80);
                    const bool tOk_j = (tpT_j > towCut && tpT_j < 80);
                    if (!rOk_j && !tOk_j) continue;

                    const double dphi  = DeltaPhi(phi_i, phi_j);
                    const int    iDphi = FindBin(dphi, dPhiBins);
                    if (iDphi < 0) continue;

                    const bool hasReco  = rOk_i && rOk_j;
                    const bool hasTruth = tOk_i && tOk_j;

                    if (hasReco && hasTruth) {
                        const double rPairW = rpT_i*rpT_j/rPTmean2;
                        const double tPairW = tpT_i*tpT_j/tPTmean2;
                        if (rPairW >= pairWeightMin && rPairW <= pairWeightMax &&
                            tPairW >= pairWeightMin && tPairW <= pairWeightMax) {
                            const int rIPw = FindBin(rPairW, pairWeightBins);
                            const int tIPw = FindBin(tPairW, pairWeightBins);
                            if (rIPw >= 0 && tIPw >= 0) {
                                const double rF3Dcen = RecoFlat3DBinCenter(RecoFlat3DIndex(iLreco,iSreco,rIPw));
                                const double tF3Dcen = TrueFlat3DBinCenter(TrueFlat3DIndex(iLtrue,iStrue,tIPw));
                                if (fillResp) {
                                    respWEEC3D[iDphi]->Fill(rF3Dcen, tF3Dcen, evtWeight);
                                    hWEEC3D_counts[iDphi]->Fill(rF3Dcen, tF3Dcen);
                                    hWEEC3D_truth[iDphi]->Fill(tF3Dcen, evtWeight);
                                    hWEEC3D_pwNum[iDphi]->Fill(tF3Dcen, evtWeight * tPairW);
                                    hWEEC3D_pwDen[iDphi]->Fill(tF3Dcen, evtWeight);
                                }
                                if (fillMeas) {
                                    hWEEC3D_meas[iDphi]->Fill(rF3Dcen, evtWeight);
                                    hWEEC3D_meas_truth[iDphi]->Fill(tF3Dcen, evtWeight);
                                }
                            }
                        }
                    } else if (hasReco) {
                        const double rPairW = rpT_i*rpT_j/rPTmean2;
                        if (rPairW >= pairWeightMin && rPairW <= pairWeightMax) {
                            const int rIPw = FindBin(rPairW, pairWeightBins);
                            if (rIPw >= 0) {
                                const double rF3Dcen = RecoFlat3DBinCenter(RecoFlat3DIndex(iLreco,iSreco,rIPw));
                                if (fillResp) {
                                    respWEEC3D[iDphi]->Fake(rF3Dcen, evtWeight);
                                    hWEEC3D_fakes[iDphi]->Fill(rF3Dcen, evtWeight);                                    
                                }
                                if (fillMeas) hWEEC3D_meas[iDphi]->Fill(rF3Dcen, evtWeight);
                            }
                        }
                    } else { // hasTruth only
                        const double tPairW = tpT_i*tpT_j/tPTmean2;
                        if (tPairW >= pairWeightMin && tPairW <= pairWeightMax) {
                            const int tIPw = FindBin(tPairW, pairWeightBins);
                            if (tIPw >= 0) {
                                const double tF3Dcen = TrueFlat3DBinCenter(TrueFlat3DIndex(iLtrue,iStrue,tIPw));
                                if (fillResp) {
                                    respWEEC3D[iDphi]->Miss(tF3Dcen, evtWeight);
                                    hWEEC3D_truth[iDphi]->Fill(tF3Dcen, evtWeight);
                                    hWEEC3D_pwNum[iDphi]->Fill(tF3Dcen, evtWeight * tPairW);
                                    hWEEC3D_pwDen[iDphi]->Fill(tF3Dcen, evtWeight);
                                    hWEEC3D_misses[iDphi]->Fill(tF3Dcen, evtWeight);
                                }
                                if (fillMeas) {
                                    hWEEC3D_meas_truth[iDphi]->Fill(tF3Dcen, evtWeight);
                                }
                            }
                        }
                    }
                }
            }
        }

        // ── fake dijet ────────────────────────────────────────────────────
        // Reco dijet with no matching truth dijet.  One reco-pair loop.
        //else if (isFake)
        if (isFake)
        {
            const double rPTmean2 = 0.25*(rLeadPT+rSublPT)*(rLeadPT+rSublPT);

            for (int i = 0; i < 24*64; ++i) {
                const int ei = i/64, pi = i%64;
                double phi_i = getPhiCenter(pi);
                const double rpT_i = recoMap[ei][pi]/cosh(getEtaCenter(ei, vtx_z));
                if (rpT_i <= towCut || rpT_i >= 80) continue;

                for (int j = i+1; j < 24*64; ++j) {
                    const int ej = j/64, pj = j%64;
                    double phi_j = getPhiCenter(pj);
                    const double rpT_j = recoMap[ej][pj]/cosh(getEtaCenter(ej, vtx_z));
                    if (rpT_j <= towCut || rpT_j >= 80) continue;

                    const double dphi   = DeltaPhi(phi_i, phi_j);
                    const double rPairW = rpT_i*rpT_j/rPTmean2;
                    if (rPairW < pairWeightMin || rPairW > pairWeightMax) continue;

                    const int iDphi = FindBin(dphi, dPhiBins);
                    const int rIPw  = FindBin(rPairW, pairWeightBins);
                    if (iDphi < 0 || rIPw < 0) continue;

                    const double rF3Dcen = RecoFlat3DBinCenter(RecoFlat3DIndex(iLreco,iSreco,rIPw));
                    if (fillResp) {
                        respWEEC3D[iDphi]->Fake(rF3Dcen, evtWeight);
                        hWEEC3D_fakes[iDphi]->Fill(rF3Dcen, evtWeight);
                    }
                    if (fillMeas) hWEEC3D_meas[iDphi]->Fill(rF3Dcen, evtWeight);
                }
            }
        }

        // ── missed dijet ──────────────────────────────────────────────────
        // Truth dijet with no matching reco dijet.  One truth-pair loop.
        //else if (isMiss)
        if (isMiss)
        {
            const double tPTmean2 = 0.25*(tLeadPT+tSublPT)*(tLeadPT+tSublPT);

            for (int i = 0; i < 24*64; ++i) {
                const int ei = i/64, pi = i%64;
                double phi_i = getPhiCenter(pi);
                const double tpT_i = truthMap[ei][pi]/cosh(getEtaCenter(ei, 0.0));
                if (tpT_i <= towCut || tpT_i >= 80) continue;

                for (int j = i+1; j < 24*64; ++j) {
                    const int ej = j/64, pj = j%64;
                    double phi_j = getPhiCenter(pj);
                    const double tpT_j = truthMap[ej][pj]/cosh(getEtaCenter(ej, 0.0));
                    if (tpT_j <= towCut || tpT_j >= 80) continue;

                    const double dphi   = DeltaPhi(phi_i, phi_j);
                    const double tPairW = tpT_i*tpT_j/tPTmean2;
                    if (tPairW < pairWeightMin || tPairW > pairWeightMax) continue;

                    const int iDphi = FindBin(dphi, dPhiBins);
                    const int tIPw  = FindBin(tPairW, pairWeightBins);
                    if (iDphi < 0 || tIPw < 0) continue;

                    const double tF3Dcen = TrueFlat3DBinCenter(TrueFlat3DIndex(iLtrue,iStrue,tIPw));
                    if (fillResp) {
                        respWEEC3D[iDphi]->Miss(tF3Dcen, evtWeight);
                        hWEEC3D_truth[iDphi]->Fill(tF3Dcen, evtWeight);
                        hWEEC3D_pwNum[iDphi]->Fill(tF3Dcen, evtWeight * tPairW);
                        hWEEC3D_pwDen[iDphi]->Fill(tF3Dcen, evtWeight);
                        hWEEC3D_misses[iDphi]->Fill(tF3Dcen, evtWeight);
                    }
                    if (fillMeas) {
                        hWEEC3D_meas_truth[iDphi]->Fill(tF3Dcen, evtWeight);
                    }
                }
            }
        }

        /*
        // ── measured-half reco only (kHalf: not matched/fake/miss on resp side) ──
        // Events where !fillResp but fillMeas and recoInRange — pure reco fill.
        else if (!fillResp && fillMeas && recoInRange)
        {
            const double rPTmean2 = 0.25*(rLeadPT+rSublPT)*(rLeadPT+rSublPT);

            for (int i = 0; i < 24*64; ++i) {
                const int ei = i/64, pi = i%64;
                if (recoMap[ei][pi] <= 0.25 || recoMap[ei][pi] >= 80) continue;
                double phi_i = getPhiCenter(pi);
                while (phi_i >  TMath::Pi()) phi_i -= 2*TMath::Pi();
                while (phi_i < -TMath::Pi()) phi_i += 2*TMath::Pi();
                const double rpT_i = recoMap[ei][pi]/cosh(getEtaCenter(ei, vtx_z));

                for (int j = i+1; j < 24*64; ++j) {
                    const int ej = j/64, pj = j%64;
                    if (recoMap[ej][pj] <= 0.25 || recoMap[ej][pj] >= 80) continue;
                    double phi_j = getPhiCenter(pj);
                    while (phi_j >  TMath::Pi()) phi_j -= 2*TMath::Pi();
                    while (phi_j < -TMath::Pi()) phi_j += 2*TMath::Pi();
                    const double rpT_j = recoMap[ej][pj]/cosh(getEtaCenter(ej, vtx_z));

                    const double dphi   = DeltaPhi(phi_i, phi_j);
                    const double rPairW = rpT_i*rpT_j/rPTmean2;
                    if (rPairW < pairWeightMin || rPairW > pairWeightMax) continue;

                    const int iDphi = FindBin(dphi, dPhiBins);
                    const int rIPw  = FindBin(rPairW, pairWeightBins);
                    if (iDphi < 0 || rIPw < 0) continue;

                    hWEEC3D_meas[iDphi]->Fill(
                        RecoFlat3DBinCenter(RecoFlat3DIndex(iLreco,iSreco,rIPw)), evtWeight);
                }
            }
        }

        // ── measured-half truth (kHalf only) ─────────────────────────────
        // Fill truth-tower pairs for odd-half events that are trueInRange
        // but were not matched/fake/miss (i.e. !fillResp).  Events that
        // were matched/fake/miss already filled hWEEC3D_meas_truth above.
        if (!fillResp && fillMeas && trueInRange)
        {
            const double tPTmean2 = 0.25*(tLeadPT+tSublPT)*(tLeadPT+tSublPT);

            for (int i = 0; i < 24*64; ++i) {
                const int ei = i/64, pi = i%64;
                if (truthMap[ei][pi] <= 0.25 || truthMap[ei][pi] >= 80) continue;
                double phi_i = getPhiCenter(pi);
                while (phi_i >  TMath::Pi()) phi_i -= 2*TMath::Pi();
                while (phi_i < -TMath::Pi()) phi_i += 2*TMath::Pi();
                const double tpT_i = truthMap[ei][pi]/cosh(getEtaCenter(ei, 0.0));

                for (int j = i+1; j < 24*64; ++j) {
                    const int ej = j/64, pj = j%64;
                    if (truthMap[ej][pj] <= 0.25 || truthMap[ej][pj] >= 80) continue;
                    double phi_j = getPhiCenter(pj);
                    while (phi_j >  TMath::Pi()) phi_j -= 2*TMath::Pi();
                    while (phi_j < -TMath::Pi()) phi_j += 2*TMath::Pi();
                    const double tpT_j = truthMap[ej][pj]/cosh(getEtaCenter(ej, 0.0));

                    const double dphi   = DeltaPhi(phi_i, phi_j);
                    const double tPairW = tpT_i*tpT_j/tPTmean2;
                    if (tPairW < pairWeightMin || tPairW > pairWeightMax) continue;

                    const int iDphi = FindBin(dphi, dPhiBins);
                    const int tIPw  = FindBin(tPairW, pairWeightBins);
                    if (iDphi < 0 || tIPw < 0) continue;

                    hWEEC3D_meas_truth[iDphi]->Fill(
                        TrueFlat3DBinCenter(TrueFlat3DIndex(iLtrue,iStrue,tIPw)), evtWeight);
                }
            }
        }
        */

    } // end event loop

    simFile->Close(); delete simFile;

    // ── write output ──────────────────────────────────────────────────────
    std::cout << "writing file " << outName << std::endl;
    std::cout.flush();

    fOut->cd();

    respJetPt->Write("response_jet_pT");
    hJetPt_truth_matched->Write();
    hJetPt_truth        ->Write();
    hJetPt_reco_matched ->Write();
    hJetPt_reco         ->Write();
    hJetPt_meas         ->Write();
    if (halfClosure) hJetPt_truth_meas->Write();

    std::cout << "done writing jet plots" << std::endl;

    for (int k = 0; k < nDphi; ++k)
        respWEEC3D[k]->Write(std::format("response_wEEC3D_{}", k).c_str());

    std::cout << "done writing response matrices" << std::endl;

    for (int k = 0; k < nDphi; ++k)
        hWEEC3D_truth[k]->Write();

    std::cout << "done writing truth plots" << std::endl;

    for (int k = 0; k < nDphi; ++k) {
        hWEEC3D_pwNum[k]->Write();
        hWEEC3D_pwDen[k]->Write();
    }

    std::cout << "done writing pair weight mean hists" << std::endl;

    for (int k = 0; k < nDphi; ++k)
        hWEEC3D_meas[k]->Write();

    std::cout << "done writing meas plots" << std::endl;

    if (halfClosure) {
        for (int k = 0; k < nDphi; ++k)
            hWEEC3D_meas_truth[k]->Write();
        std::cout << "done writing meas truth plots" << std::endl;
    }

    for (int k = 0; k < nDphi; ++k)
        hWEEC3D_counts[k]->Write();

    for (int k = 0; k < nDphi; ++k)
        hWEEC3D_misses[k]->Write();

    for (int k = 0; k < nDphi; ++k)
        hWEEC3D_fakes[k]->Write();

    std::cout << "done writing counts hists" << std::endl;

    hVtxMC->Write();

    std::cout << "done writing vtx hist" << std::endl;

    hWEEC_particle  ->Write();
    hWEEC_truthTower->Write();

    std::cout << "done writing cross-check wEEC hists" << std::endl;

    fOut->Close();
    std::cout << "Written: " << outName << "\n";
}
