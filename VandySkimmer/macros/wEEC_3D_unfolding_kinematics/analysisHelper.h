#pragma once

#include <fastjet/PseudoJet.hh>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <format>

R__LOAD_LIBRARY(libVandyClasses.so)

// ═════════════════════════════════════════════════════════════════════════════
//  Analysis mode
//    kFull — response and pseudo-data from ALL MC events
//    kHalf — response from random half, pseudo-data from other half
//    kData — real data; use fillData.C instead of fillResponse.C
//    kVtx  — vtx distribution pass only; skips all tower-pair loops
// ═════════════════════════════════════════════════════════════════════════════
enum class Mode { kFull, kHalf, kData, kVtx };

inline std::string ModeLabel(Mode mode)
{
    switch (mode) {
        case Mode::kFull: return "fullClosure";
        case Mode::kHalf: return "halfClosure";
        case Mode::kData: return "dataClosure";
        case Mode::kVtx:  return "vtx";
        default:          return "unknown";
    }
}

// ═════════════════════════════════════════════════════════════════════════════
//  Dijet pT binning — separate reco and truth spaces
//
//  Reco bins span only the reco-accessible region (lead >= recoLeadMin,
//  subl >= recoSublMin). These are the bins RooUnfold unfolds into and the
//  only bins that appear in the measured histogram.
//
//  Truth bins extend below the reco thresholds to capture the smearing region.
//  Events with truth pT in [trueLeadMin, recoLeadMin) or
//  [trueSublMin, recoSublMin) are pure misses in the response matrix —
//  they correctly account for events that smeared above the reco threshold
//  without requiring RooUnfold to invert a region with no reco data.
//
//  The unfolded output spans the full truth bin space but only bins with
//  lead pT >= recoLeadMin and subl pT >= recoSublMin are reported.
// ═════════════════════════════════════════════════════════════════════════════

// Reco-side bins: only bins accessible at reco level
const std::vector<double> recoLeadPtBins = { 15, 20, 30, 40, 50, 60, 80 };
//const std::vector<double> recoLeadPtBins = { 30, 40, 50, 60, 80 };
const std::vector<double> recoSublPtBins = {  8, 15, 20, 30, 40, 50, 60, 80 };
//const std::vector<double> recoSublPtBins = {  15, 20, 30, 40, 50, 60, 80 };
const int nRecoLead = (int)recoLeadPtBins.size() - 1;  // 6
const int nRecoSubl = (int)recoSublPtBins.size() - 1;  // 7

// Truth-side bins: extend into the smearing region below reco thresholds
const std::vector<double> trueLeadPtBins = { 15, 20, 30, 40, 50, 60, 80 };
const std::vector<double> trueSublPtBins = {  8, 15, 20, 30, 40, 50, 60, 80 };
const int nTrueLead = (int)trueLeadPtBins.size() - 1;  // 6
const int nTrueSubl = (int)trueSublPtBins.size() - 1;  // 7

// Legacy aliases used by dijet pT unfolding (reco-side only)
const std::vector<double>& leadPtBins = recoLeadPtBins;
const std::vector<double>& sublPtBins = recoSublPtBins;
const int nLead = nRecoLead;
const int nSubl = nRecoSubl;

// ═════════════════════════════════════════════════════════════════════════════
//  Kinematically-filtered flat dijet pT index
//
//  A (iL, iS) pair is valid when  sublPtBins[iS] < leadPtBins[iL+1].
//  This allows same-bin combinations like 30-40/30-40 (a 38/32 GeV dijet is
//  perfectly physical) while excluding impossible pairs like 40-50/50-60.
//
//  The index tables below are built immediately from the bin-edge vectors
//  defined above, as inline variables, so they are guaranteed consistent
//  everywhere and require no explicit initialisation call.
//
//  API — sizes:
//    nRecoFlat   number of valid reco (iL,iS) pairs
//    nTrueFlat   number of valid truth (iL,iS) pairs
//    nRecoFlat3D nRecoFlat * nPairWeight
//    nTrueFlat3D nTrueFlat * nPairWeight
//
//  API — 2D index:
//    RecoFlatIndex(iL,iS)  → compressed flat index, or -1 if forbidden
//    TrueFlatIndex(iL,iS)  → compressed flat index, or -1 if forbidden
//    RecoFlatToIJ(f,iL,iS) → recover (iL,iS) from flat index
//    TrueFlatToIJ(f,iL,iS) → recover (iL,iS) from flat index
//    RecoFlatBinCenter(f)  → bin centre for histogram Fill/Fake/Miss
//    TrueFlatBinCenter(f)  → bin centre for histogram Fill/Miss
//
//  API — 3D index (pair weight as innermost dimension, so pair-weight bins
//        for a given (iL,iS) are always contiguous):
//    RecoFlat3DIndex(iL,iS,iPw), TrueFlat3DIndex(iL,iS,iPw)
//    RecoFlat3DToIJK(iF,iL,iS,iPw), TrueFlat3DToIJK(iF,iL,iS,iPw)
//    RecoFlat3DBinCenter(iF), TrueFlat3DBinCenter(iF)
// ═════════════════════════════════════════════════════════════════════════════

// ── helper structs ────────────────────────────────────────────────────────────
struct DijetPair  { int iL, iS; };

// ── build valid pair lists (defined as inline to initialise exactly once) ─────

// Reco valid pairs
inline std::vector<DijetPair> MakeRecoPairs()
{
    std::vector<DijetPair> v;
    for (int iL = 0; iL < nRecoLead; ++iL)
        for (int iS = 0; iS < nRecoSubl; ++iS)
            if (recoSublPtBins[iS] < recoLeadPtBins[iL+1])
                v.push_back({iL, iS});
    return v;
}
inline std::map<std::pair<int,int>,int> MakeRecoFlatMap(const std::vector<DijetPair>& pairs)
{
    std::map<std::pair<int,int>,int> m;
    for (int f = 0; f < (int)pairs.size(); ++f)
        m[{pairs[f].iL, pairs[f].iS}] = f;
    return m;
}

// Truth valid pairs
inline std::vector<DijetPair> MakeTruePairs()
{
    std::vector<DijetPair> v;
    for (int iL = 0; iL < nTrueLead; ++iL)
        for (int iS = 0; iS < nTrueSubl; ++iS)
            if (trueSublPtBins[iS] < trueLeadPtBins[iL+1])
                v.push_back({iL, iS});
    return v;
}
inline std::map<std::pair<int,int>,int> MakeTrueFlatMap(const std::vector<DijetPair>& pairs)
{
    std::map<std::pair<int,int>,int> m;
    for (int f = 0; f < (int)pairs.size(); ++f)
        m[{pairs[f].iL, pairs[f].iS}] = f;
    return m;
}

// ── global index tables (initialised once at program start) ───────────────────
inline const std::vector<DijetPair>          gRecoPairs    = MakeRecoPairs();
inline const std::map<std::pair<int,int>,int> gRecoFlatMap = MakeRecoFlatMap(gRecoPairs);
inline const std::vector<DijetPair>          gTruePairs    = MakeTruePairs();
inline const std::map<std::pair<int,int>,int> gTrueFlatMap = MakeTrueFlatMap(gTruePairs);

// ── size accessors ────────────────────────────────────────────────────────────
// nPairWeight is defined later in the file, so nRecoFlat3D / nTrueFlat3D are
// deferred to after nPairWeight is defined (see below).
inline int nRecoFlat  () { return (int)gRecoPairs.size(); }
inline int nTrueFlat  () { return (int)gTruePairs.size(); }

// ── 2D index functions ────────────────────────────────────────────────────────
inline int RecoFlatIndex(int iL, int iS)
{
    auto it = gRecoFlatMap.find({iL, iS});
    return (it != gRecoFlatMap.end()) ? it->second : -1;
}
inline int TrueFlatIndex(int iL, int iS)
{
    auto it = gTrueFlatMap.find({iL, iS});
    return (it != gTrueFlatMap.end()) ? it->second : -1;
}
inline void RecoFlatToIJ(int f, int& iL, int& iS)
    { iL = gRecoPairs[f].iL; iS = gRecoPairs[f].iS; }
inline void TrueFlatToIJ(int f, int& iL, int& iS)
    { iL = gTruePairs[f].iL; iS = gTruePairs[f].iS; }
inline double RecoFlatBinCenter(int f) { return f + 0.5; }
inline double TrueFlatBinCenter(int f) { return f + 0.5; }

// Legacy aliases — dijet pT unfolding uses reco-side only
inline int    FlatIndex    (int iL, int iS) { return RecoFlatIndex(iL, iS); }
inline void   FlatToIJ     (int f, int& iL, int& iS) { RecoFlatToIJ(f, iL, iS); }
inline double FlatBinCenter(int f) { return RecoFlatBinCenter(f); }
// nFlat is used in a few legacy places — provide as a function for consistency
inline int nFlat() { return nRecoFlat(); }

// ─────────────────────────────────────────────────────────────────────────────
//  ΔΦ binning — 33 edges = 32 bins
// ─────────────────────────────────────────────────────────────────────────────
const std::vector<double> original_dPhiBins = {
    0.0000000,  0.041592654,
    0.14159265, 0.24159265, 0.34159265, 0.44159265, 0.54159265,
    0.64159265, 0.74159265, 0.84159265, 0.94159265, 1.0415927,
    1.1415927,  1.2415927,  1.3415927,  1.4415927,  1.5415927,
    1.6415927,  1.7415927,  1.8415927,  1.9415927,  2.0415927,
    2.1415927,  2.2415927,  2.3415927,  2.4415927,  2.5415927,
    2.6415927,  2.7415927,  2.8415927,  2.9415927,  3.0415927,
    3.1415927
};
const std::vector<double> new_dPhiBins = {
//    0.0000000, 0.041592654, 0.14159265,
    0.0000000, 0.14159265,
    0.34159265, 0.54159265, 0.74159265, 0.94159265,
    1.2415927, 1.5415927, 1.8415927, 2.1415927,
    2.3415927, 2.5415927, 2.7415927, 2.9415927,
    3.0415927, 3.1415927
};

//const std::vector<double> dPhiBins = original_dPhiBins;
const std::vector<double> dPhiBins = new_dPhiBins;
const int nDphi = (int)dPhiBins.size() - 1;  // 32


// ─────────────────────────────────────────────────────────────────────────────
//  Pair weight binning — 20 log-uniform bins from 1e-3 to 1
//
//  Pair weight: w_ij = pT_i * pT_j / <pT>^2
//  Log-uniform spacing ensures good resolution across the dynamic range.
//  Edges are computed at runtime to avoid floating-point literal clutter.
// ─────────────────────────────────────────────────────────────────────────────
const int    nPairWeight      = 17;
const double pairWeightMin    = 0.0;    // physical lower edge of histogram
const double pairWeightMax    = 2.0;
const double pairWeightLogMin = 5e-5;   // lower edge of the log-uniform region

// Pair weight bins:
//   Bin 0  : [0,    1e-3) — linear underflow bin for very small weights
//   Bins 1-20: log-uniform from 1e-3 to 1.0  (20 bins)
//   Bin 21 : [1.0,  2.0) — linear overflow bin
//
// The underflow bin ensures pairs with pairWeight < 1e-3 are not silently
// dropped; they are counted and can be trimmed by the count threshold in
// wEEC_doUnfolding.C if they have insufficient statistics.
inline std::vector<double> MakePairWeightBins()
{
    std::vector<double> edges(nPairWeight + 1);
    // Bin 0: [0, 1e-3)
    edges[0] = 0.0;
    // Bins 1-20: log-uniform from 1e-3 to 1.0
    const int    nLogBins = nPairWeight-2;
    const double logMin   = std::log10(pairWeightLogMin);
    const double logMax   = std::log10(0.2);
    for (int i = 0; i <= nLogBins; ++i)
        edges[i + 1] = std::pow(10.0, logMin + (logMax - logMin) * i / nLogBins);
    // Bin 21: [1.0, 2.0)
    edges[nPairWeight] = pairWeightMax;
    return edges;
}
const std::vector<double> pairWeightBins = MakePairWeightBins();

// ═════════════════════════════════════════════════════════════════════════════
//  Cross-section weights  (σ_k / N_gen)
//  Update N_gen if your full sample size differs from 9 600 000
// ═════════════════════════════════════════════════════════════════════════════
static const std::map<int,double> xsec = {
    {12,  1.4903e+06},
    {20,  6.2623e+04},
    {30,  2.5298e+03},
    {40,  1.3553e+02},
    {50,  7.3113e+00},
    {60,  3.3261e-01}
};
inline double getSampleWeight(int jetSample)
{
    return xsec.at(jetSample) / 9600000.0;
}

// ═════════════════════════════════════════════════════════════════════════════
//  pT-hat slice cuts  (applied to max truth jet pT)
// ═════════════════════════════════════════════════════════════════════════════
inline float get_jet_pTLow(int sample)
{
    switch (sample) {
        case  5: return   7;
        case 12: return  14;
        case 20: return  21;
        case 30: return  32;
        case 40: return  42;
        case 50: return  52;
        case 60: return  62;
        default: throw std::invalid_argument("Invalid jet sample");
    }
}
inline float get_jet_pTHigh(int sample)
{
    switch (sample) {
        case  5: return  14;
        case 12: return  21;
        case 20: return  32;
        case 30: return  42;
        case 40: return  52;
        case 50: return  62;
        case 60: return 1000;
        default: throw std::invalid_argument("Invalid jet sample");
    }
}

// ═════════════════════════════════════════════════════════════════════════════
//  Dijet validity cuts
// ═════════════════════════════════════════════════════════════════════════════
const double recoLeadMin = 20.0;
const double recoSublMin =  8.0;
const double trueLeadMin = 15.0;
const double trueSublMin =  8.0;

// No asymmetry cut applied — retaining all kinematically valid dijets
// to avoid biasing the pair weight and ΔΦ distributions.

// ═════════════════════════════════════════════════════════════════════════════
//  Geometric jet matching ΔR cut
// ═════════════════════════════════════════════════════════════════════════════
const double geoMatchDR = 0.3;



// ═════════════════════════════════════════════════════════════════════════════
//  Reco accessibility predicate
//
//  Returns true when the truth (iL, iS) dijet pT bin has at least partial
//  overlap with the reco-accessible region, i.e. when the upper edge of the
//  truth lead bin exceeds recoLeadMin AND the upper edge of the truth subl
//  bin exceeds recoSublMin.
//
//  Bins that fail this test are pure-miss bins at reco level.  The response
//  matrix correctly encodes them as misses, but the unfolded result in these
//  bins is unconstrained by measured data and should be excluded from all
//  projections and reported results.
// ═════════════════════════════════════════════════════════════════════════════
inline bool IsRecoAccessible(int iL, int iS)
{
    return trueLeadPtBins[iL + 1] > 30 &&
           trueSublPtBins[iS + 1] > 15;
}

// ═════════════════════════════════════════════════════════════════════════════
//  General bin-finding
// ═════════════════════════════════════════════════════════════════════════════
inline int FindBin(double val, const std::vector<double>& edges)
{
    int bin = (int)(std::upper_bound(edges.begin(), edges.end(), val)
                    - edges.begin()) - 1;
    return (bin >= 0 && bin < (int)edges.size() - 1) ? bin : -1;
}

// ═════════════════════════════════════════════════════════════════════════════
//  ΔΦ — returns value in [0, π]
// ═════════════════════════════════════════════════════════════════════════════
inline double DeltaPhi(double phi1, double phi2)
{
    double dphi = std::abs(phi1 - phi2);
    if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
    return dphi;
}

// ═════════════════════════════════════════════════════════════════════════════
//  ΔR
// ═════════════════════════════════════════════════════════════════════════════
inline double DeltaR(double eta1, double phi1, double eta2, double phi2)
{
    double deta = eta1 - eta2;
    double dphi = DeltaPhi(phi1, phi2);
    return std::sqrt(deta*deta + dphi*dphi);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Dijet geometric matching
// ═════════════════════════════════════════════════════════════════════════════
inline bool GetGeoMatch(
    std::vector<JetInfo>& recoJets,
    std::vector<JetInfo>& truthJets,
    double rLeadPT, double rSublPT,
    double tLeadPT, double tSublPT)
{
    // Find the specific jets whose pT matches the EventInfo values exactly.
    // This avoids ambiguity from sorting — EventInfo already determined which
    // jets are lead and subl, so we just locate them in the vector.
    const double pTtol = 0.01;  // GeV — tight, since EventInfo stores same value

    JetInfo* rLead = nullptr;
    JetInfo* rSubl = nullptr;
    JetInfo* tLead = nullptr;
    JetInfo* tSubl = nullptr;

    for (auto& j : recoJets) {
        if (!rLead && std::abs(j.pt() - rLeadPT) < pTtol) rLead = &j;
        else if (!rSubl && std::abs(j.pt() - rSublPT) < pTtol) rSubl = &j;
    }
    for (auto& j : truthJets) {
        if (!tLead && std::abs(j.pt() - tLeadPT) < pTtol) tLead = &j;
        else if (!tSubl && std::abs(j.pt() - tSublPT) < pTtol) tSubl = &j;
    }

    // If any of the four jets couldn't be located, can't match
    if (!rLead || !rSubl || !tLead || !tSubl) return false;

    fastjet::PseudoJet rL(rLead->px(), rLead->py(), rLead->pz(), rLead->e());
    fastjet::PseudoJet rS(rSubl->px(), rSubl->py(), rSubl->pz(), rSubl->e());
    fastjet::PseudoJet tL(tLead->px(), tLead->py(), tLead->pz(), tLead->e());
    fastjet::PseudoJet tS(tSubl->px(), tSubl->py(), tSubl->pz(), tSubl->e());

    return DeltaR(rL.pseudorapidity(), rL.phi_std(),
                  tL.pseudorapidity(), tL.phi_std()) < geoMatchDR &&
           DeltaR(rS.pseudorapidity(), rS.phi_std(),
                  tS.pseudorapidity(), tS.phi_std()) < geoMatchDR;
}

// ═════════════════════════════════════════════════════════════════════════════
//  Input file path
// ═════════════════════════════════════════════════════════════════════════════
inline std::string SimFilePath(int jetSample, int seg)
{
    return std::format(
        "/sphenix/tg/tg01/jets/sgross/VandyDSTs/Jet{0}/"
        "VandyDSTs_run2pp_ana521_2025p007_v001_Jet{0}-28-{1:06d}_to-{2:06d}.root",
        jetSample, seg, seg + 24);
}

// ═════════════════════════════════════════════════════════════════════════════
//  Calorimeter geometry
//
//  etaEdges / phiEdges: fixed (unshifted) tower boundaries
//  shiftedEtaEdges: vertex-corrected boundaries, updated by shiftEtaEdges()
//
//  Truth towers use the UNSHIFTED grid (vtx=0) so that reco and truth towers
//  share the same (iEta, iPhi) index space for matched-pair filling.
// ═════════════════════════════════════════════════════════════════════════════
const std::vector<double> etaEdges = {
    -1.1000000, -1.0083333, -0.91666667, -0.82500000, -0.73333333,
    -0.64166667, -0.55000000, -0.45833333, -0.36666667, -0.27500000,
    -0.18333333, -0.091666667, 1.1102230e-16, 0.091666667, 0.18333333,
     0.27500000,  0.36666667,  0.45833333,  0.55000000,  0.64166667,
     0.73333333,  0.82500000,  0.91666667,  1.0083333,   1.1000000 };

const std::vector<double> phiEdges = {
    -0.053619781, 0.044554989, 0.14272976, 0.24090453, 0.33907930,
     0.43725407,  0.53542884,  0.63360361, 0.73177838, 0.82995315,
     0.92812792,  1.0263027,   1.1244775,  1.2226522,  1.3208270,
     1.4190018,   1.5171765,   1.6153513,  1.7135261,  1.8117009,
     1.9098756,   2.0080504,   2.1062252,  2.2043999,  2.3025747,
     2.4007495,   2.4989242,   2.5970990,  2.6952738,  2.7934486,
     2.8916233,   2.9897981,   3.0879729,  3.1861476,  3.2843224,
     3.3824972,   3.4806720,   3.5788467,  3.6770215,  3.7751963,
     3.8733710,   3.9715458,   4.0697206,  4.1678953,  4.2660701,
     4.3642449,   4.4624197,   4.5605944,  4.6587692,  4.7569440,
     4.8551187,   4.9532935,   5.0514683,  5.1496431,  5.2478178,
     5.3459926,   5.4441674,   5.5423421,  5.6405169,  5.7386917,
     5.8368664,   5.9350412,   6.0332160,  6.1313908,  6.2295655 };

std::vector<double> shiftedEtaEdges;
const double R = 124.5;

void shiftEtaEdges(double vtx)
{
    shiftedEtaEdges.clear();
    for (auto edge : etaEdges) {
        double z = R * sinh(edge);
        shiftedEtaEdges.push_back(asinh((z - vtx) / R));
    }
}

double getEtaCenter(int etaBin, double vtx)
{
    if (vtx != 0.0)
        return 0.5 * (shiftedEtaEdges[etaBin] + shiftedEtaEdges[etaBin+1]);
    return 0.5 * (etaEdges[etaBin] + etaEdges[etaBin+1]);
}

double getPhiCenter(int phiBin)
{
    return 0.5 * (phiEdges[phiBin] + phiEdges[phiBin+1]);
}

int getEtaBin(double eta, double vtx)
{
    const std::vector<double>& edges = (vtx != 0) ? shiftedEtaEdges : etaEdges;
    return FindBin(eta, edges);
}

int getPhiBin(double phi)
{
    while (phi < phiEdges.front()) phi += 2 * TMath::Pi();
    while (phi > phiEdges.back())  phi -= 2 * TMath::Pi();
    return FindBin(phi, phiEdges);
}

// ═════════════════════════════════════════════════════════════════════════════
//  wEEC flat 2D index helpers
//
//  The (ΔΦ, pair weight) space is stored as a flattened 1D histogram with
//  nFlat2D = nDphi * nPairWeight bins, where:
//    iFlat2D = iDphi * nPairWeight + iPairWeight
//
//  This matches the layout used by the wEEC RooUnfoldResponse matrices.
//  Convert to a proper TH2D only for drawing.
// ═════════════════════════════════════════════════════════════════════════════
const int nFlat2D = nDphi * nPairWeight;  // 32 * 22 = 704

inline int   Flat2DIndex(int iDphi, int iPw) { return iDphi * nPairWeight + iPw; }
inline void  Flat2DToIJ (int iF2D, int& iDphi, int& iPw)
             { iDphi = iF2D / nPairWeight; iPw = iF2D % nPairWeight; }
inline double Flat2DBinCenter(int iF2D) { return iF2D + 0.5; }

// Convert a flat 1D wEEC histogram to a proper TH2D for plotting.
// The caller owns the returned histogram.
inline TH2D* Flat1DToTH2D(TH1D* h1D, const char* name, const char* title = "")
{
    TH2D* h2D = new TH2D(name, title,
                          nDphi,       dPhiBins.data(),
                          nPairWeight, pairWeightBins.data());
    h2D->Sumw2();
    for (int iF2D = 0; iF2D < nFlat2D; ++iF2D) {
        int iDphi, iPw; Flat2DToIJ(iF2D, iDphi, iPw);
        h2D->SetBinContent(iDphi+1, iPw+1, h1D->GetBinContent(iF2D+1));
        h2D->SetBinError  (iDphi+1, iPw+1, h1D->GetBinError  (iF2D+1));
    }
    return h2D;
}

// ═════════════════════════════════════════════════════════════════════════════
//  wEEC 3D flat index helpers — separate reco and truth spaces
//
//  Layout: pair weight is the innermost dimension, so all nPairWeight bins for
//  a given (iL,iS) are contiguous. This makes projecting over pair weight a
//  simple stride rather than a scattered gather.
//
//    iFlat3D = TrueFlatIndex(iL,iS) * nPairWeight + iPw
//
//  Only kinematically valid (iL,iS) pairs are included (see 2D index above),
//  so every bin in the 3D histogram has a nonzero response matrix row and
//  RooUnfold never receives a zero-prior truth bin.
//
//  Sizes are functions (not constants) because they depend on nPairWeight,
//  which is defined immediately above, and on nTrueFlat()/nRecoFlat().
// ═════════════════════════════════════════════════════════════════════════════
inline int nRecoFlat3D() { return nRecoFlat() * nPairWeight; }
inline int nTrueFlat3D() { return nTrueFlat() * nPairWeight; }

// Reco-side flat3D — returns -1 for forbidden (iL,iS) combinations
inline int RecoFlat3DIndex(int iL, int iS, int iPw)
{
    int f2 = RecoFlatIndex(iL, iS);
    return (f2 >= 0) ? f2 * nPairWeight + iPw : -1;
}
inline void RecoFlat3DToIJK(int iF, int& iL, int& iS, int& iPw)
{
    int f2 = iF / nPairWeight;
    iPw = iF % nPairWeight;
    RecoFlatToIJ(f2, iL, iS);
}
inline double RecoFlat3DBinCenter(int iF) { return iF + 0.5; }

// Truth-side flat3D — returns -1 for forbidden (iL,iS) combinations
inline int TrueFlat3DIndex(int iL, int iS, int iPw)
{
    int f2 = TrueFlatIndex(iL, iS);
    return (f2 >= 0) ? f2 * nPairWeight + iPw : -1;
}
inline void TrueFlat3DToIJK(int iF, int& iL, int& iS, int& iPw)
{
    int f2 = iF / nPairWeight;
    iPw = iF % nPairWeight;
    TrueFlatToIJ(f2, iL, iS);
}
inline double TrueFlat3DBinCenter(int iF) { return iF + 0.5; }

// Project a flat 3D truth histogram for one ΔΦ bin down to a single wEEC
// value, summing only over truth (iL,iS) bins that are reco-accessible.
// Bins below the reco threshold have no measured support and are excluded
// so they don't dilute the inclusive projection.
inline void ProjectWEEC3DSlice(TH1D* h1D, double& val, double& err)
{
    double sum = 0.0, sumw2 = 0.0;
    for (int ft = 0; ft < nTrueFlat(); ++ft) {
        int iL, iS; TrueFlatToIJ(ft, iL, iS);
        if (!IsRecoAccessible(iL, iS)) continue;
        for (int iPw = 0; iPw < nPairWeight; ++iPw) {
            int    iF   = ft * nPairWeight + iPw;
            double wCen = 0.5 * (pairWeightBins[iPw] + pairWeightBins[iPw+1]);
            double v    = h1D->GetBinContent(iF + 1);
            double e    = h1D->GetBinError  (iF + 1);
            sum   += wCen * v;
            sumw2 += wCen * wCen * e * e;
        }
    }
    val = sum;
    err = std::sqrt(sumw2);
}

// Project for a specific (iL, iS) truth dijet pT bin.
// Returns immediately with val=err=0 if the pair is kinematically forbidden.
// Because pair weight is the innermost dimension, the nPairWeight bins for
// this (iL,iS) are contiguous starting at ft*nPairWeight.
inline void ProjectWEEC3DSliceForBin(TH1D* h1D, int iL, int iS,
                                      double& val, double& err)
{
    val = 0; err = 0;
    int ft = TrueFlatIndex(iL, iS);
    if (ft < 0) return;  // kinematically forbidden
    double sum = 0.0, sumw2 = 0.0;
    int iF0 = ft * nPairWeight;
    for (int iPw = 0; iPw < nPairWeight; ++iPw) {
        double wCen = 0.5 * (pairWeightBins[iPw] + pairWeightBins[iPw+1]);
        double v    = h1D->GetBinContent(iF0 + iPw + 1);
        double e    = h1D->GetBinError  (iF0 + iPw + 1);
        sum   += wCen * v;
        sumw2 += wCen * wCen * e * e;
    }
    val = sum;
    err = std::sqrt(sumw2);
}

// ── Weighted-mean pair weight variants ───────────────────────────────────────
// These replace the bin-center approximation with the luminosity-weighted mean
// pair weight computed from hPWNum (sum of evtWeight*pairWeight) and hPWDen
// (sum of evtWeight) filled in fillResponse.C.  For cells where hPWDen==0
// (no events), the bin-center is used as a fallback.

inline void ProjectWEEC3DSlice(TH1D* h1D, double& val, double& err,
                                TH1D* hPWNum, TH1D* hPWDen)
{
    double sum = 0.0, sumw2 = 0.0;
    for (int ft = 0; ft < nTrueFlat(); ++ft) {
        int iL, iS; TrueFlatToIJ(ft, iL, iS);
        if (!IsRecoAccessible(iL, iS)) continue;
        int iF0 = ft * nPairWeight;
        for (int iPw = 0; iPw < nPairWeight; ++iPw) {
            int    iF   = iF0 + iPw;
            double den  = hPWDen ? hPWDen->GetBinContent(iF + 1) : 0.0;
            double wMean = (den > 0)
                ? hPWNum->GetBinContent(iF + 1) / den
                : 0.5 * (pairWeightBins[iPw] + pairWeightBins[iPw+1]);
            double v    = h1D->GetBinContent(iF + 1);
            double e    = h1D->GetBinError  (iF + 1);
            sum   += wMean * v;
            sumw2 += wMean * wMean * e * e;
        }
    }
    val = sum;
    err = std::sqrt(sumw2);
}

inline void ProjectWEEC3DSliceForBin(TH1D* h1D, int iL, int iS,
                                      double& val, double& err,
                                      TH1D* hPWNum, TH1D* hPWDen)
{
    val = 0; err = 0;
    int ft = TrueFlatIndex(iL, iS);
    if (ft < 0) return;
    double sum = 0.0, sumw2 = 0.0;
    int iF0 = ft * nPairWeight;
    for (int iPw = 0; iPw < nPairWeight; ++iPw) {
        int    iF    = iF0 + iPw;
        double den   = hPWDen ? hPWDen->GetBinContent(iF + 1) : 0.0;
        double wMean = (den > 0)
            ? hPWNum->GetBinContent(iF + 1) / den
            : 0.5 * (pairWeightBins[iPw] + pairWeightBins[iPw+1]);
        double v    = h1D->GetBinContent(iF + 1);
        double e    = h1D->GetBinError  (iF + 1);
        sum   += wMean * v;
        sumw2 += wMean * wMean * e * e;
    }
    val = sum;
    err = std::sqrt(sumw2);
}

// ── Covariance-aware projection overloads ────────────────────────────────────
// These use the full TMatrixD covariance from RooUnfold::Errunfold() so that
// off-diagonal correlations between unfolded bins are properly accounted for
// when projecting onto a subset of the flat-3D space.
//
// For a linear combination P = sum_i w_i * x_i over bins in set S:
//   Var(P) = sum_{i,j in S} w_i * w_j * C_{ij}
//
// The covariance matrix is 0-based with size nTrueFlat3D() × nTrueFlat3D().
// Pass nullptr for cov to fall back to naive diagonal (GetBinError) propagation.

inline void ProjectWEEC3DSlice(TH1D* h1D, double& val, double& err,
                                const TMatrixD* cov)
{
    double sum = 0.0, var = 0.0;
    // Collect (flat3D index, weight) pairs for all accessible bins
    std::vector<std::pair<int,double>> bins;
    for (int ft = 0; ft < nTrueFlat(); ++ft) {
        int iL, iS; TrueFlatToIJ(ft, iL, iS);
        if (!IsRecoAccessible(iL, iS)) continue;
        for (int iPw = 0; iPw < nPairWeight; ++iPw) {
            int    iF   = ft * nPairWeight + iPw;
            double wCen = 0.5 * (pairWeightBins[iPw] + pairWeightBins[iPw+1]);
            sum += wCen * h1D->GetBinContent(iF + 1);
            bins.push_back({iF, wCen});
        }
    }
    if (cov) {
        for (auto [i, wi] : bins)
        for (auto [j, wj] : bins)
            var += wi * wj * (*cov)(i, j);
    } else {
        for (auto [i, wi] : bins) {
            double e = h1D->GetBinError(i + 1);
            var += wi * wi * e * e;
        }
    }
    val = sum;
    err = std::sqrt(std::max(var, 0.0));
}

inline void ProjectWEEC3DSliceForBin(TH1D* h1D, int iL, int iS,
                                      double& val, double& err,
                                      const TMatrixD* cov)
{
    val = 0; err = 0;
    int ft = TrueFlatIndex(iL, iS);
    if (ft < 0) return;
    double sum = 0.0, var = 0.0;
    int iF0 = ft * nPairWeight;
    for (int iPw = 0; iPw < nPairWeight; ++iPw) {
        int    iF   = iF0 + iPw;
        double wCen = 0.5 * (pairWeightBins[iPw] + pairWeightBins[iPw+1]);
        sum += wCen * h1D->GetBinContent(iF + 1);
    }
    if (cov) {
        for (int iPw = 0; iPw < nPairWeight; ++iPw) {
            int    iF = iF0 + iPw;
            double wi = 0.5 * (pairWeightBins[iPw] + pairWeightBins[iPw+1]);
            for (int jPw = 0; jPw < nPairWeight; ++jPw) {
                int    jF = iF0 + jPw;
                double wj = 0.5 * (pairWeightBins[jPw] + pairWeightBins[jPw+1]);
                var += wi * wj * (*cov)(iF, jF);
            }
        }
    } else {
        for (int iPw = 0; iPw < nPairWeight; ++iPw) {
            int    iF = iF0 + iPw;
            double wCen = 0.5 * (pairWeightBins[iPw] + pairWeightBins[iPw+1]);
            double e    = h1D->GetBinError(iF + 1);
            var += wCen * wCen * e * e;
        }
    }
    val = sum;
    err = std::sqrt(std::max(var, 0.0));
}

// Covariance-aware variants with weighted mean pair weight
inline void ProjectWEEC3DSlice(TH1D* h1D, double& val, double& err,
                                TH1D* hPWNum, TH1D* hPWDen,
                                const TMatrixD* cov)
{
    double sum = 0.0, var = 0.0;
    std::vector<std::pair<int,double>> bins;
    for (int ft = 0; ft < nTrueFlat(); ++ft) {
        int iL, iS; TrueFlatToIJ(ft, iL, iS);
        if (!IsRecoAccessible(iL, iS)) continue;
        int iF0 = ft * nPairWeight;
        for (int iPw = 0; iPw < nPairWeight; ++iPw) {
            int    iF  = iF0 + iPw;
            double den = hPWDen ? hPWDen->GetBinContent(iF + 1) : 0.0;
            double w   = (den > 0)
                ? hPWNum->GetBinContent(iF + 1) / den
                : 0.5 * (pairWeightBins[iPw] + pairWeightBins[iPw+1]);
            sum += w * h1D->GetBinContent(iF + 1);
            bins.push_back({iF, w});
        }
    }
    if (cov) {
        for (auto [i, wi] : bins)
        for (auto [j, wj] : bins)
            var += wi * wj * (*cov)(i, j);
    } else {
        for (auto [i, wi] : bins) {
            double e = h1D->GetBinError(i + 1);
            var += wi * wi * e * e;
        }
    }
    val = sum;
    err = std::sqrt(std::max(var, 0.0));
}

inline void ProjectWEEC3DSliceForBin(TH1D* h1D, int iL, int iS,
                                      double& val, double& err,
                                      TH1D* hPWNum, TH1D* hPWDen,
                                      const TMatrixD* cov)
{
    val = 0; err = 0;
    int ft = TrueFlatIndex(iL, iS);
    if (ft < 0) return;
    double sum = 0.0, var = 0.0;
    int iF0 = ft * nPairWeight;
    // Build per-iPw weights
    std::vector<double> ws(nPairWeight);
    for (int iPw = 0; iPw < nPairWeight; ++iPw) {
        int    iF  = iF0 + iPw;
        double den = hPWDen ? hPWDen->GetBinContent(iF + 1) : 0.0;
        ws[iPw]    = (den > 0)
            ? hPWNum->GetBinContent(iF + 1) / den
            : 0.5 * (pairWeightBins[iPw] + pairWeightBins[iPw+1]);
        sum += ws[iPw] * h1D->GetBinContent(iF + 1);
    }
    if (cov) {
        for (int iPw = 0; iPw < nPairWeight; ++iPw)
        for (int jPw = 0; jPw < nPairWeight; ++jPw)
            var += ws[iPw] * ws[jPw] * (*cov)(iF0 + iPw, iF0 + jPw);
    } else {
        for (int iPw = 0; iPw < nPairWeight; ++iPw) {
            double e = h1D->GetBinError(iF0 + iPw + 1);
            var += ws[iPw] * ws[iPw] * e * e;
        }
    }
    val = sum;
    err = std::sqrt(std::max(var, 0.0));
}

// ═════════════════════════════════════════════════════════════════════════════
//  wEEC normalization helpers
// ═════════════════════════════════════════════════════════════════════════════

// Project a flat 1D wEEC histogram onto the ΔΦ axis by summing bin contents
// weighted by the pair weight bin center in each pair weight slice.
// This correctly accounts for the log-uniform pair weight binning — a plain
// sum would integrate w * bin_width rather than w itself.
// Returns a new TH1D owned by the caller.
inline TH1D* ProjectWEEC(TH1D* h1D, const char* name)
{
    TH1D* hProj = new TH1D(name, ";#Delta#phi;wEEC",
                            nDphi, dPhiBins.data());
    hProj->Sumw2();
    for (int iDphi = 0; iDphi < nDphi; ++iDphi) {
        double sum   = 0.0;
        double sumw2 = 0.0;
        for (int iPw = 0; iPw < nPairWeight; ++iPw) {
            int    iF2D = Flat2DIndex(iDphi, iPw);
            double wCen = 0.5 * (pairWeightBins[iPw] + pairWeightBins[iPw+1]);
            double val  = h1D->GetBinContent(iF2D + 1);
            double err  = h1D->GetBinError  (iF2D + 1);
            sum   += wCen * val;
            sumw2 += wCen * wCen * err * err;
        }
        hProj->SetBinContent(iDphi + 1, sum);
        hProj->SetBinError  (iDphi + 1, std::sqrt(sumw2));
    }
    return hProj;
}

// Normalize a 1D ΔΦ histogram by integral then bin width.
// Used only for the final projected wEEC result.
inline bool NormalizeWEEC(TH1D* h)
{
    double integral = h->Integral();
    if (integral <= 0) return false;
    h->Scale(1.0 / integral);
    h->Scale(1.0, "width");
    return true;
}

// ═════════════════════════════════════════════════════════════════════════════
//  Vertex z reweighting
//
//  The vtx_z distribution in data differs from MC.  A per-bin correction is
//  applied when filling the response matrix for data unfolding so that the
//  MC vtx_z shape matches data.
//
//  Binning: 20 bins of 1 cm width covering [-10, 10) cm, matching the vtx_z
//  acceptance cut applied in both fillResponse.C and fillData.C.
//
//  vtxWeights[] is populated by makeVtxWeights.C after it forms the ratio
//  hVtxData / hVtxMC (both normalized to unit area within [-10, 10) cm).
//  Paste the printed array contents from makeVtxWeights.C output here.
//
//  Until makeVtxWeights.C has been run, all weights default to 1.0
//  (no correction applied).
// ═════════════════════════════════════════════════════════════════════════════
constexpr double vtxZMin   = -10.0;
constexpr double vtxZMax   =  10.0;
constexpr int    nVtxBins  =  20;     // 1 cm per bin
constexpr double vtxBinWidth = (vtxZMax - vtxZMin) / nVtxBins;  // 1.0 cm

// ── vtx weight array ─────────────────────────────────────────────────────────
// Bin i covers [vtxZMin + i*vtxBinWidth, vtxZMin + (i+1)*vtxBinWidth).
// Replace this array with the output of makeVtxWeights.C once available.
constexpr double vtxWeights[nVtxBins] = {
    0.962696,  // bin  0: [-10, -9) cm
    0.987321,  // bin  1: [-9, -8) cm
    1.006085,  // bin  2: [-8, -7) cm
    1.020161,  // bin  3: [-7, -6) cm
    1.034171,  // bin  4: [-6, -5) cm
    1.037306,  // bin  5: [-5, -4) cm
    1.046184,  // bin  6: [-4, -3) cm
    1.040758,  // bin  7: [-3, -2) cm
    1.041977,  // bin  8: [-2, -1) cm
    1.042028,  // bin  9: [-1, 0) cm
    1.027416,  // bin 10: [0, 1) cm
    1.030098,  // bin 11: [1, 2) cm
    1.012770,  // bin 12: [2, 3) cm
    1.012976,  // bin 13: [3, 4) cm
    0.996425,  // bin 14: [4, 5) cm
    0.975807,  // bin 15: [5, 6) cm
    0.966945,  // bin 16: [6, 7) cm
    0.944144,  // bin 17: [7, 8) cm
    0.921764,  // bin 18: [8, 9) cm
    0.888391,  // bin 19: [9, 10) cm
};

// Return the vtx reweight factor for a given vtx_z value.
// Returns 1.0 for any vtx_z outside [-10, 10) cm (should not occur after the
// acceptance cut, but guarded for safety).
inline double GetVtxWeight(double vtx_z)
{
    if (vtx_z < vtxZMin || vtx_z >= vtxZMax) return 0.0;
    int bin = static_cast<int>((vtx_z - vtxZMin) / vtxBinWidth);
    if (bin < 0 || bin >= nVtxBins) return 0.0;
    return vtxWeights[bin];
}
