#include "HepMC3/ReaderAscii.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Units.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "/sphenix/user/hjheng/sPHENIXRepo/analysis/LightFlavorRatios/util/binning.h"

// Hardcoded MVTX global displacement (x,y,z)= (6.2955, -1.06874, -6.66289) mm, and the minimum radius of MVTX from simulation/g4simulation/g4mvtx/PHG4MvtxDefs.h
const double MVTX_global_displacement_x_cm = 6.2955 / 10.0;  // convert mm to cm
const double MVTX_global_displacement_y_cm = -1.06874 / 10.0;
const double MVTX_global_displacement_z_cm = -6.66289 / 10.0;
const double MVTX_rmin_cm = 24.610 / 10.0; // convert mm to cm, innermost layer

// Return true if this particle is the last copy in a same-PID chain
// This helps avoid double counting intermediate copies in HepMC records
bool is_last_copy(const HepMC3::ConstGenParticlePtr &p)
{
    if (!p)
        return false;

    auto end_vtx = p->end_vertex();
    if (!end_vtx)
        return true;

    for (const auto &d : end_vtx->particles_out())
    {
        if (!d)
            continue;
        if (d->pid() == p->pid())
            return false;
    }
    return true;
}

// Convert HepMC length-unit value to cm
double length_unit_to_cm(const HepMC3::GenEvent &event)
{
    const auto u = event.length_unit();

    if (u == HepMC3::Units::CM)
        return 1.0;
    if (u == HepMC3::Units::MM)
        return 0.1;

    return 1.0;
}

// decay length in cm, using production and decay vertices of the particle, option for transverse or 3D
double decay_length_cm(const HepMC3::ConstGenParticlePtr &p, const HepMC3::GenEvent &event, bool transverse = true)
{
    if (!p)
        return -999.0;

    auto prod_vtx = p->production_vertex();
    auto end_vtx = p->end_vertex();

    if (!prod_vtx || !end_vtx)
        return -999.0;

    const auto &prod = prod_vtx->position();
    const auto &dec = end_vtx->position();

    const double dx = dec.x() - prod.x();
    const double dy = dec.y() - prod.y();
    const double dz = dec.z() - prod.z();

    if (transverse)
    {
        const double L = std::sqrt(dx * dx + dy * dy);
        return L * length_unit_to_cm(event);
    }
    else
    {
        const double L = std::sqrt(dx * dx + dy * dy + dz * dz);
        return L * length_unit_to_cm(event);
    }
}

// Distance from the origin to the displaced inner MVTX boundary along azimuth phi.
double mvtx_inner_radius_at_phi_cm(const double phi)
{
    // convert this phi in HEP definition to [0, 2pi]
    double phi_local = phi;
    if (phi_local < 0)
    {
        phi_local += 2 * M_PI;
    }

    const double cos_phi = std::cos(phi_local);
    const double sin_phi = std::sin(phi_local);

    // calculate the angle for the displacement vector, make it in the range of [0, 2pi] as well
    double disp_phi = std::atan2(MVTX_global_displacement_y_cm, MVTX_global_displacement_x_cm);
    if (disp_phi < 0)
    {
        disp_phi += 2 * M_PI;
    }
    double d = std::sqrt(MVTX_global_displacement_x_cm * MVTX_global_displacement_x_cm + MVTX_global_displacement_y_cm * MVTX_global_displacement_y_cm);
    double rho = d * std::cos(phi_local - disp_phi) + std::sqrt(std::pow(MVTX_rmin_cm, 2) - std::pow(d * std::sin(phi_local - disp_phi), 2));

    return rho;
}

// Check whether this Lambda comes directly from a Xi parent at its production vertex
// 3312: from charged Xi (+/-3312)
// 3322: from neutral Xi (+/-3322)
// 0: otherwise
int lambda_parent_xi_type(const HepMC3::ConstGenParticlePtr &lambda)
{
    if (!lambda)
        return 0;

    auto prod_vtx = lambda->production_vertex();
    if (!prod_vtx)
        return 0;

    for (const auto &parent : prod_vtx->particles_in())
    {
        if (!parent)
            continue;

        const int apid = std::abs(parent->pid());

        if (apid == 3312)
            return 3312; // charged Xi
        if (apid == 3322)
            return 3322; // neutral Xi
    }

    return 0;
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " input.hepmc output.root\n";
        return 1;
    }

    const char *hepmc_file = argv[1];
    const char *root_file = argv[2];

    HepMC3::ReaderAscii reader(hepmc_file);
    HepMC3::GenEvent event;

    TFile *fout = new TFile(root_file, "RECREATE");

    // Tree for selected Lambdas
    TTree *tree = new TTree("lambda_tree", "lambda_tree");

    int event_id = -1;
    int lambda_pid = 0;
    int parent_xi_type = 0; // 0, 3312, 3322
    double lambda_pt = -999.;
    double lambda_eta = -999.;
    double lambda_y = -999.;
    double lambda_phi = -999.;
    double lambda_mass = -999.;
    double lambda_transverse_decay_length_cm = -999.;

    tree->Branch("event", &event_id, "event/I");
    tree->Branch("lambda_pid", &lambda_pid, "lambda_pid/I");
    tree->Branch("parent_xi_type", &parent_xi_type, "parent_xi_type/I");
    tree->Branch("lambda_pt", &lambda_pt, "lambda_pt/D");
    tree->Branch("lambda_eta", &lambda_eta, "lambda_eta/D");
    tree->Branch("lambda_y", &lambda_y, "lambda_y/D");
    tree->Branch("lambda_mass", &lambda_mass, "lambda_mass/D");
    tree->Branch("lambda_transverse_decay_length_cm", &lambda_transverse_decay_length_cm, "lambda_transverse_decay_length_cm/D");

    const auto &lambda_pt_bins = BinInfo::final_pt_bins.bins;
    const auto &lambda_rapidity_bins = BinInfo::final_rapidity_bins.bins;
    const auto &lambda_eta_bins = BinInfo::final_eta_bins.bins;
    const auto &lambda_phi_bins = BinInfo::final_phi_bins.bins;

    // Denominator: all selected Lambdas
    TH1D *h_lambda_pt_all = new TH1D("h_lambda_pt_all", ";#Lambda^{0}(+c.c) p_{T} [GeV];Counts", lambda_pt_bins.size() - 1, lambda_pt_bins.data());
    TH1D *h_lambda_rapidity_all = new TH1D("h_lambda_rapidity_all", ";#Lambda^{0}(+c.c) rapidity;Counts", lambda_rapidity_bins.size() - 1, lambda_rapidity_bins.data());
    TH1D *h_lambda_eta_all = new TH1D("h_lambda_eta_all", ";#Lambda^{0}(+c.c) #eta;Counts", lambda_eta_bins.size() - 1, lambda_eta_bins.data());
    TH1D *h_lambda_phi_all = new TH1D("h_lambda_phi_all", ";#Lambda^{0}(+c.c) #phi;Counts", lambda_phi_bins.size() - 1, lambda_phi_bins.data());
    // Numerators: selected Lambdas from Xi
    TH1D *h_lambda_pt_from_xi_all = new TH1D("h_lambda_pt_from_xi_all", ";#Lambda^{0}(+c.c) p_{T} [GeV];Counts", lambda_pt_bins.size() - 1, lambda_pt_bins.data());
    TH1D *h_lambda_pt_from_xi_charged = new TH1D("h_lambda_pt_from_xi_charged", ";#Lambda^{0}(+c.c) p_{T} [GeV];Counts", lambda_pt_bins.size() - 1, lambda_pt_bins.data());
    TH1D *h_lambda_pt_from_xi_neutral = new TH1D("h_lambda_pt_from_xi_neutral", ";#Lambda^{0}(+c.c) p_{T} [GeV];Counts", lambda_pt_bins.size() - 1, lambda_pt_bins.data());
    TH1D *h_lambda_rapidity_from_xi_all = new TH1D("h_lambda_rapidity_from_xi_all", ";#Lambda^{0}(+c.c) rapidity;Counts", lambda_rapidity_bins.size() - 1, lambda_rapidity_bins.data());
    TH1D *h_lambda_rapidity_from_xi_charged = new TH1D("h_lambda_rapidity_from_xi_charged", ";#Lambda^{0}(+c.c) rapidity;Counts", lambda_rapidity_bins.size() - 1, lambda_rapidity_bins.data());
    TH1D *h_lambda_rapidity_from_xi_neutral = new TH1D("h_lambda_rapidity_from_xi_neutral", ";#Lambda^{0}(+c.c) rapidity;Counts", lambda_rapidity_bins.size() - 1, lambda_rapidity_bins.data());
    TH1D *h_lambda_eta_from_xi_all = new TH1D("h_lambda_eta_from_xi_all", ";#Lambda^{0}(+c.c) #eta;Counts", lambda_eta_bins.size() - 1, lambda_eta_bins.data());
    TH1D *h_lambda_eta_from_xi_charged = new TH1D("h_lambda_eta_from_xi_charged", ";#Lambda^{0}(+c.c) #eta;Counts", lambda_eta_bins.size() - 1, lambda_eta_bins.data());
    TH1D *h_lambda_eta_from_xi_neutral = new TH1D("h_lambda_eta_from_xi_neutral", ";#Lambda^{0}(+c.c) #eta;Counts", lambda_eta_bins.size() - 1, lambda_eta_bins.data());
    TH1D *h_lambda_phi_from_xi_all = new TH1D("h_lambda_phi_from_xi_all", ";#Lambda^{0}(+c.c) #phi;Counts", lambda_phi_bins.size() - 1, lambda_phi_bins.data());
    TH1D *h_lambda_phi_from_xi_charged = new TH1D("h_lambda_phi_from_xi_charged", ";#Lambda^{0}(+c.c) #phi;Counts", lambda_phi_bins.size() - 1, lambda_phi_bins.data());
    TH1D *h_lambda_phi_from_xi_neutral = new TH1D("h_lambda_phi_from_xi_neutral", ";#Lambda^{0}(+c.c) #phi;Counts", lambda_phi_bins.size() - 1, lambda_phi_bins.data());
    // Event counter (for check)
    TH1D *h_event_counter = new TH1D("h_event_counter", ";dummy;events", 1, 0.5, 1.5);
    // Diagnostic histograms
    TH1D *h_lambda_transverse_decay_length_cm = new TH1D("h_lambda_transverse_decay_length_cm", ";#Lambda transverse decay length [cm];Counts", 100, 0, 5);
    TH2D *h_phi_decaylengthcut = new TH2D("h_phi_decaylengthcut", ";#Lambda #phi [rad];Transverse decay length cut [cm]", 128, -3.2, 3.2, 100, 0, 4);

    int iev = 0;

    while (!reader.failed())
    {
        reader.read_event(event);
        if (reader.failed())
            break;

        ++iev;

        if (iev % 5000 == 0)
        {
            std::cout << "Processing event " << iev << std::endl;
        }

        event_id = iev;
        h_event_counter->Fill(1.0);

        for (const auto &p : event.particles())
        {
            if (!p)
                continue;

            // Select Lambda0 / anti-Lambda0
            if (std::abs(p->pid()) != 3122)
                continue;

            // Avoid double counting if same-PID copies exist in a chain
            if (!is_last_copy(p))
                continue;

            // Need both production and decay vertices for the decay-length cut
            if (!p->production_vertex() || !p->end_vertex())
                continue;

            lambda_pid = p->pid();
            lambda_pt = p->momentum().pt();
            lambda_eta = p->momentum().eta();
            lambda_y = p->momentum().rap();
            lambda_phi = p->momentum().phi();
            lambda_mass = p->momentum().m();
            lambda_transverse_decay_length_cm = decay_length_cm(p, event, true); // transverse decay length
            h_lambda_transverse_decay_length_cm->Fill(lambda_transverse_decay_length_cm);

            // Kinematic selection on Lambda
            if (std::abs(lambda_eta) >= 1.3)
                continue;

            // Topological selection on transverse decay length: require the decay to occur before the displaced inner MVTX boundary at this azimuth.
            const double mvtx_transverse_cut_cm = mvtx_inner_radius_at_phi_cm(lambda_phi);
            h_phi_decaylengthcut->Fill(lambda_phi, mvtx_transverse_cut_cm);
            if (mvtx_transverse_cut_cm < 0 ||
                lambda_transverse_decay_length_cm > mvtx_transverse_cut_cm ||
                lambda_transverse_decay_length_cm < 0)
                continue;

            // Denominator
            h_lambda_pt_all->Fill(lambda_pt);
            h_lambda_rapidity_all->Fill(lambda_y);
            h_lambda_eta_all->Fill(lambda_eta);
            h_lambda_phi_all->Fill(lambda_phi);

            // Parent classification
            parent_xi_type = lambda_parent_xi_type(p);

            if (parent_xi_type == 3312)
            {
                h_lambda_pt_from_xi_all->Fill(lambda_pt);
                h_lambda_pt_from_xi_charged->Fill(lambda_pt);
                h_lambda_rapidity_from_xi_all->Fill(lambda_y);
                h_lambda_rapidity_from_xi_charged->Fill(lambda_y);
                h_lambda_eta_from_xi_all->Fill(lambda_eta);
                h_lambda_eta_from_xi_charged->Fill(lambda_eta);
                h_lambda_phi_from_xi_all->Fill(lambda_phi);
                h_lambda_phi_from_xi_charged->Fill(lambda_phi);
            }
            else if (parent_xi_type == 3322)
            {
                h_lambda_pt_from_xi_all->Fill(lambda_pt);
                h_lambda_pt_from_xi_neutral->Fill(lambda_pt);
                h_lambda_rapidity_from_xi_all->Fill(lambda_y);
                h_lambda_rapidity_from_xi_neutral->Fill(lambda_y);
                h_lambda_eta_from_xi_all->Fill(lambda_eta);
                h_lambda_eta_from_xi_neutral->Fill(lambda_eta);
                h_lambda_phi_from_xi_all->Fill(lambda_phi);
                h_lambda_phi_from_xi_neutral->Fill(lambda_phi);
            }

            tree->Fill();
        }

        event.clear();
    }

    // Feeddown fractions
    TH1D *h_feeddown_frac_xi_all = (TH1D *)h_lambda_pt_from_xi_all->Clone("h_feeddown_frac_xi_all");
    h_feeddown_frac_xi_all->SetTitle(";#Lambda^{0}(+c.c) p_{T} [GeV];#Lambda from Xi / all #Lambda");

    TH1D *h_feeddown_frac_xi_charged = (TH1D *)h_lambda_pt_from_xi_charged->Clone("h_feeddown_frac_xi_charged");
    h_feeddown_frac_xi_charged->SetTitle(";#Lambda^{0}(+c.c) p_{T} [GeV];#Lambda from charged Xi / all #Lambda");

    TH1D *h_feeddown_frac_xi_neutral = (TH1D *)h_lambda_pt_from_xi_neutral->Clone("h_feeddown_frac_xi_neutral");
    h_feeddown_frac_xi_neutral->SetTitle(";#Lambda^{0}(+c.c) p_{T} [GeV];#Lambda from neutral Xi / all #Lambda");

    // Binomial errors are appropriate since numerator is a subset of denominator.
    h_feeddown_frac_xi_all->Divide(h_lambda_pt_from_xi_all, h_lambda_pt_all, 1.0, 1.0, "B");
    h_feeddown_frac_xi_charged->Divide(h_lambda_pt_from_xi_charged, h_lambda_pt_all, 1.0, 1.0, "B");
    h_feeddown_frac_xi_neutral->Divide(h_lambda_pt_from_xi_neutral, h_lambda_pt_all, 1.0, 1.0, "B");

    fout->Write();
    fout->Close();

    std::cout << "Processed events = " << iev << "\n";
    std::cout << "Saved ROOT file  = " << root_file << "\n";

    return 0;
}