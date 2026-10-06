void draw2Dhistogram(TH2 *hist, std::string xtitle, std::string ytitle, bool profile, std::vector<std::string> addinfo, std::string filename)
{
    TCanvas *c = new TCanvas("c", "c", 800, 700);
    gPad->SetRightMargin(0.14);
    c->cd();
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    hist->GetYaxis()->SetTitle(ytitle.c_str());
    hist->Draw("COLZ");
    if (profile)
    {
        TProfile *prof = hist->ProfileX();
        prof->SetMarkerSize(0.5);
        prof->SetMarkerStyle(20);
        prof->SetMarkerColor(kRed);
        prof->SetLineColor(kRed);
        prof->SetLineWidth(2);
        prof->Draw("PE same");
    }
    // Add additional information to the histogram title using TLatex
    TLatex *latex = new TLatex();
    latex->SetTextSize(0.04);
    latex->SetTextAlign(12);
    latex->SetNDC();
    for (size_t i = 0; i < addinfo.size(); ++i)
    {
        latex->DrawLatex(0.2, 0.88 - i * 0.05, addinfo[i].c_str());
    }
    c->SaveAs(Form("%s.png", filename.c_str()));
    c->SaveAs(Form("%s.pdf", filename.c_str()));
    delete c;
}

void SeedRelEff(std::string infilename = "./test_pp.root")
{
    std::string plotdir = "./figure/figure-SeedRelEff";
    system(Form("mkdir -p %s", plotdir.c_str()));

    TH2 *hM_NRecotkl_Raw_vs_Nsiliconseeds_vtx = new TH2F("hM_NRecotkl_Raw_vs_Nsiliconseeds_vtx", "NRecotkl_Raw vs Nsiliconseeds_vtx;Nsiliconseeds_vtx;NRecotkl_Raw", 50, -0.5, 49.5, 50, -0.5, 49.5);
    TH2 *hM_NRecotkl_Raw_vs_ratio = new TH2F("hM_NRecotkl_Raw_vs_ratio", "NRecotkl_Raw vs Nsiliconseeds_vtx/NRecotkl_Raw;Nsiliconseeds_vtx/NRecotkl_Raw;NRecotkl_Raw", 50, -0.5, 49.5, 20, 0, 2);
    TH2 *hM_NRecotkl_Raw_vs_MBDchargeSum = new TH2F("hM_NRecotkl_Raw_vs_MBDchargeSum", "NRecotkl_Raw vs MBD_charge_sum;MBD_charge_sum;NRecotkl_Raw", 50, -0.5, 49.5, 50, -0.5, 49.5);
    TH2 *hM_Nsiliconseeds_vs_MBDchargeSum = new TH2F("hM_Nsiliconseeds_vs_MBDchargeSum", "Nsiliconseeds_vtx vs MBD_charge_sum;MBD_charge_sum;Nsiliconseeds_vtx", 50, -0.5, 49.5, 50, -0.5, 49.5);

    TFile *f = new TFile(infilename.c_str(), "READ");
    TTree *t = (TTree *)f->Get("minitree");
    int event;
    bool firedTrig12_vtxle10cm;
    int Nsiliconseeds_vtx, NRecotkl_Raw;
    float MBD_charge_sum;
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("firedTrig12_vtxle10cm", &firedTrig12_vtxle10cm);
    t->SetBranchAddress("Nsiliconseeds_vtx", &Nsiliconseeds_vtx);
    t->SetBranchAddress("NRecotkl_Raw", &NRecotkl_Raw);
    t->SetBranchAddress("MBD_charge_sum", &MBD_charge_sum);
    for (int ev = 0; ev < t->GetEntries(); ev++)
    {
        t->GetEntry(ev);
        // cout << "Event " << event << ": Nsiliconseeds_vtx = " << Nsiliconseeds_vtx << ", NRecotkl_Raw = " << NRecotkl_Raw << endl;
        if (!firedTrig12_vtxle10cm)
            continue;

        hM_NRecotkl_Raw_vs_Nsiliconseeds_vtx->Fill(NRecotkl_Raw, Nsiliconseeds_vtx);
        if (NRecotkl_Raw > 0)
        {
            hM_NRecotkl_Raw_vs_ratio->Fill(NRecotkl_Raw, (float)Nsiliconseeds_vtx / NRecotkl_Raw);
        }
        hM_NRecotkl_Raw_vs_MBDchargeSum->Fill(NRecotkl_Raw, MBD_charge_sum);
        hM_Nsiliconseeds_vs_MBDchargeSum->Fill(Nsiliconseeds_vtx, MBD_charge_sum);
    }

    draw2Dhistogram(hM_NRecotkl_Raw_vs_Nsiliconseeds_vtx, "N_{INTT doublets}", "N_{silicon seeds}", true, {}, Form("%s/hM_NRecotkl_Raw_vs_Nsiliconseeds_vtx", plotdir.c_str()));
    draw2Dhistogram(hM_NRecotkl_Raw_vs_ratio, "N_{INTT doublets}", "N_{silicon seeds}/N_{INTT doublets}", true, {}, Form("%s/hM_NRecotkl_Raw_vs_ratio", plotdir.c_str()));
    draw2Dhistogram(hM_NRecotkl_Raw_vs_MBDchargeSum, "N_{INTT doublets}", "MBD charge sum", false, {}, Form("%s/hM_NRecotkl_Raw_vs_MBDchargeSum", plotdir.c_str()));
    draw2Dhistogram(hM_Nsiliconseeds_vs_MBDchargeSum, "N_{silicon seeds}", "MBD charge sum", false, {}, Form("%s/hM_Nsiliconseeds_vs_MBDchargeSum", plotdir.c_str()));
}
