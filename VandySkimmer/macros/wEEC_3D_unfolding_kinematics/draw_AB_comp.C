#include "analysisHelper.h"
#include "RooUnfoldResponse.h"

static const int kColA    = kRed   + 1;
static const int kColB = kBlue  + 1;
static const int kMarA    = 21;
static const int kMarB = 20;

void setMarkers(TH1D *h, bool A)
{
    if(A)
    {
        h->SetMarkerColor(kColA);
        h->SetLineColor(kColA);
        h->SetMarkerStyle(kMarA);
        return;
    }

    h->SetMarkerColor(kColB);
    h->SetLineColor(kColB);
    h->SetMarkerStyle(kMarB);
    return;
}


void draw_AB_comp(const char* outDir = ".")
{
    TFile *respFile = new TFile(std::format("{}/response-all-halfClosure.root",outDir).c_str(),"READ");
    if(!respFile)
    {
        std::cerr << "No half closure file in directory " << outDir << ". Try a different directory" << std::endl;
        return;
    }

    TCanvas *c2 = new TCanvas("c2","",2*957,957);
    gStyle->SetOptStat(0);
    c2->SaveAs(std::format("{}/Plots/AB_comp_all.pdf[",outDir).c_str());

    for(int k=0; k<nDphi; k++)
    {
        std::cout << "working on dPhi bin " << k << std::endl;

        c2->Clear();
        c2->Divide(2,1);

        TCanvas *c1 = new TCanvas("c1","",2*957,957);
        c1->Divide(2,1);
        
        RooUnfoldResponse *resp = (RooUnfoldResponse*)respFile->Get(std::format("response_wEEC3D_{}",k).c_str());
        TH1D *hTruthA = (TH1D*)resp->Htruth();
        setMarkers(hTruthA, true);
        TH1D *hTruthB = (TH1D*)respFile->Get(std::format("hWEEC3D_meas_truth_{}",k).c_str());
        setMarkers(hTruthB, false);
        hTruthB->SetTitle(std::format("Truth AB Comparison #Delta#phi bin {}",k).c_str());

        TH1D *hRecoA = (TH1D*)resp->Hmeasured();
        setMarkers(hRecoA, true);
        TH1D *hRecoB = (TH1D*)respFile->Get(std::format("hWEEC3D_meas_{}",k).c_str());
        setMarkers(hRecoB, false);
        hRecoB->SetTitle(std::format("Reco AB Comparison #Delta#phi bin {}",k).c_str());

        std::cout << "got all hists" << std::endl;

        c1->cd(1);
        c1->cd(1)->SetLogy();

        TRatioPlot *rpTruth = new TRatioPlot(hTruthB, hTruthA, "divsym");
        rpTruth->SetH1DrawOpt("P");
        rpTruth->SetH2DrawOpt("P");

        rpTruth->Draw();
        
        rpTruth->GetLowerPad()->SetLogy();
        rpTruth->GetLowerRefYaxis()->SetTitle("\"Data\" / Matrix");

        rpTruth->GetUpperPad()->cd();
        TLegend *leg = new TLegend(0.65,0.7,0.85,0.85);
        leg->AddEntry(hTruthA,"Matrix Half","P");
        leg->AddEntry(hTruthB,"\"Data\" Half","P");
        leg->Draw();

        c1->cd(2);
        c1->cd(2)->SetLogy();

        TRatioPlot *rpReco = new TRatioPlot(hRecoB, hRecoA, "divsym");
        rpReco->SetH1DrawOpt("P");
        rpReco->SetH2DrawOpt("P");

        rpReco->Draw();

        rpReco->GetLowerPad()->SetLogy();
        rpReco->GetLowerRefYaxis()->SetTitle("\"Data\" / Matrix");
        rpReco->Draw();

        rpReco->GetUpperPad()->cd();
        leg->Draw();

        c1->SaveAs(std::format("{}/Plots/AB_comp_{}.png",outDir,k).c_str());

        c2->cd(1);
        c2->cd(1)->SetLogy();

        rpTruth->Draw();
        rpTruth->GetUpperPad()->cd();
        leg->Draw();

        c2->cd(2);
        c2->cd(2)->SetLogy();

        rpReco->Draw();
        rpReco->GetUpperPad()->cd();
        leg->Draw();

        c2->SaveAs(std::format("{}/Plots/AB_comp_all.pdf",outDir,k).c_str());

        std::cout << "done with dPhi bin " << k << std::endl;
    
        delete c1;
        delete hTruthA; delete hTruthB;
        delete resp; delete hRecoA; delete hRecoB;
        delete rpTruth; //delete rpReco;
        delete leg;
    }
    
    c2->SaveAs(std::format("{}/Plots/AB_comp_all.pdf]",outDir).c_str());

}