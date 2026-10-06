#include "analysisHelper.h"
#include "RooUnfoldResponse.h"

void draw_3DClosure(const char* outDir   = "/sphenix/tg/tg01/jets/bkimelman/VandyDSTs_wEEC_3D_unfolding_kinematics_Apr26_2026/")
{
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetPaintTextFormat(".2f");
    const std::string label   = ModeLabel(Mode::kFull);

    const std::string plotDir = std::format("{}/Plots", outDir);

    std::string wEECFile  = std::format("{}/wEEC-{}.root",   outDir, label);

    TFile* fWEEC = new TFile(wEECFile.c_str(),"READ");
    if (!fWEEC || fWEEC->IsZombie()) { std::cerr<<"Cannot open "<<wEECFile<<"\n"; return; }

    // Open response file immediately — needed for truth reference throughout
    TFile* fResp = new TFile(std::format("{}/response-all-fullClosure.root",outDir).c_str(),"READ");

    std::vector<TH1D*> unfHist(nDphi, nullptr);
    std::vector<TH1D*> truthHist(nDphi, nullptr);
    std::vector<TH1D*> measHist(nDphi, nullptr);
    std::vector<TH1D*> respTruthHist(nDphi, nullptr);
    std::vector<TH1D*> respMeasHist(nDphi, nullptr);   

    for(int i=0; i<nDphi; i++)
    {
        unfHist[i] = (TH1D*)fWEEC->Get(std::format("hWEEC3D_unfolded_{}",i).c_str());
        truthHist[i] = (TH1D*)fResp->Get(std::format("hWEEC3D_truth_{}",i).c_str());
        measHist[i] = (TH1D*)fResp->Get(std::format("hWEEC3D_meas_{}",i).c_str());
        RooUnfoldResponse *resp = (RooUnfoldResponse*)fResp->Get(std::format("response_wEEC3D_{}",i).c_str());
        respTruthHist[i] = (TH1D*)resp->Htruth();
        respMeasHist[i] = (TH1D*)resp->Hmeasured();
        delete resp;
    }

    TCanvas *c1 = new TCanvas("c1","",8*592,4*592);

    //std::pair<vector<TH1D*>, vector<TH1D*>> combos[3] = {std::make_pair(unfHist, truthHist), std::make_pair(truthHist, respTruthHist), std::make_pair(measHist, respMeasHist)};
    std::string names[3] = {"unfolded_truth", "truth_respTruth", "meas_respMeas"};

    for(int i=0; i<3; i++)
    {
        c1->Clear();
        c1->Divide(8,4);

        vector<TH1D*> rat(nDphi, nullptr);
        vector<TGraph*> gr(nDphi, nullptr);
        vector<TGraph*> grG(nDphi, nullptr);

        for(int j=0; j<nDphi; j++)
        {
            c1->cd(j+1);
            c1->cd(j+1)->SetLogy(0);
            if(i == 0)
            {
                rat[j] = (TH1D*)unfHist[j]->Clone();
                rat[j]->Divide(truthHist[j]);
            }
            else if(i == 1)
            {
                rat[j] = (TH1D*)truthHist[j]->Clone();
                rat[j]->Divide(respTruthHist[j]);
            }
             else if(i == 2)
            {
                rat[j] = (TH1D*)measHist[j]->Clone();
                rat[j]->Divide(respMeasHist[j]);
            }
            rat[j]->SetMarkerStyle(20);


            gr[j] = new TGraph();
            grG[j] = new TGraph();
            double minV = 1.0;
            double maxV = 1.0;

            for(int k=1; k<=rat[j]->GetNbinsX(); k++)
            {
                double bc = rat[j]->GetBinContent(k);
                if(bc == 0.0) continue;

                if(bc < minV) minV = bc;
                if(bc > maxV) maxV = bc;

                if(bc == 1.0)
                {
                    grG[j]->AddPoint(rat[j]->GetXaxis()->GetBinCenter(k), bc);
                }
                else
                {
                    gr[j]->AddPoint(rat[j]->GetXaxis()->GetBinCenter(k), bc);
                }
            }
            gr[j]->SetMarkerStyle(21);
            gr[j]->SetMarkerColor(kRed);

            grG[j]->SetMarkerStyle(22);
            grG[j]->SetMarkerColor(kGreen+2);

            rat[j]->GetYaxis()->SetRangeUser(0.75*minV, 1.25*maxV);
            rat[j]->Draw("P");
            if(gr[j]->GetN()>0) gr[j]->Draw("PSAME");
            if(grG[j]->GetN()>0) grG[j]->Draw("PSAME");
        }

        c1->SaveAs(std::format("{}/ratio_{}.pdf",plotDir,names[i]).c_str());
    }



}