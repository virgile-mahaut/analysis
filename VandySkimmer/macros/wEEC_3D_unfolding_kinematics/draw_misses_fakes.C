#include "analysisHelper.h"

void draw_misses_fakes(float towCut = 0.25, Mode mode = Mode::kFull)
{

    std::string towCutStr = std::format("{}",towCut);
    std::replace(towCutStr.begin(), towCutStr.end(), '.', 'p');
    std::string outDir = std::format("/sphenix/tg/tg01/jets/bkimelman/VandyDSTs_wEEC_3D_unfolding_kinematics_Apr30_2026_{}/",towCutStr);

    TFile *fResp = new TFile(std::format("{}/response-all-fullClosure.root",outDir).c_str(),"READ");
    if(!fResp)
    {
        //std::cerr << "Response file " << respFile << " does not exit" << std::endl;
        return;
    }

    //std::vector<TH1D*> hTruth(nDphi, nullptr);
    //std::vector<TH1D*> hReco(nDphi, nullptr);
    //std::vector<TH1D*> hMiss(nDphi, nullptr);
    //std::vector<TH1D*> hFake(nDphi, nullptr);

    TH1D *hMisses = new TH1D("misses","Misses;#Delta#phi;Miss Rate",nDphi, &dPhiBins[0]);
    hMisses->SetMarkerStyle(20);
    hMisses->SetMarkerColor(kRed+1);

    TH1D *hFakes = new TH1D("fakes","Fakes;#Delta#phi;Fake Rate",nDphi, &dPhiBins[0]);
    hFakes->SetMarkerStyle(21);
    hFakes->SetMarkerColor(kBlue+1);

    vector<TGraphAsymmErrors*> hTruth_pw(nDphi, nullptr);
    vector<TGraphAsymmErrors*> hReco_pw(nDphi, nullptr);
    vector<TGraphAsymmErrors*> hMisses_pw(nDphi, nullptr);
    vector<TGraphAsymmErrors*> hFakes_pw(nDphi, nullptr);
    for(int k=0; k<nDphi; ++k) {
        hTruth_pw[k] = new TGraphAsymmErrors();
        hTruth_pw[k]->SetMarkerStyle(20);
        hTruth_pw[k]->SetMarkerColor(kBlue+1);
        hTruth_pw[k]->SetLineColor(kBlue+1);

        hReco_pw[k] = new TGraphAsymmErrors();
        hReco_pw[k]->SetMarkerStyle(20);
        hReco_pw[k]->SetMarkerColor(kBlue+1);
        hReco_pw[k]->SetLineColor(kBlue+1);

        hMisses_pw[k] = new TGraphAsymmErrors();
        hMisses_pw[k]->SetMarkerStyle(21);
        hMisses_pw[k]->SetMarkerColor(kRed+1);
        hMisses_pw[k]->SetLineColor(kRed+1);

        hFakes_pw[k] = new TGraphAsymmErrors();
        hFakes_pw[k]->SetMarkerStyle(21);
        hFakes_pw[k]->SetMarkerColor(kRed+1);
        hFakes_pw[k]->SetLineColor(kRed+1);
    }

    vector<double> tMin(nDphi, 1e20);
    vector<double> tMax(nDphi, 0.0);
    vector<double> rMin(nDphi, 1e20);
    vector<double> rMax(nDphi, 0.0);

    for(int k=0; k<nDphi; ++k) {
        TH1D *hTruth = (TH1D*)fResp->Get(std::format("hWEEC3D_truth_{}",k).c_str());
        hTruth->SetDirectory(0);
        double truthIntegral = hTruth->Integral();

        TH1D *hReco = (TH1D*)fResp->Get(std::format("hWEEC3D_meas_{}",k).c_str());
        hReco->SetDirectory(0);
        double recoIntegral = hReco->Integral();

        TH1D *hMiss = (TH1D*)fResp->Get(std::format("hWEEC3D_misses_{}",k).c_str());
        hMiss->SetDirectory(0);
        double missIntegral = hMiss->Integral();

        TH1D *hFake = (TH1D*)fResp->Get(std::format("hWEEC3D_fakes_{}",k).c_str());
        hFake->SetDirectory(0);
        double fakeIntegral = hFake->Integral();

        double missRate = (truthIntegral > 0) ? missIntegral / truthIntegral : 0.0;
        hMisses->SetBinContent(k+1, missRate);

        double fakeRate = (recoIntegral > 0) ? fakeIntegral / recoIntegral : 0.0;
        hFakes->SetBinContent(k+1, fakeRate);

        int iLreco = FindBin(35, recoLeadPtBins);
        int iSreco = FindBin(35, recoSublPtBins);
        int iLtrue = FindBin(35, trueLeadPtBins);
        int iStrue = FindBin(35, trueSublPtBins);

        for(int pw=0; pw<nPairWeight; pw++) {
            int rf = RecoFlat3DIndex(iLreco, iSreco, pw);
            int tf = TrueFlat3DIndex(iLtrue, iStrue, pw);

            double truthPW = hTruth->GetBinContent(tf);
            double recoPW = hReco->GetBinContent(rf);
            double missPW = hMiss->GetBinContent(tf);
            double fakePW = hFake->GetBinContent(rf);

            double pwBinCenter = 0.5*(pairWeightBins[pw] + pairWeightBins[pw+1]);
            double pWELow = pwBinCenter - pairWeightBins[pw];
            double pWEHigh = pairWeightBins[pw+1] - pwBinCenter;

            if(truthPW > 0)
            {
                if(truthPW < tMin[k]) tMin[k] = truthPW;
                if(truthPW > tMax[k]) tMax[k] = truthPW;
                hTruth_pw[k]->AddPoint(pwBinCenter, truthPW);
                hTruth_pw[k]->SetPointError(hTruth_pw[k]->GetN()-1, pWELow, pWEHigh, 0.0, 0.0);
                hMisses_pw[k]->AddPoint(0.5*(pairWeightBins[pw] + pairWeightBins[pw+1]), missPW / truthPW);
                hMisses_pw[k]->SetPointError(hMisses_pw[k]->GetN()-1, pWELow, pWEHigh, 0.0, 0.0);
            }
            if(recoPW > 0)
            {
                if(recoPW < rMin[k]) rMin[k] = recoPW;
                if(recoPW > rMax[k]) rMax[k] = recoPW;
                hReco_pw[k]->AddPoint(pwBinCenter, recoPW);
                hReco_pw[k]->SetPointError(hReco_pw[k]->GetN()-1, pWELow, pWEHigh, 0.0, 0.0);
                hFakes_pw[k]->AddPoint(pwBinCenter, fakePW / recoPW);
                hFakes_pw[k]->SetPointError(hFakes_pw[k]->GetN()-1, pWELow, pWEHigh, 0.0, 0.0);
            }
        }

        delete hTruth; delete hReco;
        delete hMiss; delete hFake;
    }

    TCanvas *c1 = new TCanvas();
    gStyle->SetOptStat(0);

    TH2D *tmp = new TH2D("tmp","Fake and Miss Rates;#Delta#phi;Rate",1,0,TMath::Pi(),1,0,1);
    tmp->Draw();

    TLine *l25 = new TLine(0,0.25,TMath::Pi(),0.25);
    l25->SetLineColor(kBlack);
    l25->SetLineStyle(2);
    l25->Draw("same");

    TLine *l5 = new TLine(0,0.5,TMath::Pi(),0.5);
    l5->SetLineColor(kBlack);
    l5->SetLineStyle(2);
    l5->Draw("same");

    TLine *l75 = new TLine(0,0.75,TMath::Pi(),0.75);
    l75->SetLineColor(kBlack);
    l75->SetLineStyle(2);
    l75->Draw("same");

    hMisses->Draw("PSAME");
    hFakes->Draw("PSAME");

    TLegend *leg = new TLegend(0.7,0.7,0.85,0.85);
    leg->AddEntry(hMisses,"Misses","P");
    leg->AddEntry(hFakes,"Fakes","P");
    leg->Draw();
    
    c1->SaveAs(std::format("{}/Plots/misses_and_fakes-{}-{}.png",outDir,ModeLabel(mode),towCutStr).c_str());

    TCanvas *c2 = new TCanvas("c2","",2*957,957);
    //TH2D *tmp2 = new TH2D("tmp2","Fake and Miss Rates;#Delta#phi;Rate",1,1e-6,2.0,1,0,1);
    c2->SaveAs(std::format("{}/Plots/misses_and_fakes_dPhi-{}-{}.pdf[",outDir,ModeLabel(mode),towCutStr).c_str());

    for(int k=0; k<nDphi; ++k) {
        c2->Clear();
        c2->Divide(2,1);
        //c2->cd(1)->SetLogx(); 
        //c2->cd(1)->SetLogy();

        c2->cd(1);
        auto *p1L = new TPad("p1L","",0,0,1,1);
        p1L->SetLogy();
        p1L->SetLogx();
        p1L->Draw();
        p1L->cd();

        TH1F *hlT = p1L->DrawFrame(1e-6,0.1*tMin[k],2.0,10*tMax[k]);
        hlT->SetXTitle("#Delta#phi");
        hlT->SetYTitle("x-sec counts");
        hlT->SetTitle(std::format("Truth Distribution and Miss Rates p_{{T,tower}}>{} GeV. #Delta#phi bin {} {:.2f}-{:.2f} 30<p_{{T,lead}}<40, 30<p_{{T,sub}}<40",towCut,k,dPhiBins[k],dPhiBins[k+1]).c_str());
        hlT->GetYaxis()->SetAxisColor(kBlue+1);
        hlT->GetYaxis()->SetLabelColor(kBlue+1);
        hlT->GetYaxis()->SetTitleColor(kBlue+1);
        hlT->Draw();

        hTruth_pw[k]->Draw("PSAME");

        c2->cd(1);
        TPad *p1R = new TPad("p1R","",0,0,1,1);
        p1R->SetFillStyle(4000);
        p1R->SetFillColor(0);
        p1R->SetFrameFillStyle(4000);
        p1R->SetFrameLineColor(0);
        p1R->SetFrameBorderMode(0);
        p1R->SetLogx();
        p1R->Draw();
        p1R->cd();

        TH1F *hrT = p1R->DrawFrame(1e-6,0,2.0,1.0);
        hrT->GetXaxis()->SetLabelOffset(99);
        hrT->GetYaxis()->SetLabelOffset(99);
        hrT->GetXaxis()->SetAxisColor(0);
        hrT->GetXaxis()->SetTickLength(0);
        hrT->GetYaxis()->SetAxisColor(0);
        hrT->GetYaxis()->SetTickLength(0);        

        TGaxis *abT = new TGaxis(1e-6,0,2.0,0,1e-6,2,510,"G");
        abT->SetTickLength(0);
        abT->SetLabelOffset(99);
        abT->Draw();

        TGaxis *atT = new TGaxis(1e-6,1,2.0,1,1e-6,2,0,"G");
        atT->SetTickLength(0);
        atT->SetLabelOffset(99);
        atT->Draw();

        TGaxis *alT = new TGaxis(1e-6,0,1e-6,1,0.1*tMin[k],10*tMax[k],510,"G");
        alT->SetTickLength(0);
        alT->SetLabelOffset(99);
        alT->SetLineColor(kBlue+1);
        alT->Draw();

        TGaxis *arT = new TGaxis(2.0,0,2.0,1.0,0,1,510,"+L");
        arT->SetLineColor(kRed+1);
        arT->SetLabelColor(kRed+1);
        arT->SetTitleColor(kRed+1);
        arT->SetTitle("Miss Rate");
        arT->Draw();

        hMisses_pw[k]->Draw("PSAME");


        //
        c2->cd(2);
        auto *p2L = new TPad("p2L","",0,0,1,1);
        p2L->SetLogy();
        p2L->SetLogx();
        p2L->Draw();
        p2L->cd();

        TH1F *hlR = p2L->DrawFrame(1e-6,0.1*rMin[k],2.0,10*rMax[k]);
        hlR->SetXTitle("#Delta#phi");
        hlR->SetYTitle("x-sec counts");
        hlR->SetTitle(std::format("Reco Distribution and Fake Rates p_{{T,tower}}>{} GeV. #Delta#phi bin {} {:.2f}-{:.2f} 30<p_{{T,lead}}<40, 30<p_{{T,sub}}<40",towCut,k,dPhiBins[k],dPhiBins[k+1]).c_str());
        hlR->GetYaxis()->SetAxisColor(kBlue+1);
        hlR->GetYaxis()->SetLabelColor(kBlue+1);
        hlR->GetYaxis()->SetTitleColor(kBlue+1);
        hlR->Draw();

        hReco_pw[k]->Draw("PSAME");

        c2->cd(2);
        TPad *p2R = new TPad("p2R","",0,0,1,1);
        p2R->SetFillStyle(4000);
        p2R->SetFillColor(0);
        p2R->SetFrameFillStyle(4000);
        p2R->SetFrameLineColor(0);
        p2R->SetFrameBorderMode(0);        
        p2R->SetLogx();
        p2R->Draw();
        p2R->cd();

        TH1F *hrR = p2R->DrawFrame(1e-6,0,2.0,1.0);
        hrR->GetXaxis()->SetLabelOffset(99);
        hrR->GetYaxis()->SetLabelOffset(99);
        hrR->GetXaxis()->SetAxisColor(0);
        hrR->GetXaxis()->SetTickLength(0);
        hrR->GetYaxis()->SetAxisColor(0);
        hrR->GetYaxis()->SetTickLength(0);    

        TGaxis *abR = new TGaxis(1e-6,0,2.0,0,1e-6,2,510,"G");
        abR->SetTickLength(0);
        abR->SetLabelOffset(99);
        abR->Draw();

        TGaxis *atR = new TGaxis(1e-6,1,2.0,1,1e-6,2,0,"G");
        atR->SetTickLength(0);
        atR->SetLabelOffset(99);
        atR->Draw();

        TGaxis *alR = new TGaxis(1e-6,0,1e-6,1,0.1*rMin[k],10*rMax[k],510,"G");
        alR->SetTickLength(0);
        alR->SetLabelOffset(99);
        alR->SetLineColor(kBlue+1);
        alR->Draw();

        TGaxis *arR = new TGaxis(2.0,0,2.0,1.0,0,1,510,"+L");
        arR->SetLineColor(kRed+1);
        arR->SetLabelColor(kRed+1);
        arR->SetTitleColor(kRed+1);
        arR->SetTitle("Fake Rate");
        arR->Draw();


        hFakes_pw[k]->Draw("PSAME");

        c2->SaveAs(std::format("{}/Plots/misses_and_fakes_dPhi-{}-{}.pdf",outDir,ModeLabel(mode),towCutStr).c_str());

    }
    c2->SaveAs(std::format("{}/Plots/misses_and_fakes_dPhi-{}-{}.pdf]",outDir,ModeLabel(mode),towCutStr).c_str());

    

}