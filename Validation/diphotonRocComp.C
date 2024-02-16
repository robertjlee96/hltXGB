#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCut.h"
#include "TPaletteAxis.h"
#include "TMVA/Reader.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <algorithm>

void diphotonRocComp(){
    gROOT->Reset();
    gStyle->SetPalette(1);
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitle(0);
    //gROOT->ProcessLine("gErrorIgnoreLevel = 6001;");
    
    int nBinsPlot = 502; // From 0 to 1 in xgbScore with one overflow and one underflow
    int binFactor = 50;// Increase binning by factor of X. Get binFactor extra bins above and below plot
    
    
    //Put these in by hand (maybe get auto later?) to place point at specific signal/background value
    //double sigEffHLT = 56774.0/66456.0;
    //double bkgEffHLT = 334.0/24486.0;
    
    double sigEffHLT = 53814.0/61078.0;
    double bkgEffHLT = 334.0/24486.0;
    
    
    string fOut[5];
    string outStrGen = "outROC/0517/newGGH_TrainingMassCutComp_";
    fOut[0] = outStrGen + "BB";
    fOut[1] = outStrGen + "BE";
    fOut[2] = outStrGen + "EE";
    fOut[3] = outStrGen + "All";
    
    string fOutSig[5];
    fOutSig[0] = fOut[0] + "_SigScore";
    fOutSig[1] = fOut[1] + "_SigScore";
    fOutSig[2] = fOut[2] + "_SigScore";
    fOutSig[3] = fOut[3] + "_SigScore";
    
    string fOutBkg[5];
    fOutBkg[0] = fOut[0] + "_BkgScore";
    fOutBkg[1] = fOut[1] + "_BkgScore";
    fOutBkg[2] = fOut[2] + "_BkgScore";
    fOutBkg[3] = fOut[3] + "_BkgScore";
    bool zoom = true;

    string eta[4] = {"Barrel-Barrel","Barrel-Endcap","Endcap-Endcap","All"};
    //string fOutF = "idMVAPlots/GJetIDMVAF_LRSearch18_Endca_0127.rooot";
        
    int nBins = nBinsPlot*binFactor;
    double binSize = 1.0/(nBinsPlot-2); // Plotting bin size
    double plotLow = 0.0 - binSize; //Min value to have 1 bin of size binsize below 0
    double plotHigh = 1.0 + binSize; // Same for max
    double binSizeFine = (plotHigh-plotLow)/nBins;
    
    string fileName1;
    string f1Label = "MD09 LR0.1 M > 60";
    fileName1 = "validationNTuples/0516/newGGHandData_M60Model_MD9LR10.root";
    string tNameSig1 = "sigTree";
    string tNameBkg1 = "bkgTree";
    
    string fileName2;
    string f2Label = "MD09 LR0.1 M > 90";
    fileName2 = "validationNTuples/0516/newGGHandData_M90Model_MD9LR10.root";
    string tNameSig2 = "sigTree";
    string tNameBkg2 = "bkgTree";
    
    int iLow = 0;
    int iHigh = 4;
    
    TH1F *hSig1[4];
    TH1F *hBkg1[4];
    //TH1F *hSigPresel1[4];
    //TH1F *hBkgPresel1[4];
    TGraph *sigEff_vs_bkgEff1[4];
    TGraph *hltPointPlot[4];
    TGraph *sigEff_vs_bkgEffPresel1[4];
    
    TH1F *hSig2[4];
    TH1F *hBkg2[4];
    //TH1F *hSigPresel2[4];
    //TH1F *hBkgPresel2[4];
    TGraph *sigEff_vs_bkgEff2[4];
    //TGraph *sigEff_vs_bkgEffPresel2[4];
    
    for(int i = iLow; i < iHigh; i++){
        string titlePrompt1 = "hSig1" + to_string(i);
        string titleFake1 = "hBkg1" + to_string(i);
        hSig1[i] = new TH1F (titlePrompt1.c_str(),"",nBins,plotLow,plotHigh);
        hBkg1[i] = new TH1F (titleFake1.c_str(),"",nBins,plotLow,plotHigh);
        
        //string titlePromptPresel1 = "hSigPresel1" + to_string(i);
        //string titleFakePresel1 = "hBkgPresel1" + to_string(i);
        //hSigPresel1[i] = new TH1F (titlePromptPresel1.c_str(),"",nBins,plotLow,plotHigh);
        //hBkgPresel1[i] = new TH1F (titleFakePresel1.c_str(),"",nBins,plotLow,plotHigh);
        
        string titlePrompt2 = "hSig2" + to_string(i);
        string titleFake2 = "hBkg2" + to_string(i);
        hSig2[i] = new TH1F (titlePrompt2.c_str(),"",nBins,plotLow,plotHigh);
        hBkg2[i] = new TH1F (titleFake2.c_str(),"",nBins,plotLow,plotHigh);

        //string titlePromptPresel2 = "hSigPresel2" + to_string(i);
        //string titleFakePresel2 = "hBkgPresel2" + to_string(i);
        //hSigPresel2[i] = new TH1F (titlePromptPresel2.c_str(),"",nBins,plotLow,plotHigh);
        //hBkgPresel2[i] = new TH1F (titleFakePresel2.c_str(),"",nBins,plotLow,plotHigh);
    }
    
   //FOR FIRST ROC CURVE
    TFile *f1 = new TFile(fileName1.c_str());
    
    TTreeReader sigReader1(tNameSig1.c_str(), f1);
    TTreeReaderArray<Float_t> sigScores1(sigReader1, "xgbScore");
    TTreeReaderArray<Float_t> sigEta1(sigReader1, "eta");
    TTreeReaderArray<Float_t> sigMass1(sigReader1, "mass");
    
    TTreeReader bkgReader1(tNameBkg1.c_str(), f1);
    TTreeReaderArray<Float_t> bkgScores1(bkgReader1, "xgbScore");
    TTreeReaderArray<Float_t> bkgEta1(bkgReader1, "eta");
    TTreeReaderArray<Float_t> bkgMass1(bkgReader1, "mass");
    

    while (sigReader1.Next()) {
        if(abs(sigEta1[0]) < 1.4442 && abs(sigEta1[1]) < 1.4442 && sigMass1[0] > 95){
            //if(*sigLeadScores1 > *sigSubScores1)hSigPresel1[0]->Fill(*sigSubScores1,*sigSubWeights1);
            //if(*sigLeadScores1 < *sigSubScores1)hSigPresel1[0]->Fill(*sigLeadScores1,*sigLeadWeights1);
            if(sigScores1[0] > sigScores1[1])hSig1[0]->Fill(sigScores1[1]);
            if(sigScores1[0] < sigScores1[1])hSig1[0]->Fill(sigScores1[0]);
        }
        if(((abs(sigEta1[0]) < 1.4442 && abs(sigEta1[1]) > 1.556) || (abs(sigEta1[0]) > 1.556 && abs(sigEta1[1]) < 1.4442)) && sigMass1[0] > 95.0){
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 > *sigSubScores1)hSigPresel1[1]->Fill(*sigSubScores1,*sigSubWeights1);
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 < *sigSubScores1)hSigPresel1[1]->Fill(*sigLeadScores1,*sigLeadWeights1);
            if(sigScores1[0] > sigScores1[1])hSig1[1]->Fill(sigScores1[1]);
            if(sigScores1[0] < sigScores1[1])hSig1[1]->Fill(sigScores1[0]);
        }
        if(abs(sigEta1[0]) > 1.556 && abs(sigEta1[1]) > 1.556 && sigMass1[0] > 95){
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 > *sigSubScores1)hSigPresel1[2]->Fill(*sigSubScores1,*sigSubWeights1);
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 < *sigSubScores1)hSigPresel1[2]->Fill(*sigLeadScores1,*sigLeadWeights1);
            if(sigScores1[0] > sigScores1[1])hSig1[2]->Fill(sigScores1[1]);
            if(sigScores1[0] < sigScores1[1])hSig1[2]->Fill(sigScores1[0]);
        }
    }
    
    while (bkgReader1.Next()) {
        if(abs(bkgEta1[0]) < 1.4442 && abs(bkgEta1[1]) < 1.4442 && bkgMass1[0] > 95){
            //if(*sigLeadScores1 > *sigSubScores1)hSigPresel1[0]->Fill(*sigSubScores1,*sigSubWeights1);
            //if(*sigLeadScores1 < *sigSubScores1)hSigPresel1[0]->Fill(*sigLeadScores1,*sigLeadWeights1);
            if(bkgScores1[0] > bkgScores1[1])hBkg1[0]->Fill(bkgScores1[1]);
            if(bkgScores1[0] < bkgScores1[1])hBkg1[0]->Fill(bkgScores1[0]);
        }
        if(((abs(bkgEta1[0]) < 1.4442 && abs(bkgEta1[1]) > 1.556) || (abs(bkgEta1[0]) > 1.556 && abs(bkgEta1[1]) < 1.4442)) && bkgMass1[0] > 95.0){
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 > *sigSubScores1)hSigPresel1[1]->Fill(*sigSubScores1,*sigSubWeights1);
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 < *sigSubScores1)hSigPresel1[1]->Fill(*sigLeadScores1,*sigLeadWeights1);
            if(bkgScores1[0] > bkgScores1[1])hBkg1[1]->Fill(bkgScores1[1]);
            if(bkgScores1[0] < bkgScores1[1])hBkg1[1]->Fill(bkgScores1[0]);
        }
        if(abs(bkgEta1[0]) > 1.556 && abs(bkgEta1[1]) > 1.556 && bkgMass1[0] > 95){
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 > *sigSubScores1)hSigPresel1[2]->Fill(*sigSubScores1,*sigSubWeights1);
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 < *sigSubScores1)hSigPresel1[2]->Fill(*sigLeadScores1,*sigLeadWeights1);
            if(bkgScores1[0] > bkgScores1[1])hBkg1[2]->Fill(bkgScores1[1]);
            if(bkgScores1[0] < bkgScores1[1])hBkg1[2]->Fill(bkgScores1[0]);
        }
    }

    //FOR SECOND ROC CURVE
    TFile *f2 = new TFile(fileName2.c_str());
    
    TTreeReader sigReader2(tNameSig2.c_str(), f2);
    TTreeReaderArray<Float_t> sigScores2(sigReader2, "xgbScore");
    TTreeReaderArray<Float_t> sigEta2(sigReader2, "eta");
    TTreeReaderArray<Float_t> sigMass2(sigReader2, "mass");
    
    TTreeReader bkgReader2(tNameBkg2.c_str(), f2);
    TTreeReaderArray<Float_t> bkgScores2(bkgReader2, "xgbScore");
    TTreeReaderArray<Float_t> bkgEta2(bkgReader2, "eta");
    TTreeReaderArray<Float_t> bkgMass2(bkgReader2, "mass");
    
    
    while (sigReader2.Next()) {
        if(abs(sigEta2[0]) < 1.4442 && abs(sigEta2[1]) < 1.4442 && sigMass2[0] > 60){
            //if(*sigLeadScores1 > *sigSubScores1)hSigPresel1[0]->Fill(*sigSubScores1,*sigSubWeights1);
            //if(*sigLeadScores1 < *sigSubScores1)hSigPresel1[0]->Fill(*sigLeadScores1,*sigLeadWeights1);
            if(sigScores2[0] > sigScores2[1])hSig2[0]->Fill(sigScores2[1]);
            if(sigScores2[0] < sigScores2[1])hSig2[0]->Fill(sigScores2[0]);
        }
        if(((abs(sigEta2[0]) < 1.4442 && abs(sigEta2[1]) > 1.556) || (abs(sigEta2[0]) > 1.556 && abs(sigEta2[1]) < 1.4442)) && sigMass2[0] > 95.0){
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 > *sigSubScores1)hSigPresel1[1]->Fill(*sigSubScores1,*sigSubWeights1);
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 < *sigSubScores1)hSigPresel1[1]->Fill(*sigLeadScores1,*sigLeadWeights1);
            if(sigScores2[0] > sigScores2[1])hSig2[1]->Fill(sigScores2[1]);
            if(sigScores2[0] < sigScores2[1])hSig2[1]->Fill(sigScores2[0]);
        }
        if(abs(sigEta2[0]) > 1.556 && abs(sigEta2[1]) > 1.556 && sigMass2[0] > 60){
            //if(sigLeadVars2[21] == 1 && *sigLeadScores1 > *sigSubScores1)hSigPresel1[2]->Fill(*sigSubScores1,*sigSubWeights1);
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 < *sigSubScores1)hSigPresel1[2]->Fill(*sigLeadScores1,*sigLeadWeights1);
            if(sigScores2[0] > sigScores2[1])hSig2[2]->Fill(sigScores2[1]);
            if(sigScores2[0] < sigScores2[1])hSig2[2]->Fill(sigScores2[0]);
        }
    }
    
    while (bkgReader2.Next()) {
        if(abs(bkgEta2[0]) < 1.4442 && abs(bkgEta2[1]) < 1.4442 && bkgMass2[0] > 60){
            //if(*sigLeadScores1 > *sigSubScores1)hSigPresel1[0]->Fill(*sigSubScores1,*sigSubWeights1);
            //if(*sigLeadScores1 < *sigSubScores1)hSigPresel1[0]->Fill(*sigLeadScores1,*sigLeadWeights1);
            if(bkgScores2[0] > bkgScores2[1])hBkg2[0]->Fill(bkgScores2[1]);
            if(bkgScores2[0] < bkgScores2[1])hBkg2[0]->Fill(bkgScores2[0]);
        }
        if(((abs(bkgEta2[0]) < 1.4442 && abs(bkgEta2[1]) > 1.556) || (abs(bkgEta2[0]) > 1.556 && abs(bkgEta2[1]) < 1.4442)) && bkgMass2[0] > 95.0){
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 > *sigSubScores1)hSigPresel1[1]->Fill(*sigSubScores1,*sigSubWeights1);
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 < *sigSubScores1)hSigPresel1[1]->Fill(*sigLeadScores1,*sigLeadWeights1);
            if(bkgScores2[0] > bkgScores2[1])hBkg2[1]->Fill(bkgScores2[1]);
            if(bkgScores2[0] < bkgScores2[1])hBkg2[1]->Fill(bkgScores2[0]);
        }
        if(abs(bkgEta2[0]) > 1.556 && abs(bkgEta2[1]) > 1.556 && bkgMass2[0] > 60){
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 > *sigSubScores1)hSigPresel1[2]->Fill(*sigSubScores1,*sigSubWeights1);
            //if(sigLeadVars1[21] == 1 && *sigLeadScores1 < *sigSubScores1)hSigPresel1[2]->Fill(*sigLeadScores1,*sigLeadWeights1);
            if(bkgScores2[0] > bkgScores2[1])hBkg2[2]->Fill(bkgScores2[1]);
            if(bkgScores2[0] < bkgScores2[1])hBkg2[2]->Fill(bkgScores2[0]);
        }
    }
    
    for(int i = iLow; i < iHigh; i++){


        double maxVal = 0.0;
        int maxInt = -1;

        string title = "RoC Curves with Data Bkg " + eta[i];
        string canTitle = "canR" + to_string(i);
        
        string titlePromptCan = "xgbScores for Signal " + eta[i];
        string canTitlePrompt = "canS" + to_string(i);
        string titleFakeCan = "xgbScores for Background " + eta[i];
        string canTitleFake = "canB" + to_string(i);

        //TCanvas *canF = new TCanvas ("canF","canF",10,10,1600,900);
        //canF->SetGrid();
        TCanvas *can = new TCanvas (canTitle.c_str(),canTitle.c_str(),10,10,1600,900);
        can->SetGrid();
        
//        TCanvas *canS = new TCanvas (canTitlePrompt.c_str(),canTitlePrompt.c_str(),10,10,1600,900);
//        canS->SetGrid();
//        TCanvas *canB = new TCanvas (canTitleFake.c_str(),canTitleFake.c_str(),10,10,1600,900);
//        canB->SetGrid();

        // TLegend *legendF = new TLegend(0.35,0.6,0.65,0.9,"","brNDC");
//        TLegend *legend = new TLegend(0.25,0.1,0.55,0.4,"","brNDC");
        TLegend *legend = new TLegend(0.35,0.1,0.65,0.4,"","brNDC");
        
        TLegend *legendS = new TLegend(0.35,0.1,0.65,0.4,"","brNDC");
        TLegend *legendB = new TLegend(0.35,0.1,0.65,0.4,"","brNDC");

        if (i == 3){
            for (int p = 0; p < 3; p++){
                hSig1[i]->Add(hSig1[p]);
                //hSigPresel1[i]->Add(hSigPresel1[p]);

                hBkg1[i]->Add(hBkg1[p]);
                //hBkgPresel1[i]->Add(hBkgPresel1[p]);

                hSig2[i]->Add(hSig2[p]);
                //hSigPresel2[i]->Add(hSigPresel2[p]);

                hBkg2[i]->Add(hBkg2[p]);
                //hBkgPresel2[i]->Add(hBkgPresel2[p]);

            }
        }

        double nEventsP1 = hSig1[i]->GetSumOfWeights();
        double nEventsF1 = hBkg1[i]->GetSumOfWeights();

        double nEventsP2 = hSig2[i]->GetSumOfWeights();
        double nEventsF2 = hBkg2[i]->GetSumOfWeights();

        int nCuts = nBins;
        double stepVal = binSizeFine;
        double cutVal = 0.0 - stepVal;
    //
        float Nsig1[nCuts+1], Nbkg1[nCuts+1];
        float sigEff1[nCuts+1], bkgEff1[nCuts+1], bkgEff1Scaled[nCuts+1];
        //float NsigPresel1[nCuts+1], NbkgPresel1[nCuts+1];
        //float sigEffPresel1[nCuts+1], bkgEffPresel1[nCuts+1];

        float Nsig2[nCuts+1], Nbkg2[nCuts+1];
        float sigEff2[nCuts+1], bkgEff2[nCuts+1], bkgEff2Scaled[nCuts+1];
        //float NsigPresel2[nCuts+1], NbkgPresel2[nCuts+1];
        //float sigEffPresel2[nCuts+1], bkgEffPresel2[nCuts+1];

        float cutsVal[nCuts+1];
        float mvaResCutVal = 0.0 - stepVal;

        int mvaSMaxBin = hSig1[i]->GetXaxis()->FindBin(1.0) + 1;
        int mvaBMaxBin = hBkg1[i]->GetXaxis()->FindBin(1.0) + 1;

        for(int k = 0; k < nCuts+1; k++){
            mvaResCutVal+= stepVal;
            cutsVal[k] = mvaResCutVal;

            int mvaBinSig = hSig1[i]->GetXaxis()->FindBin(mvaResCutVal);
            Nsig1[k] = hSig1[i]->Integral(mvaBinSig,mvaSMaxBin);
            //NsigPresel1[k] = hSigPresel1[i]->Integral(mvaBinSig,mvaSMaxBin);

            Nsig2[k] = hSig2[i]->Integral(mvaBinSig,mvaSMaxBin);
            //NsigPresel2[k] = hSigPresel2[i]->Integral(mvaBinSig,mvaSMaxBin);

            int mvaBinBkg = hBkg1[i]->GetXaxis()->FindBin(mvaResCutVal);
            Nbkg1[k] = hBkg1[i]->Integral(mvaBinBkg,mvaBMaxBin);
            //NbkgPresel1[k] = hBkgPresel1[i]->Integral(mvaBinBkg,mvaBMaxBin);

            Nbkg2[k] = hBkg2[i]->Integral(mvaBinBkg,mvaBMaxBin);
            //NbkgPresel2[k] = hBkgPresel2[i]->Integral(mvaBinBkg,mvaBMaxBin);

//            sigEff1[k] = Nsig1[k]/NsigPresel1[0];
            //bkgEff1Scaled[k] = Nbkg1[k]/NbkgPresel1[0];
            sigEff1[k] = Nsig1[k]/Nsig1[0];
            bkgEff1[k] = Nbkg1[k]/Nbkg1[0];
//
            //sigEffPresel1[k] = NsigPresel1[k]/NsigPresel1[0];
            //bkgEffPresel1[k] = NbkgPresel1[k]/NbkgPresel1[0];

            //sigEff2[k] = Nsig2[k]/NsigPresel2[0];
            //bkgEff2Scaled[k] = Nbkg2[k]/NbkgPresel2[0];
            sigEff2[k] = Nsig2[k]/Nsig2[0];
            bkgEff2[k] = Nbkg2[k]/Nbkg2[0];
            
            //sigEffPresel2[k] = NsigPresel2[k]/NsigPresel2[0];
            //bkgEffPresel2[k] = NbkgPresel2[k]/NbkgPresel2[0];
        }

        
        sigEff_vs_bkgEff1[i] = new TGraph (nCuts, bkgEff1, sigEff1);
        //sigEff_vs_bkgEffPresel1[i] = new TGraph (nCuts, bkgEffPresel1, sigEffPresel1);
        //sigEff_vs_bkgEffPresel1[i]->SetLineStyle(9);

        sigEff_vs_bkgEff2[i] = new TGraph (nCuts, bkgEff2, sigEff2);
        //sigEff_vs_bkgEffPresel2[i] = new TGraph (nCuts, bkgEffPresel2, sigEffPresel2);
        sigEff_vs_bkgEff2[i]->SetLineColor(2);
        sigEff_vs_bkgEff2[i]->SetMarkerColor(2);
        //sigEff_vs_bkgEffPresel2[i]->SetLineColor(2);
        //sigEff_vs_bkgEffPresel2[i]->SetLineStyle(9);
//
//        if (sigEff1[0] >= sigEff2[0]){
//            sigEff_vs_bkgEff1[i]->SetTitle(title.c_str());
//            sigEff_vs_bkgEff1[i]->GetXaxis()->SetTitle("Background Efficiency");
//            sigEff_vs_bkgEff1[i]->GetYaxis()->SetTitle("Signal Efficiency");
//            if (zoom == true)sigEff_vs_bkgEff1[i]->GetXaxis()->SetRangeUser(0.0,1.5);
//        }
//        if (sigEff2[0] > sigEff1[0]){
//            sigEff_vs_bkgEff2[i]->SetTitle(title.c_str());
//            sigEff_vs_bkgEff2[i]->GetXaxis()->SetTitle("Background Efficiency");
//            sigEff_vs_bkgEff2[i]->GetYaxis()->SetTitle("Signal Efficiency");
//            if (zoom == true)sigEff_vs_bkgEff2[i]->GetXaxis()->SetRangeUser(0.0,1.5);
//        }
        sigEff_vs_bkgEff1[i]->SetTitle(title.c_str());
        sigEff_vs_bkgEff1[i]->GetXaxis()->SetTitle("Background Efficiency");
        sigEff_vs_bkgEff1[i]->GetYaxis()->SetTitle("Signal Efficiency");
////        if (zoom == true)sigEff_vs_bkgEff1[i]->GetXaxis()->SetRangeUser(0.0,1.5);

        TF1 f1("f1",[&](double *bkgEff1, double *sigEff1){ return sigEff_vs_bkgEff1[i]->Eval(bkgEff1[0]); },0,1,0); //Gets values of TGraph
        double AUC1 = f1.Integral(0,1);
        stringstream stream1;
        stream1 << fixed << setprecision(5) << AUC1;
        string AUCStr1 = f1Label + " All, AUC = " + stream1.str();
        cout<<AUCStr1.c_str()<<endl;

//        TF1 f2("f2",[&](double *bkgEffPresel1, double *sigEffPresel1){ return sigEff_vs_bkgEffPresel1[i]->Eval(bkgEffPresel1[0]); },0,1,0); //Gets values of TGraph
//        double AUCPresel1 = f2.Integral(0,1);
//        stringstream streamPresel1;
//        streamPresel1 << fixed << setprecision(5) << AUCPresel1;
//        string AUCStrPresel1 = f1Label + " W/ Presel, AUC = " + streamPresel1.str();
        //cout<<AUCStrPresel1<<endl;

        TF1 f3("f3",[&](double *bkgEff2, double *sigEff2){ return sigEff_vs_bkgEff2[i]->Eval(bkgEff2[0]); },0,1,0); //Gets values of TGraph
        double AUC2 = f3.Integral(0,1);
        stringstream stream2;
        stream2 << fixed << setprecision(5) << AUC2;
        string AUCStr2 = f2Label + " All, AUC = " + stream2.str();
        cout<<AUCStr2.c_str()<<endl;
        
//        TF1 f4("f4",[&](double *bkgEffPresel2, double *sigEffPresel2){ return sigEff_vs_bkgEffPresel2[i]->Eval(bkgEffPresel2[0]); },0,1,0); //Gets values of TGraph
//        double AUCPresel2 = f4.Integral(0,1);
//        stringstream streamPresel2;
//        streamPresel2 << fixed << setprecision(5) << AUCPresel2;
//        string AUCStrPresel2 = f2Label + " W/ Presel, AUC = " + streamPresel2.str();
        //cout<<AUCStrPresel2<<endl;

        can->cd();
        string label1 = AUCStr1.c_str();
        //string label1 = f1Label.c_str();
        //string label2 = f2Label.c_str();
        //string labelPresel1 = AUCStrPresel1.c_str();
        string label2 = AUCStr2.c_str();
        //string labelPresel2 = AUCStrPresel2.c_str();

        legend->AddEntry(sigEff_vs_bkgEff1[i],label1.c_str());
        //legend->AddEntry(sigEff_vs_bkgEffPresel1[i],labelPresel1.c_str());

        legend->AddEntry(sigEff_vs_bkgEff2[i],label2.c_str());
        //legend->AddEntry(sigEff_vs_bkgEffPresel2[i],labelPresel2.c_str());
//
        double sigEffHLTList[1] = {sigEffHLT};
        double bkgEffHLTList[1] = {bkgEffHLT};
//
        //hltPointPlot[i] = new TGraph (1, bkgEffHLTList, sigEffHLTList);
        hltPointPlot[i] = new TGraph ();
        //hltPointPlot[i]->SetPoint(0,-1,-1);
        hltPointPlot[i]->SetPoint(0,bkgEffHLT,sigEffHLT);
        
//        hltPointPlot[i]->SetLineColor(3);
//        hltPointPlot[i]->SetMarkerColor(3);
        
//        auto hltPointPlot = new TGraph();
//        hltPointPlot->AddPoint(bkgEffHLT,sigEffHLT);
        hltPointPlot[i]->SetMarkerColor(kBlue);
        hltPointPlot[i]->SetMarkerStyle(kFullCircle);

        
//        if (sigEff2[0] >= sigEff1[0]){
//            sigEff_vs_bkgEff2[i]->Draw();
//            sigEff_vs_bkgEff1[i]->Draw("same");
//        }
//        if (sigEff1[0] > sigEff2[0]){
//            sigEff_vs_bkgEff1[i]->Draw();
//            sigEff_vs_bkgEff2[i]->Draw("same");
//        }
        
        sigEff_vs_bkgEff1[i]->Draw();
        sigEff_vs_bkgEff2[i]->Draw("same");
        hltPointPlot[i]->Draw("sameP");
        //sigEff_vs_bkgEffPresel1[i]->Draw("same");
        //sigEff_vs_bkgEffPresel2[i]->Draw("same");

//        sigEff_vs_bkgEffPresel1[i]->Draw();
//        sigEff_vs_bkgEffPresel2[i]->Draw("same");
//        sigEff_vs_bkgEff1[i]->Draw("same");
//        sigEff_vs_bkgEff2[i]->Draw("same");
        
        //hltPointPlot[i]->Draw("same");

        
        
        //        sigEff_vs_bkgEff1[i]->Draw();

////        sigEff_vs_bkgEffPresel2[i]->Draw("same");

        legend->Draw();
        string outRoot = fOut[i] + ".root";
        string outPng = fOut[i] + ".png";
        can->Print(outRoot.c_str());
        can->Print(outPng.c_str());
        
//        canS->cd();
//        canS->SetLogy();
//        hSig1[i]->Draw();
//        hSigPresel1[i]->SetLineColor(2);
//        hSigPresel1[i]->Draw("same");
//        legendS->AddEntry(hSig1[i],("Without Presel(" + to_string(hSig1[i]->Integral()) + ")" ).c_str());
//        legendS->AddEntry(hSigPresel1[i],("With Presel(" + to_string(hSigPresel1[i]->Integral()) + ")" ).c_str());
//        legendS->Draw("same");
//        string outRootS = fOutSig[i] + ".root";
//        string outPngS = fOutSig[i] + ".png";
//        canS->Print(outRootS.c_str());
//        canS->Print(outPngS.c_str());
//
//        canB->cd();
//        canB->SetLogy();
//        hBkg1[i]->Draw();
//        hBkgPresel1[i]->SetLineColor(2);
//        hBkgPresel1[i]->Draw("same");
//        legendB->AddEntry(hBkg1[i],("Without Presel(" + to_string(hBkg1[i]->Integral()) + ")" ).c_str());
//        legendB->AddEntry(hBkgPresel1[i],("With Presel(" + to_string(hBkgPresel1[i]->Integral()) + ")" ).c_str());
//        legendB->Draw("same");
//        string outRootB = fOutBkg[i] + ".root";
//        string outPngB = fOutBkg[i] + ".png";
//        canB->Print(outRootB.c_str());
//        canB->Print(outPngB.c_str());
                    
    }
}
