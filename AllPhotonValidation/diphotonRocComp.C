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
#include "./pairSelection.h"

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
    string outStrGen = "outROC/0814/MD9LRComp";
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
    
    int nCurves = 3;

    //string fileNames[3];
    //string fLabels[3];
    string tNameSig = "sigTree";
    string tNameBkg = "bkgTree";
    
    //string fLabel1 = "MD 9 LR 0.15, No Scale Pos Weight";
    //string fileName1 = "validationNTuples/0808/M9LR15_NoTrkIso_M60_GGH13andDataD_CombinedAllTrain_AllCombinedValidated.root";
    
    string fLabel1 = "MD 9 LR 0.15";
    string fileName1 = "validationNTuples/0814/M9L15_GGH13andDataD_CombinedAll_NoTrkIso_M60_0814_AllCombValidation_Validated.root";
    
    string fLabel2 = "MD 9 LR 0.25";
    string fileName2 = "validationNTuples/0814/M9L25_GGH13andDataD_CombinedAll_NoTrkIso_M60_0814_AllCombValidation_Validated.root";

    string fLabel3 = "MD 9 LR 0.35";
    string fileName3 = "validationNTuples/0814/M9L35_GGH13andDataD_CombinedAll_NoTrkIso_M60_0814_AllCombValidation_Validated.root";

    string fileNames[3] = {fileName1,fileName2,fileName3};
    string fLabels[3] = {fLabel1,fLabel2,fLabel3};
    
    int iLow = 0;
    int iHigh = 4;
    
    TH1F *hSig1[4];
    TH1F *hSig2[4];
    TH1F *hSig3[4];

    TH1F *hBkg1[4];
    TH1F *hBkg2[4];
    TH1F *hBkg3[4];

    TGraph *sigEff_vs_bkgEff1[4];
    TGraph *sigEff_vs_bkgEff2[4];
    TGraph *sigEff_vs_bkgEff3[4];
    TGraph *hltPointPlot[4];

    for(int i = iLow; i < iHigh; i++){
        string titlePrompt1 = "hSig1" + to_string(i);
        string titleFake1 = "hBkg1" + to_string(i);
        hSig1[i] = new TH1F (titlePrompt1.c_str(),"",nBins,plotLow,plotHigh);
        hBkg1[i] = new TH1F (titleFake1.c_str(),"",nBins,plotLow,plotHigh);
        
        string titlePrompt2 = "hSig2" + to_string(i);
        string titleFake2 = "hBkg2" + to_string(i);
        hSig2[i] = new TH1F (titlePrompt2.c_str(),"",nBins,plotLow,plotHigh);
        hBkg2[i] = new TH1F (titleFake2.c_str(),"",nBins,plotLow,plotHigh);
        
        string titlePrompt3 = "hSig3" + to_string(i);
        string titleFake3 = "hBkg13" + to_string(i);
        hSig3[i] = new TH1F (titlePrompt3.c_str(),"",nBins,plotLow,plotHigh);
        hBkg3[i] = new TH1F (titleFake3.c_str(),"",nBins,plotLow,plotHigh);
    }
    
   //FOR FIRST ROC CURVE
    TFile *f1 = new TFile(fileName1.c_str());

    TTreeReader inReaderSig1(tNameSig.c_str(), f1);
    TTreeReaderArray<Float_t> inEtSig1(inReaderSig1,"et");
    TTreeReaderArray<Float_t> inEtaSig1(inReaderSig1, "eta");
    TTreeReaderArray<Float_t> inPhiSig1(inReaderSig1, "phi");
    TTreeReaderArray<Float_t> inXGBSig1(inReaderSig1, "xgbScore");
    TTreeReaderValue<int> inNEgsSig1(inReaderSig1,"nEgs");

    TTreeReader inReaderBkg1(tNameBkg.c_str(), f1);
    TTreeReaderArray<Float_t> inEtBkg1(inReaderBkg1,"et");
    TTreeReaderArray<Float_t> inEtaBkg1(inReaderBkg1, "eta");
    TTreeReaderArray<Float_t> inPhiBkg1(inReaderBkg1, "phi");
    TTreeReaderArray<Float_t> inXGBBkg1(inReaderBkg1, "xgbScore");
    TTreeReaderValue<int> inNEgsBkg1(inReaderBkg1,"nEgs");


    while (inReaderSig1.Next()) {
        int iOutXGBSig,jOutXGBSig;
        double massCalcXGBSig = calcMassXGB(&inXGBSig1,0.25,0.25,&inEtSig1,&inEtaSig1,&inPhiSig1,*inNEgsSig1,&iOutXGBSig,&jOutXGBSig);
        if(abs(inEtaSig1[iOutXGBSig]) < 1.4442 && abs(inEtaSig1[jOutXGBSig]) < 1.4442 && massCalcXGBSig > 95){
            if(inXGBSig1[iOutXGBSig] > inXGBSig1[jOutXGBSig])hSig1[0]->Fill(inXGBSig1[jOutXGBSig]);
            if(inXGBSig1[iOutXGBSig] < inXGBSig1[jOutXGBSig])hSig1[0]->Fill(inXGBSig1[iOutXGBSig]);
        }
        if(((abs(inEtaSig1[iOutXGBSig]) < 1.4442 && abs(inEtaSig1[jOutXGBSig]) > 1.556) || (abs(inEtaSig1[iOutXGBSig]) > 1.556 && abs(inEtaSig1[jOutXGBSig]) < 1.4442)) && massCalcXGBSig > 95){
            if(inXGBSig1[iOutXGBSig] > inXGBSig1[jOutXGBSig])hSig1[1]->Fill(inXGBSig1[jOutXGBSig]);
            if(inXGBSig1[iOutXGBSig] < inXGBSig1[jOutXGBSig])hSig1[1]->Fill(inXGBSig1[iOutXGBSig]);
        }
        if(abs(inEtaSig1[iOutXGBSig]) > 1.556 && abs(inEtaSig1[jOutXGBSig]) > 1.556 && massCalcXGBSig > 95){
            if(inXGBSig1[iOutXGBSig] > inXGBSig1[jOutXGBSig])hSig1[2]->Fill(inXGBSig1[jOutXGBSig]);
            if(inXGBSig1[iOutXGBSig] < inXGBSig1[jOutXGBSig])hSig1[2]->Fill(inXGBSig1[iOutXGBSig]);
        }
    }

    while (inReaderBkg1.Next()) {
        int iOutXGBBkg,jOutXGBBkg;
        double massCalcXGBBkg = calcMassXGB(&inXGBBkg1,0.25,0.25,&inEtBkg1,&inEtaBkg1,&inPhiBkg1,*inNEgsBkg1,&iOutXGBBkg,&jOutXGBBkg);
        if(abs(inEtaBkg1[iOutXGBBkg]) < 1.4442 && abs(inEtaBkg1[jOutXGBBkg]) < 1.4442 && massCalcXGBBkg > 95){
            if(inXGBBkg1[iOutXGBBkg] > inXGBBkg1[jOutXGBBkg])hBkg1[0]->Fill(inXGBBkg1[jOutXGBBkg]);
            if(inXGBBkg1[iOutXGBBkg] < inXGBBkg1[jOutXGBBkg])hBkg1[0]->Fill(inXGBBkg1[iOutXGBBkg]);
        }
        if(((abs(inEtaBkg1[iOutXGBBkg]) < 1.4442 && abs(inEtaBkg1[jOutXGBBkg]) > 1.556) || (abs(inEtaBkg1[iOutXGBBkg]) > 1.556 && abs(inEtaBkg1[jOutXGBBkg]) < 1.4442)) && massCalcXGBBkg > 95.0){
            if(inXGBBkg1[iOutXGBBkg] > inXGBBkg1[jOutXGBBkg])hBkg1[1]->Fill(inXGBBkg1[jOutXGBBkg]);
            if(inXGBBkg1[iOutXGBBkg] < inXGBBkg1[jOutXGBBkg])hBkg1[1]->Fill(inXGBBkg1[iOutXGBBkg]);
        }
        if(abs(inEtaBkg1[iOutXGBBkg]) > 1.556 && abs(inEtaBkg1[jOutXGBBkg]) > 1.556 && massCalcXGBBkg > 95){
            if(inXGBBkg1[iOutXGBBkg] > inXGBBkg1[jOutXGBBkg])hBkg1[2]->Fill(inXGBBkg1[jOutXGBBkg]);
            if(inXGBBkg1[iOutXGBBkg] < inXGBBkg1[jOutXGBBkg])hBkg1[2]->Fill(inXGBBkg1[iOutXGBBkg]);
        }
    }

    //FOR SECOND ROC CURVE
    TFile *f2 = new TFile(fileName2.c_str());

    TTreeReader inReaderSig2(tNameSig.c_str(), f2);
    TTreeReaderArray<Float_t> inEtSig2(inReaderSig2,"et");
    TTreeReaderArray<Float_t> inEtaSig2(inReaderSig2, "eta");
    TTreeReaderArray<Float_t> inPhiSig2(inReaderSig2, "phi");
    TTreeReaderArray<Float_t> inXGBSig2(inReaderSig2, "xgbScore");
    TTreeReaderValue<int> inNEgsSig2(inReaderSig2,"nEgs");

    TTreeReader inReaderBkg2(tNameBkg.c_str(), f2);
    TTreeReaderArray<Float_t> inEtBkg2(inReaderBkg2,"et");
    TTreeReaderArray<Float_t> inEtaBkg2(inReaderBkg2, "eta");
    TTreeReaderArray<Float_t> inPhiBkg2(inReaderBkg2, "phi");
    TTreeReaderArray<Float_t> inXGBBkg2(inReaderBkg2, "xgbScore");
    TTreeReaderValue<int> inNEgsBkg2(inReaderBkg2,"nEgs");


    while (inReaderSig2.Next()) {
        int iOutXGBSig,jOutXGBSig;
        double massCalcXGBSig = calcMassXGB(&inXGBSig2,0.25,0.25,&inEtSig2,&inEtaSig2,&inPhiSig2,*inNEgsSig2,&iOutXGBSig,&jOutXGBSig);
        if(abs(inEtaSig2[iOutXGBSig]) < 1.4442 && abs(inEtaSig2[jOutXGBSig]) < 1.4442 && massCalcXGBSig > 95){
            if(inXGBSig2[iOutXGBSig] > inXGBSig2[jOutXGBSig])hSig2[0]->Fill(inXGBSig2[jOutXGBSig]);
            if(inXGBSig2[iOutXGBSig] < inXGBSig2[jOutXGBSig])hSig2[0]->Fill(inXGBSig2[iOutXGBSig]);
        }
        if(((abs(inEtaSig2[iOutXGBSig]) < 1.4442 && abs(inEtaSig2[jOutXGBSig]) > 1.556) || (abs(inEtaSig2[iOutXGBSig]) > 1.556 && abs(inEtaSig2[jOutXGBSig]) < 1.4442)) && massCalcXGBSig > 95){
            if(inXGBSig2[iOutXGBSig] > inXGBSig2[jOutXGBSig])hSig2[1]->Fill(inXGBSig2[jOutXGBSig]);
            if(inXGBSig2[iOutXGBSig] < inXGBSig2[jOutXGBSig])hSig2[1]->Fill(inXGBSig2[iOutXGBSig]);
        }
        if(abs(inEtaSig2[iOutXGBSig]) > 1.556 && abs(inEtaSig2[jOutXGBSig]) > 1.556 && massCalcXGBSig > 95){
            if(inXGBSig2[iOutXGBSig] > inXGBSig2[jOutXGBSig])hSig2[2]->Fill(inXGBSig2[jOutXGBSig]);
            if(inXGBSig2[iOutXGBSig] < inXGBSig2[jOutXGBSig])hSig2[2]->Fill(inXGBSig2[iOutXGBSig]);
        }
    }

    while (inReaderBkg2.Next()) {
        int iOutXGBBkg,jOutXGBBkg;
        double massCalcXGBBkg = calcMassXGB(&inXGBBkg2,0.25,0.25,&inEtBkg2,&inEtaBkg2,&inPhiBkg2,*inNEgsBkg2,&iOutXGBBkg,&jOutXGBBkg);
        if(abs(inEtaBkg2[iOutXGBBkg]) < 1.4442 && abs(inEtaBkg2[jOutXGBBkg]) < 1.4442 && massCalcXGBBkg > 95){
            if(inXGBBkg2[iOutXGBBkg] > inXGBBkg2[jOutXGBBkg])hBkg2[0]->Fill(inXGBBkg2[jOutXGBBkg]);
            if(inXGBBkg2[iOutXGBBkg] < inXGBBkg2[jOutXGBBkg])hBkg2[0]->Fill(inXGBBkg2[iOutXGBBkg]);
        }
        if(((abs(inEtaBkg2[iOutXGBBkg]) < 1.4442 && abs(inEtaBkg2[jOutXGBBkg]) > 1.556) || (abs(inEtaBkg2[iOutXGBBkg]) > 1.556 && abs(inEtaBkg2[jOutXGBBkg]) < 1.4442)) && massCalcXGBBkg > 95.0){
            if(inXGBBkg2[iOutXGBBkg] > inXGBBkg2[jOutXGBBkg])hBkg2[1]->Fill(inXGBBkg2[jOutXGBBkg]);
            if(inXGBBkg2[iOutXGBBkg] < inXGBBkg2[jOutXGBBkg])hBkg2[1]->Fill(inXGBBkg2[iOutXGBBkg]);
        }
        if(abs(inEtaBkg2[iOutXGBBkg]) > 1.556 && abs(inEtaBkg2[jOutXGBBkg]) > 1.556 && massCalcXGBBkg > 95){
            if(inXGBBkg2[iOutXGBBkg] > inXGBBkg2[jOutXGBBkg])hBkg2[2]->Fill(inXGBBkg2[jOutXGBBkg]);
            if(inXGBBkg2[iOutXGBBkg] < inXGBBkg2[jOutXGBBkg])hBkg2[2]->Fill(inXGBBkg2[iOutXGBBkg]);
        }
    }
    
    //FOR THIRD ROC CURVE
    TFile *f3 = new TFile(fileName3.c_str());
    
    TTreeReader inReaderSig3(tNameSig.c_str(), f3);
    TTreeReaderArray<Float_t> inEtSig3(inReaderSig3,"et");
    TTreeReaderArray<Float_t> inEtaSig3(inReaderSig3, "eta");
    TTreeReaderArray<Float_t> inPhiSig3(inReaderSig3, "phi");
    TTreeReaderArray<Float_t> inXGBSig3(inReaderSig3, "xgbScore");
    TTreeReaderValue<int> inNEgsSig3(inReaderSig3,"nEgs");
    
    TTreeReader inReaderBkg3(tNameBkg.c_str(), f3);
    TTreeReaderArray<Float_t> inEtBkg3(inReaderBkg3,"et");
    TTreeReaderArray<Float_t> inEtaBkg3(inReaderBkg3, "eta");
    TTreeReaderArray<Float_t> inPhiBkg3(inReaderBkg3, "phi");
    TTreeReaderArray<Float_t> inXGBBkg3(inReaderBkg3, "xgbScore");
    TTreeReaderValue<int> inNEgsBkg3(inReaderBkg3,"nEgs");
    
    
    while (inReaderSig3.Next()) {
        int iOutXGBSig,jOutXGBSig;
        double massCalcXGBSig = calcMassXGB(&inXGBSig3,0.25,0.25,&inEtSig3,&inEtaSig3,&inPhiSig3,*inNEgsSig3,&iOutXGBSig,&jOutXGBSig);
        if(abs(inEtaSig3[iOutXGBSig]) < 1.4442 && abs(inEtaSig3[jOutXGBSig]) < 1.4442 && massCalcXGBSig > 95){
            if(inXGBSig3[iOutXGBSig] > inXGBSig3[jOutXGBSig])hSig3[0]->Fill(inXGBSig3[jOutXGBSig]);
            if(inXGBSig3[iOutXGBSig] < inXGBSig3[jOutXGBSig])hSig3[0]->Fill(inXGBSig3[iOutXGBSig]);
        }
        if(((abs(inEtaSig3[iOutXGBSig]) < 1.4442 && abs(inEtaSig3[jOutXGBSig]) > 1.556) || (abs(inEtaSig3[iOutXGBSig]) > 1.556 && abs(inEtaSig3[jOutXGBSig]) < 1.4442)) && massCalcXGBSig > 95){
            if(inXGBSig3[iOutXGBSig] > inXGBSig3[jOutXGBSig])hSig3[1]->Fill(inXGBSig3[jOutXGBSig]);
            if(inXGBSig3[iOutXGBSig] < inXGBSig3[jOutXGBSig])hSig3[1]->Fill(inXGBSig3[iOutXGBSig]);
        }
        if(abs(inEtaSig3[iOutXGBSig]) > 1.556 && abs(inEtaSig3[jOutXGBSig]) > 1.556 && massCalcXGBSig > 95){
            if(inXGBSig3[iOutXGBSig] > inXGBSig3[jOutXGBSig])hSig3[2]->Fill(inXGBSig3[jOutXGBSig]);
            if(inXGBSig3[iOutXGBSig] < inXGBSig3[jOutXGBSig])hSig3[2]->Fill(inXGBSig3[iOutXGBSig]);
        }
    }
    
    while (inReaderBkg3.Next()) {
        int iOutXGBBkg,jOutXGBBkg;
        double massCalcXGBBkg = calcMassXGB(&inXGBBkg3,0.25,0.25,&inEtBkg3,&inEtaBkg3,&inPhiBkg3,*inNEgsBkg3,&iOutXGBBkg,&jOutXGBBkg);
        if(abs(inEtaBkg3[iOutXGBBkg]) < 1.4442 && abs(inEtaBkg3[jOutXGBBkg]) < 1.4442 && massCalcXGBBkg > 95){
            if(inXGBBkg3[iOutXGBBkg] > inXGBBkg3[jOutXGBBkg])hBkg3[0]->Fill(inXGBBkg3[jOutXGBBkg]);
            if(inXGBBkg3[iOutXGBBkg] < inXGBBkg3[jOutXGBBkg])hBkg3[0]->Fill(inXGBBkg3[iOutXGBBkg]);
        }
        if(((abs(inEtaBkg3[iOutXGBBkg]) < 1.4442 && abs(inEtaBkg3[jOutXGBBkg]) > 1.556) || (abs(inEtaBkg3[iOutXGBBkg]) > 1.556 && abs(inEtaBkg3[jOutXGBBkg]) < 1.44422)) && massCalcXGBBkg > 95.0){
            if(inXGBBkg3[iOutXGBBkg] > inXGBBkg3[jOutXGBBkg])hBkg3[1]->Fill(inXGBBkg3[jOutXGBBkg]);
            if(inXGBBkg3[iOutXGBBkg] < inXGBBkg3[jOutXGBBkg])hBkg3[1]->Fill(inXGBBkg3[iOutXGBBkg]);
        }
        if(abs(inEtaBkg3[iOutXGBBkg]) > 1.556 && abs(inEtaBkg3[jOutXGBBkg]) > 1.556 && massCalcXGBBkg > 95){
            if(inXGBBkg3[iOutXGBBkg] > inXGBBkg3[jOutXGBBkg])hBkg3[2]->Fill(inXGBBkg3[jOutXGBBkg]);
            if(inXGBBkg3[iOutXGBBkg] < inXGBBkg3[jOutXGBBkg])hBkg3[2]->Fill(inXGBBkg3[iOutXGBBkg]);
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

        TCanvas *can = new TCanvas (canTitle.c_str(),canTitle.c_str(),10,10,1600,900);
        can->SetGrid();


        TLegend *legend = new TLegend(0.35,0.1,0.65,0.4,"","brNDC");

        TLegend *legendS = new TLegend(0.35,0.1,0.65,0.4,"","brNDC");
        TLegend *legendB = new TLegend(0.35,0.1,0.65,0.4,"","brNDC");

        if (i == 3){
            for (int p = 0; p < 3; p++){
                hSig1[i]->Add(hSig1[p]);
                hBkg1[i]->Add(hBkg1[p]);

                hSig2[i]->Add(hSig2[p]);
                hBkg2[i]->Add(hBkg2[p]);
                
                hSig3[i]->Add(hSig3[p]);
                hBkg3[i]->Add(hBkg3[p]);

            }
        }

        double nEventsP1 = hSig1[i]->GetSumOfWeights();
        double nEventsF1 = hBkg1[i]->GetSumOfWeights();

        double nEventsP2 = hSig2[i]->GetSumOfWeights();
        double nEventsF2 = hBkg2[i]->GetSumOfWeights();
        
        double nEventsP3 = hSig3[i]->GetSumOfWeights();
        double nEventsF3 = hBkg3[i]->GetSumOfWeights();

        int nCuts = nBins;
        double stepVal = binSizeFine;
        double cutVal = 0.0 - stepVal;
    //
        float Nsig1[nCuts+1], Nbkg1[nCuts+1];
        float sigEff1[nCuts+1], bkgEff1[nCuts+1];

        float Nsig2[nCuts+1], Nbkg2[nCuts+1];
        float sigEff2[nCuts+1], bkgEff2[nCuts+1];

        float Nsig3[nCuts+1], Nbkg3[nCuts+1];
        float sigEff3[nCuts+1], bkgEff3[nCuts+1];

        
        float cutsVal[nCuts+1];
        float mvaResCutVal = 0.0 - stepVal;

        int mvaSMaxBin = hSig1[i]->GetXaxis()->FindBin(1.0) + 1;
        int mvaBMaxBin = hBkg1[i]->GetXaxis()->FindBin(1.0) + 1;

        for(int k = 0; k < nCuts+1; k++){
            mvaResCutVal+= stepVal;
            cutsVal[k] = mvaResCutVal;

            int mvaBinSig = hSig1[i]->GetXaxis()->FindBin(mvaResCutVal);
            Nsig1[k] = hSig1[i]->Integral(mvaBinSig,mvaSMaxBin);
            Nsig2[k] = hSig2[i]->Integral(mvaBinSig,mvaSMaxBin);
            Nsig3[k] = hSig3[i]->Integral(mvaBinSig,mvaSMaxBin);

            int mvaBinBkg = hBkg1[i]->GetXaxis()->FindBin(mvaResCutVal);
            Nbkg1[k] = hBkg1[i]->Integral(mvaBinBkg,mvaBMaxBin);
            Nbkg2[k] = hBkg2[i]->Integral(mvaBinBkg,mvaBMaxBin);
            Nbkg3[k] = hBkg3[i]->Integral(mvaBinBkg,mvaBMaxBin);

            sigEff1[k] = Nsig1[k]/Nsig1[0];
            bkgEff1[k] = Nbkg1[k]/Nbkg1[0];

            sigEff2[k] = Nsig2[k]/Nsig2[0];
            bkgEff2[k] = Nbkg2[k]/Nbkg2[0];
            
            sigEff3[k] = Nsig3[k]/Nsig3[0];
            bkgEff3[k] = Nbkg3[k]/Nbkg3[0];

        }


        sigEff_vs_bkgEff1[i] = new TGraph (nCuts, bkgEff1, sigEff1);

        sigEff_vs_bkgEff2[i] = new TGraph (nCuts, bkgEff2, sigEff2);
        sigEff_vs_bkgEff2[i]->SetLineColor(2);
        sigEff_vs_bkgEff2[i]->SetMarkerColor(2);
        
        sigEff_vs_bkgEff3[i] = new TGraph (nCuts, bkgEff3, sigEff3);
        sigEff_vs_bkgEff3[i]->SetLineColor(3);
        sigEff_vs_bkgEff3[i]->SetMarkerColor(3);
        
        sigEff_vs_bkgEff1[i]->SetTitle(title.c_str());
        sigEff_vs_bkgEff1[i]->GetXaxis()->SetTitle("Background Efficiency");
        sigEff_vs_bkgEff1[i]->GetYaxis()->SetTitle("Signal Efficiency");

        TF1 f1("f1",[&](double *bkgEff1, double *sigEff1){ return sigEff_vs_bkgEff1[i]->Eval(bkgEff1[0]); },0,1,0); //Gets values of TGraph
        double AUC1 = f1.Integral(0,1);
        stringstream stream1;
        stream1 << fixed << setprecision(5) << AUC1;
        string AUCStr1 = fLabel1 + " All, AUC = " + stream1.str();
        cout<<AUCStr1.c_str()<<endl;


        TF1 f2("f2",[&](double *bkgEff2, double *sigEff2){ return sigEff_vs_bkgEff2[i]->Eval(bkgEff2[0]); },0,1,0); //Gets values of TGraph
        double AUC2 = f2.Integral(0,1);
        stringstream stream2;
        stream2 << fixed << setprecision(5) << AUC2;
        string AUCStr2 = fLabel2 + " All, AUC = " + stream2.str();
        cout<<AUCStr2.c_str()<<endl;

        TF1 f3("f3",[&](double *bkgEff3, double *sigEff3){ return sigEff_vs_bkgEff3[i]->Eval(bkgEff3[0]); },0,1,0); //Gets values of TGraph
        double AUC3 = f3.Integral(0,1);
        stringstream stream3;
        stream3 << fixed << setprecision(5) << AUC3;
        string AUCStr3 = fLabel3 + " All, AUC = " + stream3.str();
        cout<<AUCStr3.c_str()<<endl;

        can->cd();
        string label1 = AUCStr1.c_str();
        legend->AddEntry(sigEff_vs_bkgEff1[i],label1.c_str());

        string label2 = AUCStr2.c_str();
        legend->AddEntry(sigEff_vs_bkgEff2[i],label2.c_str());
        
        string label3 = AUCStr3.c_str();
        legend->AddEntry(sigEff_vs_bkgEff3[i],label3.c_str());

        hltPointPlot[i] = new TGraph ();
        hltPointPlot[i]->SetPoint(0,bkgEffHLT,sigEffHLT);


        hltPointPlot[i]->SetMarkerColor(kBlue);
        hltPointPlot[i]->SetMarkerStyle(kFullCircle);

        sigEff_vs_bkgEff1[i]->Draw();
        sigEff_vs_bkgEff2[i]->Draw("same");
        sigEff_vs_bkgEff3[i]->Draw("same");
        hltPointPlot[i]->Draw("sameP");


        legend->Draw();
        string outRoot = fOut[i] + ".root";
        string outPng = fOut[i] + ".png";
        can->Print(outRoot.c_str());
        can->Print(outPng.c_str());


    }
}


//FAILED ATTEMPT TO GENERALIZE TO LOOP OVER CURVES
//for (int i = 0; i < nCurves; i++){
//    TFile *f = new TFile(fileNames[i].c_str());
//
//    TTreeReader inReaderSig(tNameSig.c_str(), f);
//    TTreeReaderArray<Float_t> inEtSig(inReaderSig,"et");
//    TTreeReaderArray<Float_t> inEtaSig(inReaderSig, "eta");
//    TTreeReaderArray<Float_t> inPhiSig(inReaderSig, "phi");
//    TTreeReaderArray<Float_t> inXGBSig(inReaderSig, "xgbScore");
//    TTreeReaderValue<int> inNEgsSig(inReaderSig,"nEgs");
//
//    TTreeReader inReaderBkg(tNameBkg.c_str(), f);
//    TTreeReaderArray<Float_t> inEtBkg(inReaderBkg,"et");
//    TTreeReaderArray<Float_t> inEtaBkg(inReaderBkg, "eta");
//    TTreeReaderArray<Float_t> inPhiBkg(inReaderBkg, "phi");
//    TTreeReaderArray<Float_t> inXGBBkg(inReaderBkg, "xgbScore");
//    TTreeReaderValue<int> inNEgsBkg(inReaderBkg,"nEgs");
//
//
//    while (inReaderSig.Next()) {
//        int iOutXGBSig,jOutXGBSig;
//        double massCalcXGBSig = calcMassXGB(&inXGBSig,0.25,0.25,&inEtSig,&inEtaSig,&inPhiSig,*inNEgsSig,&iOutXGBSig,&jOutXGBSig);
//        if(abs(inEtaSig[iOutXGBSig]) < 1.4442 && abs(inEtaSig[jOutXGBSig]) < 1.4442 && massCalcXGBSig > 95){
//            if(inXGBSig[iOutXGBSig] > inXGBSig[jOutXGBSig])hSig[i][0]->Fill(inXGBSig[jOutXGBSig]);
//            if(inXGBSig[iOutXGBSig] < inXGBSig[jOutXGBSig])hSig[i][0]->Fill(inXGBSig[iOutXGBSig]);
//        }
//        if(((abs(inEtaSig[iOutXGBSig]) < 1.4442 && abs(inEtaSig[jOutXGBSig]) > 1.556) || (abs(inEtaSig[iOutXGBSig]) > 1.556 && abs(inEtaSig[jOutXGBSig]) < 1.4442)) && massCalcXGBSig > 95){
//            if(inXGBSig[iOutXGBSig] > inXGBSig[jOutXGBSig])hSig[i][1]->Fill(inXGBSig[jOutXGBSig]);
//            if(inXGBSig[iOutXGBSig] < inXGBSig[jOutXGBSig])hSig[i][1]->Fill(inXGBSig[iOutXGBSig]);
//        }
//        if(abs(inEtaSig[iOutXGBSig]) > 1.556 && abs(inEtaSig[jOutXGBSig]) > 1.556 && massCalcXGBSig > 95){
//            if(inXGBSig[iOutXGBSig] > inXGBSig[jOutXGBSig])hSig[i][2]->Fill(inXGBSig[jOutXGBSig]);
//            if(inXGBSig[iOutXGBSig] < inXGBSig[jOutXGBSig])hSig[i][2]->Fill(inXGBSig[iOutXGBSig]);
//        }
//    }
//
//    while (inReaderBkg.Next()) {
//        int iOutXGBBkg,jOutXGBBkg;
//        double massCalcXGBBkg = calcMassXGB(&inXGBBkg,0.25,0.25,&inEtBkg,&inEtaBkg,&inPhiBkg,*inNEgsBkg,&iOutXGBBkg,&jOutXGBBkg);
//        if(abs(inEtaBkg[iOutXGBBkg]) < 1.4442 && abs(inEtaBkg[jOutXGBBkg]) < 1.4442 && massCalcXGBBkg > 95){
//            if(inXGBBkg[iOutXGBBkg] > inXGBBkg[jOutXGBBkg])hBkg[i][0]->Fill(inXGBBkg[jOutXGBBkg]);
//            if(inXGBBkg[iOutXGBBkg] < inXGBBkg[jOutXGBBkg])hBkg[i][0]->Fill(inXGBBkg[iOutXGBBkg]);
//        }
//        if(((abs(inEtaBkg[iOutXGBBkg]) < 1.4442 && abs(inEtaBkg[jOutXGBBkg]) > 1.556) || (abs(inEtaBkg[iOutXGBBkg]) > 1.556 && abs(inEtaBkg[jOutXGBBkg]) < 1.4442)) && massCalcXGBBkg > 95){
//            if(inXGBBkg[iOutXGBBkg] > inXGBBkg[jOutXGBBkg])hBkg[i][1]->Fill(inXGBBkg[jOutXGBBkg]);
//            if(inXGBBkg[iOutXGBBkg] < inXGBBkg[jOutXGBBkg])hBkg[i][1]->Fill(inXGBBkg[iOutXGBBkg]);
//        }
//        if(abs(inEtaBkg[iOutXGBBkg]) > 1.556 && abs(inEtaBkg[jOutXGBBkg]) > 1.556 && massCalcXGBBkg > 95){
//            if(inXGBBkg[iOutXGBBkg] > inXGBBkg[jOutXGBBkg])hBkg[i][2]->Fill(inXGBBkg[jOutXGBBkg]);
//            if(inXGBBkg[iOutXGBBkg] < inXGBBkg[jOutXGBBkg])hBkg[i][2]->Fill(inXGBBkg[iOutXGBBkg]);
//        }
//    }
//
//
//    for(int i = iLow; i < iHigh; i++){
//
//
//        double maxVal = 0.0;
//        int maxInt = -1;
//
//        string title = "RoC Curves with Data Bkg " + eta[i];
//        string canTitle = "canR" + to_string(i);
//
//        string titlePromptCan = "xgbScores for Signal " + eta[i];
//        string canTitlePrompt = "canS" + to_string(i);
//        string titleFakeCan = "xgbScores for Background " + eta[i];
//        string canTitleFake = "canB" + to_string(i);
//
//        TCanvas *can = new TCanvas (canTitle.c_str(),canTitle.c_str(),10,10,1600,900);
//        can->SetGrid();
//
//        TLegend *legend = new TLegend(0.35,0.1,0.65,0.4,"","brNDC");
//
//        TLegend *legendS = new TLegend(0.35,0.1,0.65,0.4,"","brNDC");
//        TLegend *legendB = new TLegend(0.35,0.1,0.65,0.4,"","brNDC");
//
//        if (i == 3){
//            for (int j = 0; j < nCurves; j++){
//                for (int p = 0; p < 3; p++){
//                    hSig[j][i]->Add(hSig[j][p]);
//                    hBkg[j][i]->Add(hBkg[j][p]);
//                }
//            }
//        }
//
//        int nCuts = nBins;
//        double stepVal = binSizeFine;
//        double cutVal = 0.0 - stepVal;
//
//        double cutsVal[nCuts+1];
//
//        //Fill array with all cut value locations
//        for(int k = 0; k < nCuts+1; k++){
//            cutVal+= stepVal;
//            cutsVal[k] = cutVal;
//        }
//
//        double Nsig[3][nCuts+1], Nbkg[3][nCuts+1];
//        double sigEff[3][nCuts+1], bkgEff[3][nCuts+1];
//
//
//        for (int j = 0; j < nCurves; j++){
//            double nEventsP = hSig[j][i]->GetSumOfWeights();
//            double nEventsF = hBkg[j][i]->GetSumOfWeights();
//            int mvaSMaxBin = hSig[j][i]->GetXaxis()->FindBin(1.0) + 1;
//            int mvaBMaxBin = hBkg[j][i]->GetXaxis()->FindBin(1.0) + 1;
//
//            for(int k = 0; k < nCuts+1; k++){
//                double mvaResCutVal = cutsVal[k];
//                int mvaBinSig = hSig[j][i]->GetXaxis()->FindBin(mvaResCutVal);
//                Nsig[j][k] = hSig[j][i]->Integral(mvaBinSig,mvaSMaxBin);
//                cout<<"Nsig["<<j<<"]["<<k<<" = "<<Nsig[j][k]<<endl;
//
//                int mvaBinBkg = hBkg[j][i]->GetXaxis()->FindBin(mvaResCutVal);
//                Nbkg[j][k] = hBkg[j][i]->Integral(mvaBinBkg,mvaBMaxBin);
//                cout<<"Nbkg["<<j<<"]["<<k<<" = "<<Nbkg[j][k]<<endl;
//
//                sigEff[j][k] = Nsig[j][k]/Nsig[j][0];
//                bkgEff[j][k] = Nbkg[j][k]/Nbkg[j][0];
//
//                cout<<"sigEff["<<j<<"]["<<k<<" = "<<sigEff[j][k]<<endl;
//                cout<<"bkgEff["<<j<<"]["<<k<<" = "<<bkgEff[j][k]<<endl;
//
//            }
//
//            sigEff_vs_bkgEff[j][i] = new TGraph (nCuts, bkgEff[j], sigEff[j]);
//
//            if(j == 0){
//                sigEff_vs_bkgEff[j][i]->SetTitle(title.c_str());
//                sigEff_vs_bkgEff[j][i]->GetXaxis()->SetTitle("Background Efficiency");
//                sigEff_vs_bkgEff[j][i]->GetYaxis()->SetTitle("Signal Efficiency");
//            }
//
//            if (j == 1){
//                sigEff_vs_bkgEff[j][i]->SetLineColor(2);
//                sigEff_vs_bkgEff[j][i]->SetMarkerColor(2);
//            }
//            if (j == 2){
//                sigEff_vs_bkgEff[j][i]->SetLineColor(3);
//                sigEff_vs_bkgEff[j][i]->SetMarkerColor(3);
//            }
//            if (j == 3){
//                sigEff_vs_bkgEff[j][i]->SetLineColor(4);
//                sigEff_vs_bkgEff[j][i]->SetMarkerColor(4);
//            }
//
//
//            //                TF1 fit("fit",[&](double *(bkgEff[j]), double *(sigEff[j])){ return sigEff_vs_bkgEff[j][i]->Eval(bkgEff[j][0]); },0,1,0); //Gets values of TGraph
//            //                double AUC = fit.Integral(0,1);
//            double AUC = sigEff_vs_bkgEff[j][i]->Integral(0,1);
//            stringstream stream;
//            stream << fixed << setprecision(5) << AUC;
//            string AUCStr = fLabels[j] + " All, AUC = " + stream.str();
//            cout<<AUCStr.c_str()<<endl;
//
//
//            can->cd();
//            string label = AUCStr.c_str();
//
//            legend->AddEntry(sigEff_vs_bkgEff[j][i],label.c_str());
//
//            if(j==0)sigEff_vs_bkgEff[j][i]->Draw();
//            else sigEff_vs_bkgEff[j][i]->Draw("same");
//
//
//        }
//        hltPointPlot[i] = new TGraph ();
//        hltPointPlot[i]->SetPoint(0,bkgEffHLT,sigEffHLT);
//
//        hltPointPlot[i]->SetMarkerColor(kBlue);
//        hltPointPlot[i]->SetMarkerStyle(kFullCircle);
//
//        hltPointPlot[i]->Draw("sameP");
//
//        legend->Draw();
//
//        string outRoot = fOut[i] + ".root";
//        string outPng = fOut[i] + ".png";
//        can->Print(outRoot.c_str());
//        can->Print(outPng.c_str());
//
//
//    }
//}
