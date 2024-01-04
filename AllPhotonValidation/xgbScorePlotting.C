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

void xgbScorePlotting(){
    gROOT->Reset();
    gStyle->SetPalette(1);
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitle(0);
    
    string fileName = "validationNTuples/1117/M7L25_GGH13andDataD_NoTrkIso_M60_PdgIDTrain_1117_PdgIDVali_1117.root";
    string genTitleStringSignal = "GGH Signal ";
    string genTitleStringBkg = "Data ";
    string dirStr ="sigBkgVarPlots/1117/MD7L25_GGH13andDataD_NoTrkIso_M60_pdgIDCutValiANDTrain_Added/";
    string plotType = "Added";
    string out2D = dirStr + "highandLowScore_XGBCut_M95_1117";
    string outLead = dirStr + "highScore_XGBCut_M95_1117";
    string outSub = dirStr + "lowScore_XGBCut_M95_1117";
    //string plotType = "Separate";
    //string dirStr ="varPlots/0130/DataRelaxedUnseededChoose2/";
    
    int nEta;
    
    nEta = 4;
    string etaLabels[4] = {"Barrel-Barrel","Barrel-Endcap","Endcap-Endcap","All #eta"};
    string etaFLabels[4] = {"_BB","_BE","_EE","_All"};
    
    //New cuts (11/09) chosen by hand (just Barrel and Endcap)
    double leadCuts1[2] = {0.85,0.90};
    double subCuts1[2] = {0.0,0.0};//If First score above leadCuts1, second score must be above subCuts1
    double leadCuts2[2] = {0.75,0.8};
    double subCuts2[2] = {0.015,0.02};//If First score between leadCuts1 and leadCuts2, second score must be above subCuts2
    double subCuts3[2] = {0.075,0.075};//If First score below leadCuts2, second score must be above subCuts3
    
    double nBins = 500;
    
    double limLow = 0.0;
    double limHigh = 1.0;
    
    string outNameGen = dirStr;
    
    TH1F *hSigAllLead[4];
    TH1F *hSigAllSub[4];
    TH1F *hBkgAllLead[4];
    TH1F *hBkgAllSub[4];
    
    TH1F *hSigPassLooseXGBLead[4];
    TH1F *hSigPassLooseXGBSub[4];
    TH1F *hBkgPassLooseXGBLead[4];
    TH1F *hBkgPassLooseXGBSub[4];
    
    TH1F *hSigPassTightXGBLead[4];
    TH1F *hSigPassTightXGBSub[4];
    TH1F *hBkgPassTightXGBLead[4];
    TH1F *hBkgPassTightXGBSub[4];
    
    TH1F *hSigPassTightXGBTightMassLead[4];
    TH1F *hSigPassTightXGBTightMassSub[4];
    TH1F *hBkgPassTightXGBTightMassLead[4];
    TH1F *hBkgPassTightXGBTightMassSub[4];
    
    TH1F *hSigPassStdLead[4];
    TH1F *hSigPassStdSub[4];
    TH1F *hBkgPassStdLead[4];
    TH1F *hBkgPassStdSub[4];
    
    TH2F *hSig2DScores[4];
    TH2F *hSig2DScoresPass[4];
    TH2F *hBkg2DScores[4];
    TH2F *hBkg2DScoresPass[4];
    
    
    TCanvas *can = new TCanvas ("can","can",10,10,1600,900);
    
    for (int e = 0; e < nEta; e++){
        hSigAllLead[e] = new TH1F ("hSigAllLead","",nBins,limLow,limHigh);
        hSigAllSub[e] = new TH1F ("hSigAllSub","",nBins,limLow,limHigh);
        hBkgAllLead[e] = new TH1F ("hBkgAllLead","",nBins,limLow,limHigh);
        hBkgAllSub[e] = new TH1F ("hBkgAllSub","",nBins,limLow,limHigh);
        
        hSigPassLooseXGBLead[e] = new TH1F ("hSigPassLooseXGBLead","",nBins,limLow,limHigh);
        hSigPassLooseXGBSub[e] = new TH1F ("hSigPassLooseXGBSub","",nBins,limLow,limHigh);
        hBkgPassLooseXGBLead[e] = new TH1F ("hBkgPassLooseXGBLead","",nBins,limLow,limHigh);
        hBkgPassLooseXGBSub[e] = new TH1F ("hBkgPassLooseXGBSub","",nBins,limLow,limHigh);
        
        hSigPassTightXGBLead[e] = new TH1F ("hSigPassTightXGBLead","",nBins,limLow,limHigh);
        hSigPassTightXGBSub[e] = new TH1F ("hSigPassTightXGBSub","",nBins,limLow,limHigh);
        hBkgPassTightXGBLead[e] = new TH1F ("hBkgPassTightXGBLead","",nBins,limLow,limHigh);
        hBkgPassTightXGBSub[e] = new TH1F ("hBkgPassTightXGBSub","",nBins,limLow,limHigh);
        
        hSigPassTightXGBTightMassLead[e] = new TH1F ("hSigPassTightXGBTightMassLead","",nBins,limLow,limHigh);
        hSigPassTightXGBTightMassSub[e] = new TH1F ("hSigPassTightXGBTightMassSub","",nBins,limLow,limHigh);
        hBkgPassTightXGBTightMassLead[e] = new TH1F ("hBkgPassTightXGBTightMassLead","",nBins,limLow,limHigh);
        hBkgPassTightXGBTightMassSub[e] = new TH1F ("hBkgPassTightXGBTightMassSub","",nBins,limLow,limHigh);
        
        hSigPassStdLead[e] = new TH1F ("hSigPassStdLead","",nBins,limLow,limHigh);
        hSigPassStdSub[e] = new TH1F ("hSigPassStdSub","",nBins,limLow,limHigh);
        hBkgPassStdLead[e] = new TH1F ("hBkgPassStdLead","",nBins,limLow,limHigh);
        hBkgPassStdSub[e] = new TH1F ("hBkgPassStdSub","",nBins,limLow,limHigh);
        
        hSig2DScores[e] = new TH2F ("hSig2DScores","",nBins,limLow,limHigh,nBins,limLow,limHigh);
        hSig2DScoresPass[e] = new TH2F ("hSig2DScoresPass","",nBins,limLow,limHigh,nBins,limLow,limHigh);
        hBkg2DScores[e] = new TH2F ("hBkg2DScores","",nBins,limLow,limHigh,nBins,limLow,limHigh);
        hBkg2DScoresPass[e] = new TH2F ("hBkg2DScoresPass","",nBins,limLow,limHigh,nBins,limLow,limHigh);
        
    }
    TFile *f = new TFile(fileName.c_str());

    //FIRST handle Signal
    TTreeReader inReaderSig("sigTree", f);
    inReaderSig.Restart();
    
    TTreeReaderArray<Float_t> inESig(inReaderSig,"energy");
    TTreeReaderArray<Float_t> inEtSig(inReaderSig,"et");
    TTreeReaderArray<Float_t> inEtaSig(inReaderSig, "eta");
    TTreeReaderArray<Float_t> inPhiSig(inReaderSig, "phi");
    TTreeReaderArray<Float_t> inR9Sig(inReaderSig,"r9HLT");
    TTreeReaderArray<Float_t> inHovrESig(inReaderSig,"hOvrE");
    TTreeReaderArray<Float_t> inSigIEtaIEtaSig(inReaderSig,"sigmaIEtaIEta");
    TTreeReaderArray<Float_t> inPhoIsoSig(inReaderSig,"trkIsoPho");
    TTreeReaderArray<Float_t> inEcalIsoSig(inReaderSig,"ecalPFIso");
    TTreeReaderArray<Float_t> inXGBSig(inReaderSig, "xgbScore");
    TTreeReaderArray<int> inPdgIDSig(inReaderSig,"pdgID");

    TTreeReaderValue<int> inNEgsSig(inReaderSig,"nEgs");
    TTreeReaderValue<int> inTriggerBitsSig(inReaderSig,"triggerBits");
    TTreeReaderValue<int> inPassFailStdSig(inReaderSig,"passFailStd");
    TTreeReaderValue<int> inPassFailDoubleL1Sig(inReaderSig,"passFailL1Double");
    TTreeReaderValue<int> inPassFailSingleL1Sig(inReaderSig,"passFailL1Single");
    
    //NORMAL ETA THRESHOLDS ARE 1.444 and 1.556
    while (inReaderSig.Next()) {
        int iOutXGBSig,jOutXGBSig;
        double massCalcXGBSig = calcMassXGB(&inXGBSig,0.25,0.25,&inEtSig,&inEtaSig,&inPhiSig,*inNEgsSig,&iOutXGBSig,&jOutXGBSig);
        double drXGBSig = calcDR(inEtaSig[iOutXGBSig], inEtaSig[jOutXGBSig], inPhiSig[iOutXGBSig], inPhiSig[jOutXGBSig]);
        //if(inEtSig[iOutXGBSig] > 14.25 && inEtSig[jOutXGBSig] > 14.25 && abs(inEtaSig[iOutXGBSig]) < 2.5 && abs(inEtaSig[jOutXGBSig]) < 2.5 && iOutXGBSig != -1 && jOutXGBSig != -1 && massCalcXGBSig > 95.0){
        if(inPdgIDSig[iOutXGBSig] == 22 && inPdgIDSig[jOutXGBSig] == 22 && inEtSig[iOutXGBSig] > 14.25 && inEtSig[jOutXGBSig] > 14.25 && abs(inEtaSig[iOutXGBSig]) < 2.5 && abs(inEtaSig[jOutXGBSig]) < 2.5 && iOutXGBSig != -1 && jOutXGBSig != -1 && massCalcXGBSig > 95.0){
            double indexHigh, indexLow;
            double outOne,outTwo;
            //outOne = inXGBSig[iOutXGBSig];
            //outTwo = inXGBSig[jOutXGBSig];
            if(inXGBSig[iOutXGBSig] > inXGBSig[jOutXGBSig]){
                outOne = inXGBSig[iOutXGBSig];
                indexHigh = iOutXGBSig;
                outTwo = inXGBSig[jOutXGBSig];
                indexLow = jOutXGBSig;
            }
            if(inXGBSig[iOutXGBSig] < inXGBSig[jOutXGBSig]){
                outOne = inXGBSig[jOutXGBSig];
                indexHigh = jOutXGBSig;
                outTwo = inXGBSig[iOutXGBSig];
                indexLow = iOutXGBSig;
            }

            int e;//Combined eta index -- 0=BB, 1=BE, 2=EE
            int e1,e2;//Individual eta index -- 0=Barrel, 1=Endcap
            if (abs(inEtaSig[indexHigh]) < 1.5 && abs(inEtaSig[indexLow]) < 1.5){
                e = 0;
                e1 = 0;//Barrel
                e2 = 0;//Barrel
            }
            if (abs(inEtaSig[indexHigh]) < 1.5 && abs(inEtaSig[indexLow]) > 1.5){
                e = 1;
                e1 = 0;//Barrel
                e2 = 1;//Endcap
            }
            if(abs(inEtaSig[indexHigh]) > 1.5 && abs(inEtaSig[indexLow]) < 1.5){
                e = 1;
                e1 = 1;//Endcap
                e2 = 0;//Barrel
            }
            if (abs(inEtaSig[indexHigh]) > 1.5 && abs(inEtaSig[indexLow]) > 1.5){
                e = 2;
                e1 = 1;//Endcap
                e2 = 1;//Endcap
            }
            
            hSigAllLead[e]->Fill(outOne);
            hSigAllSub[e]->Fill(outTwo);
            
            hSig2DScores[e]->Fill(outOne,outTwo);
            if((inXGBSig[indexHigh] > leadCuts1[e1] && inXGBSig[indexLow] > subCuts1[e2])
               || (inXGBSig[indexHigh] > leadCuts2[e1] && inXGBSig[indexLow] > subCuts2[e2])
               || (inXGBSig[indexHigh] < leadCuts2[e1] && inXGBSig[indexLow] > subCuts3[e2])
               ){
                hSig2DScoresPass[e]->Fill(outOne,outTwo);
                if(massCalcXGBSig > 50.0){
                    hSigPassTightXGBLead[e]->Fill(outOne);
                    hSigPassTightXGBSub[e]->Fill(outTwo);
                }
                if(massCalcXGBSig > 95.0){
                    hSigPassTightXGBTightMassLead[e]->Fill(outOne);
                    hSigPassTightXGBTightMassSub[e]->Fill(outTwo);
                }
            }
            else{
                hSigPassLooseXGBLead[e]->Fill(outOne);
                hSigPassLooseXGBSub[e]->Fill(outTwo);
            }
        }//XGB Eta selection
        
        int iOutHLTSig,jOutHLTSig;
        double massCalcHLTSig = calcMassHLT(&inEtSig, &inEtaSig, &inPhiSig, &inR9Sig, &inHovrESig, &inSigIEtaIEtaSig, &inPhoIsoSig, &inEcalIsoSig, *inNEgsSig, &iOutHLTSig, &jOutHLTSig);
        double drHLTSig = calcDR(inEtaSig[iOutHLTSig], inEtaSig[jOutHLTSig], inPhiSig[iOutHLTSig], inPhiSig[jOutHLTSig]);
        //if(inEtSig[iOutHLTSig] > 14.25 && inEtSig[jOutHLTSig] > 14.25 && abs(inEtaSig[iOutHLTSig]) < 2.5 && abs(inEtaSig[jOutHLTSig]) < 2.5 && iOutHLTSig != -1 && jOutHLTSig != -1 && massCalcHLTSig > 95.0){
        if(inPdgIDSig[iOutHLTSig] == 22 && inPdgIDSig[jOutHLTSig] == 22 && inEtSig[iOutHLTSig] > 14.25 && inEtSig[jOutHLTSig] > 14.25 && abs(inEtaSig[iOutHLTSig]) < 2.5 && abs(inEtaSig[jOutHLTSig]) < 2.5 && iOutHLTSig != -1 && jOutHLTSig != -1 && massCalcHLTSig > 95.0){
            double outOne,outTwo;
            if(inXGBSig[iOutHLTSig] > inXGBSig[jOutHLTSig]){
                outOne = inXGBSig[iOutHLTSig];
                outTwo = inXGBSig[jOutHLTSig];
            }
            if(inXGBSig[iOutHLTSig] < inXGBSig[jOutHLTSig]){
                outOne = inXGBSig[jOutHLTSig];
                outTwo = inXGBSig[iOutHLTSig];
            }
            
            int e;
            if (abs(inEtaSig[iOutHLTSig]) < 1.5 && abs(inEtaSig[jOutHLTSig]) < 1.5)e = 0;
            if ((abs(inEtaSig[iOutHLTSig]) > 1.5 && abs(inEtaSig[jOutHLTSig]) < 1.5)
                || (abs(inEtaSig[iOutHLTSig]) < 1.5 && abs(inEtaSig[jOutHLTSig]) > 1.5))e = 1;
            if (abs(inEtaSig[iOutHLTSig]) > 1.5 && abs(inEtaSig[jOutHLTSig]) > 1.5)e = 2;
            
            //if(*inPassFailStdSig == 1)hSig2DScoresPass[e]->Fill(outOne,outTwo);
            
            if(*inPassFailStdSig == 1 && (*inPassFailSingleL1Sig == 1 || *inPassFailDoubleL1Sig == 1)){
                hSigPassStdLead[e]->Fill(outOne);
                hSigPassStdSub[e]->Fill(outTwo);
            }
        }//XGB Eta selection
    }
    
    //THEN handle Backgrouund
    TTreeReader inReaderBkg("bkgTree", f);
    inReaderBkg.Restart();
    
    TTreeReaderArray<Float_t> inEBkg(inReaderBkg,"energy");
    TTreeReaderArray<Float_t> inEtBkg(inReaderBkg,"et");
    TTreeReaderArray<Float_t> inEtaBkg(inReaderBkg, "eta");
    TTreeReaderArray<Float_t> inPhiBkg(inReaderBkg, "phi");
    TTreeReaderArray<Float_t> inR9Bkg(inReaderBkg,"r9HLT");
    TTreeReaderArray<Float_t> inHovrEBkg(inReaderBkg,"hOvrE");
    TTreeReaderArray<Float_t> inSigIEtaIEtaBkg(inReaderBkg,"sigmaIEtaIEta");
    TTreeReaderArray<Float_t> inPhoIsoBkg(inReaderBkg,"trkIsoPho");
    TTreeReaderArray<Float_t> inEcalIsoBkg(inReaderBkg,"ecalPFIso");
    TTreeReaderArray<Float_t> inXGBBkg(inReaderBkg, "xgbScore");
    
    TTreeReaderValue<int> inNEgsBkg(inReaderBkg,"nEgs");
    TTreeReaderValue<int> inTriggerBitsBkg(inReaderBkg,"triggerBits");
    TTreeReaderValue<int> inPassFailStdBkg(inReaderBkg,"passFailStd");
    TTreeReaderValue<int> inPassFailDoubleL1Bkg(inReaderBkg,"passFailL1Double");
    TTreeReaderValue<int> inPassFailSingleL1Bkg(inReaderBkg,"passFailL1Single");
    
    //NORMAL ETA THRESHOLDS ARE 1.444 and 1.556
    while (inReaderBkg.Next()) {
        int iOutXGBBkg,jOutXGBBkg;
        double massCalcXGBBkg = calcMassXGB(&inXGBBkg,0.25,0.25,&inEtBkg,&inEtaBkg,&inPhiBkg,*inNEgsBkg,&iOutXGBBkg,&jOutXGBBkg);
        double drXGBBkg = calcDR(inEtaBkg[iOutXGBBkg], inEtaBkg[jOutXGBBkg], inPhiBkg[iOutXGBBkg], inPhiBkg[jOutXGBBkg]);
        if(inEtBkg[iOutXGBBkg] > 14.25 && inEtBkg[jOutXGBBkg] > 14.25 && abs(inEtaBkg[iOutXGBBkg]) < 2.5 && abs(inEtaBkg[jOutXGBBkg]) < 2.5 && iOutXGBBkg != -1 && jOutXGBBkg != -1 && massCalcXGBBkg > 95.0){
            //cout<<"For XGB, (i,j) = ("<<iOutXGBBkg<<","<<jOutXGBBkg<<"), mass = "<<massCalcXGBBkg<<endl;
            double outOne,outTwo;
            double indexHigh,indexLow;
            //outOne = inXGBBkg[iOutXGBBkg];
            //outTwo = inXGBBlg[jOutXGBBkg];
            if(inXGBBkg[iOutXGBBkg] > inXGBBkg[jOutXGBBkg]){
                outOne = inXGBBkg[iOutXGBBkg];
                indexHigh = iOutXGBBkg;
                outTwo = inXGBBkg[jOutXGBBkg];
                indexLow = jOutXGBBkg;
            }
            if(inXGBBkg[iOutXGBBkg] < inXGBBkg[jOutXGBBkg]){
                outOne = inXGBBkg[jOutXGBBkg];
                indexHigh = iOutXGBBkg;
                outTwo = inXGBBkg[iOutXGBBkg];
                indexLow = jOutXGBBkg;
            }
            int e;//Combined eta index -- 0=BB, 1=BE, 2=EE
            int e1,e2;//Individual eta index -- 0=Barrel, 1=Endcap
            if (abs(inEtaBkg[indexHigh]) < 1.5 && abs(inEtaBkg[indexLow]) < 1.5){
                e = 0;
                e1 = 0;//Barrel
                e2 = 0;//Barrel
            }
            if (abs(inEtaBkg[indexHigh]) < 1.5 && abs(inEtaBkg[indexLow]) > 1.5){
                e = 1;
                e1 = 0;//Barrel
                e2 = 1;//Endcap
            }
            if(abs(inEtaBkg[indexHigh]) > 1.5 && abs(inEtaBkg[indexLow]) < 1.5){
                e = 1;
                e1 = 1;//Endcap
                e2 = 0;//Barrel
            }
            if (abs(inEtaBkg[indexHigh]) > 1.5 && abs(inEtaBkg[indexLow]) > 1.5){
                e = 2;
                e1 = 1;//Endcap
                e2 = 1;//Endcap
            }
            
            hBkgAllLead[e]->Fill(outOne);
            hBkgAllSub[e]->Fill(outTwo);
            
            hBkg2DScores[e]->Fill(outOne,outTwo);
            if((inXGBBkg[indexHigh] > leadCuts1[e1] && inXGBBkg[indexLow] > subCuts1[e2])
               || (inXGBBkg[indexHigh] > leadCuts2[e1] && inXGBBkg[indexLow] > subCuts2[e2])
               || (inXGBBkg[indexHigh] < leadCuts2[e1] && inXGBBkg[indexLow] > subCuts3[e2])
               ){
                hBkg2DScoresPass[e]->Fill(outOne,outTwo);
                if(massCalcXGBBkg > 50.0){
                    hBkgPassTightXGBLead[e]->Fill(outOne);
                    hBkgPassTightXGBSub[e]->Fill(outTwo);
                }
                if(massCalcXGBBkg > 95.0){
                    hBkgPassTightXGBTightMassLead[e]->Fill(outOne);
                    hBkgPassTightXGBTightMassSub[e]->Fill(outTwo);
                }
            }
            else{
                hBkgPassLooseXGBLead[e]->Fill(outOne);
                hBkgPassLooseXGBSub[e]->Fill(outTwo);
            }
        }//XGB Eta selection
        
        int iOutHLTBkg,jOutHLTBkg;
        double massCalcHLTBkg = calcMassHLT(&inEtBkg, &inEtaBkg, &inPhiBkg, &inR9Bkg, &inHovrEBkg, &inSigIEtaIEtaBkg, &inPhoIsoBkg, &inEcalIsoBkg, *inNEgsBkg, &iOutHLTBkg, &jOutHLTBkg);
        double drHLTBkg = calcDR(inEtaBkg[iOutHLTBkg], inEtaBkg[jOutHLTBkg], inPhiBkg[iOutHLTBkg], inPhiBkg[jOutHLTBkg]);
        if(inEtBkg[iOutHLTBkg] > 14.25 && inEtBkg[jOutHLTBkg] > 14.25 && abs(inEtaBkg[iOutHLTBkg]) < 2.5 && abs(inEtaBkg[jOutHLTBkg]) < 2.5 && iOutHLTBkg != -1 && jOutHLTBkg != -1 && massCalcHLTBkg > 95.0){
            double outOne,outTwo;

            if(inXGBBkg[iOutHLTBkg] > inXGBBkg[jOutHLTBkg]){
                outOne = inXGBBkg[iOutHLTBkg];
                outTwo = inXGBBkg[jOutHLTBkg];
            }
            if(inXGBBkg[iOutHLTBkg] < inXGBBkg[jOutHLTBkg]){
                outOne = inXGBBkg[jOutHLTBkg];
                outTwo = inXGBBkg[iOutHLTBkg];
            }
            
            int e;
            if (abs(inEtaBkg[iOutHLTBkg]) < 1.5 && abs(inEtaBkg[jOutHLTBkg]) < 1.5)e = 0;
            if ((abs(inEtaBkg[iOutHLTBkg]) > 1.5 && abs(inEtaBkg[jOutHLTBkg]) < 1.5)
                || (abs(inEtaBkg[iOutHLTBkg]) < 1.5 && abs(inEtaBkg[jOutHLTBkg]) > 1.5))e = 1;
            if (abs(inEtaBkg[iOutHLTBkg]) > 1.5 && abs(inEtaBkg[jOutHLTBkg]) > 1.5)e = 2;
            
            //if(*inPassFailStdBkg == 1)hBkg2DScoresPass[e]->Fill(outOne,outTwo);
            
            if(*inPassFailStdBkg == 1 && (*inPassFailSingleL1Bkg == 1 || *inPassFailDoubleL1Bkg == 1)){
                hBkgPassStdLead[e]->Fill(outOne);
                hBkgPassStdSub[e]->Fill(outTwo);
            }
        }//XGB Eta selection
    }
    
    for (int p = 0; p < 3; p++){
        
        hSigAllLead[3]->Add(hSigAllLead[p]);
        hSigAllSub[3]->Add(hSigAllSub[p]);
        hBkgAllLead[3]->Add(hBkgAllLead[p]);
        hBkgAllSub[3]->Add(hBkgAllSub[p]);
        
        hSigPassLooseXGBLead[3]->Add(hSigPassLooseXGBLead[p]);
        hSigPassLooseXGBSub[3]->Add(hSigPassLooseXGBSub[p]);
        hBkgPassLooseXGBLead[3]->Add(hBkgPassLooseXGBLead[p]);
        hBkgPassLooseXGBSub[3]->Add(hBkgPassLooseXGBSub[p]);
        
        hSigPassTightXGBLead[3]->Add(hSigPassTightXGBLead[p]);
        hSigPassTightXGBSub[3]->Add(hSigPassTightXGBSub[p]);
        hBkgPassTightXGBLead[3]->Add(hBkgPassTightXGBLead[p]);
        hBkgPassTightXGBSub[3]->Add(hBkgPassTightXGBSub[p]);
        
        hSigPassTightXGBTightMassLead[3]->Add(hSigPassTightXGBTightMassLead[p]);
        hSigPassTightXGBTightMassSub[3]->Add(hSigPassTightXGBTightMassSub[p]);
        hBkgPassTightXGBTightMassLead[3]->Add(hBkgPassTightXGBTightMassLead[p]);
        hBkgPassTightXGBTightMassSub[3]->Add(hBkgPassTightXGBTightMassSub[p]);
        
        hSigPassStdLead[3]->Add(hSigPassStdLead[p]);
        hSigPassStdSub[3]->Add(hSigPassStdSub[p]);
        hBkgPassStdLead[3]->Add(hBkgPassStdLead[p]);
        hBkgPassStdSub[3]->Add(hBkgPassStdSub[p]);
        
        hSig2DScores[3]->Add(hSig2DScores[p]);
        hSig2DScoresPass[3]->Add(hSig2DScoresPass[p]);
        hBkg2DScores[3]->Add(hBkg2DScores[p]);
        hBkg2DScoresPass[3]->Add(hBkg2DScoresPass[p]);
        
    }
    
    for (int e = 0; e < nEta; e++){
        string outNameLead = outNameGen + "xgbScoresLead" + etaLabels[e];
        string outNameSub = outNameGen + "xgbScoresSub" + etaLabels[e];
        string outName2D = outNameGen + "2DxgbPlot" + etaLabels[e];
        
        string label1A = "All";
        string label1PLoose = "Failing XGBScore Cut";
        string label1PTight = "Passing XGBScore Cut";
        string label1PTightTightMass = label1PTight + " , M > 95";
        string label1PStd = "Passing Std. HLT";
        string label1PStdManual = "Passing Std. HLT (Manual)";
        
        string label2A = "All";

        string label2PLoose = "Failing XGBScore Cut";
        string label2PTight = "Passing XGBScore Cut";
        string label2PTightTightMass = label1PTight + " , M > 95";
        string label2PStd = "Passing Std. HLT";
        string label2PStdManual = "Passing Std. HLT (Manual)";
        
        char buff1[100], buff2[100], buff3[100], buff4[100], buff5[100], buff6[100], buff7[100], buff8[100], buff9[100], buff10[100], buff11[100], buff12[100];
        
        snprintf(buff1, sizeof(buff1), "(%0.1f)", hSigAllLead[e]->GetSumOfWeights());
        string nEvents1 = buff1;
        label1A += nEvents1;
        
        snprintf(buff2, sizeof(buff2), "(%0.1f)", hSigPassLooseXGBLead[e]->GetSumOfWeights());
        string nEvents1PLoose = buff2;
        label1PLoose += nEvents1PLoose;
        
        snprintf(buff3, sizeof(buff3), "(%0.1f)", hSigPassTightXGBLead[e]->GetSumOfWeights());
        string nEvents1PTight = buff3;
        label1PTight += nEvents1PTight;
        
        snprintf(buff4, sizeof(buff4), "(%0.1f)", hSigPassTightXGBTightMassLead[e]->GetSumOfWeights());
        string nEvents1PTightTightMass = buff4;
        label1PTightTightMass += nEvents1PTightTightMass;
        
        snprintf(buff5, sizeof(buff5), "(%0.1f)", hSigPassStdLead[e]->GetSumOfWeights());
        string nEvents1PStd = buff5;
        label1PStd += nEvents1PStd;
        
        snprintf(buff6, sizeof(buff6), "(%0.1f)", hBkgAllLead[e]->GetSumOfWeights());
        string nEvents2 = buff6;
        label2A += nEvents2;
        
        snprintf(buff7, sizeof(buff7), "(%0.1f)", hBkgPassLooseXGBLead[e]->GetSumOfWeights());
        string nEvents2PLoose = buff7;
        label2PLoose += nEvents2PLoose;
        
        snprintf(buff8, sizeof(buff8), "(%0.1f)", hBkgPassTightXGBLead[e]->GetSumOfWeights());
        string nEvents2PTight = buff8;
        label2PTight += nEvents2PTight;
        
        snprintf(buff9, sizeof(buff9), "(%0.1f)", hBkgPassTightXGBTightMassLead[e]->GetSumOfWeights());
        string nEvents2PTightTightMass = buff9;
        label2PTightTightMass += nEvents2PTightTightMass;
        
        snprintf(buff10, sizeof(buff10), "(%0.1f)", hBkgPassStdLead[e]->GetSumOfWeights());
        string nEvents2PStd = buff10;
        label2PStd += nEvents2PStd;
        
        string plotTitle2DSig = genTitleStringSignal + etaLabels[e];
        
        string plotTitleLeadSig = genTitleStringSignal + etaLabels[e] + ";xgbLead";
        THStack *hStackLeadSig = new THStack("hStackLeadSig",plotTitleLeadSig.c_str());
        
        string plotTitleSubSig = genTitleStringSignal + etaLabels[e] + ";xgbSub";
        THStack *hStackSubSig = new THStack("hStackSubSig",plotTitleSubSig.c_str());
        
        string plotTitleLeadBkg = genTitleStringBkg + etaLabels[e] + ";xgbLead";
        THStack *hStackLeadBkg = new THStack("hStackLeadBkg",plotTitleLeadBkg.c_str());
        
        string plotTitleSubBkg = genTitleStringBkg + etaLabels[e] + ";xgbSub";
        THStack *hStackSubBkg = new THStack("hStackSubBkg",plotTitleSubBkg.c_str());
        
        TLegend *legendLeadSig = new TLegend(0.50,0.65,0.90,0.9,"","brNDC");
        TLegend *legendSubSig = new TLegend(0.50,0.65,0.90,0.9,"","brNDC");
        
        TLegend *legendLeadBkg = new TLegend(0.50,0.65,0.90,0.9,"","brNDC");
        TLegend *legendSubBkg = new TLegend(0.50,0.65,0.90,0.9,"","brNDC");
        
        //START WITH LEAD SIGNAL
        
        TH1F *hSigAllLeadDraw = DrawOverflow(hSigAllLead[e]);
        TH1F *hStackHistoSigLead = (TH1F*)hSigAllLeadDraw->Clone();
        hStackHistoSigLead->Reset();
        hStackHistoSigLead->GetXaxis()->SetRange(0,hStackHistoSigLead->GetNbinsX());
        hStackLeadSig->SetHistogram(hStackHistoSigLead);
        
        hStackLeadSig->GetHistogram()->GetXaxis()->SetTitle("xgbLead");
        hStackLeadSig->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
        hStackLeadSig->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
        hStackLeadSig->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
        
        hStackLeadSig->Add(hSigAllLeadDraw);
        legendLeadSig->AddEntry(hSigAllLeadDraw,label1A.c_str(),"pl");
        
        TH1F *hSigPassLooseLeadDraw = DrawOverflow(hSigPassLooseXGBLead[e]);
        hSigPassLooseLeadDraw->SetLineColor(2);
        hStackLeadSig->Add(hSigPassLooseLeadDraw);
        legendLeadSig->AddEntry(hSigPassLooseLeadDraw,label1PLoose.c_str(),"pl");
        
        TH1F *hSigPassStdLeadDraw = DrawOverflow(hSigPassStdLead[e]);
        hSigPassStdLeadDraw->SetLineColor(3);
        hStackLeadSig->Add(hSigPassStdLeadDraw);
        legendLeadSig->AddEntry(hSigPassStdLeadDraw,label1PStd.c_str(),"pl");
        
        TH1F *hSigPassTightLeadDraw = DrawOverflow(hSigPassTightXGBLead[e]);
        hSigPassTightLeadDraw->SetLineColor(4);
        hStackLeadSig->Add(hSigPassTightLeadDraw);
        legendLeadSig->AddEntry(hSigPassTightLeadDraw,label1PTight.c_str(),"pl");
        
        TH1F *hSigPassTightTightMassLeadDraw = DrawOverflow(hSigPassTightXGBTightMassLead[e]);
        hSigPassTightTightMassLeadDraw->SetLineColor(4);
        hStackLeadSig->Add(hSigPassTightTightMassLeadDraw);
        legendLeadSig->AddEntry(hSigPassTightTightMassLeadDraw,label1PTightTightMass.c_str(),"pl");
        
        //THEN SUB SIGNAL
        
        TH1F *hSigAllSubDraw = DrawOverflow(hSigAllSub[e]);
        TH1F *hStackHistoSigSub = (TH1F*)hSigAllSubDraw->Clone();
        hStackHistoSigSub->Reset();
        hStackHistoSigSub->GetXaxis()->SetRange(0,hStackHistoSigSub->GetNbinsX());
        hStackSubSig->SetHistogram(hStackHistoSigSub);
        
        hStackSubSig->GetHistogram()->GetXaxis()->SetTitle("xgbSub");
        hStackSubSig->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
        hStackSubSig->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
        hStackSubSig->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
        
        hStackSubSig->Add(hSigAllSubDraw);
        legendSubSig->AddEntry(hSigAllSubDraw,label1A.c_str(),"pl");
        
        TH1F *hSigPassLooseSubDraw = DrawOverflow(hSigPassLooseXGBSub[e]);
        hSigPassLooseSubDraw->SetLineColor(2);
        hStackSubSig->Add(hSigPassLooseSubDraw);
        legendSubSig->AddEntry(hSigPassLooseSubDraw,label1PLoose.c_str(),"pl");
        
        TH1F *hSigPassStdSubDraw = DrawOverflow(hSigPassStdSub[e]);
        hSigPassStdSubDraw->SetLineColor(3);
        hStackSubSig->Add(hSigPassStdSubDraw);
        legendSubSig->AddEntry(hSigPassStdSubDraw,label1PStd.c_str(),"pl");
        
        TH1F *hSigPassTightSubDraw = DrawOverflow(hSigPassTightXGBSub[e]);
        hSigPassTightSubDraw->SetLineColor(4);
        hStackSubSig->Add(hSigPassTightSubDraw);
        legendSubSig->AddEntry(hSigPassTightSubDraw,label1PTight.c_str(),"pl");
        
        TH1F *hSigPassTightTightMassSubDraw = DrawOverflow(hSigPassTightXGBTightMassSub[e]);
        hSigPassTightTightMassSubDraw->SetLineColor(4);
        hStackSubSig->Add(hSigPassTightTightMassSubDraw);
        legendSubSig->AddEntry(hSigPassTightTightMassSubDraw,label1PTightTightMass.c_str(),"pl");
        
        //THEN LEAD BKG
        
        TH1F *hBkgAllLeadDraw = DrawOverflow(hBkgAllLead[e]);
        TH1F *hStackHistoBkgLead = (TH1F*)hBkgAllLeadDraw->Clone();
        hStackHistoBkgLead->Reset();
        hStackHistoBkgLead->GetXaxis()->SetRange(0,hStackHistoBkgLead->GetNbinsX());
        hStackLeadBkg->SetHistogram(hStackHistoBkgLead);
        
        hStackLeadBkg->GetHistogram()->GetXaxis()->SetTitle("xgbLead");
        hStackLeadBkg->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
        hStackLeadBkg->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
        hStackLeadBkg->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
        
        hStackLeadBkg->Add(hBkgAllLeadDraw);
        legendLeadBkg->AddEntry(hBkgAllLeadDraw,label2A.c_str(),"pl");
        
        TH1F *hBkgPassLooseLeadDraw = DrawOverflow(hBkgPassLooseXGBLead[e]);
        hBkgPassLooseLeadDraw->SetLineColor(2);
        hStackLeadBkg->Add(hBkgPassLooseLeadDraw);
        legendLeadBkg->AddEntry(hBkgPassLooseLeadDraw,label2PLoose.c_str(),"pl");
        
        TH1F *hBkgPassStdLeadDraw = DrawOverflow(hBkgPassStdLead[e]);
        hBkgPassStdLeadDraw->SetLineColor(3);
        hStackLeadBkg->Add(hBkgPassStdLeadDraw);
        legendLeadBkg->AddEntry(hBkgPassStdLeadDraw,label2PStd.c_str(),"pl");
        
        TH1F *hBkgPassTightLeadDraw = DrawOverflow(hBkgPassTightXGBLead[e]);
        hBkgPassTightLeadDraw->SetLineColor(4);
        hStackLeadBkg->Add(hBkgPassTightLeadDraw);
        legendLeadBkg->AddEntry(hBkgPassTightLeadDraw,label2PTight.c_str(),"pl");
        
        TH1F *hBkgPassTightTightMassLeadDraw = DrawOverflow(hBkgPassTightXGBTightMassLead[e]);
        hBkgPassTightTightMassLeadDraw->SetLineColor(4);
        hStackLeadBkg->Add(hBkgPassTightTightMassLeadDraw);
        legendLeadBkg->AddEntry(hBkgPassTightTightMassLeadDraw,label2PTightTightMass.c_str(),"pl");
        
        //THEN SUB BKG
        
        TH1F *hBkgAllSubDraw = DrawOverflow(hBkgAllSub[e]);
        TH1F *hStackHistoBkgSub = (TH1F*)hBkgAllSubDraw->Clone();
        hStackHistoBkgSub->Reset();
        hStackHistoBkgSub->GetXaxis()->SetRange(0,hStackHistoBkgSub->GetNbinsX());
        hStackSubBkg->SetHistogram(hStackHistoBkgSub);
        
        hStackSubBkg->GetHistogram()->GetXaxis()->SetTitle("xgbSub");
        hStackSubBkg->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
        hStackSubBkg->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
        hStackSubBkg->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
        
        hStackSubBkg->Add(hBkgAllSubDraw);
        legendSubBkg->AddEntry(hBkgAllSubDraw,label2A.c_str(),"pl");
        
        TH1F *hBkgPassLooseSubDraw = DrawOverflow(hBkgPassLooseXGBSub[e]);
        hBkgPassLooseSubDraw->SetLineColor(2);
        hStackSubBkg->Add(hBkgPassLooseSubDraw);
        legendSubBkg->AddEntry(hBkgPassLooseSubDraw,label2PLoose.c_str(),"pl");
        
        TH1F *hBkgPassStdSubDraw = DrawOverflow(hBkgPassStdSub[e]);
        hBkgPassStdSubDraw->SetLineColor(3);
        hStackSubBkg->Add(hBkgPassStdSubDraw);
        legendSubBkg->AddEntry(hBkgPassStdSubDraw,label2PStd.c_str(),"pl");
        
        TH1F *hBkgPassTightSubDraw = DrawOverflow(hBkgPassTightXGBSub[e]);
        hBkgPassTightSubDraw->SetLineColor(4);
        hStackSubBkg->Add(hBkgPassTightSubDraw);
        legendSubBkg->AddEntry(hBkgPassTightSubDraw,label2PTight.c_str(),"pl");
        
        TH1F *hBkgPassTightTightMassSubDraw = DrawOverflow(hBkgPassTightXGBTightMassSub[e]);
        hBkgPassTightTightMassSubDraw->SetLineColor(4);
        hStackSubBkg->Add(hBkgPassTightTightMassSubDraw);
        legendSubBkg->AddEntry(hBkgPassTightTightMassSubDraw,label2PTightTightMass.c_str(),"pl");
        
        can->Clear();
        can->Divide(2,1);
        can->cd(1);
        gPad->SetGrid();
        hStackLeadSig->Draw("nostackhist");
        legendLeadSig->Draw("same");
        can->cd(2);
        gPad->SetGrid();
        hStackLeadBkg->Draw("nostackhist");
        legendLeadBkg->Draw("same");
        can->SaveAs((outLead + etaFLabels[e] +".png").c_str());
        can->SaveAs((outLead + etaFLabels[e] +".root").c_str());
        
        can->Clear();
        can->Divide(2,1);
        can->cd(1);
        gPad->SetGrid();
        hStackSubSig->Draw("nostackhist");
        legendSubSig->Draw("same");
        can->cd(2);
        gPad->SetGrid();
        hStackSubBkg->Draw("nostackhist");
        legendSubBkg->Draw("same");
        can->SaveAs((outSub + etaFLabels[e] +".png").c_str());
        can->SaveAs((outSub + etaFLabels[e] +".root").c_str());
    
        can->Clear();
        can->Divide(2,1);
        can->cd(1);
        gPad->SetGrid();
        hSig2DScores[e]->SetTitle((genTitleStringSignal + etaLabels[e]).c_str());
        hSig2DScores[e]->GetXaxis()->SetTitle("xgbHigh");
        hSig2DScores[e]->GetYaxis()->SetTitle("xgbLow");
        hSig2DScores[e]->Draw("box");
        hSig2DScoresPass[e]->SetLineColor(2);
        hSig2DScoresPass[e]->Draw("boxsame");
        can->cd(2);
        gPad->SetGrid();
        hBkg2DScores[e]->SetTitle((genTitleStringBkg + etaLabels[e]).c_str());
        hBkg2DScores[e]->GetXaxis()->SetTitle("xgbHigh");
        hBkg2DScores[e]->GetYaxis()->SetTitle("xgbLow");
        hBkg2DScores[e]->Draw("box");
        hBkg2DScoresPass[e]->SetLineColor(2);
        hBkg2DScoresPass[e]->Draw("boxsame");
        can->SaveAs((out2D + etaFLabels[e] +".png").c_str());
        can->SaveAs((out2D + etaFLabels[e] +".root").c_str());
    }
}
