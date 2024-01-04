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
#include <vector>
#include <math.h>
string ToStringWithPrecision(double inVal, const int n = 2);
TH1F *DrawOverflow(TH1F* h);
double calcDR(float eta1, float eta2, float phi1, float phi2);
double calcMassXGB(TTreeReaderArray<float>* xgbScore,double xgbScoreB, double xgbScoreE,TTreeReaderArray<float>* et, TTreeReaderArray<float>* eta, TTreeReaderArray<float>* phi, int nEgsIn, int* pho1, int* pho2);
double calcMassHLT(TTreeReaderArray<float>* et, TTreeReaderArray<float>* eta, TTreeReaderArray<float>* phi, TTreeReaderArray<float>* r9, TTreeReaderArray<float>* hOvrE, TTreeReaderArray<float>* sIeIe, TTreeReaderArray<float>* phoIso, TTreeReaderArray<float>* ecalIso, int nEgsIn, int* pho1, int* pho2);
void oldDiphotonSigBkgComp(){
    gROOT->Reset();
    gStyle->SetPalette(1);
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitle(0);
    
    string fileName = "validationNTuples/0606/GGH13andData_AllPhotonsValidated.root";
    string genTitleStringSignal = "GGH Signal ";
    string genTitleStringBkg = "Data ";
    string dirStr ="sigBkgVarPlots/0621/TESTOLDCODE/";
    string plotType = "Added";
    //string plotType = "Separate";
    //string dirStr ="varPlots/0130/DataRelaxedUnseededChoose2/";
        
    int nEta;
    
    nEta = 4;
    string etaLabels[4] = {"Barrel-Barrel","Barrel-Endcap","Endcap-Endcap","All #eta"};
    string etaFLabels[4] = {"_BB","_BE","_EE","_All"};
    
    double looseCut = 0.1;
    double tightCutBB = 0.25;
    double tightCutBE = 0.09182;
    double tightCutEE = 0.12375;
    double tightCuts[4] = {tightCutBB,tightCutBE,tightCutEE,-1.0};
    
    string looseCutString = ToStringWithPrecision(looseCut,2);
    string tightCutStringBB = ToStringWithPrecision(tightCutBB,2);
    string tightCutStringBE = ToStringWithPrecision(tightCutBE,4);
    string tightCutStringEE = ToStringWithPrecision(tightCutEE,4);
    string tightCutStrings[4] = {tightCutStringBB,tightCutStringBE,tightCutStringEE,"Mixed"};
    
    string outNameGen = dirStr;
    
    int nVars = 18;
    string varNames[18] = {"nEgs","triggerBits","mass","dR","et","energy","rawEnergy","eta","etaWidth","phi","phiWidth","sigmaIEtaIEta","ecalPFIso","r9HLT","s4","hOvrE","trkIsoPho","xgbScore"};
    double limsLow[18] = {1,0.0,0.0,0.0,0.0,0.0,0.0,-3,0.0,-3.15,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    //double limsHigh[16] = {3,4,500,300,600,600,3.15,0.1,3.15,0.2,0.2,20.0,1.2,1.2,1.0,8.0};
    double limsHigh[18] = {10,4,500,6,300,600,600,3,0.1,3.15,0.2,0.1,40.0,1.2,1.2,1.0,40.0,1.0};
    //double limsHigh[16] = {3,4,500,300,300,600,3.15,0.1,3.15,0.2,2.0,200.0,5.0,1.2,2.0,200.0};
    
    int nBins[18] = {9,4,500,600,600,600,600,200,400,400,400,400,400,400,400,400,400,400};
    
    TCanvas *can = new TCanvas ("can","can",10,10,1600,900);
    //for(int i = 10; i < 11; i++){
    for(int i = 0; i < nVars; i++){
    //for(int i = 2; i < 4; i++){
    //for(int i = 0; i < 1; i++){
        // for(int i = 14; i < nVars; i++){
        
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
        
        TH1F *hSigPassStdManualLead[4];
        TH1F *hSigPassStdManualSub[4];
        TH1F *hBkgPassStdManualLead[4];
        TH1F *hBkgPassStdManualSub[4];
        
        for (int e = 0; e < nEta; e++){
            hSigAllLead[e] = new TH1F ("hSigAllLead","",nBins[i],limsLow[i],limsHigh[i]);
            hSigAllSub[e] = new TH1F ("hSigAllSub","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgAllLead[e] = new TH1F ("hBkgAllLead","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgAllSub[e] = new TH1F ("hBkgAllSub","",nBins[i],limsLow[i],limsHigh[i]);

            hSigPassLooseXGBLead[e] = new TH1F ("hSigPassLooseXGBLead","",nBins[i],limsLow[i],limsHigh[i]);
            hSigPassLooseXGBSub[e] = new TH1F ("hSigPassLooseXGBSub","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassLooseXGBLead[e] = new TH1F ("hBkgPassLooseXGBLead","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassLooseXGBSub[e] = new TH1F ("hBkgPassLooseXGBSub","",nBins[i],limsLow[i],limsHigh[i]);

            hSigPassTightXGBLead[e] = new TH1F ("hSigPassTightXGBLead","",nBins[i],limsLow[i],limsHigh[i]);
            hSigPassTightXGBSub[e] = new TH1F ("hSigPassTightXGBSub","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassTightXGBLead[e] = new TH1F ("hBkgPassTightXGBLead","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassTightXGBSub[e] = new TH1F ("hBkgPassTightXGBSub","",nBins[i],limsLow[i],limsHigh[i]);
            
            hSigPassTightXGBTightMassLead[e] = new TH1F ("hSigPassTightXGBTightMassLead","",nBins[i],limsLow[i],limsHigh[i]);
            hSigPassTightXGBTightMassSub[e] = new TH1F ("hSigPassTightXGBTightMassSub","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassTightXGBTightMassLead[e] = new TH1F ("hBkgPassTightXGBTightMassLead","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassTightXGBTightMassSub[e] = new TH1F ("hBkgPassTightXGBTightMassSub","",nBins[i],limsLow[i],limsHigh[i]);
           
            hSigPassStdLead[e] = new TH1F ("hSigPassStdLead","",nBins[i],limsLow[i],limsHigh[i]);
            hSigPassStdSub[e] = new TH1F ("hSigPassStdSub","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassStdLead[e] = new TH1F ("hBkgPassStdLead","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassStdSub[e] = new TH1F ("hBkgPassStdSub","",nBins[i],limsLow[i],limsHigh[i]);

            hSigPassStdManualLead[e] = new TH1F ("hSigPassStdManualLead","",nBins[i],limsLow[i],limsHigh[i]);
            hSigPassStdManualSub[e] = new TH1F ("hSigPassStdManualSub","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassStdManualLead[e] = new TH1F ("hBkgPassStdManualLead","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassStdManualSub[e] = new TH1F ("hBkgPassStdManualSub","",nBins[i],limsLow[i],limsHigh[i]);

        }
        
        int nPassEta = 0;
        int nFailEta = 0;
        int nPassEtaSig = 0;
        int nFailEtaSig = 0;
        
        for (int e = 0; e < nEta; e++){
            
            string outName = outNameGen + varNames[i] + etaFLabels[e];
            
            TFile *f = new TFile(fileName.c_str());
            
            string varNameTmp;
            if (e != 3){
                varNameTmp = varNames[i];
                string varNameTmpHere = varNameTmp;
                if (varNameTmp == "nEgs" || varNameTmp == "triggerBits" || varNameTmp == "mass" || varNameTmp == "dR")varNameTmpHere = "eta";
                
                //FIRST handle Signal
                TTreeReader inReaderSig("sigTree", f);
                inReaderSig.Restart();
                
                TTreeReaderArray<Float_t> inVarSig(inReaderSig, varNameTmpHere.c_str());
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

                TTreeReaderValue<int> inNEgsSig(inReaderSig,"nEgs");
                TTreeReaderValue<int> inTriggerBitsSig(inReaderSig,"triggerBits");
                TTreeReaderValue<int> inPassFailStdSig(inReaderSig,"passFailStd");
                TTreeReaderValue<int> inPassFailDoubleL1Sig(inReaderSig,"passFailL1Double");
                TTreeReaderValue<int> inPassFailSingleL1Sig(inReaderSig,"passFailL1Single");
                
                //NORMAL ETA THRESHOLDS ARE 1.444 and 1.556
                int nTotalSig = 0;
                int nDisagreeSig = 0;
                int nNoneChosenSig = 0;
                while (inReaderSig.Next()) {
                    int iOutXGBSig,jOutXGBSig;
                    double massCalcXGBSig = calcMassXGB(&inXGBSig,0.25,0.25,&inEtSig,&inEtaSig,&inPhiSig,*inNEgsSig,&iOutXGBSig,&jOutXGBSig);
                    double drXGBSig = calcDR(inEtaSig[iOutXGBSig], inEtaSig[jOutXGBSig], inPhiSig[iOutXGBSig], inPhiSig[jOutXGBSig]);
                    if (e == 0 && (iOutXGBSig == -1 || jOutXGBSig == -1))nNoneChosenSig+=1;
                    if(inEtSig[iOutXGBSig] > 14.25 && inEtSig[jOutXGBSig] > 14.25 && abs(inEtaSig[iOutXGBSig]) < 2.5 && abs(inEtaSig[jOutXGBSig]) < 2.5 && iOutXGBSig != -1 && jOutXGBSig != -1 && ((e == 0 && abs(inEtaSig[iOutXGBSig]) < 1.5 && abs(inEtaSig[jOutXGBSig]) < 1.5)
                       ||(e == 1 && ((abs(inEtaSig[iOutXGBSig]) < 1.5 && abs(inEtaSig[jOutXGBSig]) > 1.5) || (abs(inEtaSig[iOutXGBSig]) > 1.5 && abs(inEtaSig[jOutXGBSig]) < 1.5) ))
                       ||(e == 2 && abs(inEtaSig[iOutXGBSig]) > 1.5 && abs(inEtaSig[jOutXGBSig]) > 1.5))){
                        //cout<<"For XGB, (i,j) = ("<<iOutXGBSig<<","<<jOutXGBSig<<"), mass = "<<massCalcXGBSig<<endl;
                        double outOne,outTwo;
                        if (varNameTmp != "mass" && varNameTmp != "nEgs" && varNameTmp != "triggerBits" && varNameTmp != "dR"){
                            outOne = inVarSig[iOutXGBSig];
                            outTwo = inVarSig[jOutXGBSig];
                        }
                        else if (varNameTmp == "mass"){
                            outOne = massCalcXGBSig;
                            outTwo = massCalcXGBSig;
                        }
                        else if (varNameTmp == "nEgs"){
                            outOne = (double)(*inNEgsSig);
                            outTwo = (double)(*inNEgsSig);
                            
                        }
                        else if (varNameTmp == "triggerBits"){
                            outOne = (double)(*inTriggerBitsSig);
                            outTwo = (double)(*inTriggerBitsSig);
                            
                        }
                        else if (varNameTmp == "dR"){
                            outOne = drXGBSig;
                            outTwo = drXGBSig;
                            
                        }
                        hSigAllLead[e]->Fill(outOne);
                        hSigAllSub[e]->Fill(outTwo);
                        
                        //REPURPOSE THE LOOSE XGB CUT PLOT FOR FAILING THE TIGHT CUT. TEMPORARY
                        //if(inXGBSig[iOutXGBSig] > looseCut && inXGBSig[jOutXGBSig] > looseCut){
                        if(inXGBSig[iOutXGBSig] < tightCuts[e] || inXGBSig[jOutXGBSig] < tightCuts[e]){
                            hSigPassLooseXGBLead[e]->Fill(outOne);
                            hSigPassLooseXGBSub[e]->Fill(outTwo);
                        }
                        
                        if(inXGBSig[iOutXGBSig] > tightCuts[e] && inXGBSig[jOutXGBSig] > tightCuts[e]){
                            hSigPassTightXGBLead[e]->Fill(outOne);
                            hSigPassTightXGBSub[e]->Fill(outTwo);
                            if (massCalcXGBSig > 95.0){
                                hSigPassTightXGBTightMassLead[e]->Fill(outOne);
                                hSigPassTightXGBTightMassSub[e]->Fill(outTwo);
                            }
                        }
                    }//XGB Eta selection
                    
                    int iOutHLTSig,jOutHLTSig;
                    double massCalcHLTSig = calcMassHLT(&inEtSig, &inEtaSig, &inPhiSig, &inR9Sig, &inHovrESig, &inSigIEtaIEtaSig, &inPhoIsoSig, &inEcalIsoSig, *inNEgsSig, &iOutHLTSig, &jOutHLTSig);
                    double drHLTSig = calcDR(inEtaSig[iOutHLTSig], inEtaSig[jOutHLTSig], inPhiSig[iOutHLTSig], inPhiSig[jOutHLTSig]);
                    if(iOutHLTSig != -1 && jOutHLTSig != -1 && ((e == 0 && abs(inEtaSig[iOutHLTSig]) < 1.5 && abs(inEtaSig[jOutHLTSig]) < 1.5)||(e == 1 && ((abs(inEtaSig[iOutHLTSig]) < 1.5 && abs(inEtaSig[jOutHLTSig]) > 1.5) || (abs(inEtaSig[iOutHLTSig]) > 1.5 && abs(inEtaSig[jOutHLTSig]) < 1.5)))||(e == 2 && abs(inEtaSig[iOutHLTSig]) > 1.5 && abs(inEtaSig[jOutHLTSig]) > 1.5))){
                        //cout<<"For HLT, (i,j) = ("<<iOutHLTSig<<","<<jOutHLTSig<<"), mass = "<<massCalcHLTSig<<endl;
                        double outOne,outTwo;
                        if (varNameTmp != "mass" && varNameTmp != "nEgs" && varNameTmp != "triggerBits" && varNameTmp != "dR"){
                            outOne = inVarSig[iOutHLTSig];
                            outTwo = inVarSig[jOutHLTSig];
                        }
                        else if (varNameTmp == "mass"){
                            outOne = massCalcHLTSig;
                            outTwo = massCalcHLTSig;
                        }
                        else if (varNameTmp == "nEgs"){
                            outOne = (double)(*inNEgsSig);
                            outTwo = (double)(*inNEgsSig);
                            
                        }
                        else if (varNameTmp == "triggerBits"){
                            outOne = (double)(*inTriggerBitsSig);
                            outTwo = (double)(*inTriggerBitsSig);
                            
                        }
                        else if (varNameTmp == "dR"){
                            outOne = drHLTSig;
                            outTwo = drHLTSig;
                        }
                        if(*inPassFailStdSig == 1 && (*inPassFailSingleL1Sig == 1 || *inPassFailDoubleL1Sig == 1)){
                            hSigPassStdLead[e]->Fill(outOne);
                            hSigPassStdSub[e]->Fill(outTwo);
                        }
                        if(inEtSig[iOutHLTSig] > 30 && inEtSig[jOutHLTSig] > 22 && massCalcHLTSig > 95.0){
                        //if(inEtSig[iOutXGBSig] > 30 && inEtSig[jOutXGBSig] > 22){
                            if ((inEtaSig[iOutHLTSig] < 1.5 && inR9Sig[iOutHLTSig] > 0.5 &&  inHovrESig[iOutHLTSig] < 0.12 && !(inR9Sig[iOutHLTSig] < 0.85 && (inSigIEtaIEtaSig[iOutHLTSig] > 0.015 || inPhoIsoSig[iOutHLTSig] > 6)))
                            || (inEtaSig[iOutHLTSig] > 1.5 && inR9Sig[iOutHLTSig] > 0.8 &&  inHovrESig[iOutHLTSig] < 0.10 && !(inR9Sig[iOutHLTSig] < 0.90 && (inSigIEtaIEtaSig[iOutHLTSig] > 0.035 || inPhoIsoSig[iOutHLTSig] > 6)))){

                                if ((inEtaSig[jOutHLTSig] < 1.5 && inR9Sig[jOutHLTSig] > 0.5 &&  inHovrESig[jOutHLTSig] < 0.12 && !(inR9Sig[jOutHLTSig] < 0.85 && (inSigIEtaIEtaSig[jOutHLTSig] > 0.015 || inPhoIsoSig[jOutHLTSig] > 6)))
                                || (inEtaSig[jOutHLTSig] > 1.5 && inR9Sig[jOutHLTSig] > 0.8 &&  inHovrESig[jOutHLTSig] < 0.10 && !(inR9Sig[jOutHLTSig] < 0.90 && (inSigIEtaIEtaSig[jOutHLTSig] > 0.035 || inPhoIsoSig[jOutHLTSig] > 6)))){
                                    hSigPassStdManualLead[e]->Fill(outOne);
                                    hSigPassStdManualSub[e]->Fill(outTwo);
                                }//second photon cuts
                            }//First photon cuts
                        }//Et and Mass Cuts
                    }//HLT Eta Selection if
                    nTotalSig+=1;
                    if (iOutHLTSig != iOutXGBSig || jOutHLTSig != jOutXGBSig)nDisagreeSig+=1;
                }//While inReaderSig.Next()
                
                
                //THEN handle Bkg
                TTreeReader inReaderBkg("bkgTree", f);
                inReaderBkg.Restart();
                
                TTreeReaderArray<Float_t> inVarBkg(inReaderBkg, varNameTmpHere.c_str());
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
                
                int nTotalBkg = 0;
                int nDisagreeBkg = 0;
                int nNoneChosenBkg = 0;
                while (inReaderBkg.Next()) {
                    int iOutXGBBkg,jOutXGBBkg;
                    double massCalcXGBBkg = calcMassXGB(&inXGBBkg,0.25,0.25,&inEtBkg,&inEtaBkg,&inPhiBkg,*inNEgsBkg,&iOutXGBBkg,&jOutXGBBkg);
                    double drXGBBkg = calcDR(inEtaBkg[iOutXGBBkg], inEtaBkg[jOutXGBBkg], inPhiBkg[iOutXGBBkg], inPhiBkg[jOutXGBBkg]);
                    if (e == 0 && (iOutXGBBkg == -1 || jOutXGBBkg == -1))nNoneChosenBkg+=1;
                    if(inEtBkg[iOutXGBBkg] > 14.25 && inEtBkg[jOutXGBBkg] > 14.25 && abs(inEtaBkg[iOutXGBBkg]) < 2.5 && abs(inEtaBkg[jOutXGBBkg]) < 2.5 && iOutXGBBkg != -1 && jOutXGBBkg != -1 && ((e == 0 && abs(inEtaBkg[iOutXGBBkg]) < 1.5 && abs(inEtaBkg[jOutXGBBkg]) < 1.5)||(e == 1 && ((abs(inEtaBkg[iOutXGBBkg]) < 1.5 && abs(inEtaBkg[jOutXGBBkg]) > 1.5) || (abs(inEtaBkg[iOutXGBBkg]) > 1.5 && abs(inEtaBkg[jOutXGBBkg]) < 1.5) ))||(e == 2 && abs(inEtaBkg[iOutXGBBkg]) > 1.5 && abs(inEtaBkg[jOutXGBBkg]) > 1.5))){
                        //cout<<"i,j = "<<iOutXGBBkg<<" , "<<jOutXGBBkg<<endl;
                        double outOne,outTwo;
                        if (varNameTmp != "mass" && varNameTmp != "nEgs" && varNameTmp != "triggerBits" && varNameTmp != "dR"){
                            outOne = inVarBkg[iOutXGBBkg];
                            outTwo = inVarBkg[jOutXGBBkg];
                        }
                        else if (varNameTmp == "mass"){
                            outOne = massCalcXGBBkg;
                            outTwo = massCalcXGBBkg;
                        }
                        else if (varNameTmp == "nEgs"){
                            outOne = (double)(*inNEgsBkg);
                            outTwo = (double)(*inNEgsBkg);
                            
                        }
                        else if (varNameTmp == "triggerBits"){
                            outOne = (double)(*inTriggerBitsBkg);
                            outTwo = (double)(*inTriggerBitsBkg);
                        }
                        else if (varNameTmp == "dR"){
                            outOne = drXGBBkg;
                            outTwo = drXGBBkg;
                        }
                        hBkgAllLead[e]->Fill(outOne);
                        hBkgAllSub[e]->Fill(outTwo);
                        
                        //REPURPOSE THE LOOSE XGB CUT PLOT FOR FAILING THE TIGHT CUT. TEMPORARY
                        //if(inXGBBkg[iOutXGBBkg] > looseCut && inXGBBkg[jOutXGBBkg] > looseCut){
                        if(inXGBBkg[iOutXGBBkg] < tightCuts[e] || inXGBBkg[jOutXGBBkg] < tightCuts[e]){
                            hBkgPassLooseXGBLead[e]->Fill(outOne);
                            hBkgPassLooseXGBSub[e]->Fill(outTwo);
                        }
                        
                        if(inXGBBkg[iOutXGBBkg] > tightCuts[e] && inXGBBkg[jOutXGBBkg] > tightCuts[e]){
                            hBkgPassTightXGBLead[e]->Fill(outOne);
                            hBkgPassTightXGBSub[e]->Fill(outTwo);
                            if (massCalcXGBBkg > 95.0){
                                hBkgPassTightXGBTightMassLead[e]->Fill(outOne);
                                hBkgPassTightXGBTightMassSub[e]->Fill(outTwo);
                            }
                        }
                    }//XGB Eta selection
                    int iOutHLTBkg,jOutHLTBkg;
                    double massCalcHLTBkg = calcMassHLT(&inEtBkg, &inEtaBkg, &inPhiBkg, &inR9Bkg, &inHovrEBkg, &inSigIEtaIEtaBkg, &inPhoIsoBkg, &inEcalIsoBkg, *inNEgsBkg, &iOutHLTBkg, &jOutHLTBkg);
                    double drHLTBkg = calcDR(inEtaBkg[iOutHLTBkg], inEtaBkg[jOutHLTBkg], inPhiBkg[iOutHLTBkg], inPhiBkg[jOutHLTBkg]);
                    if(iOutHLTBkg != -1 && jOutHLTBkg != -1 && ((e == 0 && abs(inEtaBkg[iOutHLTBkg]) < 1.5 && abs(inEtaBkg[jOutHLTBkg]) < 1.5)||(e == 1 && ((abs(inEtaBkg[iOutHLTBkg]) < 1.5 && abs(inEtaBkg[jOutHLTBkg]) > 1.5) || (abs(inEtaBkg[iOutHLTBkg]) > 1.5 && abs(inEtaBkg[jOutHLTBkg]) < 1.5)))||(e == 2 && abs(inEtaBkg[iOutHLTBkg]) > 1.5 && abs(inEtaBkg[jOutHLTBkg]) > 1.5))){
                        double outOne,outTwo;
                        if (varNameTmp != "mass" && varNameTmp != "nEgs" && varNameTmp != "triggerBits" && varNameTmp != "dR"){
                            outOne = inVarBkg[iOutHLTBkg];
                            outTwo = inVarBkg[jOutHLTBkg];
                        }
                        else if (varNameTmp == "mass"){
                            outOne = massCalcHLTBkg;
                            outTwo = massCalcHLTBkg;
                        }

                        else if (varNameTmp == "nEgs"){
                            outOne = (double)(*inNEgsBkg);
                            outTwo = (double)(*inNEgsBkg);
                        }

                        else if (varNameTmp == "triggerBits"){
                            outOne = (double)(*inTriggerBitsBkg);
                            outTwo = (double)(*inTriggerBitsBkg);
                        }
                        else if (varNameTmp == "dR"){
                            outOne = drHLTBkg;
                            outTwo = drHLTBkg;
                            
                        }
                        if(*inPassFailStdBkg == 1 && (*inPassFailSingleL1Bkg == 1 || *inPassFailDoubleL1Bkg == 1)){
                            hBkgPassStdLead[e]->Fill(outOne);
                            hBkgPassStdSub[e]->Fill(outTwo);
                        }

                        if(inEtBkg[iOutHLTBkg] > 30 && inEtBkg[jOutHLTBkg] > 22 && massCalcHLTBkg > 95.0){
                            //if(inEtBkg[iOutXGBBkg] > 30 && inEtBkg[jOutXGBBkg] > 22){
                            if ((inEtaBkg[iOutHLTBkg] < 1.5 && inR9Bkg[iOutHLTBkg] > 0.5 &&  inHovrEBkg[iOutHLTBkg] < 0.12 && !(inR9Bkg[iOutHLTBkg] < 0.85 && (inSigIEtaIEtaBkg[iOutHLTBkg] > 0.015 || inPhoIsoBkg[iOutHLTBkg] > 6)))
                                || (inEtaBkg[iOutHLTBkg] > 1.5 && inR9Bkg[iOutHLTBkg] > 0.8 &&  inHovrEBkg[iOutHLTBkg] < 0.10 && !(inR9Bkg[iOutHLTBkg] < 0.90 && (inSigIEtaIEtaBkg[iOutHLTBkg] > 0.035 || inPhoIsoBkg[iOutHLTBkg] > 6)))){
                                
                                if ((inEtaBkg[jOutHLTBkg] < 1.5 && inR9Bkg[jOutHLTBkg] > 0.5 &&  inHovrEBkg[jOutHLTBkg] < 0.12 && !(inR9Bkg[jOutHLTBkg] < 0.85 && (inSigIEtaIEtaBkg[jOutHLTBkg] > 0.015 || inPhoIsoBkg[jOutHLTBkg] > 6)))
                                    || (inEtaBkg[jOutHLTBkg] > 1.5 && inR9Bkg[jOutHLTBkg] > 0.8 &&  inHovrEBkg[jOutHLTBkg] < 0.10 && !(inR9Bkg[jOutHLTBkg] < 0.90 && (inSigIEtaIEtaBkg[jOutHLTBkg] > 0.035 || inPhoIsoBkg[jOutHLTBkg] > 6)))){
                                    hBkgPassStdManualLead[e]->Fill(outOne);
                                    hBkgPassStdManualSub[e]->Fill(outTwo);
                                }//second photon cuts
                            }//First photon cuts
                        }//Et and Mass Cuts
                    }//HLT Eta Selection if
                    nTotalBkg+=1;
                    if (iOutHLTBkg != iOutXGBBkg || jOutHLTBkg != jOutXGBBkg)nDisagreeBkg+=1;
                }//While inReaderBkg.Next()
                cout<<"e = "<<e<<endl;
                cout<<"nTotalSig = "<<nTotalSig<<", nDisagreeSig = "<<nDisagreeSig<<endl;
                cout<<"nTotalBkg = "<<nTotalBkg<<", nDisagreeBkg = "<<nDisagreeBkg<<endl;
                cout<<"sig N where NONE chosen = "<<nNoneChosenSig<<endl;
                cout<<"bkg N where NONE chosen = "<<nNoneChosenBkg<<endl;
            }//if e != 3
            if (e == 3){
                for (int p = 0; p < 3; p++){
                   
                    hSigAllLead[e]->Add(hSigAllLead[p]);
                    hSigAllSub[e]->Add(hSigAllSub[p]);
                    hBkgAllLead[e]->Add(hBkgAllLead[p]);
                    hBkgAllSub[e]->Add(hBkgAllSub[p]);
                    
                    hSigPassLooseXGBLead[e]->Add(hSigPassLooseXGBLead[p]);
                    hSigPassLooseXGBSub[e]->Add(hSigPassLooseXGBSub[p]);
                    hBkgPassLooseXGBLead[e]->Add(hBkgPassLooseXGBLead[p]);
                    hBkgPassLooseXGBSub[e]->Add(hBkgPassLooseXGBSub[p]);
                        
                    hSigPassTightXGBLead[e]->Add(hSigPassTightXGBLead[p]);
                    hSigPassTightXGBSub[e]->Add(hSigPassTightXGBSub[p]);
                    hBkgPassTightXGBLead[e]->Add(hBkgPassTightXGBLead[p]);
                    hBkgPassTightXGBSub[e]->Add(hBkgPassTightXGBSub[p]);
                    
                    hSigPassTightXGBTightMassLead[e]->Add(hSigPassTightXGBTightMassLead[p]);
                    hSigPassTightXGBTightMassSub[e]->Add(hSigPassTightXGBTightMassSub[p]);
                    hBkgPassTightXGBTightMassLead[e]->Add(hBkgPassTightXGBTightMassLead[p]);
                    hBkgPassTightXGBTightMassSub[e]->Add(hBkgPassTightXGBTightMassSub[p]);

                    hSigPassStdLead[e]->Add(hSigPassStdLead[p]);
                    hSigPassStdSub[e]->Add(hSigPassStdSub[p]);
                    hBkgPassStdLead[e]->Add(hBkgPassStdLead[p]);
                    hBkgPassStdSub[e]->Add(hBkgPassStdSub[p]);

                    hSigPassStdManualLead[e]->Add(hSigPassStdManualLead[p]);
                    hSigPassStdManualSub[e]->Add(hSigPassStdManualSub[p]);
                    hBkgPassStdManualLead[e]->Add(hBkgPassStdManualLead[p]);
                    hBkgPassStdManualSub[e]->Add(hBkgPassStdManualSub[p]);
                }
            }//if e == 3
            string varNameTmp2;
            if (varNames[i] == "mass")varNameTmp2 = "M_{#gamma#gamma}";
            else if (varNames[i] == "r9HLT")varNameTmp2 = "R_{9}";
            else if (varNames[i] == "sigmaIEtaIEta")varNameTmp2 = "#sigma_{i#eta i#eta}";
            else if (varNames[i] == "etaWidth")varNameTmp2 = "#eta Width";
            else if (varNames[i] == "phiWidth")varNameTmp2 = "#phi Width";
            else if (varNames[i] == "s4")varNameTmp2 = "S_{4}";
            else if (varNames[i] == "eta")varNameTmp2 = "#eta_{SC}";
            else if (varNames[i] == "et")varNameTmp2 = "E_{T}";
            else if (varNames[i] == "dR")varNameTmp2 = "#Delta r";
            else varNameTmp2 = varNames[i];
            
            can->Clear();
                
            can->Divide(2,1);
            
            string plotTitle1 = genTitleStringSignal + etaLabels[e] + ";" + varNameTmp2;
            THStack *hStack1 = new THStack("hStack1",plotTitle1.c_str());
                
            string plotTitle2 = genTitleStringBkg + etaLabels[e] + ";" + varNameTmp2;
            THStack *hStack2 = new THStack("hStack2",plotTitle2.c_str());
                
            TLegend *legend1;
            if (varNames[i] == "r9HLT"|| varNames[i] == "s4") legend1 = new TLegend(0.1,0.65,0.5,0.9,"","brNDC");
            else legend1 = new TLegend(0.50,0.65,0.90,0.9,"","brNDC");
                
            TLegend *legend2;
            if (varNames[i] == "r9HLT"|| varNames[i] == "s4") legend2 = new TLegend(0.1,0.65,0.5,0.9,"","brNDC");
            else legend2 = new TLegend(0.50,0.65,0.90,0.9,"","brNDC");
            
            string label1A = "All";
            //string label1PLoose = "XGBScore > " + looseCutString;
            string label1PLoose = "XGBScore < " + tightCutStrings[e];
            string label1PTight = "XGBScore > " + tightCutStrings[e];
            string label1PTightTightMass = "XGBScore > " + tightCutStrings[e] + " , M > 95";
            string label1PStd = "Passing Std. HLT";
            string label1PStdManual = "Passing Std. HLT (Manual)";
            
            string label2A = "All";
//            string label2PLoose = "XGBScore > " + looseCutString;
            string label2PLoose = "XGBScore < " + tightCutStrings[e];
            string label2PTight = "XGBScore > " + tightCutStrings[e];
            string label2PTightTightMass = "XGBScore > " + tightCutStrings[e] + " , M > 95";
            string label2PStd = "Passing Std. HLT";
            string label2PStdManual = "Passing Std. HLT (Manual)";
                        
            if (plotType == "Added"){
                            
                TH1F *hSigAll = (TH1F*)hSigAllLead[e]->Clone();
                hSigAll->Add(hSigAllSub[e]);
                TH1F *hBkgAll = (TH1F*)hBkgAllLead[e]->Clone();
                hBkgAll->Add(hBkgAllSub[e]);

                TH1F *hSigPassLooseXGB = (TH1F*)hSigPassLooseXGBLead[e]->Clone();
                hSigPassLooseXGB->Add(hSigPassLooseXGBSub[e]);
                TH1F *hBkgPassLooseXGB = (TH1F*)hBkgPassLooseXGBLead[e]->Clone();
                hBkgPassLooseXGB->Add(hBkgPassLooseXGBSub[e]);

                TH1F *hSigPassTightXGB = (TH1F*)hSigPassTightXGBLead[e]->Clone();
                hSigPassTightXGB->Add(hSigPassTightXGBSub[e]);
                TH1F *hBkgPassTightXGB = (TH1F*)hBkgPassTightXGBLead[e]->Clone();
                hBkgPassTightXGB->Add(hBkgPassTightXGBSub[e]);
                
                TH1F *hSigPassTightXGBTightMass = (TH1F*)hSigPassTightXGBTightMassLead[e]->Clone();
                hSigPassTightXGBTightMass->Add(hSigPassTightXGBTightMassSub[e]);
                TH1F *hBkgPassTightXGBTightMass = (TH1F*)hBkgPassTightXGBTightMassLead[e]->Clone();
                hBkgPassTightXGBTightMass->Add(hBkgPassTightXGBTightMassSub[e]);

                TH1F *hSigPassStd = (TH1F*)hSigPassStdLead[e]->Clone();
                hSigPassStd->Add(hSigPassStdSub[e]);
                TH1F *hBkgPassStd = (TH1F*)hBkgPassStdLead[e]->Clone();
                hBkgPassStd->Add(hBkgPassStdSub[e]);
                
                TH1F *hSigPassStdManual = (TH1F*)hSigPassStdManualLead[e]->Clone();
                hSigPassStdManual->Add(hSigPassStdManualSub[e]);
                TH1F *hBkgPassStdManual = (TH1F*)hBkgPassStdManualLead[e]->Clone();
                hBkgPassStdManual->Add(hBkgPassStdManualSub[e]);
            
                char buff1[100], buff2[100], buff3[100], buff4[100], buff5[100], buff6[100], buff7[100], buff8[100], buff9[100], buff10[100], buff11[100], buff12[100];
                
                snprintf(buff1, sizeof(buff1), "(%0.1f)", hSigAll->GetSumOfWeights());
                string nEvents1 = buff1;
                label1A += nEvents1;
                
                snprintf(buff2, sizeof(buff2), "(%0.1f)", hSigPassLooseXGB->GetSumOfWeights());
                string nEvents1PLoose = buff2;
                label1PLoose += nEvents1PLoose;
                
                snprintf(buff3, sizeof(buff3), "(%0.1f)", hSigPassTightXGB->GetSumOfWeights());
                string nEvents1PTight = buff3;
                label1PTight += nEvents1PTight;
                
                snprintf(buff9, sizeof(buff9), "(%0.1f)", hSigPassTightXGBTightMass->GetSumOfWeights());
                string nEvents1PTightTightMass = buff9;
                label1PTightTightMass += nEvents1PTightTightMass;
                
                snprintf(buff4, sizeof(buff4), "(%0.1f)", hSigPassStd->GetSumOfWeights());
                string nEvents1PStd = buff4;
                label1PStd += nEvents1PStd;

                snprintf(buff11, sizeof(buff11), "(%0.1f)", hSigPassStdManual->GetSumOfWeights());
                string nEvents1PStdManual = buff11;
                label1PStdManual += nEvents1PStdManual;

                
                TH1F *hSigAllDraw = DrawOverflow(hSigAll);
                TH1F *hStackHistoSig = (TH1F*)hSigAllDraw->Clone();
                hStackHistoSig->Reset();
                hStackHistoSig->GetXaxis()->SetRange(0,hStackHistoSig->GetNbinsX());
                hStack1->SetHistogram(hStackHistoSig);
                
                hStack1->GetHistogram()->GetXaxis()->SetTitle(varNames[i].c_str());
                hStack1->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
                hStack1->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
                hStack1->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
                
                hStack1->Add(hSigAllDraw);
                legend1->AddEntry(hSigAllDraw,label1A.c_str(),"pl");
                
                TH1F *hSigPassLooseDraw = DrawOverflow(hSigPassLooseXGB);
                hSigPassLooseDraw->SetLineColor(2);
                hStack1->Add(hSigPassLooseDraw);
                legend1->AddEntry(hSigPassLooseDraw,label1PLoose.c_str(),"pl");
                
                TH1F *hSigPassStdDraw = DrawOverflow(hSigPassStd);
                hSigPassStdDraw->SetLineColor(3);
                hStack1->Add(hSigPassStdDraw);
                legend1->AddEntry(hSigPassStdDraw,label1PStd.c_str(),"pl");

                TH1F *hSigPassStdManualDraw = DrawOverflow(hSigPassStdManual);
                hSigPassStdManualDraw->SetLineColor(3);
                hSigPassStdManualDraw->SetLineStyle(9);
                hStack1->Add(hSigPassStdManualDraw);
                legend1->AddEntry(hSigPassStdManualDraw,label1PStdManual.c_str(),"pl");
                
                TH1F *hSigPassTightDraw = DrawOverflow(hSigPassTightXGB);
                hSigPassTightDraw->SetLineColor(4);
                hStack1->Add(hSigPassTightDraw);
                legend1->AddEntry(hSigPassTightDraw,label1PTight.c_str(),"pl");
                
                TH1F *hSigPassTightTightMassDraw = DrawOverflow(hSigPassTightXGBTightMass);
                hSigPassTightTightMassDraw->SetLineColor(4);
                hSigPassTightTightMassDraw->SetLineStyle(9);
                hStack1->Add(hSigPassTightTightMassDraw);
                legend1->AddEntry(hSigPassTightTightMassDraw,label1PTightTightMass.c_str(),"pl");
                
                can->cd(1);
                gPad->SetGrid();
                if (varNames[i] == "ecalPFIso" || varNames[i] == "trkIsoPho")gPad->SetLogy();
                else gPad->SetLogy(0);
                
                hStack1->Draw("nostackhist");
                legend1->Draw("same");
                                
                snprintf(buff5, sizeof(buff5), "(%0.1f)", hBkgAll->GetSumOfWeights());
                string nEvents2 = buff5;
                label2A += nEvents2;

                snprintf(buff6, sizeof(buff6), "(%0.1f)", hBkgPassLooseXGB->GetSumOfWeights());
                string nEvents2PLoose = buff6;
                label2PLoose += nEvents2PLoose;
                
                snprintf(buff7, sizeof(buff7), "(%0.1f)", hBkgPassTightXGB->GetSumOfWeights());
                string nEvents2PTight = buff7;
                label2PTight += nEvents2PTight;

                snprintf(buff10, sizeof(buff10), "(%0.1f)", hBkgPassTightXGBTightMass->GetSumOfWeights());
                string nEvents2PTightTightMass = buff10;
                label2PTightTightMass += nEvents2PTightTightMass;
                
                snprintf(buff8, sizeof(buff8), "(%0.1f)", hBkgPassStd->GetSumOfWeights());
                string nEvents2PStd = buff8;
                label2PStd += nEvents2PStd;
                
                snprintf(buff12, sizeof(buff12), "(%0.1f)", hBkgPassStdManual->GetSumOfWeights());
                string nEvents2PStdManual = buff12;
                label2PStdManual += nEvents2PStdManual;
                
                TH1F *hBkgAllDraw = DrawOverflow(hBkgAll);
                TH1F *hStackHistoBkg = (TH1F*)hBkgAllDraw->Clone();
                hStackHistoBkg->Reset();
                hStackHistoBkg->GetXaxis()->SetRange(0,hStackHistoBkg->GetNbinsX());
                hStack2->SetHistogram(hStackHistoBkg);
                
                hStack2->GetHistogram()->GetXaxis()->SetTitle(varNames[i].c_str());
                hStack2->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
                hStack2->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
                hStack2->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
                
                hStack2->Add(hBkgAllDraw);
                legend2->AddEntry(hBkgAllDraw,label2A.c_str(),"pl");
                
                TH1F *hBkgPassLooseXGBDraw = DrawOverflow(hBkgPassLooseXGB);
                hBkgPassLooseXGBDraw->SetLineColor(2);
                legend2->AddEntry(hBkgPassLooseXGBDraw,label2PLoose.c_str(),"pl");
                hStack2->Add(hBkgPassLooseXGBDraw);
                
                TH1F *hBkgPassStdDraw = DrawOverflow(hBkgPassStd);
                hBkgPassStdDraw->SetLineColor(3);
                hStack2->Add(hBkgPassStdDraw);
                legend2->AddEntry(hBkgPassStdDraw,label2PStd.c_str(),"pl");
                
                TH1F *hBkgPassStdManualDraw = DrawOverflow(hBkgPassStdManual);
                hBkgPassStdManualDraw->SetLineColor(3);
                hBkgPassStdManualDraw->SetLineStyle(9);
                hStack2->Add(hBkgPassStdManualDraw);
                legend2->AddEntry(hBkgPassStdManualDraw,label2PStdManual.c_str(),"pl");
                
                TH1F *hBkgPassTightXGBDraw = DrawOverflow(hBkgPassTightXGB);
                hBkgPassTightXGBDraw->SetLineColor(4);
                hStack2->Add(hBkgPassTightXGBDraw);
                legend2->AddEntry(hBkgPassTightXGBDraw,label2PTight.c_str(),"pl");
                
                TH1F *hBkgPassTightXGBTightMassDraw = DrawOverflow(hBkgPassTightXGBTightMass);
                hBkgPassTightXGBTightMassDraw->SetLineColor(4);
                hBkgPassTightXGBTightMassDraw->SetLineStyle(9);
                hStack2->Add(hBkgPassTightXGBTightMassDraw);
                legend2->AddEntry(hBkgPassTightXGBTightMassDraw,label2PTightTightMass.c_str(),"pl");
                
                can->cd(2);
                gPad->SetGrid();
                if (varNames[i] == "ecalPFIso" || varNames[i] == "trkIsoPho")gPad->SetLogy();
                else gPad->SetLogy(0);
                gPad->SetGrid();
                
                hStack2->Draw("nostackhist");
                legend2->Draw("same");
                
                
                can->SaveAs((outName+".png").c_str());
                can->SaveAs((outName+".root").c_str());
                
                can->Clear();
                            
                f->Close();
                f->Delete();
            }//For "Added" Plotting
            
            
            if (plotType == "Separate"){
                
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
                
                snprintf(buff6, sizeof(buff6), "(%0.1f)", hSigPassStdManualLead[e]->GetSumOfWeights());
                string nEvents1PStdManual = buff6;
                label1PStdManual += nEvents1PStdManual;
                
                TH1F *hSigAllLeadDraw = DrawOverflow(hSigAllLead[e]);
                TH1F *hStackHistoSig = (TH1F*)hSigAllLeadDraw->Clone();
                hStackHistoSig->Reset();
                hStackHistoSig->GetXaxis()->SetRange(0,hStackHistoSig->GetNbinsX());
                hStack1->SetHistogram(hStackHistoSig);
                
                hStack1->GetHistogram()->GetXaxis()->SetTitle(varNames[i].c_str());
                hStack1->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
                hStack1->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
                hStack1->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
                
                hStack1->Add(hSigAllLeadDraw);
                legend1->AddEntry(hSigAllLeadDraw,label1A.c_str(),"pl");
                TH1F *hSigAllSubDraw = DrawOverflow(hSigAllSub[e]);
                hSigAllSubDraw->SetLineStyle(3);
                hStack1->Add(hSigAllSubDraw);
            
                TH1F *hSigPassLooseLeadDraw = DrawOverflow(hSigPassLooseXGBLead[e]);
                hSigPassLooseLeadDraw->SetLineColor(2);
                hStack1->Add(hSigPassLooseLeadDraw);
                legend1->AddEntry(hSigPassLooseLeadDraw,label1PLoose.c_str(),"pl");
                TH1F *hSigPassLooseSubDraw = DrawOverflow(hSigPassLooseXGBSub[e]);
                hSigPassLooseSubDraw->SetLineColor(2);
                hSigPassLooseSubDraw->SetLineStyle(3);
                hStack1->Add(hSigPassLooseSubDraw);
               
                TH1F *hSigPassStdLeadDraw = DrawOverflow(hSigPassStdLead[e]);
                hSigPassStdLeadDraw->SetLineColor(3);
                hStack1->Add(hSigPassStdLeadDraw);
                legend1->AddEntry(hSigPassStdLeadDraw,label1PStd.c_str(),"pl");
                TH1F *hSigPassStdSubDraw = DrawOverflow(hSigPassStdSub[e]);
                hSigPassStdSubDraw->SetLineColor(3);
                hSigPassStdSubDraw->SetLineStyle(3);
                hStack1->Add(hSigPassStdSubDraw);
                
                TH1F *hSigPassStdManualLeadDraw = DrawOverflow(hSigPassStdManualLead[e]);
                hSigPassStdManualLeadDraw->SetLineColor(5);
                hStack1->Add(hSigPassStdManualLeadDraw);
                legend1->AddEntry(hSigPassStdManualLeadDraw,label1PStdManual.c_str(),"pl");
                TH1F *hSigPassStdManualSubDraw = DrawOverflow(hSigPassStdManualSub[e]);
                hSigPassStdManualSubDraw->SetLineColor(5);
                hSigPassStdManualSubDraw->SetLineStyle(3);
                hStack1->Add(hSigPassStdManualSubDraw);
                
                TH1F *hSigPassTightLeadDraw = DrawOverflow(hSigPassTightXGBLead[e]);
                hSigPassTightLeadDraw->SetLineColor(4);
                hStack1->Add(hSigPassTightLeadDraw);
                legend1->AddEntry(hSigPassTightLeadDraw,label1PTight.c_str(),"pl");
                TH1F *hSigPassTightSubDraw = DrawOverflow(hSigPassTightXGBSub[e]);
                hSigPassTightSubDraw->SetLineColor(4);
                hSigPassTightSubDraw->SetLineStyle(3);
                hStack1->Add(hSigPassTightSubDraw);
                
//
//                TH1F *hSigPassTightTightMassLeadDraw = DrawOverflow(hSigPassTightXGBTightMassLead[e]);
//                hSigPassTightTightMassLeadDraw->SetLineColor(4);
//                hSigPassTightTightMassLeadDraw->SetLineStyle(9);
//                hStack1->Add(hSigPassTightTightMassLeadDraw);
//                legend1->AddEntry(hSigPassTightTightMassDraw,label1PTightTightMass.c_str(),"pl");
//
                can->cd(1);
                gPad->SetGrid();
                if (varNames[i] == "ecalPFIso" || varNames[i] == "trkIsoPho")gPad->SetLogy();
                else gPad->SetLogy(0);
                
                hStack1->Draw("nostackhist");
                legend1->Draw("same");
                
                snprintf(buff7, sizeof(buff7), "(%0.1f)", hBkgAllLead[e]->GetSumOfWeights());
                string nEvents2 = buff7;
                label2A += nEvents2;
                
                snprintf(buff8, sizeof(buff8), "(%0.1f)", hBkgPassLooseXGBLead[e]->GetSumOfWeights());
                string nEvents2PLoose = buff8;
                label2PLoose += nEvents2PLoose;
                
                snprintf(buff9, sizeof(buff9), "(%0.1f)", hBkgPassTightXGBLead[e]->GetSumOfWeights());
                string nEvents2PTight = buff9;
                label2PTight += nEvents2PTight;
                
                snprintf(buff10, sizeof(buff10), "(%0.1f)", hBkgPassTightXGBTightMassLead[e]->GetSumOfWeights());
                string nEvents2PTightTightMass = buff10;
                label2PTightTightMass += nEvents2PTightTightMass;
                
                snprintf(buff11, sizeof(buff11), "(%0.1f)", hBkgPassStdLead[e]->GetSumOfWeights());
                string nEvents2PStd = buff11;
                label2PStd += nEvents2PStd;
                
                snprintf(buff12, sizeof(buff12), "(%0.1f)", hBkgPassStdManualLead[e]->GetSumOfWeights());
                string nEvents2PStdManual = buff12;
                label2PStdManual += nEvents2PStdManual;
                
                TH1F *hBkgAllLeadDraw = DrawOverflow(hBkgAllLead[e]);
                TH1F *hStackHistoBkg = (TH1F*)hBkgAllLeadDraw->Clone();
                hStackHistoBkg->Reset();
                hStackHistoBkg->GetXaxis()->SetRange(0,hStackHistoBkg->GetNbinsX());
                hStack2->SetHistogram(hStackHistoBkg);
                
                hStack2->GetHistogram()->GetXaxis()->SetTitle(varNames[i].c_str());
                hStack2->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
                hStack2->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
                hStack2->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
                
                hStack2->Add(hBkgAllLeadDraw);
                legend2->AddEntry(hBkgAllLeadDraw,label2A.c_str(),"pl");
                TH1F *hBkgAllSubDraw = DrawOverflow(hBkgAllSub[e]);
                hBkgAllSubDraw->SetLineStyle(3);
                hStack2->Add(hBkgAllSubDraw);
                
                TH1F *hBkgPassLooseLeadDraw = DrawOverflow(hBkgPassLooseXGBLead[e]);
                hBkgPassLooseLeadDraw->SetLineColor(2);
                hStack2->Add(hBkgPassLooseLeadDraw);
                legend2->AddEntry(hBkgPassLooseLeadDraw,label2PLoose.c_str(),"pl");
                TH1F *hBkgPassLooseSubDraw = DrawOverflow(hBkgPassLooseXGBSub[e]);
                hBkgPassLooseSubDraw->SetLineColor(2);
                hBkgPassLooseSubDraw->SetLineStyle(3);
                hStack2->Add(hBkgPassLooseSubDraw);
                
                TH1F *hBkgPassStdLeadDraw = DrawOverflow(hBkgPassStdLead[e]);
                hBkgPassStdLeadDraw->SetLineColor(3);
                hStack2->Add(hBkgPassStdLeadDraw);
                legend2->AddEntry(hBkgPassStdLeadDraw,label2PStd.c_str(),"pl");
                TH1F *hBkgPassStdSubDraw = DrawOverflow(hBkgPassStdSub[e]);
                hBkgPassStdSubDraw->SetLineColor(3);
                hBkgPassStdSubDraw->SetLineStyle(3);
                hStack2->Add(hBkgPassStdSubDraw);
                
                TH1F *hBkgPassStdManualLeadDraw = DrawOverflow(hBkgPassStdManualLead[e]);
                hBkgPassStdManualLeadDraw->SetLineColor(5);
                hStack2->Add(hBkgPassStdManualLeadDraw);
                legend2->AddEntry(hBkgPassStdManualLeadDraw,label2PStdManual.c_str(),"pl");
                TH1F *hBkgPassStdManualSubDraw = DrawOverflow(hBkgPassStdManualSub[e]);
                hBkgPassStdManualSubDraw->SetLineColor(5);
                hBkgPassStdManualSubDraw->SetLineStyle(3);
                hStack2->Add(hBkgPassStdManualSubDraw);
                
                TH1F *hBkgPassTightLeadDraw = DrawOverflow(hBkgPassTightXGBLead[e]);
                hBkgPassTightLeadDraw->SetLineColor(4);
                hStack2->Add(hBkgPassTightLeadDraw);
                legend2->AddEntry(hBkgPassTightLeadDraw,label2PTight.c_str(),"pl");
                TH1F *hBkgPassTightSubDraw = DrawOverflow(hBkgPassTightXGBSub[e]);
                hBkgPassTightSubDraw->SetLineColor(4);
                hBkgPassTightSubDraw->SetLineStyle(3);
                hStack2->Add(hBkgPassTightSubDraw);
               
                can->cd(2);
                gPad->SetGrid();
                if (varNames[i] == "ecalPFIso" || varNames[i] == "trkIsoPho")gPad->SetLogy();
                else gPad->SetLogy(0);
                gPad->SetGrid();
                
                hStack2->Draw("nostackhist");
                legend2->Draw("same");
                
                
                can->SaveAs((outName+".png").c_str());
                can->SaveAs((outName+".root").c_str());
                
                can->Clear();
                
                f->Close();
                f->Delete();
            }//For "Separate" Plotting

            
        }//eta loop
        for (int p = 0; p < 4; p++){
            hSigAllLead[p]->Delete();
            hSigAllSub[p]->Delete();
            hSigPassLooseXGBLead[p]->Delete();
            hSigPassLooseXGBSub[p]->Delete();
            hSigPassTightXGBLead[p]->Delete();
            hSigPassTightXGBSub[p]->Delete();
            hSigPassTightXGBTightMassLead[p]->Delete();
            hSigPassTightXGBTightMassSub[p]->Delete();
            hSigPassStdLead[p]->Delete();
            hSigPassStdSub[p]->Delete();
            hSigPassStdManualLead[p]->Delete();
            hSigPassStdManualSub[p]->Delete();
            
            hBkgAllLead[p]->Delete();
            hBkgAllSub[p]->Delete();
            hBkgPassLooseXGBLead[p]->Delete();
            hBkgPassLooseXGBSub[p]->Delete();
            hBkgPassTightXGBLead[p]->Delete();
            hBkgPassTightXGBSub[p]->Delete();
            hBkgPassTightXGBTightMassLead[p]->Delete();
            hBkgPassTightXGBTightMassSub[p]->Delete();
            hBkgPassStdLead[p]->Delete();
            hBkgPassStdSub[p]->Delete();
            hBkgPassStdManualLead[p]->Delete();
            hBkgPassStdManualSub[p]->Delete();
        }
    }// variable loop (i)
}//Function definition


TH1F *DrawOverflow(TH1F* h){
    int nBins    = h->GetNbinsX()+1;
    Double_t *xbins= new Double_t[nBins+1];
    for (UInt_t i=0;i<nBins;i++)
        xbins[i]=h->GetBinLowEdge(i+1);
    xbins[nBins]=xbins[nBins-1]+h->GetBinWidth(nBins);
    //book a temporary histogram having extra bins for overflows
    TH1F *htmp = new TH1F(h->GetName(), h->GetTitle(), nBins, xbins);
    htmp->Sumw2();
    //fill the new histogram including the overflows
    for (UInt_t i=1; i<=nBins; i++) {
        htmp->SetBinContent(htmp->FindBin(htmp->GetBinCenter(i)),h->GetBinContent(i));
    }
    //htmp->SetBinContent(htmp->FindBin(h->GetBinLowEdge(1)-1), h->GetBinContent(0));
    htmp->SetBinContent(0, h->GetBinContent(0));
    htmp->GetXaxis()->SetRange(0,nBins);
    return htmp;
}

string ToStringWithPrecision(double inVal, const int n = 2){
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << inVal;
    return std::move(out).str();
}
double calcMassXGB(TTreeReaderArray<float>* xgbScore,double xgbCutB, double xgbCutE, TTreeReaderArray<float>* et, TTreeReaderArray<float>* eta, TTreeReaderArray<float>* phi, int nEgsIn, int* pho1, int* pho2){
    double mass = -999.0;
    int iOut = -1;
    int jOut = -1;
    int bothPass = 0;
    int someSelectionMade = 0;
    if(nEgsIn > 1){
        for(int i = 0; i < nEgsIn;i++){
            //if ((*et)[i] > 14.25 && fabs((*eta)[i]) < 2.5){
            if (fabs((*eta)[i]) < 2.5 && (*et)[i] > 14.25 && ((abs((*eta)[i]) < 1.5 && (*xgbScore)[i] > xgbCutB) ||(abs((*eta)[i]) > 1.5 && (*xgbScore)[i] > xgbCutE))){
                ROOT::Math::PtEtaPhiMVector v1((*et)[i],(*eta)[i],(*phi)[i],0.0);
                for(int j = i+1; j < nEgsIn; j ++){
                    //if(fabs((*eta)[j]) < 2.5 && (*et)[j] > 14.25){
                    if (fabs((*eta)[j]) < 2.5 && (*et)[j] > 14.25 && ((abs((*eta)[j]) < 1.5 && (*xgbScore)[j] > xgbCutB) ||(abs((*eta)[j]) > 1.5 && (*xgbScore)[j] > xgbCutE))){
                        ROOT::Math::PtEtaPhiMVector v2((*et)[j],(*eta)[j],(*phi)[j],0.0);
                        auto v3 = v1+v2;
                        double massTmp = v3.M();
                        double dR = calcDR((*eta)[i],(*eta)[j],(*phi)[i],(*phi)[j]);
                        //double dRTmp = calcDR((*eta)[i], (*eta)[j], (*phi)[i], (*phi)[j])
                        if(massTmp > mass){
                            bothPass = 1;
                            mass = massTmp;
                            if (massTmp > 60 && dR > 0.5){
                                iOut = i;
                                jOut = j;
                                someSelectionMade = 1;
                            }//Mass cut for photon selection. If mass is less than 60, don't use this pair.
                        }//massTmp > mass
                    }//If et and eta good [j]
                }//j Loop
            }//If et and eta good [i]
        }//i Loop for passing both
        if(bothPass != 1){
            double scoreLead = -1;
            double scoreSub = -1;
            for(int i = 0; i < nEgsIn;i++){
                if(fabs((*eta)[i]) < 2.5 && (*et)[i] > 14.25){
                    double scoreLeadTmp = (*xgbScore)[i];
                    if(scoreLeadTmp > scoreLead){
                        scoreLead = scoreLeadTmp;
                        for(int j = i+1; j < nEgsIn; j ++){
                            if(fabs((*eta)[j]) < 2.5 && (*et)[j] > 14.25){
                                double scoreSubTmp = (*xgbScore)[j];
                                if(scoreSubTmp > scoreSub){
                                    scoreLead = scoreSubTmp;
                                    double dR = calcDR((*eta)[i],(*eta)[j],(*phi)[i],(*phi)[j]);
                                    if (dR > 0.5){
                                        iOut = i;
                                        jOut = j;
                                        someSelectionMade = 1;
                                    }
                                }//If score oof j is larger than highest sub score
                            }//Et and Eta cuts on second photon
                        }//j Loop
                    }//If score of i is larger than highest lead score
                }//Et and Eta cuts on first photon
            }//i Loop for not passing both
            ROOT::Math::PtEtaPhiMVector v1((*et)[iOut],(*eta)[iOut],(*phi)[iOut],0.0);
            ROOT::Math::PtEtaPhiMVector v2((*et)[jOut],(*eta)[jOut],(*phi)[jOut],0.0);
            auto v3 = v1+v2;
            mass = v3.M();
        }//If bothPass != 1
    }//nEgs
    if(someSelectionMade != 1 || nEgsIn < 2){
        iOut = -1;
        jOut = -1;
        mass = -999;
    }
    *pho1 = iOut;
    *pho2 = jOut;
    return mass;
}
double calcMassHLT(TTreeReaderArray<float>* et, TTreeReaderArray<float>* eta, TTreeReaderArray<float>* phi, TTreeReaderArray<float>* r9, TTreeReaderArray<float>* hOvrE, TTreeReaderArray<float>* sIeIe, TTreeReaderArray<float>* phoIso, TTreeReaderArray<float>* ecalIso, int nEgsIn, int* pho1, int* pho2){
    double mass = -999.0;
    int iOut = -1;
    int jOut = -1;
    int bothPass = 0;
    int hltPassFail = 0;
    if(nEgsIn > 1){
        for(int i = 0; i < nEgsIn;i++){
            //if ((*et)[i] > 14.25 && fabs((*eta)[i]) < 2.5){
            if ((*et)[i] > 30.0 && fabs((*eta)[i]) < 2.55 && ((fabs((*eta)[i]) < 1.479 && (*r9)[i] > 0.5 && (*hOvrE)[i] < 0.12 && !((*r9)[i] < 0.85 && ((*sIeIe)[i] > 0.015 || (*phoIso)[i] > 6))) || (fabs((*eta)[i]) > 1.479 && (*r9)[i] > 0.8 && (*hOvrE)[i] < 0.10 && !((*r9)[i] < 0.90 && ((*sIeIe)[i] > 0.035 || (*phoIso)[i] > 6))))){
                ROOT::Math::PtEtaPhiMVector v1((*et)[i],(*eta)[i],(*phi)[i],0.0);
                for(int j = i+1; j < nEgsIn; j ++){
                    //if(fabs((*eta)[j]) < 2.5 && (*et)[j] > 14.25){
                    if ((*et)[j] > 22.0 && fabs((*eta)[j]) < 2.55 && ((fabs((*eta)[j]) < 1.479 && (*r9)[j] > 0.5 && (*hOvrE)[j] < 0.12 && !((*r9)[j] < 0.85 && ((*sIeIe)[j] > 0.015 || (*phoIso)[j] > 6))) || (fabs((*eta)[j]) > 1.479 && (*r9)[j] > 0.8 && (*hOvrE)[j] < 0.10 && !((*r9)[j] < 0.90 && ((*sIeIe)[j] > 0.035 || (*phoIso)[j] > 6))))){
                        ROOT::Math::PtEtaPhiMVector v2((*et)[j],(*eta)[j],(*phi)[j],0.0);
                        auto v3 = v1+v2;
                        double massTmp = v3.M();
                        //double dRTmp = calcDR((*eta)[i], (*eta)[j], (*phi)[i], (*phi)[j])
                        if(massTmp > mass){
                            bothPass = 1;
                            iOut = i;
                            jOut = j;
                        }//massTmp > mass
                    }//If et and eta good [j]
                }//j Loop
            }//If et and eta good [i]
        }//i Loop for passing both
        if (bothPass == 1){
            for(int i = 0; i < nEgsIn; i ++){
                if ((*et)[i] > 22.0 && fabs((*eta)[i]) < 2.55 && ((fabs((*eta)[i]) < 1.5 && (*r9)[i] > 0.5 && (*hOvrE)[i] < 0.12) || (fabs((*eta)[i]) > 1.5 && (*r9)[i] > 0.8 && (*hOvrE)[i] < 0.10))){
                    ROOT::Math::PtEtaPhiMVector v1((*et)[i],(*eta)[i],(*phi)[i],0.0);
                    for(int j = i+1; j < nEgsIn; j ++){
                        if ((*et)[j] > 22.0 && fabs((*eta)[j]) < 2.55 && ((fabs((*eta)[j]) < 1.5 && (*r9)[j] > 0.5 && (*hOvrE)[j] < 0.12 && !((*r9)[j] < 0.85 && ((*sIeIe)[j] > 0.015 || (*phoIso)[j] > 6 || (*ecalIso)[j] > 6))) || (fabs((*eta)[j]) > 1.5 && (*r9)[j] > 0.8 && (*hOvrE)[j] < 0.10 && !((*r9)[j] < 0.90 && ((*sIeIe)[j] > 0.035 || (*phoIso)[j] > 6 || (*ecalIso)[j] > 6))))){
                            ROOT::Math::PtEtaPhiMVector v2((*et)[j],(*eta)[j],(*phi)[j],0.0);
                            auto v3 = v1+v2;
                            double massTmp = v3.M();
                            if(massTmp > mass){
                                mass = massTmp;
                            }
                        }
                    }
                }
            }
        }
    }//nEgs
    if(bothPass != 1){
        mass = -999;
        iOut = -1;
        jOut = -1;
    }
    *pho1 = iOut;
    *pho2 = jOut;
    return mass;
}
double calcDR(float eta1, float eta2, float phi1, float phi2){
    double DR = 0.0;
    double dE = eta1-eta2;
    double dP = phi1-phi2;
    if (dP > M_PI) dP -= 2*M_PI;
    if (dP < -M_PI) dP += 2*M_PI;
    
    DR = sqrt(dE*dE + dP*dP);
    
    return DR;
}
