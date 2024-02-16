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

TH1F *DrawOverflow(TH1F* h);
void singlePhotonSigBkgComp(){
    gROOT->Reset();
    gStyle->SetPalette(1);
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitle(0);
    
    string fileName = "validationNTuples/0517/SinglePhoton_GGHandData_MD9LR10_M60.root";
    string genTitleStringSignal = "GGH Signal ";
    string genTitleStringBkg = "Data ";
    string dirStr ="varPlots/0519/MD9LR10_M60/";
    //string dirStr ="varPlots/0130/DataRelaxedUnseededChoose2/";
    
    int nEta;
    
    nEta = 3;
    string etaLabels[3] = {"Barrel","Endcap","All #eta"};
    string etaFLabels[3] = {"_B","_E","_All"};
    
    string outNameGen = dirStr;
    
    int nVars = 17;
    string varNames[17] = {"nEgs","triggerBits","mass","et","energy","rawEnergy","eta","etaWidth","phi","phiWidth","sigmaIEtaIEta","ecalPFIso","r9HLT","s4","hOvrE","trkIsoPho","xgbScore"};
    double limsLow[17] = {1,0.0,0.0,0.0,0.0,0.0,-3,0.0,-3.15,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    //double limsHigh[16] = {3,4,500,300,600,600,3.15,0.1,3.15,0.2,0.2,20.0,1.2,1.2,1.0,8.0};
    double limsHigh[17] = {3,4,500,300,600,600,3,0.1,3.15,0.2,0.1,40.0,1.2,1.2,1.0,40.0,1.0};
    //double limsHigh[16] = {3,4,500,300,300,600,3.15,0.1,3.15,0.2,2.0,200.0,5.0,1.2,2.0,200.0};
    
    int nBins[17] = {2,4,500,600,600,600,200,400,400,400,400,400,400,400,400,400,400};
    
    TCanvas *can = new TCanvas ("can","can",10,10,1600,900);
    //for(int i = 10; i < 11; i++){
    for(int i = 0; i < nVars; i++){
        // for(int i = 14; i < nVars; i++){
        
        TH1F *hSigAll[3];
        TH1F *hBkgAll[3];
        
        TH1F *hSigPass[3];
        TH1F *hBkgPass[3];
        
        TH1F *hSigPassTightMass[3];
        TH1F *hBkgPassTightMass[3];
        
        TH1F *hSigPassXGB[3];
        TH1F *hBkgPassXGB[3];
        
        TH1F *hSigPassL1TwoPhotons[3];
        TH1F *hBkgPassL1TwoPhotons[3];
        
        TH1F *hSigPassStd[3];
        TH1F *hBkgPassStd[3];
        
        TH1F *hSigPassStdManual[3];
        TH1F *hBkgPassStdManual[3];
        
        for (int e = 0; e < nEta; e++){
            hSigAll[e] = new TH1F ("hSigAll","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgAll[e] = new TH1F ("hBkgAll","",nBins[i],limsLow[i],limsHigh[i]);
            
            hSigPass[e] = new TH1F ("hSigPass","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPass[e] = new TH1F ("hBkgPass","",nBins[i],limsLow[i],limsHigh[i]);
            
            hSigPassTightMass[e] = new TH1F ("hSigTightPass","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassTightMass[e] = new TH1F ("hBkgTightPass","",nBins[i],limsLow[i],limsHigh[i]);
            
            hSigPassXGB[e] = new TH1F ("hSigPassXGB","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassXGB[e] = new TH1F ("hBkgPassXGB","",nBins[i],limsLow[i],limsHigh[i]);
            
            hSigPassL1TwoPhotons[e] = new TH1F ("hSigPassL1TwoPhotons","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassL1TwoPhotons[e] = new TH1F ("hBkgPassL1TwoPhotons","",nBins[i],limsLow[i],limsHigh[i]);
            
            hSigPassStd[e] = new TH1F ("hSigPassStd","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassStd[e] = new TH1F ("hBkgPassStd","",nBins[i],limsLow[i],limsHigh[i]);
            
            hSigPassStdManual[e] = new TH1F ("hSigPassStdManual","",nBins[i],limsLow[i],limsHigh[i]);
            hBkgPassStdManual[e] = new TH1F ("hBkgPassStdManual","",nBins[i],limsLow[i],limsHigh[i]);
            
        }
        for (int e = 0; e < nEta; e++){
            
            string outName = outNameGen + varNames[i] + etaFLabels[e];
            
            TFile *f = new TFile(fileName.c_str());
            
            string varNameTmp;
            if (e != 2){
                varNameTmp = varNames[i];
                
                //FIRST handle Signal
                TTreeReader inReaderSig("sigTree", f);
                inReaderSig.Restart();
                
                TTreeReaderValue<Float_t> inVarSig(inReaderSig, varNameTmp.c_str());
                TTreeReaderValue<Float_t> inESig(inReaderSig,"energy");
                TTreeReaderValue<Float_t> inMSig(inReaderSig,"mass");
                TTreeReaderValue<Float_t> inNEgsSig(inReaderSig,"nEgs");
                TTreeReaderValue<Float_t> inEtSig(inReaderSig,"et");
                TTreeReaderValue<Float_t> inEtaSig(inReaderSig, "eta");
                TTreeReaderValue<Float_t> inR9Sig(inReaderSig,"r9HLT");
                TTreeReaderValue<Float_t> inHovrESig(inReaderSig,"hOvrE");
                TTreeReaderValue<Float_t> inSigIEtaIEtaSig(inReaderSig,"sigmaIEtaIEta");
                TTreeReaderValue<Float_t> inPhoIsoSig(inReaderSig,"trkIsoPho");
                TTreeReaderValue<Float_t> inXGBSig(inReaderSig,"xgbScore");
                
                
                TTreeReaderValue<int> inPassFailStdSig(inReaderSig,"passFailStd");
                TTreeReaderValue<int> inPassFailDoubleL1Sig(inReaderSig,"passFailL1Double");
                TTreeReaderValue<int> inPassFailSingleL1Sig(inReaderSig,"passFailL1Single");
                                
                while (inReaderSig.Next()) {
                    if((*inPassFailSingleL1Sig == 1 || *inPassFailDoubleL1Sig == 1) && ((e == 0 && abs(*inEtaSig) < 1.5) || (e == 1 && abs(*inEtaSig) > 1.5))){
                        hSigAll[e]->Fill(*inVarSig);
                        
                        if(*inNEgsSig > 1 && (*inPassFailSingleL1Sig == 1 || *inPassFailDoubleL1Sig == 1) && !(!(varNames[i] =="mass") && *inMSig < 60)){
                            hSigPass[e]->Fill(*inVarSig);
                            if(!(!(varNames[i] =="mass") && *inMSig < 95))hSigPassTightMass[e]->Fill(*inVarSig);
                        }
                        if(*inNEgsSig > 1){
                            hSigPassL1TwoPhotons[e]->Fill(*inVarSig);
                        }
                        if(*inXGBSig > 0.25){
                            hSigPassXGB[e]->Fill(*inVarSig);
                        }
                        if(*inPassFailStdSig == 1){
                            hSigPassStd[e]->Fill(*inVarSig);
                        }
                        if(e == 0){
                            if(*inEtSig > 22 && *inHovrESig < 0.12 && *inR9Sig > 0.5 && !(*inR9Sig < 0.85 && (*inSigIEtaIEtaSig > 0.015 || *inPhoIsoSig < 6)) && *inMSig > 95.0){
                                hSigPassStdManual[e]->Fill(*inVarSig);
                            }
                        }
                        if(e == 1){
                            if(*inEtSig > 22 && *inHovrESig < 0.1 && *inR9Sig > 0.8 && !(*inR9Sig < 0.9 && (*inSigIEtaIEtaSig > 0.035 || *inPhoIsoSig < 6)) && *inMSig > 95.0){
                                hSigPassStdManual[e]->Fill(*inVarSig);
                            }
                        }
                    }//Eta Selection if
                }//While inReaderSig.Next()
                                
                //THEN handle Bkg
                TTreeReader inReaderBkg("bkgTree", f);
                inReaderBkg.Restart();
                
                TTreeReaderValue<Float_t> inVarBkg(inReaderBkg, varNameTmp.c_str());
                TTreeReaderValue<Float_t> inEBkg(inReaderBkg,"energy");
                TTreeReaderValue<Float_t> inMBkg(inReaderBkg,"mass");
                TTreeReaderValue<Float_t> inNEgsBkg(inReaderBkg,"nEgs");
                TTreeReaderValue<Float_t> inEtBkg(inReaderBkg,"et");
                TTreeReaderValue<Float_t> inEtaBkg(inReaderBkg, "eta");
                TTreeReaderValue<Float_t> inR9Bkg(inReaderBkg,"r9HLT");
                TTreeReaderValue<Float_t> inHovrEBkg(inReaderBkg,"hOvrE");
                TTreeReaderValue<Float_t> inSigIEtaIEtaBkg(inReaderBkg,"sigmaIEtaIEta");
                TTreeReaderValue<Float_t> inPhoIsoBkg(inReaderBkg,"trkIsoPho");
                TTreeReaderValue<Float_t> inXGBBkg(inReaderBkg,"xgbScore");
                
                TTreeReaderValue<int> inPassFailStdBkg(inReaderBkg,"passFailStd");
                TTreeReaderValue<int> inPassFailDoubleL1Bkg(inReaderBkg,"passFailL1Double");
                TTreeReaderValue<int> inPassFailSingleL1Bkg(inReaderBkg,"passFailL1Single");
                
                while (inReaderBkg.Next()) {
                    if((*inPassFailSingleL1Bkg == 1 || *inPassFailDoubleL1Bkg == 1) && ((e == 0 && abs(*inEtaBkg) < 1.5) || (e == 1 && abs(*inEtaBkg) > 1.5))){
                        hBkgAll[e]->Fill(*inVarBkg);
                        
                        if(*inNEgsBkg > 1 && (*inPassFailSingleL1Bkg == 1 || *inPassFailDoubleL1Bkg == 1) && !(!(varNames[i] =="mass") && *inMBkg < 60)){
                            hBkgPass[e]->Fill(*inVarBkg);
                            if(!(!(varNames[i] =="mass") && *inMBkg < 95))hBkgPassTightMass[e]->Fill(*inVarBkg);
                        }
                        
                        if(*inNEgsBkg > 1){
                            hBkgPassL1TwoPhotons[e]->Fill(*inVarBkg);
                        }
                        if(*inPassFailStdBkg == 1){
                            hBkgPassStd[e]->Fill(*inVarBkg);
                        }
                        if(*inXGBBkg > 0.25){
                            hBkgPassXGB[e]->Fill(*inVarBkg);
                        }
                        if(e == 0){
                            if(*inEtBkg > 22 && *inHovrEBkg < 0.12 && *inR9Bkg > 0.5 && !(*inR9Bkg < 0.85 && (*inSigIEtaIEtaBkg > 0.015 || *inPhoIsoBkg < 6)) && *inMBkg > 95.0){
                                hBkgPassStdManual[e]->Fill(*inVarBkg);
                            }
                        }
                        if(e == 1){
                            if(*inEtBkg > 22 && *inHovrEBkg < 0.1 && *inR9Bkg > 0.8 && !(*inR9Bkg < 0.9 && (*inSigIEtaIEtaBkg > 0.035 || *inPhoIsoBkg < 6)) && *inMBkg > 95.0){
                                hBkgPassStdManual[e]->Fill(*inVarBkg);
                            }
                        }
                    }//Eta Selection if
                }//While inReaderBkg.Next()
            }//if e != 2
                        if (e == 2){
                for (int p = 0; p < 2; p++){
                  
                    
                    hSigAll[e]->Add(hSigAll[p]);
                    hBkgAll[e]->Add(hBkgAll[p]);
                        
                    hSigPass[e]->Add(hSigPass[p]);
                    hBkgPass[e]->Add(hBkgPass[p]);
                    
                    hSigPassTightMass[e]->Add(hSigPassTightMass[p]);
                    hBkgPassTightMass[e]->Add(hBkgPassTightMass[p]);
                        
                    hSigPassXGB[e]->Add(hSigPassXGB[p]);
                    hBkgPassXGB[e]->Add(hBkgPassXGB[p]);
                        
                    hSigPassL1TwoPhotons[e]->Add(hSigPassL1TwoPhotons[p]);
                    hBkgPassL1TwoPhotons[e]->Add(hBkgPassL1TwoPhotons[p]);
                        
                    hSigPassStd[e]->Add(hSigPassStd[p]);
                    hBkgPassStd[e]->Add(hBkgPassStd[p]);
                        
                    hSigPassStdManual[e]->Add(hSigPassStdManual[p]);
                    hBkgPassStdManual[e]->Add(hBkgPassStdManual[p]);
                        
                }
            }//if e == 2
            string varNameTmp2;
            if (varNames[i] == "mass")varNameTmp2 = "M_{#gamma#gamma}";
            else if (varNames[i] == "r9HLT")varNameTmp2 = "R_{9}";
            else if (varNames[i] == "sigmaIEtaIEta")varNameTmp2 = "#sigma_{i#eta i#eta}";
            else if (varNames[i] == "etaWidth")varNameTmp2 = "#eta Width";
            else if (varNames[i] == "phiWidth")varNameTmp2 = "#phi Width";
            else if (varNames[i] == "s4")varNameTmp2 = "S_{4}";
            else if (varNames[i] == "eta")varNameTmp2 = "#eta_{SC}";
            else if (varNames[i] == "et")varNameTmp2 = "E_{T}";
            else varNameTmp2 = varNames[i];
            
            can->Clear();
                
            can->Divide(2,1);
            
            
            string plotTitle1 = genTitleStringSignal + etaLabels[e] + ";" + varNameTmp2;
            THStack *hStack1 = new THStack("hStack1",plotTitle1.c_str());
            //hStack1->SetHistogram(new TH1F("hstot1","",nBins[i],limsLow[i],limsHigh[i]));
                
            string plotTitle2 = genTitleStringBkg + etaLabels[e] + ";" + varNameTmp2;
            THStack *hStack2 = new THStack("hStack2",plotTitle2.c_str());
            //hStack2->SetHistogram(new TH1F("hstot2","",nBins[i],limsLow[i],limsHigh[i]));
                
            TLegend *legend1;
            if (varNames[i] == "r9HLT"|| varNames[i] == "s4") legend1 = new TLegend(0.1,0.65,0.5,0.9,"","brNDC");
            else legend1 = new TLegend(0.50,0.65,0.90,0.9,"","brNDC");
                
            TLegend *legend2;
            if (varNames[i] == "r9HLT"|| varNames[i] == "s4") legend2 = new TLegend(0.1,0.65,0.5,0.9,"","brNDC");
            else legend2 = new TLegend(0.50,0.65,0.90,0.9,"","brNDC");
            
            string label1A = "All";
            string label12Pho = "2 Photons";
            string label1P = "Passing Relaxed Diphoton Cuts";
            string label1PTightMass = "Passing Relaxed Cuts & M>95";
            string label1XGB = "Passing xgbScore > 0.25";
            string label1PStd = "Passing Std. HLT";
            string label1PStdM = "Passing Std. HLT (Manual)";
            
            string label2A = "All";
            string label22Pho = "2 Photons";
            string label2P = "Passing Relaxed Diphoton Cuts";
            string label2PTightMass = "Passing Relaxed Cuts & M>95";
            string label2XGB = "Passing xgbScore > 0.25";
            string label2PStd = "Passing Std. HLT";
            string label2PStdM = "Passing Std. HLT (Manual)";

            
            char buff1[100], buff2[100], buff3[100], buff4[100], buff5[100], buff6[100], buff7[100], buff8[100], buff9[100], buff10[100], buff11[100], buff12[100], buff13[100], buff14[100];
                
            snprintf(buff1, sizeof(buff1), "(%0.1f)", hSigAll[e]->GetSumOfWeights());
            string nEvents1 = buff1;
            label1A += nEvents1;
            
            snprintf(buff2, sizeof(buff2), "(%0.1f)", hSigPass[e]->GetSumOfWeights());
            string nEvents1P = buff2;
            label1P += nEvents1P;

            snprintf(buff3, sizeof(buff3), "(%0.1f)", hSigPassXGB[e]->GetSumOfWeights());
            string nEvents1XGB = buff3;
            label1XGB += nEvents1XGB;
            
            snprintf(buff4, sizeof(buff4), "(%0.1f)", hSigPassL1TwoPhotons[e]->GetSumOfWeights());
            string nEvents1PL12 = buff4;
            label12Pho += nEvents1PL12;
                        
            snprintf(buff5, sizeof(buff5), "(%0.1f)", hSigPassTightMass[e]->GetSumOfWeights());
            string nEvents1PTight = buff5;
            label1PTightMass += nEvents1PTight;
            
            snprintf(buff6, sizeof(buff6), "(%0.1f)", hSigPassStd[e]->GetSumOfWeights());
            string nEvents1PStd = buff6;
            label1PStd += nEvents1PStd;
            
            snprintf(buff7, sizeof(buff7), "(%0.1f)", hSigPassStdManual[e]->GetSumOfWeights());
            string nEvents1PStdM = buff7;
            label1PStdM += nEvents1PStdM;
            
            TH1F *hSigAllDraw = DrawOverflow(hSigAll[e]);
            TH1F *hStackHistoSig = (TH1F*)hSigAllDraw->Clone();
            hStackHistoSig->Reset();
            hStackHistoSig->GetXaxis()->SetRange(0,hStackHistoSig->GetNbinsX());
            hStack1->SetHistogram(hStackHistoSig);
            
            hStack1->GetHistogram()->GetXaxis()->SetTitle(varNames[i].c_str());
            hStack1->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
            hStack1->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
            hStack1->GetHistogram()->GetYaxis()->SetLabelSize(0.025);

            TH1F *hSigPassL1TwoPhotonsDraw = DrawOverflow(hSigPassL1TwoPhotons[e]);
            hSigPassL1TwoPhotonsDraw->SetLineColor(6);
            hStack1->Add(hSigPassL1TwoPhotonsDraw);
            legend1->AddEntry(hSigPassL1TwoPhotonsDraw,label12Pho.c_str(),"pl");
            
            TH1F *hSigPassDraw = DrawOverflow(hSigPass[e]);
            hSigPassDraw->SetLineColor(3);
            hStack1->Add(hSigPassDraw);
            legend1->AddEntry(hSigPassDraw,label1P.c_str(),"pl");
            
            TH1F *hSigPassXGBDraw = DrawOverflow(hSigPassXGB[e]);
            hSigPassXGBDraw->SetLineColor(5);
            hStack1->Add(hSigPassXGBDraw);
            legend1->AddEntry(hSigPassXGBDraw,label1XGB.c_str(),"pl");
           
            TH1F *hSigPassTightMassDraw = DrawOverflow(hSigPassTightMass[e]);
            hSigPassTightMassDraw->SetLineColor(3);
            hSigPassTightMassDraw->SetLineStyle(10);
            legend1->AddEntry(hSigPassTightMassDraw,label1PTightMass.c_str(),"pl");
            hStack1->Add(hSigPassTightMassDraw);
            
            TH1F *hSigPassStdDraw = DrawOverflow(hSigPassStd[e]);
            hSigPassStdDraw->SetLineColor(4);
            hStack1->Add(hSigPassStdDraw);
            legend1->AddEntry(hSigPassStdDraw,label1PStd.c_str(),"pl");
            
            TH1F *hSigPassStdManualDraw = DrawOverflow(hSigPassStdManual[e]);
            hSigPassStdManualDraw->SetLineColor(4);
            hSigPassStdManualDraw->SetLineStyle(9);
            hStack1->Add(hSigPassStdManualDraw);
            legend1->AddEntry(hSigPassStdManualDraw,label1PStdM.c_str(),"pl");
            
            can->cd(1);
            gPad->SetGrid();
            if (varNames[i] == "ecalPFIso" || varNames[i] == "trkIsoPho")gPad->SetLogy();
            else gPad->SetLogy(0);
            
            hStack1->Draw("nostackhist");
            legend1->Draw("same");
            
            snprintf(buff8, sizeof(buff8), "(%0.1f)", hBkgAll[e]->GetSumOfWeights());
            string nEvents2 = buff8;
            label2A += nEvents1;
            
            snprintf(buff9, sizeof(buff9), "(%0.1f)", hBkgPass[e]->GetSumOfWeights());
            string nEvents2P = buff9;
            label2P += nEvents2P;
            
            snprintf(buff10, sizeof(buff10), "(%0.1f)", hBkgPassXGB[e]->GetSumOfWeights());
            string nEvents2XGB = buff10;
            label2XGB += nEvents2XGB;
            
            snprintf(buff11, sizeof(buff11), "(%0.1f)", hBkgPassL1TwoPhotons[e]->GetSumOfWeights());
            string nEvents2PL12 = buff11;
            label22Pho += nEvents2PL12;
            
            snprintf(buff12, sizeof(buff12), "(%0.1f)", hBkgPassTightMass[e]->GetSumOfWeights());
            string nEvents2PTight = buff12;
            label2PTightMass += nEvents2PTight;
            
            snprintf(buff13, sizeof(buff13), "(%0.1f)", hBkgPassStd[e]->GetSumOfWeights());
            string nEvents2PStd = buff13;
            label2PStd += nEvents2PStd;
            
            snprintf(buff14, sizeof(buff14), "(%0.1f)", hBkgPassStdManual[e]->GetSumOfWeights());
            string nEvents2PStdM = buff14;
            label2PStdM += nEvents2PStdM;
            
            TH1F *hBkgAllDraw = DrawOverflow(hBkgAll[e]);
            TH1F *hStackHistoBkg = (TH1F*)hBkgAllDraw->Clone();
            hStackHistoBkg->Reset();
            hStackHistoBkg->GetXaxis()->SetRange(0,hStackHistoBkg->GetNbinsX());
            hStack2->SetHistogram(hStackHistoBkg);
            
            hStack2->GetHistogram()->GetXaxis()->SetTitle(varNames[i].c_str());
            hStack2->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
            hStack2->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
            hStack2->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
            

            TH1F *hBkgPassL1TwoPhotonsDraw = DrawOverflow(hBkgPassL1TwoPhotons[e]);
            hBkgPassL1TwoPhotonsDraw->SetLineColor(6);
            hStack2->Add(hBkgPassL1TwoPhotonsDraw);
            legend2->AddEntry(hBkgPassL1TwoPhotonsDraw,label22Pho.c_str(),"pl");
            
            TH1F *hBkgPassDraw = DrawOverflow(hBkgPass[e]);
            hBkgPassDraw->SetLineColor(3);
            hStack2->Add(hBkgPassDraw);
            legend2->AddEntry(hBkgPassDraw,label2P.c_str(),"pl");
            
            TH1F *hBkgPassXGBDraw = DrawOverflow(hBkgPassXGB[e]);
            hBkgPassXGBDraw->SetLineColor(5);
            hStack2->Add(hBkgPassXGBDraw);
            legend2->AddEntry(hBkgPassXGBDraw,label2XGB.c_str(),"pl");
            
            TH1F *hBkgPassTightMassDraw = DrawOverflow(hBkgPassTightMass[e]);
            hBkgPassTightMassDraw->SetLineColor(3);
            hBkgPassTightMassDraw->SetLineStyle(10);
            legend2->AddEntry(hBkgPassTightMassDraw,label2PTightMass.c_str(),"pl");
            hStack2->Add(hBkgPassTightMassDraw);
            
            TH1F *hBkgPassStdDraw = DrawOverflow(hBkgPassStd[e]);
            hBkgPassStdDraw->SetLineColor(4);
            hStack2->Add(hBkgPassStdDraw);
            legend2->AddEntry(hBkgPassStdDraw,label2PStd.c_str(),"pl");
            
            TH1F *hBkgPassStdManualDraw = DrawOverflow(hBkgPassStdManual[e]);
            hBkgPassStdManualDraw->SetLineColor(4);
            hBkgPassStdManualDraw->SetLineStyle(9);
            hStack2->Add(hBkgPassStdManualDraw);
            legend2->AddEntry(hBkgPassStdManualDraw,label2PStdM.c_str(),"pl");
            
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
            
        }//eta loop
        for (int p = 0; p < 3; p++){
            hSigAll[p]->Delete();
            hSigPass[p]->Delete();
            hSigPassXGB[p]->Delete();
            hSigPassTightMass[p]->Delete();
            hSigPassL1TwoPhotons[p]->Delete();
            hSigPassStd[p]->Delete();
            hSigPassStdManual[p]->Delete();
            
            hBkgAll[p]->Delete();
            hBkgPass[p]->Delete();
            hBkgPassXGB[p]->Delete();
            hBkgPassTightMass[p]->Delete();
            hBkgPassL1TwoPhotons[p]->Delete();
            hBkgPassStd[p]->Delete();
            hBkgPassStdManual[p]->Delete();
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

//
//TH1F *DrawOverUnderflow(TH1F* h){
//    //function to paint the histogram h with an extra bin for overflows
//    UInt_t nx    = h->GetNbinsX()+1;
//    Double_t *xbins= new Double_t[nx+2];
//    for (UInt_t i=1;i<nx;i++)
//        xbins[i]=h->GetBinLowEdge(i);
//    xbins[0]=xbins[1]-h->GetBinWidth(nx);
//    xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
//    //book a temporary histogram having extra bin for overflow, and another for underflow
//    TH1F *htmp = new TH1F(h->GetName(), h->GetTitle(), nx, xbins);
//    htmp->Sumw2();
//    //fill the new histogram including the overflows
//    for (UInt_t i=0; i<=nx+1; i++) {
//        htmp->SetBinContent(htmp->FindBin(htmp->GetBinCenter(i)),h->GetBinContent(i));
//    }
//    htmp->SetBinContent(htmp->FindBin(h->GetBinLowEdge(1)-1), h->GetBinContent(0));
//    // Restore the number of entries
//    //htmp->SetEntries(h->GetEffectiveEntries());
//    return htmp;
//}

//TH1F *DrawOverUnderflowTest(TH1F* h, double limLow, double limHigh){
//
//    int nBins = h->GetNbinsX();
//
//    double overFlowVal = h->GetBinContent(h->GetNbinsX()+1);
//    double underFlowVal = h->GetBinContent(0);
//
//    double binWidth = h->GetBinWidth();
//    double limLowNew = limLow - binWidth;
//
//    TH1F* hNew =
//
//
//
//}






