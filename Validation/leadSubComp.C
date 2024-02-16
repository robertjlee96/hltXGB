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
void leadSubComp(){
    gROOT->Reset();
    gStyle->SetPalette(1);
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitle(0);
  
    string fileName = "validationNTuples/0425/GGHandData_MD13LR07.root";
    string genTitleString = "Data ";
    string treeString = "bkgTree";
    string dirStr ="diphotonVarPlots/0425/MD13LR07_Data_RelaxedV24/";
    //string dirStr ="varPlots/0130/DataRelaxedUnseededChoose2/";
  
    int nEta;
 
    nEta = 4;
    string etaLabels[4] = {"Barrel-Barrel","Barrel-Endcap","Endcap-Endcap","All #eta"};
    string etaFLabels[4] = {"_BB","_BE","_EE","_All"};
    
    string outNameGen = dirStr;
    
    int nVars = 20;
    string varNames[20] = {"xgbScore","nEgs","triggerBits","mass","et","energy","rawEnergy","eta","etaWidth","phi","phiWidth","sigmaIEtaIEta","ecalPFIso","hcalPFIso","r9HLT","r9Recalc","s4","hOvrE","trkIso","trkIsoPho"};
    double limsLow[20] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,-3.15,0.0,-3.15,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double limsHigh[20] = {1.0,5,4,500,300,600,600,3.15,0.1,3.15,0.2,0.2,20.0,20.0,1.2,1.2,1.2,1.0,20.0,20.0};
   
    int nBins[20] = {400,5,4,500,600,600,600,400,400,400,400,400,400,400,400,400,400,500,400,400};
    
    TCanvas *can = new TCanvas ("can","can",10,10,1600,900);
    //for(int i = 0; i < 1; i++){
    for(int i = 0; i < nVars; i++){
    //for(int i = 16; i < nVars; i++){
        TH1F *hFirstAll[4];
        TH1F *hSecondAll[4];

        TH1F *hFirstPass[4];
        TH1F *hSecondPass[4];
       
        TH1F *hFirstPassLooseCut[4];
        TH1F *hSecondPassLooseCut[4];

        TH1F *hFirstPassTightCut[4];
        TH1F *hSecondPassTightCut[4];
	
        TH1F *hFirstPassStd[4];
        TH1F *hSecondPassStd[4];

        TH1F *hFirstFail[4];
        TH1F *hSecondFail[4];

        for (int e = 0; e < nEta; e++){
            hFirstAll[e] = new TH1F ("hFirstAll","",nBins[i],limsLow[i],limsHigh[i]);
            hSecondAll[e] = new TH1F ("hSecondAll","",nBins[i],limsLow[i],limsHigh[i]);

            hFirstPass[e] = new TH1F ("hFirstPass","",nBins[i],limsLow[i],limsHigh[i]);
            hSecondPass[e] = new TH1F ("hSecondPass","",nBins[i],limsLow[i],limsHigh[i]);

            hFirstPassLooseCut[e] = new TH1F ("hFirstPassLooseCut","",nBins[i],limsLow[i],limsHigh[i]);
            hSecondPassLooseCut[e] = new TH1F ("hSecondPassLooseCut","",nBins[i],limsLow[i],limsHigh[i]);

            hFirstPassTightCut[e] = new TH1F ("hFirstPassTightCut","",nBins[i],limsLow[i],limsHigh[i]);
            hSecondPassTightCut[e] = new TH1F ("hSecondPassTightCut","",nBins[i],limsLow[i],limsHigh[i]);

            hFirstPassStd[e] = new TH1F ("hFirstPassStd","",nBins[i],limsLow[i],limsHigh[i]);
            hSecondPassStd[e] = new TH1F ("hSecondPassStd","",nBins[i],limsLow[i],limsHigh[i]);

            hFirstFail[e] = new TH1F ("hFirstFail","",nBins[i],limsLow[i],limsHigh[i]);
            hSecondFail[e] = new TH1F ("hSecondFail","",nBins[i],limsLow[i],limsHigh[i]);
	   
        }
        for (int e = 0; e < nEta; e++){
                                
            string outName = outNameGen + varNames[i] + etaFLabels[e];
            
            TFile *f = new TFile(fileName.c_str());
    
            string varNameTmp;
            if (e != 3){
                TTreeReader inReader(treeString.c_str(), f);
                inReader.Restart();
                varNameTmp = varNames[i];
        
                TTreeReaderArray<Float_t> inVar(inReader, varNameTmp.c_str());
                TTreeReaderArray<Float_t> inXGB(inReader,"xgbScore");
                TTreeReaderArray<Float_t> inE(inReader,"energy");
                TTreeReaderArray<Float_t> inM(inReader,"mass");
                TTreeReaderValue<Float_t> inNEgs(inReader,"nEgs");
                TTreeReaderArray<Float_t> inEt(inReader,"et");
		
                TTreeReaderArray<Float_t> inEta(inReader, "eta");

                TTreeReaderValue<int> inPassFailStd(inReader,"passFailStd");
                TTreeReaderValue<int> inPassFailL1Double(inReader,"passFailL1Double");
                TTreeReaderValue<int> inPassFailL1Single(inReader,"passFailL1Single");
                  
                //FIRST handle standard variables
                int q = 0;
                while (inReader.Next()) {
                    if(*inNEgs > 1 && ((e == 0 && abs(inEta[0]) < 1.5 && abs(inEta[1]) < 1.5)
                        || (e == 1 && ((abs(inEta[0]) > 1.5 && abs(inEta[1]) < 1.5) || (abs(inEta[0]) < 1.5 && abs(inEta[1]) > 1.5)))
				       || (e == 2 && abs(inEta[0]) > 1.5 && abs(inEta[1]) > 1.5))){
		    //if(varNameTmp == "hOvrE" || varNameTmp == "ecalPFIso" || varNameTmp == "trkIsoPho"){
                           if(varNameTmp == "blahblahblah"){
                               double varMax;
                               if(inVar[0] < inVar[1]) varMax = inVar[0];
                               else varMax = inVar[1];
					       
                               if ((varNameTmp != "mass") || (varNameTmp == "mass")){
                                   hFirstAll[e]->Fill(varMax);
                                   hSecondAll[e]->Fill(varMax);
			
                                   if(inM[0] > 60){
                                       hFirstPass[e]->Fill(varMax);
                                       hSecondPass[e]->Fill(varMax);
                                   }
                                   else{
                                       hFirstFail[e]->Fill(varMax);
                                       hSecondFail[e]->Fill(varMax);
                                   }

                                   if(inXGB[0] > 0.1 && inXGB[1] > 0.1){
                                       hFirstPassLooseCut[e]->Fill(varMax);
                                       hSecondPassLooseCut[e]->Fill(varMax);
                                   }

                                   if(inXGB[0] > 0.5 && inXGB[1] > 0.5){
                                       hFirstPassTightCut[e]->Fill(varMax);
                                       hSecondPassTightCut[e]->Fill(varMax);
                                   }

                                   if(*inPassFailStd == 1){
                                       hFirstPassStd[e]->Fill(varMax);
                                       hSecondPassStd[e]->Fill(varMax);
                                   }
	        
                               }//All vars besides mass have a 95 mass cut
                           }//Plot maximum value for some vars
                           else{
                               if ((varNameTmp != "mass") || (varNameTmp == "mass")){
                                   hFirstAll[e]->Fill(inVar[0]);
                                   hSecondAll[e]->Fill(inVar[1]);
			
                                   if(inM[0] > 60){
                                       hFirstPass[e]->Fill(inVar[0]);
                                       hSecondPass[e]->Fill(inVar[1]);
                                   }
                                   else{
                                       hFirstFail[e]->Fill(inVar[0]);
                                       hSecondFail[e]->Fill(inVar[1]);
                                   }
			
                                   if(inXGB[0] > 0.1 && inXGB[1] > 0.1){
                                       hFirstPassLooseCut[e]->Fill(inVar[0]);
                                       hSecondPassLooseCut[e]->Fill(inVar[1]);
                                   }

                                   if(inXGB[0] > 0.5 && inXGB[1] > 0.5){
                                       hFirstPassTightCut[e]->Fill(inVar[0]);
                                       hSecondPassTightCut[e]->Fill(inVar[1]);
                                   }

                                   if(*inPassFailStd == 1){
                                       hFirstPassStd[e]->Fill(inVar[0]);
                                       hSecondPassStd[e]->Fill(inVar[1]);
                                   }
		
                               }//All vars besides mass have a 95 mass cut
                           }//All other vars plotted normaly
                       }//Eta Selection if
                }//While inReader.Next()
            }
            if (e == 3){
                for (int p = 0; p < 3; p++){
                    hFirstAll[e]->Add(hFirstAll[p]);
                    hSecondAll[e]->Add(hSecondAll[p]);

                    hFirstPass[e]->Add(hFirstPass[p]);
                    hSecondPass[e]->Add(hSecondPass[p]);
                    
                    hFirstPassLooseCut[e]->Add(hFirstPassLooseCut[p]);
                    hSecondPassLooseCut[e]->Add(hSecondPassLooseCut[p]);

                    hFirstPassTightCut[e]->Add(hFirstPassTightCut[p]);
                    hSecondPassTightCut[e]->Add(hSecondPassTightCut[p]);
        
                    hFirstPassStd[e]->Add(hFirstPassStd[p]);
                    hSecondPassStd[e]->Add(hSecondPassStd[p]);

                    hFirstFail[e]->Add(hFirstFail[p]);
                    hSecondFail[e]->Add(hSecondFail[p]);

                }
            }
            can->Clear();

            can->Divide(2,1);

            string plotTitle1 = genTitleString + "First Photon for " + etaLabels[e] + ";" + varNames[i];
            THStack *hStack1 = new THStack("hStack1",plotTitle1.c_str());
            hStack1->SetHistogram(new TH1F("hstot1","",nBins[i],limsLow[i],limsHigh[i]));

            string plotTitle2 = genTitleString + "Second Photon for " + etaLabels[e] + ";" + varNames[i];
            THStack *hStack2 = new THStack("hStack2",plotTitle2.c_str());
            hStack2->SetHistogram(new TH1F("hstot2","",nBins[i],limsLow[i],limsHigh[i]));
	    
            TLegend *legend1;
            if (varNames[i] == "r9HLT"|| varNames[i] == "r9Recalc" || varNames[i] == "s4") legend1 = new TLegend(0.01,0.65,0.4,0.9,"","brNDC");
            else legend1 = new TLegend(0.50,0.65,0.89,0.9,"","brNDC");

            TLegend *legend2;
            if (varNames[i] == "r9HLT"|| varNames[i] == "r9Recalc" || varNames[i] == "s4") legend2 = new TLegend(0.01,0.65,0.4,0.9,"","brNDC");
            else legend2 = new TLegend(0.50,0.65,0.89,0.9,"","brNDC");

            string label1A = "All";
            string label1PL1 = "xgbScore > 0.1";
            string label1PL12 = "xgbScore > 0.5";
            string label1P = "Passing L1 & Relaxed Diphoton Cuts";
            string label1PStd = "Passing L1 & Std. HLT";
            string label1F = "Failing L1 or Relaxed Cuts";

            string label2A = "All";
            string label2PL1 = "xgbScore > 0.1";
            string label2PL12 = "xgbScore > 0.5";
            string label2P = "Passing L1 & Relaxed Diphoton Cuts";
            string label2PStd = "Passing L1 & Std. HLT";
            string label2F = "Failing L1 or Relaxed Cuts";
                            
            char buff1[100], buff2[100], buff3[100], buff4[100], buff5[100], buff6[100], buff7[100], buff8[100], buff9[100], buff10[100], buff11[100], buff12[100];
            
            snprintf(buff1, sizeof(buff1), "(%0.1f)", hFirstAll[e]->GetSumOfWeights());
            string nEvents1 = buff1;
            label1A += nEvents1;

            snprintf(buff2, sizeof(buff2), "(%0.1f)", hFirstPassLooseCut[e]->GetSumOfWeights());
            string nEvents1PL1 = buff2;
            label1PL1 += nEvents1PL1;

            snprintf(buff3, sizeof(buff3), "(%0.1f)", hFirstPassTightCut[e]->GetSumOfWeights());
            string nEvents1PL12 = buff3;
            label1PL12 += nEvents1PL12;

            snprintf(buff4, sizeof(buff4), "(%0.1f)", hFirstPass[e]->GetSumOfWeights());
            string nEvents1P = buff4;
            label1P += nEvents1P;

            snprintf(buff5, sizeof(buff5), "(%0.1f)", hFirstPassStd[e]->GetSumOfWeights());
            string nEvents1PStd = buff5;
            label1PStd += nEvents1PStd;

            snprintf(buff6, sizeof(buff6), "(%0.1f)", hFirstFail[e]->GetSumOfWeights());
            string nEvents1F = buff6;
            label1F += nEvents1F;

            TH1F *hFirstAllDraw = DrawOverflow(hFirstAll[e]);
            TH1F *hStackHistoFirst = (TH1F*)hFirstAllDraw->Clone();
            hStackHistoFirst->Reset();
            hStackHistoFirst->GetXaxis()->SetRange(0,hStackHistoFirst->GetNbinsX());
            hStack1->SetHistogram(hStackHistoFirst);
            
            hStack1->GetHistogram()->GetXaxis()->SetTitle(varNames[i].c_str());
            //hStack->GetHistogram()->GetXaxis()->SetNdivisions(-520);
            hStack1->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
            hStack1->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
            
            hStack1->Add(hFirstAllDraw);
            legend1->AddEntry(hFirstAllDraw,label1A.c_str(),"pl");

            TH1F *hFirstPassLooseCutDraw = DrawOverflow(hFirstPassLooseCut[e]);
            hFirstPassLooseCutDraw->SetLineColor(7);
            hStack1->Add(hFirstPassLooseCutDraw);
            legend1->AddEntry(hFirstPassLooseCutDraw,label1PL1.c_str(),"pl");

	    TH1F *hFirstPassTightCutDraw = DrawOverflow(hFirstPassTightCut[e]);
	    hFirstPassTightCutDraw->SetLineColor(6);
	    hStack1->Add(hFirstPassTightCutDraw);
	    legend1->AddEntry(hFirstPassTightCutDraw,label1PL12.c_str(),"pl");

	    //TH1F *hFirstPassDraw = DrawOverflow(hFirstPass[e]);
	    //hFirstPassDraw->SetLineColor(3);
	    //hStack1->Add(hFirstPassDraw);
	    //legend1->AddEntry(hFirstPassDraw,label1P.c_str(),"pl");

	    TH1F *hFirstPassStdDraw = DrawOverflow(hFirstPassStd[e]);
	    hFirstPassStdDraw->SetLineColor(4);
	    hStack1->Add(hFirstPassStdDraw);
	    legend1->AddEntry(hFirstPassStdDraw,label1PStd.c_str(),"pl");

	    //TH1F *hFirstFailDraw = DrawOverflow(hFirstFail[e]);
	    //hFirstFailDraw->SetLineColor(2);
	    //hStack1->Add(hFirstFailDraw);
	    //legend1->AddEntry(hFirstFailDraw,label1F.c_str(),"pl");

	    can->cd(1);

	    hStack1->Draw("nostackhist");
	    legend1->Draw("same");

	    snprintf(buff7, sizeof(buff7), "(%0.1f)", hSecondAll[e]->GetSumOfWeights());
            string nEvents2 = buff7;
            label2A += nEvents2;

	    snprintf(buff8, sizeof(buff8), "(%0.1f)", hSecondPassLooseCut[e]->GetSumOfWeights());
            string nEvents2PL1 = buff8;
            label2PL1 += nEvents2PL1;

	    snprintf(buff9, sizeof(buff9), "(%0.1f)", hSecondPassTightCut[e]->GetSumOfWeights());
            string nEvents2PL12 = buff9;
            label2PL12 += nEvents2PL12;

	    snprintf(buff10, sizeof(buff10), "(%0.1f)", hSecondPass[e]->GetSumOfWeights());
            string nEvents2P = buff10;
            label2P += nEvents2P;

	    snprintf(buff11, sizeof(buff11), "(%0.1f)", hSecondPassStd[e]->GetSumOfWeights());
            string nEvents2PStd = buff11;
            label2PStd += nEvents2PStd;

	    snprintf(buff12, sizeof(buff12), "(%0.1f)", hSecondFail[e]->GetSumOfWeights());
            string nEvents2F = buff12;
            label2F += nEvents2F;

	    TH1F *hSecondAllDraw = DrawOverflow(hSecondAll[e]);
            TH1F *hStackHistoSecond = (TH1F*)hSecondAllDraw->Clone();
            hStackHistoSecond->Reset();
            hStackHistoSecond->GetXaxis()->SetRange(0,hStackHistoSecond->GetNbinsX());
            hStack2->SetHistogram(hStackHistoSecond);
            
            hStack2->GetHistogram()->GetXaxis()->SetTitle(varNames[i].c_str());
            //hStack->GetHistogram()->GetXaxis()->SetNdivisions(-520);
            hStack2->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
            hStack2->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
            
            hStack2->Add(hSecondAllDraw);
            legend2->AddEntry(hSecondAllDraw,label2A.c_str(),"pl");

	    TH1F *hSecondPassLooseCutDraw = DrawOverflow(hSecondPassLooseCut[e]);
	    hSecondPassLooseCutDraw->SetLineColor(7);
	    hStack2->Add(hSecondPassLooseCutDraw);
	    legend2->AddEntry(hSecondPassLooseCutDraw,label2PL1.c_str(),"pl");
	    
	    TH1F *hSecondPassTightCutDraw = DrawOverflow(hSecondPassTightCut[e]);
	    hSecondPassTightCutDraw->SetLineColor(6);
	    hStack2->Add(hSecondPassTightCutDraw);
	    legend2->AddEntry(hSecondPassTightCutDraw,label2PL12.c_str(),"pl");

	    //TH1F *hSecondPassDraw = DrawOverflow(hSecondPass[e]);
	    //hSecondPassDraw->SetLineColor(3);
	    //hStack2->Add(hSecondPassDraw);
	    //legend2->AddEntry(hSecondPassDraw,label2P.c_str(),"pl");

	    TH1F *hSecondPassStdDraw = DrawOverflow(hSecondPassStd[e]);
	    hSecondPassStdDraw->SetLineColor(4);
	    hStack2->Add(hSecondPassStdDraw);
	    legend2->AddEntry(hSecondPassStdDraw,label2PStd.c_str(),"pl");

	    //TH1F *hSecondFailDraw = DrawOverflow(hSecondFail[e]);
	    //hSecondFailDraw->SetLineColor(2);
	    //hStack2->Add(hSecondFailDraw);
	    //legend2->AddEntry(hSecondFailDraw,label2F.c_str(),"pl");

	    can->cd(2);

	    hStack2->Draw("nostackhist");
	    legend2->Draw("same");

            can->SaveAs((outName+".png").c_str());
            can->SaveAs((outName+".root").c_str());
            
            can->Clear();

            f->Close();
            f->Delete();
	    /*
	    hFirstAllDraw->Delete();
	    hFirstPassL1Draw->Delete();
	    hFirstPassTightCutDraw->Delete();
	    hFirstPassDraw->Delete();
	    hFirstPassStdDraw->Delete();
	    hFirstFailDraw->Delete();

	    hSecondAllDraw->Delete();
	    hSecondPassL1Draw->Delete();
	    hSecondPassTightCutDraw->Delete();
	    hSecondPassDraw->Delete();
	    hSecondPassStdDraw->Delete();
	    hSecondFailDraw->Delete();
	    */
        }
        
        for (int p = 0; p < 4; p++){
            hFirstAll[p]->Delete();
            hFirstPass[p]->Delete();
            hFirstPassLooseCut[p]->Delete();
            hFirstPassTightCut[p]->Delete();
            hFirstPassStd[p]->Delete();
            hFirstFail[p]->Delete();

            hSecondAll[p]->Delete();
            hSecondPass[p]->Delete();            
            hSecondPassLooseCut[p]->Delete();
            hSecondPassTightCut[p]->Delete();
            hSecondPassStd[p]->Delete();
            hSecondFail[p]->Delete();
        }
    }

}

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
    htmp->SetBinContent(htmp->FindBin(h->GetBinLowEdge(1)-1), h->GetBinContent(0));
    return htmp;
}

