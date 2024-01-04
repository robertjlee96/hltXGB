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

double calcMassXGB(TTreeReaderArray<float>* xgbScore,double xgbScoreB, double xgbScoreE,TTreeReaderArray<float>* et, TTreeReaderArray<float>* eta, TTreeReaderArray<float>* phi, int nEgsIn, int* pho1, int* pho2);
void diphotonEtaSplitFindXGBCut(){
    gROOT->Reset();
    gStyle->SetPalette(1);
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitle(0);
    
    int nBinsPlot = 1002; // From 0 to 1 in xgbScore with one overflow and one underflow
    double binFactor = 1;// Increase binning by factor of X. Get binFactor extra bins above and below plot
    
    //double prop = 0.5*0.163225;
    //double propPresel = 0.0624773;
    
    int nBins = nBinsPlot*binFactor;
    double binSize = 1.0/(nBinsPlot-2); // Plotting bin size
    double plotLow = 0.0 - binSize; //Min value to have 1 bin of size binsize below 0
    double plotHigh = 1.0 + binSize; // Same for max
    double binSizeFine = (plotHigh-plotLow)/nBins;
    
    string fileName = "validationNTuples/0606/GGH13andData_AllPhotonsValidated.root";
    
    TH1F *hSigAll = new TH1F ("hSigAll","",nBins,plotLow,plotHigh);
    TH1F *hSigHLTAll = new TH1F ("hSigHLTAll","",nBins,plotLow,plotHigh);
    TH1F *hBkgAll = new TH1F ("hBkgAll","",nBins,plotLow,plotHigh);
    TH1F *hBkgHLTAll = new TH1F ("hBkgHLTAll","",nBins,plotLow,plotHigh);
    
    TH1F *hSigBB = new TH1F ("hSigBB","",nBins,plotLow,plotHigh);
    TH1F *hSigHLTBB = new TH1F ("hPrompHLTBB","",nBins,plotLow,plotHigh);
    TH1F *hBkgBB = new TH1F ("hBkgBB","",nBins,plotLow,plotHigh);
    TH1F *hBkgHLTBB = new TH1F ("hBkgHLTBB","",nBins,plotLow,plotHigh);
    
    TH1F *hSigBE= new TH1F ("hSigBE","",nBins,plotLow,plotHigh);
    TH1F *hSigHLTBE = new TH1F ("hPrompHLTBE","",nBins,plotLow,plotHigh);
    TH1F *hBkgBE = new TH1F ("hBkgBE","",nBins,plotLow,plotHigh);
    TH1F *hBkgHLTBE = new TH1F ("hBkgHLTBE","",nBins,plotLow,plotHigh);
    
    TH1F *hSigEE= new TH1F ("hSigEE","",nBins,plotLow,plotHigh);
    TH1F *hSigHLTEE = new TH1F ("hPrompHLTEE","",nBins,plotLow,plotHigh);
    TH1F *hBkgEE = new TH1F ("hBkgEE","",nBins,plotLow,plotHigh);
    TH1F *hBkgHLTEE = new TH1F ("hBkgHLTEE","",nBins,plotLow,plotHigh);

    
    TFile *f = new TFile(fileName.c_str());
    
    double maxVal = 0.0;
    int maxInt = -1;
    
    TTreeReader inReaderSig("sigTree", f);
    inReaderSig.Restart();
    
    TTreeReaderArray<Float_t> inEtSig(inReaderSig,"et");
    TTreeReaderArray<Float_t> inEtaSig(inReaderSig, "eta");
    TTreeReaderArray<Float_t> inPhiSig(inReaderSig, "phi");
    TTreeReaderArray<Float_t> inXGBSig(inReaderSig, "xgbScore");
    
    TTreeReaderValue<int> inNEgsSig(inReaderSig,"nEgs");
    TTreeReaderValue<int> inPassFailStdSig(inReaderSig,"passFailStd");
    //Start with all signal photons
    while (inReaderSig.Next()) {
        int iOutXGBSig,jOutXGBSig;
        double massCalcXGBSig = calcMassXGB(&inXGBSig,0.25,0.25,&inEtSig,&inEtaSig,&inPhiSig,*inNEgsSig,&iOutXGBSig,&jOutXGBSig);
        //First for All Eta
        int indexHigh,indexLow;
        if (inXGBSig[iOutXGBSig] > inXGBSig[jOutXGBSig]){
            indexHigh = iOutXGBSig;
            indexLow = jOutXGBSig;
        }
        else if (inXGBSig[iOutXGBSig] < inXGBSig[jOutXGBSig]){
            indexHigh = jOutXGBSig;
            indexLow = iOutXGBSig;
        }
        
        if(massCalcXGBSig > 120.0 && massCalcXGBSig < 128.5){
            hSigAll->Fill(inXGBSig[indexLow]);
            hSigAll->Fill(inXGBSig[indexLow]);
            if (*inPassFailStdSig == 1){
                hSigHLTAll->Fill(inXGBSig[indexLow]);
                hSigHLTAll->Fill(inXGBSig[indexLow]);
            }
            
            if(abs(inEtaSig[iOutXGBSig]) < 1.4442 && abs(inEtaSig[jOutXGBSig]) < 1.4442){
                hSigBB->Fill(inXGBSig[indexLow]);
                hSigBB->Fill(inXGBSig[indexLow]);
                if (*inPassFailStdSig == 1){
                    hSigHLTBB->Fill(inXGBSig[indexLow]);
                    hSigHLTBB->Fill(inXGBSig[indexLow]);
                }
            }
            if((abs(inEtaSig[iOutXGBSig]) > 1.556 && abs(inEtaSig[jOutXGBSig]) < 1.4442) || (abs(inEtaSig[iOutXGBSig]) < 1.4442 && abs(inEtaSig[jOutXGBSig]) > 1.556)){
                hSigBE->Fill(inXGBSig[indexLow]);
                hSigBE->Fill(inXGBSig[indexLow]);
                if (*inPassFailStdSig == 1){
                    hSigHLTBE->Fill(inXGBSig[indexLow]);
                    hSigHLTBE->Fill(inXGBSig[indexLow]);
                }
            }
            if(abs(inEtaSig[iOutXGBSig]) > 1.556 && abs(inEtaSig[jOutXGBSig]) > 1.556){
                hSigEE->Fill(inXGBSig[indexLow]);
                hSigEE->Fill(inXGBSig[indexLow]);
                if (*inPassFailStdSig == 1){
                    hSigHLTEE->Fill(inXGBSig[indexLow]);
                    hSigHLTEE->Fill(inXGBSig[indexLow]);
                }
            }
        }
    }
    TTreeReader inReaderBkg("bkgTree", f);
    inReaderBkg.Restart();
    
    TTreeReaderArray<Float_t> inEtBkg(inReaderBkg,"et");
    TTreeReaderArray<Float_t> inEtaBkg(inReaderBkg, "eta");
    TTreeReaderArray<Float_t> inPhiBkg(inReaderBkg, "phi");
    TTreeReaderArray<Float_t> inXGBBkg(inReaderBkg, "xgbScore");
    
    TTreeReaderValue<int> inNEgsBkg(inReaderBkg,"nEgs");
    TTreeReaderValue<int> inPassFailStdBkg(inReaderBkg,"passFailStd");
    
    while (inReaderBkg.Next()) {
        int iOutXGBBkg,jOutXGBBkg;
        double massCalcXGBBkg = calcMassXGB(&inXGBBkg,0.25,0.25,&inEtBkg,&inEtaBkg,&inPhiBkg,*inNEgsBkg,&iOutXGBBkg,&jOutXGBBkg);
        //First for All Eta
        int indexHigh,indexLow;
        if (inXGBBkg[iOutXGBBkg] > inXGBBkg[jOutXGBBkg]){
            indexHigh = iOutXGBBkg;
            indexLow = jOutXGBBkg;
        }
        else if (inXGBBkg[iOutXGBBkg] < inXGBBkg[jOutXGBBkg]){
            indexHigh = jOutXGBBkg;
            indexLow = iOutXGBBkg;
        }
        
        if(massCalcXGBBkg > 120.0 && massCalcXGBBkg < 128.5){
            hBkgAll->Fill(inXGBBkg[indexLow]);
            hBkgAll->Fill(inXGBBkg[indexLow]);
            if (*inPassFailStdBkg == 1){
                hBkgHLTAll->Fill(inXGBBkg[indexLow]);
                hBkgHLTAll->Fill(inXGBBkg[indexLow]);
            }
            
            if(abs(inEtaBkg[iOutXGBBkg]) < 1.4442 && abs(inEtaBkg[jOutXGBBkg]) < 1.4442){
                hBkgBB->Fill(inXGBBkg[indexLow]);
                hBkgBB->Fill(inXGBBkg[indexLow]);
                if (*inPassFailStdBkg == 1){
                    hBkgHLTBB->Fill(inXGBBkg[indexLow]);
                    hBkgHLTBB->Fill(inXGBBkg[indexLow]);
                }
            }
            if((abs(inEtaBkg[iOutXGBBkg]) > 1.556 && abs(inEtaBkg[jOutXGBBkg]) < 1.4442) || (abs(inEtaBkg[iOutXGBBkg]) < 1.4442 && abs(inEtaBkg[jOutXGBBkg]) > 1.556)){
                hBkgBE->Fill(inXGBBkg[indexLow]);
                hBkgBE->Fill(inXGBBkg[indexLow]);
                if (*inPassFailStdBkg == 1){
                    hBkgHLTBE->Fill(inXGBBkg[indexLow]);
                    hBkgHLTBE->Fill(inXGBBkg[indexLow]);
                }
            }
            if(abs(inEtaBkg[iOutXGBBkg]) > 1.556 && abs(inEtaBkg[jOutXGBBkg]) > 1.556){
                hBkgEE->Fill(inXGBBkg[indexLow]);
                hBkgEE->Fill(inXGBBkg[indexLow]);
                if (*inPassFailStdBkg == 1){
                    hBkgHLTEE->Fill(inXGBBkg[indexLow]);
                    hBkgHLTEE->Fill(inXGBBkg[indexLow]);
                }
            }
        }
    }
    
  
    
    //First, determine cut to match presel for ALL eta.
    int nCuts = nBins;
    double cutStep = 2.0/nCuts;
    double cutVal = 0.0 - cutStep;
    
    double propIDMVACut;
    double equivalentIDMVACut;
    
    int mvaMaxBin = hSigAll->GetXaxis()->FindBin(1) + 1;
    
    double totalIntBkg = hBkgAll->Integral();
    double totalIntSig = hSigAll->Integral();
    
    //double propHLT = 0.5*(hBkgHLTAll->Integral()/totalIntBkg); //FOR BKG PRESEL
    double propHLT = hSigHLTAll->Integral()/totalIntSig; //FOR SIG PRESEL
    cout<<"propHLT = "<<propHLT<<endl;
    
    for(int i = 0; i < nCuts; i++){
        cutVal+= cutStep;
        
        int mvaBin = hBkgAll->GetXaxis()->FindBin(cutVal);
        //double tmpIDMVAInt = hBkgAll->Integral(mvaBin,mvaMaxBin); //FOR BKG EFFICIENCY
	double tmpIDMVAInt = hSigAll->Integral(mvaBin,mvaMaxBin); //FOR SIG EFFICIENCY
        
	//propIDMVACut = tmpIDMVAInt/totalIntBkg; //FOR BKG EFFICIENCY
        propIDMVACut = tmpIDMVAInt/totalIntSig; //FOR SIG EFFICIENCY
        if(propIDMVACut <= propHLT)break;
    }
    cout<<"propIDMVACut = "<<propIDMVACut<<endl;
    int binsOnEachSide = 8;
    
    //Let's attempt to doubly-constrain the numbers here. We will simultaneously constrain the B & E differential sig/bkg efficiencies to be close together (delta) and constrain the total bkg efficiency to be close to the preselection efficiency (epsilon)
    double delta = 0.01; //B & E differential sig/bkg are within this proportion of each other.
    double epsilon = 0.01; //Total bkg efficiency is within the proportion of preselection efficiency.
    
    double cutVals[10002], bbDiffSigOvrBkg[10002], bbTotalEff[10002], beDiffSigOvrBkg[10002], beTotalEff[10002], eeDiffSigOvrBkg[10002], eeTotalEff[10002];
    double bbCut, beCut,eeCut;
    bbCut = 0.0 - cutStep;
    beCut = 0.0 - cutStep;
    eeCut = 0.0 - cutStep;
    
    //Let's first make arrays of the Differential sig/bkg efficiencies, cuts, and total efficiencies.
    for(int i = 0; i < nCuts; i++){
        bbCut += cutStep;
        cutVals[i] = bbCut;
        int mvaBin = hBkgBB->GetXaxis()->FindBin(bbCut);
        bbDiffSigOvrBkg[i] = hSigBB->Integral(mvaBin-binsOnEachSide,mvaBin+binsOnEachSide)/hBkgBB->Integral(mvaBin-binsOnEachSide,mvaBin+binsOnEachSide);
        //barrelTotalEff[i] = hBkgB->Integral(mvaBin,mvaMaxBin)/hBkgB->Integral(); //FOR BKG EFFICIENCY
	bbTotalEff[i] = hSigBB->Integral(mvaBin,mvaMaxBin)/hSigBB->Integral(); //FOR SIG EFFICIENCY
    }
    for(int i = 0; i < nCuts; i++){
        endcapCut+= cutStep;
        int mvaBin = hBkgE->GetXaxis()->FindBin(endcapCut);
        endcapDiffSigOvrBkg[i] = hSigE->Integral(mvaBin-binsOnEachSide,mvaBin+binsOnEachSide)/hBkgE->Integral(mvaBin-binsOnEachSide,mvaBin+binsOnEachSide);
        //endcapTotalEff[i] = hBkgE->Integral(mvaBin,mvaMaxBin)/hBkgE->Integral(); //FOR BKG EFFICIENCY
        endcapTotalEff[i] = hSigE->Integral(mvaBin,mvaMaxBin)/hSigE->Integral(); //FOR SIG EFFICIENCY
    }
    
    //Now let's do a double-loop. Go through each Barrel cut, then each endcap cut for each barrel cut.
    double finalCutB, finalCutE;
    
    double targetTotal = propHLT;
    
    for(int i = 0; i < nCuts; i++){
        double tmpCutB = cutVals[i];
        double tmpDiffEffB = barrelDiffSigOvrBkg[i];
        double tmpTotalEffB = barrelTotalEff[i];
        for(int j = 0; j < nCuts; j++){
            double tmpCutE = cutVals[j];
            double tmpDiffEffE = endcapDiffSigOvrBkg[j];
            double tmpTotalEffE = endcapTotalEff[j];
	    //double totalEffAllEta = (tmpTotalEffB*hBkgB->Integral() + tmpTotalEffE*hBkgE->Integral())/hBkgAll->Integral(); //FOR BKG EFFICIENCY
            double totalEffAllEta = (tmpTotalEffB*hSigB->Integral() + tmpTotalEffE*hSigE->Integral())/hSigAll->Integral(); //FOR SIG EFFICIENCY
//            cout<<"totalEffAllEta = "<<totalEffAllEta<<endl;
            if (abs(totalEffAllEta-targetTotal)/targetTotal < epsilon && abs(tmpDiffEffB-tmpDiffEffE)/tmpDiffEffB < delta){
                cout<<"For cutB = "<<tmpCutB<<", & cutE = "<<tmpCutE<<endl;
                cout<<"tmpDiffEffB = "<<tmpDiffEffB<<", tmpTotalEffB = "<<tmpTotalEffB<<endl;
                cout<<"tmpDiffEffE = "<<tmpDiffEffE<<", tmpTotalEffE = "<<tmpTotalEffE<<endl;
                cout<<"Total Eff = "<<totalEffAllEta<<endl;
                finalCutB = tmpCutB;
                finalCutE = tmpCutE;
                break;
            }
        }
    }
    cout<<"Final Cuts: Barrel = "<<finalCutB<<", Endcap = "<<finalCutE<<endl;
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
                        double dR = calcDR((*eta)[i],(*eta)[j],(*phi)[i],(*phi)[j])
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
                                    double dR = calcDR((*eta)[i],(*eta)[j],(*phi)[i],(*phi)[j])
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
