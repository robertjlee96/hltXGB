#!/usr/bin/env python -W ignore::DeprecationWarning
import numpy as np
import uproot
from sklearn.model_selection import train_test_split
from sklearn.utils.fixes import _Iterable as Iterable, _Sized as Sized
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import accuracy_score
from ROOT import TFile, TTree, TGraph
from ROOT import gROOT, AddressOf
import xgboost as xgb
from xgboost import XGBClassifier
import matplotlib 
matplotlib.use('agg') #Sketchy way of using alternate backend for problematic Tkinter package (FIX?)
import matplotlib.pyplot as plt
import diphotonUtils
#import oldXGBoost2tmva
import time,pickle
from tqdm import tqdm
import sys
from tempfile import TemporaryFile
import time
from array import array
#import joblib
#import gc

# for sklearn, see
np.random.seed(1337)

inFileName ="../NTuples/outNTuples/GGH13andData_HLTSelectionTest_DiphotonValidation_0602.root"
modelNamesB = []
mNamesGen =  []
mNamesGen.append("0516/M9LR10_GGH13andData_M60_0516")

modelNamesB.append("../Training/barrelOut/" + mNamesGen[0] + "_Barrel.model")
#modelNamesB.append("../Training/barrelOut/Model_NewGGH_MD13LR07_M90_Barrel_0504.Model")
#modelNamesB.append("../Training/barrelOut/Model_NewGGH_MD17LR05_M90_Barrel_0504.Model")

modelNamesE = []
modelNamesE.append("../Training/endcapOut/" + mNamesGen[0] + "_Endcap.model")
#modelNamesE.append("../Training/endcapOut/Model_NewGGH_MD9LR10_M90_Endcap_0504.model")
#modelNamesE.append("../Training/endcapOut/Model_NewGGH_MD9LR10_M90_Endcap_0504.model")

fNames = []
fNames.append("validationNTuples/0602/GGH13andData_HLTSelectionTes_DiphotonValidation")
#fNames.append("validationNTuples/0504/newGGHandData_MD13LR07")
#fNames.append("validationNTuples/0504/newGGHandData_MD17LR05")

fin = uproot.open(inFileName)
prompt = fin['sigTree']
fake = fin['bkgTree']

## for barrel
geometry_selection = lambda tree: np.abs(tree.array('eta')) < 2.5
ptCut = lambda tree: tree.array('et') > 0.0
#massCut = lambda tree: tree.array('mass') > 60.0
massCut = lambda tree: tree.array('mass') > 0.0

#diphotonUtils.load_file(fin,geometry_selection,ptCut,massCut)

inputValuesSigLead, inputValuesSigSub, inputValuesBkgLead, inputValuesBkgSub, targetValuesSigLead, targetValuesSigSub, targetValuesBkgLead,targetValuesBkgSub, origWeightsSigLead, origWeightsSigSub,origWeightsBkgLead,origWeightsBkgSub, varValuesSigLead, varValuesSigSub, varValuesBkgLead, varValuesBkgSub, inputVarNames, varNames = diphotonUtils.load_file(fin,geometry_selection,ptCut,massCut)


print len(inputValuesBkgLead)
print len(inputValuesBkgLead[0])


for i in range(0,len(fNames)):
    fNameRoot = fNames[i] + ".root"
    outFileRoot = TFile(fNameRoot,'recreate')
   
    barrelSigLead = inputValuesSigLead.T[np.abs(varValuesSigLead[7,:]) < 1.5]
    barrelSigLeadVars = varValuesSigLead.T[np.abs(varValuesSigLead[7,:]) < 1.5]
    barrelSigSub = inputValuesSigSub.T[np.abs(varValuesSigSub[7,:]) < 1.5]
    barrelSigSubVars = varValuesSigSub.T[np.abs(varValuesSigSub[7,:]) < 1.5]

    barrelBkgLead = inputValuesBkgLead.T[np.abs(varValuesBkgLead[7,:]) < 1.5]
    barrelBkgLeadVars = varValuesBkgLead.T[np.abs(varValuesBkgLead[7,:]) < 1.5]
    barrelBkgSub = inputValuesBkgSub.T[np.abs(varValuesBkgSub[7,:]) < 1.5]
    barrelBkgSubVars = varValuesBkgSub.T[np.abs(varValuesBkgSub[7,:]) < 1.5]
    
    modelFileB = open(modelNamesB[i],'rb')
    modelB = pickle.load(modelFileB)

    endcapSigLead = inputValuesSigLead.T[np.abs(varValuesSigLead[7,:]) > 1.5]
    endcapSigLeadVars = varValuesSigLead.T[np.abs(varValuesSigLead[7,:]) > 1.5]
    endcapSigSub = inputValuesSigSub.T[np.abs(varValuesSigSub[7,:]) > 1.5]
    endcapSigSubVars = varValuesSigSub.T[np.abs(varValuesSigSub[7,:]) > 1.5]
    
    endcapBkgLead = inputValuesBkgLead.T[np.abs(varValuesBkgLead[7,:]) > 1.5]
    endcapBkgLeadVars = varValuesBkgLead.T[np.abs(varValuesBkgLead[7,:]) > 1.5]
    endcapBkgSub = inputValuesBkgSub.T[np.abs(varValuesBkgSub[7,:]) > 1.5]
    endcapBkgSubVars = varValuesBkgSub.T[np.abs(varValuesBkgSub[7,:]) > 1.5]
   
    #barrelBkgSubVars = varValuesBkgSub.T[np.abs(varValuesBkgSub[7,:]) < 1.4442]
    #endcapBkgSubVars = varValuesBkgSub.T[np.abs(varValuesBkgSub[7,:]) > 1.556]

    modelFileE = open(modelNamesE[i],'rb')
    modelE = pickle.load(modelFileE)

    startTime = time.time()
    
    #18 is the event index determined before splitting. Allows Lead and Sub to be recombined after evaluation
    #So XGB results are an array with the score result and the event index to be recombined into diphoton event
    yPredSigLeadB = np.stack((modelB.predict_proba(barrelSigLead)[:,1],barrelSigLeadVars[:,19]),axis=1)
    yPredSigSubB = np.stack((modelB.predict_proba(barrelSigSub)[:,1],barrelSigSubVars[:,19]),axis=1)
    yPredBkgLeadB = np.stack((modelB.predict_proba(barrelBkgLead)[:,1],barrelBkgLeadVars[:,19]),axis=1)
    yPredBkgSubB = np.stack((modelB.predict_proba(barrelBkgSub)[:,1],barrelBkgSubVars[:,19]),axis=1)

    yPredSigLeadE = np.stack((modelE.predict_proba(endcapSigLead)[:,1],endcapSigLeadVars[:,19]),axis=1)
    yPredSigSubE = np.stack((modelE.predict_proba(endcapSigSub)[:,1],endcapSigSubVars[:,19]),axis=1)
    yPredBkgLeadE = np.stack((modelE.predict_proba(endcapBkgLead)[:,1],endcapBkgLeadVars[:,19]),axis=1)
    yPredBkgSubE = np.stack((modelE.predict_proba(endcapBkgSub)[:,1],endcapBkgSubVars[:,19]),axis=1)

    yPredSigLeadTotal = np.row_stack((yPredSigLeadB,yPredSigLeadE))
    yPredSigLeadTotal[np.argsort(yPredSigLeadTotal[:,1])]
    yPredSigSubTotal = np.row_stack((yPredSigSubB,yPredSigSubE))
    yPredSigSubTotal[np.argsort(yPredSigSubTotal[:,1])]
    yPredBkgLeadTotal = np.row_stack((yPredBkgLeadB,yPredBkgLeadE))
    yPredBkgLeadTotal[np.argsort(yPredBkgLeadTotal[:,1])]
    yPredBkgSubTotal = np.row_stack((yPredBkgSubB,yPredBkgSubE))
    yPredBkgSubTotal[np.argsort(yPredBkgSubTotal[:,1])]


    endTime = time.time()
    timeSpent = endTime - startTime
    timeSpentMin = timeSpent/60.0
    print 'Run Time(min) =  {0:.4f}'.format(timeSpentMin)

    xgbScoreSig = array('f',[0.0,0.0])
    rawEnergySig = array('f',[0.0,0.0])
    r9HLTSig = array('f',[0.0,0.0])
    sigmaIEtaIEtaSig = array('f',[0.0,0.0])
    etaWidthSig = array('f',[0.0,0.0])
    phiWidthSig = array('f',[0.0,0.0])
    s4Sig = array('f',[0.0,0.0])
    trkIsoPhoSig = array('f',[0.0,0.0])
    etaSig = array('f',[0.0,0.0])
    hOvrESig = array('f',[0.0,0.0])
    ecalPFIsoSig = array('f',[0.0,0.0])
    etSig = array('f',[0.0,0.0])
    energySig = array('f',[0.0,0.0])
    dRSig = array('f',[0.0,0.0])
    massSig = array('f',[0.0,0.0])
    nEgsSig = array('f',[0.0,0.0])
    triggerBitsSig = array('f',[0.0,0.0])
    passFailStdSig = array('i',[0])
    passFailL1SingleSig = array('i',[0])
    passFailL1DoubleSig = array('i',[0])

    tSig = TTree("sigTree", "sigTree")
    tSig.Branch( 'xgbScore', xgbScoreSig, 'xgbScore[2]/F' )
    tSig.Branch( 'rawEnergy', rawEnergySig, 'rawEnergy[2]/F' )
    tSig.Branch( 'r9HLT', r9HLTSig, 'r9HLT[2]/F' )
    tSig.Branch( 'sigmaIEtaIEta', sigmaIEtaIEtaSig, 'sigmaIEtaIEta[2]/F' )
    tSig.Branch( 'etaWidth', etaWidthSig, 'etaWidth[2]/F' )
    tSig.Branch( 'phiWidth', phiWidthSig, 'phiWidth[2]/F' )
    tSig.Branch( 's4', s4Sig, 's4[2]/F' )
    tSig.Branch( 'trkIsoPho', trkIsoPhoSig, 'trkIsoPho[2]/F' )
    tSig.Branch( 'eta', etaSig, 'eta[2]/F' )
    tSig.Branch( 'hOvrE', hOvrESig, 'hOvrE[2]/F' )
    tSig.Branch( 'ecalPFIso', ecalPFIsoSig, 'ecalPFIso[2]/F' )
    tSig.Branch( 'et', etSig, 'et[2]/F' )
    tSig.Branch( 'energy', energySig, 'energy[2]/F' )
    tSig.Branch( 'dR', dRSig, 'dR[2]/F')
    tSig.Branch( 'mass', massSig, 'mass[2]/F' )
    tSig.Branch( 'nEgs', nEgsSig, 'nEgs[2]/F' )
    tSig.Branch( 'triggerBits', triggerBitsSig, 'triggerBits[2]/F' )
    tSig.Branch( 'passFailStd', passFailStdSig, 'passFailStd/I' )
    tSig.Branch( 'passFailL1Single', passFailL1SingleSig, 'passFailL1Single/I' )
    tSig.Branch( 'passFailL1Double', passFailL1DoubleSig, 'passFailL1Double/I' )

    sigLen = len(yPredSigLeadTotal)
    
    for n in range(0,sigLen):
        xgbScoreSig[0] = yPredSigLeadTotal[n][0]
        xgbScoreSig[1] = yPredSigSubTotal[n][0]
        rawEnergySig[0] = varValuesSigLead[0][n]
        rawEnergySig[1] = varValuesSigSub[0][n]
        r9HLTSig[0] = varValuesSigLead[1][n]
        r9HLTSig[1] = varValuesSigSub[1][n]
        sigmaIEtaIEtaSig[0] = varValuesSigLead[2][n]
        sigmaIEtaIEtaSig[1] = varValuesSigSub[2][n]
        etaWidthSig[0] = varValuesSigLead[3][n]
        etaWidthSig[1] = varValuesSigSub[3][n]
        phiWidthSig[0] = varValuesSigLead[4][n]
        phiWidthSig[1] = varValuesSigSub[4][n]
        s4Sig[0] = varValuesSigLead[5][n]
        s4Sig[1] = varValuesSigSub[5][n]
        trkIsoPhoSig[0] = varValuesSigLead[6][n]
        trkIsoPhoSig[1] = varValuesSigSub[6][n]
        etaSig[0] = varValuesSigLead[7][n]
        etaSig[1] = varValuesSigSub[7][n]
        hOvrESig[0] = varValuesSigLead[8][n]
        hOvrESig[1] = varValuesSigSub[8][n]
        ecalPFIsoSig[0] = varValuesSigLead[9][n]
        ecalPFIsoSig[1] = varValuesSigSub[9][n]
        etSig[0] = varValuesSigLead[10][n]
        etSig[1] = varValuesSigSub[10][n]
        energySig[0] = varValuesSigLead[11][n]
        energySig[1] = varValuesSigSub[11][n]
        dRSig[0] = varValuesSigLead[12][n]
        dRSig[1] = varValuesSigSub[12][n]
        massSig[0] = varValuesSigLead[13][n]
        massSig[1] = varValuesSigSub[13][n]
        nEgsSig[0] = int(varValuesSigLead[14][n])
        nEgsSig[1] = int(varValuesSigLead[14][n])
        passFailStdSig[0] = int(varValuesSigLead[15][n])
        passFailL1SingleSig[0] = int(varValuesSigLead[16][n])
        passFailL1DoubleSig[0] = int(varValuesSigLead[17][n])
        triggerBitsSig[0] = int(varValuesSigLead[18][n])
        triggerBitsSig[1] = int(varValuesSigLead[18][n])
        tSig.Fill()
        
    xgbScoreBkg = array('f',[0.0,0.0])
    rawEnergyBkg = array('f',[0.0,0.0])
    r9HLTBkg = array('f',[0.0,0.0])
    sigmaIEtaIEtaBkg = array('f',[0.0,0.0])
    etaWidthBkg = array('f',[0.0,0.0])
    phiWidthBkg = array('f',[0.0,0.0])
    s4Bkg = array('f',[0.0,0.0])
    trkIsoPhoBkg = array('f',[0.0,0.0])
    etaBkg = array('f',[0.0,0.0])
    hOvrEBkg = array('f',[0.0,0.0])
    ecalPFIsoBkg = array('f',[0.0,0.0])
    etBkg = array('f',[0.0,0.0])
    energyBkg = array('f',[0.0,0.0])
    dRBkg = array('f',[0.0,0.0])
    massBkg = array('f',[0.0,0.0])
    nEgsBkg = array('f',[0.0,0.0])
    triggerBitsBkg = array('f',[0.0,0.0])
    passFailStdBkg = array('i',[0])
    passFailL1SingleBkg = array('i',[0])
    passFailL1DoubleBkg = array('i',[0])

    tBkg = TTree("bkgTree", "bkgTree")
    tBkg.Branch( 'xgbScore', xgbScoreBkg, 'xgbScore[2]/F' )
    tBkg.Branch( 'rawEnergy', rawEnergyBkg, 'rawEnergy[2]/F' )
    tBkg.Branch( 'r9HLT', r9HLTBkg, 'r9HLT[2]/F' )
    tBkg.Branch( 'sigmaIEtaIEta', sigmaIEtaIEtaBkg, 'sigmaIEtaIEta[2]/F' )
    tBkg.Branch( 'etaWidth', etaWidthBkg, 'etaWidth[2]/F' )
    tBkg.Branch( 'phiWidth', phiWidthBkg, 'phiWidth[2]/F' )
    tBkg.Branch( 's4', s4Bkg, 's4[2]/F' )
    tBkg.Branch( 'trkIsoPho', trkIsoPhoBkg, 'trkIsoPho[2]/F' )
    tBkg.Branch( 'eta', etaBkg, 'eta[2]/F' )
    tBkg.Branch( 'hOvrE', hOvrEBkg, 'hOvrE[2]/F' )
    tBkg.Branch( 'ecalPFIso', ecalPFIsoBkg, 'ecalPFIso[2]/F' )
    tBkg.Branch( 'et', etBkg, 'et[2]/F' )
    tBkg.Branch( 'energy', energyBkg, 'energy[2]/F' )
    tBkg.Branch( 'dR', dRBkg, 'dR[2]/F' )
    tBkg.Branch( 'mass', massBkg, 'mass[2]/F' )
    tBkg.Branch( 'nEgs', nEgsBkg, 'nEgs[2]/F' )
    tBkg.Branch( 'triggerBits', triggerBitsBkg, 'triggerBits[2]/F' )
    tBkg.Branch( 'passFailStd', passFailStdBkg, 'passFailStd/I' )
    tBkg.Branch( 'passFailL1Single', passFailL1SingleBkg, 'passFailL1Single/I' )
    tBkg.Branch( 'passFailL1Double', passFailL1DoubleBkg, 'passFailL1Double/I' )

    bkgLen = len(yPredBkgLeadTotal)
    
    for n in range(0,bkgLen):
        xgbScoreBkg[0] = yPredBkgLeadTotal[n][0]
        xgbScoreBkg[1] = yPredBkgSubTotal[n][0]
        rawEnergyBkg[0] = varValuesBkgLead[0][n]
        rawEnergyBkg[1] = varValuesBkgSub[0][n]
        r9HLTBkg[0] = varValuesBkgLead[1][n]
        r9HLTBkg[1] = varValuesBkgSub[1][n]
        sigmaIEtaIEtaBkg[0] = varValuesBkgLead[2][n]
        sigmaIEtaIEtaBkg[1] = varValuesBkgSub[2][n]
        etaWidthBkg[0] = varValuesBkgLead[3][n]
        etaWidthBkg[1] = varValuesBkgSub[3][n]
        phiWidthBkg[0] = varValuesBkgLead[4][n]
        phiWidthBkg[1] = varValuesBkgSub[4][n]
        s4Bkg[0] = varValuesBkgLead[5][n]
        s4Bkg[1] = varValuesBkgSub[5][n]
        trkIsoPhoBkg[0] = varValuesBkgLead[6][n]
        trkIsoPhoBkg[1] = varValuesBkgSub[6][n]
        etaBkg[0] = varValuesBkgLead[7][n]
        etaBkg[1] = varValuesBkgSub[7][n]
        hOvrEBkg[0] = varValuesBkgLead[8][n]
        hOvrEBkg[1] = varValuesBkgSub[8][n]
        ecalPFIsoBkg[0] = varValuesBkgLead[9][n]
        ecalPFIsoBkg[1] = varValuesBkgSub[9][n]
        etBkg[0] = varValuesBkgLead[10][n]
        etBkg[1] = varValuesBkgSub[10][n]
        energyBkg[0] = varValuesBkgLead[11][n]
        energyBkg[1] = varValuesBkgSub[11][n]
        dRBkg[0] = varValuesBkgLead[12][n]
        dRBkg[1] = varValuesBkgSub[12][n]
        massBkg[0] = varValuesBkgLead[13][n]
        massBkg[1] = varValuesBkgSub[13][n]
        nEgsBkg[0] = int(varValuesBkgLead[14][n])
        nEgsBkg[1] = int(varValuesBkgLead[14][n])
        passFailStdBkg[0] = int(varValuesBkgLead[15][n])
        passFailL1SingleBkg[0] = int(varValuesBkgLead[16][n])
        passFailL1DoubleBkg[0] = int(varValuesBkgLead[17][n])
        triggerBitsBkg[0] = int(varValuesBkgLead[18][n])
        triggerBitsBkg[1] = int(varValuesBkgLead[18][n])
        tBkg.Fill()
        

    outFileRoot.Write()
    outFileRoot.Close()
