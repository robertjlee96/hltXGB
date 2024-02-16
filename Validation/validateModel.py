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
import validationUtils
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

inFileName ="../NTuples/outNTuples/GGH13andData_RelaxedV24_TrainingCuts_SinglePhotonValidation_0517.root"
modelNamesB = []
modelNamesB.append("../Training/barrelOut/0516/M9LR10_GGH13andData_M60_0516_Barrel.model")
#modelNamesB.append("../Training/barrelOut/Model_MD13LR07_M90_Barrel_0418.Model")
#modelNamesB.append("../Training/barrelOut/Model_MD17LR05_M90_Barrel_0418.Model")


modelNamesE = []
modelNamesE.append("../Training/endcapOut/0516/M9LR10_GGH13andData_M60_0516_Endcap.model")
#modelNamesE.append("../Training/endcapOut/Model_MD13LR07_M90_Endcap_0418.model")
#modelNamesE.append("../Training/endcapOut/Model_MD17LR05_M90_Endcap_0418.model")

fNames = []
fNames.append("validationNTuples/0517/SinglePhoton_GGHandData_MD9LR10_M60")
#fNames.append("validationNTuples/0517/GGHandData_MD13LR07")
#fNames.append("validationNTuples/0419/GGHandData_MD17LR05")


fin = uproot.open(inFileName)
prompt = fin['sigTree']
fake = fin['bkgTree']

## for barrel
geometry_selection = lambda tree: np.abs(tree.array('eta')) < 2.5
ptCut = lambda tree: tree.array('et') > 0.0
massCut = lambda tree: tree.array('mass') > 60.0

inputValuesSig, inputValuesBkg, targetValuesSig, targetValuesBkg, origWeightsSig, origWeightsBkg, varValuesSig, varValuesBkg, inputVarNames, varNames = validationUtils.load_file(fin,geometry_selection,ptCut,massCut)

for i in range(0,len(fNames)):
    fNameRoot = fNames[i] + ".root"
    outFileRoot = TFile(fNameRoot,'recreate')
   
#    FOR MODELS WITH H/E VARIABLES
    barrelSig = (((inputValuesSig[np.abs(varValuesSig[:,7]) < 1.4442]).T)[0:10]).T
    barrelSigVars = varValuesSig[np.abs(varValuesSig[:,7]) < 1.4442]
    
    barrelBkg = (((inputValuesBkg[np.abs(varValuesBkg[:,7]) < 1.4442]).T)[0:10]).T
    barrelBkgVars = varValuesBkg[np.abs(varValuesBkg[:,7]) < 1.4442]

    endcapSig = (((inputValuesSig[np.abs(varValuesSig[:,7]) > 1.556]).T)[0:10]).T
    endcapSigVars = varValuesSig[np.abs(varValuesSig[:,7]) > 1.556]
    
    endcapBkg = (((inputValuesBkg[np.abs(varValuesBkg[:,7]) > 1.556]).T)[0:10]).T
    endcapBkgVars = varValuesBkg[np.abs(varValuesBkg[:,7]) > 1.556]

    startTime = time.time()
   
    modelFileB = open(modelNamesB[i],'rb')
    modelB = pickle.load(modelFileB)
    yPredSigB = modelB.predict_proba(barrelSig)[:,1]
    yPredBkgB = modelB.predict_proba(barrelBkg)[:,1]
    
    
    modelFileE = open(modelNamesE[i],'rb')
    modelE = pickle.load(modelFileE)
    yPredSigE = modelE.predict_proba(endcapSig)[:,1]
    yPredBkgE = modelE.predict_proba(endcapBkg)[:,1]

    #yPredSigTotal = yPredSigB + yPredSigE
    

    endTime = time.time()
    timeSpent = endTime - startTime
    timeSpentMin = timeSpent/60.0
    print 'Run Time(min) =  {0:.4f}'.format(timeSpentMin)

    xgbScoreSig = array('f',[0.0])
    rawEnergySig = array('f',[0.0])
    r9HLTSig = array('f',[0.0])
    sigmaIEtaIEtaSig = array('f',[0.0])
    etaWidthSig = array('f',[0.0])
    phiWidthSig = array('f',[0.0])
    s4Sig = array('f',[0.0])
    trkIsoPhoSig = array('f',[0.0])
    etaSig = array('f',[0.0])
    hOvrESig = array('f',[0.0])
    ecalPFIsoSig = array('f',[0.0])
    nEgsSig = array('f',[0.0])
    etSig = array('f',[0.0])
    energySig = array('f',[0.0])
    massSig = array('f',[0.0])
    passFailStdSig = array('i',[0])
    passFailL1SingleSig = array('i',[0])
    passFailL1DoubleSig = array('i',[0])

    tSig = TTree("sigTree", "sigTree")
    tSig.Branch( 'xgbScore', xgbScoreSig, 'xgbScore/F' )
    tSig.Branch( 'rawEnergy', rawEnergySig, 'rawEnergy/F' )
    tSig.Branch( 'r9HLT', r9HLTSig, 'r9HLT/F' )
    tSig.Branch( 'sigmaIEtaIEta', sigmaIEtaIEtaSig, 'sigmaIEtaIEta/F' )
    tSig.Branch( 'etaWidth', etaWidthSig, 'etaWidth/F' )
    tSig.Branch( 'phiWidth', phiWidthSig, 'phiWidth/F' )
    tSig.Branch( 's4', s4Sig, 's4/F' )
    tSig.Branch( 'trkIsoPho', trkIsoPhoSig, 'trkIsoPho/F' )
    tSig.Branch( 'eta', etaSig, 'eta/F' )
    tSig.Branch( 'hOvrE', hOvrESig, 'hOvrE/F' )
    tSig.Branch( 'ecalPFIso', ecalPFIsoSig, 'ecalPFIso/F' )
    tSig.Branch( 'nEgs', nEgsSig, 'nEgs/F' )
    tSig.Branch( 'et', etSig, 'et/F' )
    tSig.Branch( 'energy', energySig, 'energy/F' )
    tSig.Branch( 'mass', massSig, 'mass/F' )
    tSig.Branch( 'passFailStd', passFailStdSig, 'passFailStd/I' )
    tSig.Branch( 'passFailL1Single', passFailL1SingleSig, 'passFailL1Single/I' )
    tSig.Branch( 'passFailL1Double', passFailL1DoubleSig, 'passFailL1Double/I' )
    
    sigLenB = len(yPredSigB)
    for n in range(0,sigLenB):
        xgbScoreSig[0] = yPredSigB[n]
        rawEnergySig[0] = barrelSigVars[n][0]
        r9HLTSig[0] = barrelSigVars[n][1]
        sigmaIEtaIEtaSig[0] = barrelSigVars[n][2]
        etaWidthSig[0] = barrelSigVars[n][3]
        phiWidthSig[0] = barrelSigVars[n][4]
        s4Sig[0] = barrelSigVars[n][5]
        trkIsoPhoSig[0] = barrelSigVars[n][6]
        etaSig[0] = barrelSigVars[n][7]
        hOvrESig[0] = barrelSigVars[n][8]
        ecalPFIsoSig[0] = barrelSigVars[n][9]
        nEgsSig[0] = barrelSigVars[n][10]
        etSig[0] = barrelSigVars[n][11]
        energySig[0] = barrelSigVars[n][12]
        massSig[0] = barrelSigVars[n][13]
        passFailStdSig[0] = int(barrelSigVars[n][14])
        passFailL1SingleSig[0] = int(barrelSigVars[n][15])
        passFailL1DoubleSig[0] = int(barrelSigVars[n][16])
        tSig.Fill()
        
    sigLenE = len(yPredSigE)
    for n in range(0,sigLenE):
        xgbScoreSig[0] = yPredSigE[n]
        rawEnergySig[0] = endcapSigVars[n][0]
        r9HLTSig[0] = endcapSigVars[n][1]
        sigmaIEtaIEtaSig[0] = endcapSigVars[n][2]
        etaWidthSig[0] = endcapSigVars[n][3]
        phiWidthSig[0] = endcapSigVars[n][4]
        s4Sig[0] = endcapSigVars[n][5]
        trkIsoPhoSig[0] = endcapSigVars[n][6]
        etaSig[0] = endcapSigVars[n][7]
        hOvrESig[0] = endcapSigVars[n][8]
        ecalPFIsoSig[0] = endcapSigVars[n][9]
        nEgsSig[0] = endcapSigVars[n][10]
        etSig[0] = endcapSigVars[n][11]
        energySig[0] = endcapSigVars[n][12]
        massSig[0] = endcapSigVars[n][13]
        passFailStdSig[0] = int(endcapSigVars[n][14])
        passFailL1SingleSig[0] = int(endcapSigVars[n][15])
        passFailL1DoubleSig[0] = int(endcapSigVars[n][16])
        tSig.Fill()
        
        
    xgbScoreBkg = array('f',[0.0])
    rawEnergyBkg = array('f',[0.0])
    r9HLTBkg = array('f',[0.0])
    sigmaIEtaIEtaBkg = array('f',[0.0])
    etaWidthBkg = array('f',[0.0])
    phiWidthBkg = array('f',[0.0])
    s4Bkg = array('f',[0.0])
    trkIsoPhoBkg = array('f',[0.0])
    etaBkg = array('f',[0.0])
    hOvrEBkg = array('f',[0.0])
    ecalPFIsoBkg = array('f',[0.0])
    nEgsBkg = array('f',[0.0])
    etBkg = array('f',[0.0])
    energyBkg = array('f',[0.0])
    massBkg = array('f',[0.0])
    passFailStdBkg = array('i',[0])
    passFailL1SingleBkg = array('i',[0])
    passFailL1DoubleBkg = array('i',[0])

    tBkg = TTree("bkgTree", "bkgTree")
    tBkg.Branch( 'xgbScore', xgbScoreBkg, 'xgbScore/F' )
    tBkg.Branch( 'rawEnergy', rawEnergyBkg, 'rawEnergy/F' )
    tBkg.Branch( 'r9HLT', r9HLTBkg, 'r9HLT/F' )
    tBkg.Branch( 'sigmaIEtaIEta', sigmaIEtaIEtaBkg, 'sigmaIEtaIEta/F' )
    tBkg.Branch( 'etaWidth', etaWidthBkg, 'etaWidth/F' )
    tBkg.Branch( 'phiWidth', phiWidthBkg, 'phiWidth/F' )
    tBkg.Branch( 's4', s4Bkg, 's4/F' )
    tBkg.Branch( 'trkIsoPho', trkIsoPhoBkg, 'trkIsoPho/F' )
    tBkg.Branch( 'eta', etaBkg, 'eta/F' )
    tBkg.Branch( 'hOvrE', hOvrEBkg, 'hOvrE/F' )
    tBkg.Branch( 'ecalPFIso', ecalPFIsoBkg, 'ecalPFIso/F' )
    tBkg.Branch( 'nEgs', nEgsBkg, 'nEgs/F' )
    tBkg.Branch( 'et', etBkg, 'et/F' )
    tBkg.Branch( 'energy', energyBkg, 'energy/F' )
    tBkg.Branch( 'mass', massBkg, 'mass/F' )
    tBkg.Branch( 'passFailStd', passFailStdBkg, 'passFailStd/I' )
    tBkg.Branch( 'passFailL1Single', passFailL1SingleBkg, 'passFailL1Single/I' )
    tBkg.Branch( 'passFailL1Double', passFailL1DoubleBkg, 'passFailL1Double/I' )
    
    bkgLenB = len(yPredBkgB)
    for n in range(0,bkgLenB):
        xgbScoreBkg[0] = yPredBkgB[n]
        rawEnergyBkg[0] = barrelBkgVars[n][0]
        r9HLTBkg[0] = barrelBkgVars[n][1]
        sigmaIEtaIEtaBkg[0] = barrelBkgVars[n][2]
        etaWidthBkg[0] = barrelBkgVars[n][3]
        phiWidthBkg[0] = barrelBkgVars[n][4]
        s4Bkg[0] = barrelBkgVars[n][5]
        trkIsoPhoBkg[0] = barrelBkgVars[n][6]
        etaBkg[0] = barrelBkgVars[n][7]
        hOvrEBkg[0] = barrelBkgVars[n][8]
        ecalPFIsoBkg[0] = barrelBkgVars[n][9]
        nEgsBkg[0] = barrelBkgVars[n][10]
        etBkg[0] = barrelBkgVars[n][11]
        energyBkg[0] = barrelBkgVars[n][12]
        massBkg[0] = barrelBkgVars[n][13]
        passFailStdBkg[0] = int(barrelBkgVars[n][14])
        passFailL1SingleBkg[0] = int(barrelBkgVars[n][15])
        passFailL1DoubleBkg[0] = int(barrelBkgVars[n][16])
        tBkg.Fill()
        
    bkgLenE = len(yPredBkgE)
    for n in range(0,bkgLenE):
        xgbScoreBkg[0] = yPredBkgE[n]
        rawEnergyBkg[0] = endcapBkgVars[n][0]
        r9HLTBkg[0] = endcapBkgVars[n][1]
        sigmaIEtaIEtaBkg[0] = endcapBkgVars[n][2]
        etaWidthBkg[0] = endcapBkgVars[n][3]
        phiWidthBkg[0] = endcapBkgVars[n][4]
        s4Bkg[0] = endcapBkgVars[n][5]
        trkIsoPhoBkg[0] = endcapBkgVars[n][6]
        etaBkg[0] = endcapBkgVars[n][7]
        hOvrEBkg[0] = endcapBkgVars[n][8]
        ecalPFIsoBkg[0] = endcapBkgVars[n][9]
        nEgsBkg[0] = endcapBkgVars[n][10]
        etBkg[0] = endcapBkgVars[n][11]
        energyBkg[0] = endcapBkgVars[n][12]
        massBkg[0] = endcapBkgVars[n][13]
        passFailStdBkg[0] = int(endcapBkgVars[n][14])
        passFailL1SingleBkg[0] = int(endcapBkgVars[n][15])
        passFailL1DoubleBkg[0] = int(endcapBkgVars[n][16])
        tBkg.Fill()
    
    outFileRoot.Write()
    outFileRoot.Close()
