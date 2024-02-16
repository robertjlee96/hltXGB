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
import allPhotonsUtils
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

inFileName ="../NTuples/outNTuples/GGH13andData_RelaxedV24_TrainingCuts_0502_DiphotonValidation.root"
modelNamesB = []
mNamesGen =  []
mNamesGen.append("0516/M9LR10_GGH13andData_M90_0516")

modelNamesB.append("../Training/barrelOut/" + mNamesGen[0] + "_Barrel.model")
#modelNamesB.append("../Training/barrelOut/Model_NewGGH_MD13LR07_M90_Barrel_0504.Model")
#modelNamesB.append("../Training/barrelOut/Model_NewGGH_MD17LR05_M90_Barrel_0504.Model")

modelNamesE = []
modelNamesE.append("../Training/endcapOut/" + mNamesGen[0] + "_Endcap.model")
#modelNamesE.append("../Training/endcapOut/Model_NewGGH_MD9LR10_M90_Endcap_0504.model")
#modelNamesE.append("../Training/endcapOut/Model_NewGGH_MD9LR10_M90_Endcap_0504.model")

fNames = []
fNames.append("validationNTuples/0516/newGGHandData_M90Model_MD9LR10")
#fNames.append("validationNTuples/0504/newGGHandData_MD13LR07")
#fNames.append("validationNTuples/0504/newGGHandData_MD17LR05")

fin = uproot.open(inFileName)
prompt = fin['sigTree']
fake = fin['bkgTree']

## for barrel
geometry_selection = lambda tree: np.abs(tree.array('eta')) < 2.5
ptCut = lambda tree: tree.array('et') > 0.0
massCut = lambda tree: tree.array('mass') > 60

#diphotonUtils.load_file(fin,geometry_selection,ptCut,massCut)

inputValuesSigLead, inputValuesSigSub, inputValuesBkgLead, inputValuesBkgSub, targetValuesSigLead, targetValuesSigSub, targetValuesBkgLead,targetValuesBkgSub, origWeightsSigLead, origWeightsSigSub,origWeightsBkgLead,origWeightsBkgSub, varValuesSigLead, varValuesSigSub, varValuesBkgLead, varValuesBkgSub, inputVarNames, varNames = diphotonUtils.load_file(fin,geometry_selection,ptCut,massCut)


print len(inputValuesSigLead)
print len(inputValuesSigLead[0])


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
    yPredSigLeadB = np.stack((modelB.predict_proba(barrelSigLead)[:,1],barrelSigLeadVars[:,18]),axis=1)
    yPredSigSubB = np.stack((modelB.predict_proba(barrelSigSub)[:,1],barrelSigSubVars[:,18]),axis=1)
    yPredBkgLeadB = np.stack((modelB.predict_proba(barrelBkgLead)[:,1],barrelBkgLeadVars[:,18]),axis=1)
    yPredBkgSubB = np.stack((modelB.predict_proba(barrelBkgSub)[:,1],barrelBkgSubVars[:,18]),axis=1)

    yPredSigLeadE = np.stack((modelE.predict_proba(endcapSigLead)[:,1],endcapSigLeadVars[:,18]),axis=1)
    yPredSigSubE = np.stack((modelE.predict_proba(endcapSigSub)[:,1],endcapSigSubVars[:,18]),axis=1)
    yPredBkgLeadE = np.stack((modelE.predict_proba(endcapBkgLead)[:,1],endcapBkgLeadVars[:,18]),axis=1)
    yPredBkgSubE = np.stack((modelE.predict_proba(endcapBkgSub)[:,1],endcapBkgSubVars[:,18]),axis=1)

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
        massSig[0] = varValuesSigLead[12][n]
        massSig[1] = varValuesSigSub[12][n]
        nEgsSig[0] = int(varValuesSigLead[13][n])
        nEgsSig[1] = int(varValuesSigLead[13][n])
        passFailStdSig[0] = int(varValuesSigLead[14][n])
        passFailL1SingleSig[0] = int(varValuesSigLead[15][n])
        passFailL1DoubleSig[0] = int(varValuesSigLead[16][n])
        triggerBitsSig[0] = int(varValuesSigLead[17][n])
        triggerBitsSig[1] = int(varValuesSigLead[17][n])
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
        massBkg[0] = varValuesBkgLead[12][n]
        massBkg[1] = varValuesBkgSub[12][n]
        nEgsBkg[0] = int(varValuesBkgLead[13][n])
        nEgsBkg[1] = int(varValuesBkgLead[13][n])
        passFailStdBkg[0] = int(varValuesBkgLead[14][n])
        passFailL1SingleBkg[0] = int(varValuesBkgLead[15][n])
        passFailL1DoubleBkg[0] = int(varValuesBkgLead[16][n])
        triggerBitsBkg[0] = int(varValuesBkgLead[17][n])
        triggerBitsBkg[1] = int(varValuesBkgLead[17][n])
        tBkg.Fill()
        

    outFileRoot.Write()
    outFileRoot.Close()
    
    
#
#    xgbScoreSigLead = array('f',[0.0])
#    weightSigLead = array('f',[0.0])
#    varValsSigLead = array('f',23*[0.0])
#
#    xgbScoreSigSub = array('f',[0.0])
#    weightSigSub = array('f',[0.0])
#    varValsSigSub = array('f',23*[0.0])
#
#    xgbScorePrompt = array('f',[0.0])
#    weightPrompt = array('f',[0.0])
#    varValsPrompt = array('f',23*[0.0])
#
#    xgbScoreFake = array('f',[0.0])
#    weightFake = array('f',[0.0])
#    varValsFake = array('f',23*[0.0])
#
#    tSigLead = TTree("sigLead", "Signal Lead Photons")
#    tSigLead.Branch( 'xgbScore', xgbScoreSigLead, 'xgbScore/F' )
#    tSigLead.Branch( 'weight', weightSigLead, 'weight/F' )
#    tSigLead.Branch( 'varVals', varValsSigLead, 'varVals[23]/F' )
#
#    tSigSub = TTree("sigSub", "Signal SubLead Photons")
#    tSigSub.Branch( 'xgbScore', xgbScoreSigSub, 'xgbScore/F' )
#    tSigSub.Branch( 'weight', weightSigSub, 'weight/F' )
#    tSigSub.Branch( 'varVals', varValsSigSub, 'varVals[23]/F' )
#
#    tPrompt = TTree("promptPhotons", "GJet Prompt")
#    tPrompt.Branch( 'xgbScore', xgbScorePrompt, 'xgbScore/F' )
#    tPrompt.Branch( 'weight', weightPrompt, 'weight/F' )
#    tPrompt.Branch( 'varVals', varValsPrompt, 'varVals[23]/F' )
#
#    tFake = TTree("fakePhotons", "GJet Fake")
#    tFake.Branch( 'xgbScore', xgbScoreFake, 'xgbScore/F' )
#    tFake.Branch( 'weight', weightFake, 'weight/F' )
#    tFake.Branch( 'varVals', varValsFake, 'varVals[23]/F' )
#
#    xgbScoreSigLeadArray = []
#    weightSigLeadArray = []
#    xgbScoreSigSubArray = []
#    weightSigSubArray = []
#
#    xgbScorePromptArray = []
#    weightPromptArray = []
#
#    xgbScoreFakeArray = []
#    weightFakeArray = []
#
##    sigSubLen = len(yPredSigLead)
##    for n in range (0,sigSubLen):
##        if targetValuesSigSub[n] == 1:
##            xgbScoreSigArray.append(y_pred[n])
##            weightSigArray.append(orig_weights[n])
##        if target_values[n] == 0:
##            xgbScoreBkgArray.append(y_pred[n])
##            weightBkgArray.append(orig_weights[n])
##
#    sigLeadLen = len(yPredSigLeadTotal)
#    for n in range(0,sigLeadLen):
#        xgbScoreSigLead[0] = yPredSigLeadTotal[n][0]
#        weightSigLead[0] = origWeightsSigLead[0][n]
#        for i in range(0,23):
#            varValsSigLead[i] = varValuesSigLead[n][i]
#        tSigLead.Fill()
#
#    sigSubLen = len(yPredSigSubTotal)
#    for n in range(0,sigSubLen):
#        xgbScoreSigSub[0] = yPredSigSubTotal[n][0]
#        weightSigSub[0] = origWeightsSigSub[0][n]
#        for i in range(0,23):
#            varValsSigSub[i] = varValuesSigSub[n][i]
#        tSigSub.Fill()
#
#    promptLen = len(yPredPromptTotal)
#    for n in range(0,promptLen):
#        xgbScorePrompt[0] = yPredPromptTotal[n][0]
#        weightPrompt[0] = origWeightsPrompt[0][n]
#        for i in range(0,23):
#            varValsPrompt[i] = varValuesPrompt[n][i]
#        tPrompt.Fill()
#
#    fakeLen = len(yPredFakeTotal)
#
#    for n in range(0,fakeLen):
#        xgbScoreFake[0] = yPredFakeTotal[n][0]
#        weightFake[0] = origWeightsFake[n]
#        for i in range(0,23):
#            varValsFake[i] = varValuesFake[n][i]
#        tFake.Fill()
##
##    for n in range(0,bkgLen):
##        xgbScoreBkg[0] = xgbScoreBkgArray[n]
##        weightBkg[0] = weightBkgArray[n]
##        for i in range(0,22):
##            varValsBkg[i] = varValsBkgIn[n][i]
##        tFake.Fill()
##
#    outFileRoot.Write()
#    outFileRoot.Close()
##













    #    xgBoostParams = []
#    aucTest = []
#    aucTrain = []
#    lossTest = []
#    lossTrain = []
#    mvaValsTest = []
#    mvaValsTrain = []
#    bestVals = []
#    fprVals = []
#    tprVals = []
#    timeVals = []

#    boostParamsTmp = []
#    boostParamsTmp.append(model.get_xgb_params()["max_depth"])
#    boostParamsTmp.append(model.get_xgb_params()["learning_rate"])
#    boostParamsTmp.append(model.get_xgb_params()["gamma"])
#    boostParamsTmp.append(model.get_xgb_params()["min_child_weight"])
#    boostParamsTmp.append(model.get_xgb_params()["max_delta_step"])
#    boostParamsTmp.append(model.get_xgb_params()["colsample_bytree"])
#    boostParamsTmp.append(model.get_xgb_params()["subsample"])
#    boostParamsTmp.append(model.get_xgb_params()["reg_alpha"])
#    boostParamsTmp.append(model.get_xgb_params()["reg_lambda"])
#    xgBoostParams.append(boostParamsTmp)
    
#     fpr, tpr, thresholds = roc_curve(target_values,y_pred,sample_weight=train_weights)
#    fprVals.append(fpr)
#    tprVals.append(tpr)

    
#    results = model.evals_result() ##NO TEST REWEIGHTING
    
#    aucTest.append(auc_val_test)
#    aucTrain.append(results['validation_0']['auc'])
#    lossTest.append(results['validation_1']['logloss'])
#    lossTrain.append(results['validation_0']['logloss'])
#
#    bestRound = model.best_iteration
#
#    bestValsTmp = [bestRound,results['validation_0']['auc'][bestRound],results['validation_1']['logloss'][bestRound],results['validation_0']['logloss'][bestRound]]
#
#    bestVals.append(bestValsTmp)
#
#
#    #MoreVarsTmvaUtils.variable_importance(model,input_variables,model_name)
#
#    del results
#    del model
#
#    np.savez(outFile, xgBoostParams = xgBoostParams, aucTest = aucTest, aucTrain = aucTrain, lossTest = lossTest, lossTrain = lossTrain, mvaValsTest = mvaValsTest, mvaValsTrain = mvaValsTrain, bestVals = bestVals, fprVals = fprVals, tprVals = tprVals, timeVals = timeVals)
