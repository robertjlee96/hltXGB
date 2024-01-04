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
import gc
from operator import itemgetter

# for sklearn, see
np.random.seed(1337)

inFileName = "../NTuples/outNTuples/GGH13andData_RelaxedV24_TrainingCuts_DiphotonTrain_0510.root"

outNamesGen = []
outNamesGen.append("Model_M9LR10_M90_0510")

modelNamesB = []
modelNamesB.append("barrelOut/" + outNamesGen[0] + "Barrel.model")
#modelNamesB.append("../Training/barrelOut/Model_MD13LR07_M90_Barrel_0418.Model")
#modelNamesB.append("../Training/barrelOut/Model_MD17LR05_M90_Barrel_0418.Model")
#
modelNamesE = []
modelNamesE.append("endcapOut/" + outNamesGen[0] + "Endcap.model")
#modelNamesE.append("../Training/endcapOut/Model_MD13LR07_M90_Endcap_0418.model")
#modelNamesE.append("../Training/endcapOut/Model_MD17LR05_M90_Endcap_0418.model")
#
#fNames = []
#fNames.append("validationNTuples/0425/GGHandData_MD9LR10")
#fNames.append("validationNTuples/0425/GGHandData_MD13LR07")
#fNames.append("validationNTuples/0425/GGHandData_MD17LR05")

fin = uproot.open(inFileName)
prompt = fin['sigTree']
fake = fin['bkgTree']

## for barrel
geometry_selection = lambda tree: np.abs(tree.array('eta')) < 2.5
ptCut = lambda tree: tree.array('et') > 14.25
massCut = lambda tree: tree.array('mass') > 90.0

#diphotonUtils.load_file(fin,geometry_selection,ptCut,massCut)

inputValuesSig, inputValuesBkg, targetValuesSig, targetValuesBkg, origWeightsSig, origWeightsBkg, varValuesSig, varValuesBkg, inputVarNames, varNames = diphotonUtils.load_file(fin,geometry_selection,ptCut,massCut)

fNameRoot = outNamesGen[0] + ".root"
outFileRoot = TFile(fNameRoot,'recreate')


barrelSig = inputValuesSig.T[np.abs(varValuesSig[7,:]) < 1.5]
barrelSigVars = varValuesSig.T[np.abs(varValuesSig[7,:]) < 1.5]
barrelSigTarget = np.ones(len(barrelSig))
barrelSigWeights = np.ones(len(barrelSig))
#barrelSigTarget = targetValuesSig[np.abs(varValuesSig[7,:]) < 1.5]
#barrelSigWeights = origWeightsSig[np.abs(varValuesSig[7,:]) < 1.5]

barrelBkg = inputValuesBkg.T[np.abs(varValuesBkg[7,:]) < 1.5]
barrelBkgVars = varValuesBkg.T[np.abs(varValuesBkg[7,:]) < 1.5]
barrelBkgTarget = np.zeros(len(barrelBkg))
barrelBkgWeights = np.ones(len(barrelBkg))
#barrelBkgTarget = targetValuesBkg[np.abs(varValuesBkg[7,:]) < 1.5]
#barrelBkgWeights = origWeightsBkg[np.abs(varValuesBkg[7,:]) < 1.5]

endcapSig = inputValuesSig.T[np.abs(varValuesSig[7,:]) > 1.5]
endcapSigVars = varValuesSig.T[np.abs(varValuesSig[7,:]) > 1.5]
endcapSigTarget = np.ones(len(endcapSig))
endcapSigWeights = np.ones(len(endcapSig))
#endcapSigTarget = targetValuesSig[np.abs(varValuesSig[7,:]) > 1.5]
#endcapSigWeights = origWeightsSig[np.abs(varValuesSig[7,:]) > 1.5]

endcapBkg = inputValuesBkg.T[np.abs(varValuesBkg[7,:]) > 1.5]
endcapBkgVars = varValuesBkg.T[np.abs(varValuesBkg[7,:]) > 1.5]
endcapBkgTarget = np.zeros(len(endcapBkg))
endcapBkgWeights = np.ones(len(endcapBkg))
#endcapBkgTarget = targetValuesBkg[np.abs(varValuesBkg[7,:]) > 1.5]
#endcapBkgWeights = origWeightsBkg[np.abs(varValuesBkg[7,:]) > 1.5]

#REMOVED .T from above arrays (varValuesBkg.T[np.abs(varValuesBkg[7,:]) > 1.5]), instead cat the arrays THEN tranpose
barrelInputTotal = np.concatenate((barrelSig,barrelBkg),axis = 0)
barrelVarsTotal = np.concatenate((barrelSigVars,barrelBkgVars),axis = 0)
barrelTargetTotal = np.concatenate((barrelSigTarget,barrelBkgTarget),axis = 0)
barrelWeightsTotal = np.concatenate((barrelSigWeights,barrelBkgWeights),axis = 0)

endcapInputTotal = np.concatenate((endcapSig,endcapBkg),axis = 0)
endcapVarsTotal = np.concatenate((endcapSigVars,endcapBkgVars),axis = 0)
endcapTargetTotal = np.concatenate((endcapSigTarget,endcapBkgTarget),axis = 0)
endcapWeightsTotal = np.concatenate((endcapSigWeights,endcapBkgWeights),axis = 0)

MD = 13
LR = 0.1
nEst = 1500

print

xTrainB, xTestB, wTrainB, dummy, yTrainB, yTestB, dummy, wTestB, varValuesTrainB, varValuesTestB = train_test_split(barrelInputTotal,barrelWeightsTotal,barrelTargetTotal,barrelWeightsTotal,barrelVarsTotal,test_size=0.25)

xTrainE, xTestE, wTrainE, dummy, yTrainE, yTestE, dummy, wTestE, varValuesTrainE, varValuesTestE = train_test_split(endcapInputTotal,endcapWeightsTotal,endcapTargetTotal,endcapWeightsTotal,endcapVarsTotal,test_size=0.25)

startTime = time.time()

print len(xTrainB)
print len(xTrainE)

sys.exit()
modelB = XGBClassifier(max_depth = MD, learning_rate = LR, subsample = 0.9, max_delta_step = 0.0, min_child_weight = 0.0, reg_alpha = 0.0, reg_lambda = 4.0, colsample_bytree = 0.65, verbosity = 1, n_estimators=nEst, nthread = 16, feature_names = inputVarNames)
evalSetB = [(xTrainB, yTrainB), (xTestB, yTestB)]
modelB.fit(xTrainB, yTrainB, sample_weight = wTrainB, eval_metric=["auc", "logloss"], eval_set=evalSetB, early_stopping_rounds=50, verbose=True)


modelE = XGBClassifier(max_depth = MD, learning_rate = LR, subsample = 0.9, max_delta_step = 0.0, min_child_weight = 0.0, reg_alpha = 0.0, reg_lambda = 4.0, colsample_bytree = 0.65, verbosity = 1, n_estimators=nEst, nthread = 16, feature_names = inputVarNames)
evalSetE = [(xTrainE, yTrainE), (xTestE, yTestE)]
modelE.fit(xTrainE, yTrainE, sample_weight = wTrainE, eval_metric=["auc", "logloss"], eval_set=evalSetE, early_stopping_rounds=50, verbose=True)


gc.collect()
pickle.dump(modelB, open(modelNamesB[0], "wb"))
pickle.dump(modelE, open(modelNamesE[0], "wb"))

##18 is the event index determined before splitting. Allows Lead and Sub to be recombined after evaluation
##So XGB results are an array with the score result and the event index to be recombined into diphoton event
yPredTrainB = np.stack((modelB.predict_proba(xTrainB)[:,1],varValuesTrainB[:,18]),axis=1)
yPredTestB = np.stack((modelB.predict_proba(xTestB)[:,1],varValuesTestB[:,18]),axis=1)

yPredTrainE = np.stack((modelE.predict_proba(xTrainE)[:,1],varValuesTrainE[:,18]),axis=1)
yPredTestE = np.stack((modelE.predict_proba(xTestE)[:,1],varValuesTestE[:,18]),axis=1)

#Split into sig/bkg arrays by label
yPredTrainSigB = yPredTrainB[yTrainB[:] == 1]
yPredTrainBkgB = yPredTrainB[yTrainB[:] != 1]
yPredTrainSigE = yPredTrainE[yTrainE[:] == 1]
yPredTrainBkgE = yPredTrainE[yTrainE[:] != 1]

varsTrainSigB = varValuesTrainB[yTrainB[:] == 1]
varsTrainBkgB = varValuesTrainB[yTrainB[:] != 1]
varsTrainSigE = varValuesTrainE[yTrainE[:] == 1]
varsTrainBkgE = varValuesTrainE[yTrainE[:] != 1]

print len(yPredTrainB)
print len(yPredTrainE)
print len(yPredTrainSigB)
print len(yPredTrainSigE)

#Add barrel and endcap for sig and bkg seperately
yPredTrainSigTotal = np.concatenate((yPredTrainSigB,yPredTrainSigE),axis = 0)
yPredTrainBkgTotal = np.concatenate((yPredTrainBkgE,yPredTrainBkgE),axis = 0)

print len(yPredTrainSigTotal)

varsTrainSigTotal = np.concatenate((varsTrainSigB,varsTrainSigE), axis = 0)
varsTrainBkgTotal = np.concatenate((varsTrainBkgB,varsTrainBkgE), axis = 0)

#Sort the sig and bkg by event index, so events with same index are next to each other
#yPredTrainSigTotal = np.argsort(yPredTrainSigTotal[:,1])
yPredTrainSigTotal = sorted(yPredTrainSigTotal, key=itemgetter(1))
yPredTrainBkgTotal = sorted(yPredTrainBkgTotal, key=itemgetter(1))

varsTrainSigTotal = sorted(varsTrainSigTotal, key=itemgetter(18))
varsTrainBkgTotal = sorted(varsTrainBkgTotal, key=itemgetter(18))


print len(yPredTrainSigTotal)

endTime = time.time()
timeSpent = endTime - startTime
timeSpentMin = timeSpent/60.0
print 'Run Time(min) =  {0:.4f}'.format(timeSpentMin)
#
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

sigLen = len(varsTrainSigTotal)
bkgLen = len(varsTrainBkgTotal)
print sigLen
print bkgLen

#for n in range(0,sigLen):
    

#    xgbScoreSig[0] = yPredSigLeadTotal[n][0]
#    xgbScoreSig[1] = yPredSigSubTotal[n][0]
#    rawEnergySig[0] = varValuesSigLead[0][n]
#    rawEnergySig[1] = varValuesSigSub[0][n]
#    r9HLTSig[0] = varValuesSigLead[1][n]
#    r9HLTSig[1] = varValuesSigSub[1][n]
#    sigmaIEtaIEtaSig[0] = varValuesSigLead[2][n]
#    sigmaIEtaIEtaSig[1] = varValuesSigSub[2][n]
#    etaWidthSig[0] = varValuesSigLead[3][n]
#    etaWidthSig[1] = varValuesSigSub[3][n]
#    phiWidthSig[0] = varValuesSigLead[4][n]
#    phiWidthSig[1] = varValuesSigSub[4][n]
#    s4Sig[0] = varValuesSigLead[5][n]
#    s4Sig[1] = varValuesSigSub[5][n]
#    trkIsoPhoSig[0] = varValuesSigLead[6][n]
#    trkIsoPhoSig[1] = varValuesSigSub[6][n]
#    etaSig[0] = varValuesSigLead[7][n]
#    etaSig[1] = varValuesSigSub[7][n]
#    hOvrESig[0] = varValuesSigLead[8][n]
#    hOvrESig[1] = varValuesSigSub[8][n]
#    ecalPFIsoSig[0] = varValuesSigLead[9][n]
#    ecalPFIsoSig[1] = varValuesSigSub[9][n]
#    etSig[0] = varValuesSigLead[10][n]
#    etSig[1] = varValuesSigSub[10][n]
#    energySig[0] = varValuesSigLead[11][n]
#    energySig[1] = varValuesSigSub[11][n]
#    massSig[0] = varValuesSigLead[12][n]
#    massSig[1] = varValuesSigSub[12][n]
#    nEgsSig[0] = int(varValuesSigLead[13][n])
#    nEgsSig[1] = int(varValuesSigLead[13][n])
#    passFailStdSig[0] = int(varValuesSigLead[14][n])
#    passFailL1SingleSig[0] = int(varValuesSigLead[15][n])
#    passFailL1DoubleSig[0] = int(varValuesSigLead[16][n])
#    triggerBitsSig[0] = int(varValuesSigLead[17][n])
#    triggerBitsSig[1] = int(varValuesSigLead[17][n])
#    tSig.Fill()

#xgbScoreBkg = array('f',[0.0,0.0])
#rawEnergyBkg = array('f',[0.0,0.0])
#r9HLTBkg = array('f',[0.0,0.0])
#sigmaIEtaIEtaBkg = array('f',[0.0,0.0])
#etaWidthBkg = array('f',[0.0,0.0])
#phiWidthBkg = array('f',[0.0,0.0])
#s4Bkg = array('f',[0.0,0.0])
#trkIsoPhoBkg = array('f',[0.0,0.0])
#etaBkg = array('f',[0.0,0.0])
#hOvrEBkg = array('f',[0.0,0.0])
#ecalPFIsoBkg = array('f',[0.0,0.0])
#etBkg = array('f',[0.0,0.0])
#energyBkg = array('f',[0.0,0.0])
#massBkg = array('f',[0.0,0.0])
#nEgsBkg = array('f',[0.0,0.0])
#triggerBitsBkg = array('f',[0.0,0.0])
#passFailStdBkg = array('i',[0])
#passFailL1SingleBkg = array('i',[0])
#passFailL1DoubleBkg = array('i',[0])
#
#tBkg = TTree("bkgTree", "bkgTree")
#tBkg.Branch( 'xgbScore', xgbScoreBkg, 'xgbScore[2]/F' )
#tBkg.Branch( 'rawEnergy', rawEnergyBkg, 'rawEnergy[2]/F' )
#tBkg.Branch( 'r9HLT', r9HLTBkg, 'r9HLT[2]/F' )
#tBkg.Branch( 'sigmaIEtaIEta', sigmaIEtaIEtaBkg, 'sigmaIEtaIEta[2]/F' )
#tBkg.Branch( 'etaWidth', etaWidthBkg, 'etaWidth[2]/F' )
#tBkg.Branch( 'phiWidth', phiWidthBkg, 'phiWidth[2]/F' )
#tBkg.Branch( 's4', s4Bkg, 's4[2]/F' )
#tBkg.Branch( 'trkIsoPho', trkIsoPhoBkg, 'trkIsoPho[2]/F' )
#tBkg.Branch( 'eta', etaBkg, 'eta[2]/F' )
#tBkg.Branch( 'hOvrE', hOvrEBkg, 'hOvrE[2]/F' )
#tBkg.Branch( 'ecalPFIso', ecalPFIsoBkg, 'ecalPFIso[2]/F' )
#tBkg.Branch( 'et', etBkg, 'et[2]/F' )
#tBkg.Branch( 'energy', energyBkg, 'energy[2]/F' )
#tBkg.Branch( 'mass', massBkg, 'mass[2]/F' )
#tBkg.Branch( 'nEgs', nEgsBkg, 'nEgs[2]/F' )
#tBkg.Branch( 'triggerBits', triggerBitsBkg, 'triggerBits[2]/F' )
#tBkg.Branch( 'passFailStd', passFailStdBkg, 'passFailStd/I' )
#tBkg.Branch( 'passFailL1Single', passFailL1SingleBkg, 'passFailL1Single/I' )
#tBkg.Branch( 'passFailL1Double', passFailL1DoubleBkg, 'passFailL1Double/I' )
#
#bkgLen = len(yPredBkgLeadTotal)
#
#for n in range(0,bkgLen):
#    xgbScoreBkg[0] = yPredBkgLeadTotal[n][0]
#    xgbScoreBkg[1] = yPredBkgSubTotal[n][0]
#    rawEnergyBkg[0] = varValuesBkgLead[0][n]
#    rawEnergyBkg[1] = varValuesBkgSub[0][n]
#    r9HLTBkg[0] = varValuesBkgLead[1][n]
#    r9HLTBkg[1] = varValuesBkgSub[1][n]
#    sigmaIEtaIEtaBkg[0] = varValuesBkgLead[2][n]
#    sigmaIEtaIEtaBkg[1] = varValuesBkgSub[2][n]
#    etaWidthBkg[0] = varValuesBkgLead[3][n]
#    etaWidthBkg[1] = varValuesBkgSub[3][n]
#    phiWidthBkg[0] = varValuesBkgLead[4][n]
#    phiWidthBkg[1] = varValuesBkgSub[4][n]
#    s4Bkg[0] = varValuesBkgLead[5][n]
#    s4Bkg[1] = varValuesBkgSub[5][n]
#    trkIsoPhoBkg[0] = varValuesBkgLead[6][n]
#    trkIsoPhoBkg[1] = varValuesBkgSub[6][n]
#    etaBkg[0] = varValuesBkgLead[7][n]
#    etaBkg[1] = varValuesBkgSub[7][n]
#    hOvrEBkg[0] = varValuesBkgLead[8][n]
#    hOvrEBkg[1] = varValuesBkgSub[8][n]
#    ecalPFIsoBkg[0] = varValuesBkgLead[9][n]
#    ecalPFIsoBkg[1] = varValuesBkgSub[9][n]
#    etBkg[0] = varValuesBkgLead[10][n]
#    etBkg[1] = varValuesBkgSub[10][n]
#    energyBkg[0] = varValuesBkgLead[11][n]
#    energyBkg[1] = varValuesBkgSub[11][n]
#    massBkg[0] = varValuesBkgLead[12][n]
#    massBkg[1] = varValuesBkgSub[12][n]
#    nEgsBkg[0] = int(varValuesBkgLead[13][n])
#    nEgsBkg[1] = int(varValuesBkgLead[13][n])
#    passFailStdBkg[0] = int(varValuesBkgLead[14][n])
#    passFailL1SingleBkg[0] = int(varValuesBkgLead[15][n])
#    passFailL1DoubleBkg[0] = int(varValuesBkgLead[16][n])
#    triggerBitsBkg[0] = int(varValuesBkgLead[17][n])
#    triggerBitsBkg[1] = int(varValuesBkgLead[17][n])
#    tBkg.Fill()
#
#
#outFileRoot.Write()
#outFileRoot.Close()
    
