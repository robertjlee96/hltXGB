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

inFileName ="../NTuples/SinglePhoton_GGHandGJet_StdPresel_Validation_0602.root"
modelNamesB = []
#modelNamesB.append("barrel/Model_MD10_LR01_M95PTM15_UpperBoundMass_NoPTReweight_UL2017_0525.model")
modelNamesB.append("barrel/Model_MD20LR045_M95PTM15_UL2017_0606.model")
modelNamesB.append("barrel/Model_MD19LR05_M95PTM15_UL2017_0606.model")
modelNamesB.append("barrel/Model_MD17LR07_M95PTM15_UL2017_0606.model")
modelNamesB.append("barrel/Model_MD15LR09_M95PTM15_UL2017_0606.model")
modelNamesB.append("barrel/Model_MD13LR1_M95PTM15_UL2017_0606.model")




#modelNamesB.append("barrel/XGB_Model_Barrel_SA_phoID_UL2017_woCorr.pkl")
modelNamesE = []
modelNamesE.append("endcap/Model_MD20LR045_M95PTM15_UL2017_0606_Endcap.model")
modelNamesE.append("endcap/Model_MD19LR05_M95PTM15_UL2017_0606_Endcap.model")
modelNamesE.append("endcap/Model_MD17LR07_M95PTM15_UL2017_0606_Endcap.model")
modelNamesE.append("endcap/Model_MD15LR09_M95PTM15_UL2017_0606_Endcap.model")
modelNamesE.append("endcap/Model_MD13LR1_M95PTM15_UL2017_0606_Endcap.model")

#modelNamesE.append("endcap/Model_MD10_LR01_M95PTM15_UpperBoundMass_NoPTReweight_UL2017_0525_Endcap.model")
#modelNamesE.append("endcap/XGB_Model_Endcap_SA_phoID_UL2017_woCorr.pkl")

fNames = []
fNames.append("validationNTuples/0607/GGHandGJet_UL2017_MD20LR045")
fNames.append("validationNTuples/0607/GGHandGJet_UL2017_MD19LR05")
fNames.append("validationNTuples/0607/GGHandGJet_UL2017_MD17LR07")
fNames.append("validationNTuples/0607/GGHandGJet_UL2017_MD15LR09")
fNames.append("validationNTuples/0607/GGHandGJet_UL2017_MD13LR1")


fin = uproot.open(inFileName)
prompt = fin['promptPhotons']
#fake = fin['fakePhotons']
fake = fin['promptPhotons']
#fake = fin['diphotonBkg']
## for barrel
geometry_selection = lambda tree: np.abs(tree.array('scEta')) < 1.4442
dptCut = lambda tree: abs((tree.array('genPt')-tree.array('pt'))/tree.array('genPt')) < 0.1
ptCut = lambda tree: tree.array('pt') > 0.0

inputValuesSigLead, inputValuesSigSub, inputValuesPrompt, inputValuesFake, targetValuesSigLead, targetValuesSigSub, targetValuesPrompt,targetValuesFake, origWeightsSigLead, origWeightsSigSub,origWeightsPrompt,origWeightsFake, varValuesSigLead, varValuesSigSub, varValuesPrompt, varValuesFake, inputVarNames, varNames = diphotonUtils.load_file(fin,geometry_selection,ptCut)


for i in range(0,len(fNames)):
    fNameRoot = fNames[i] + ".root"
    outFileRoot = TFile(fNameRoot,'recreate')
   
#    FOR MODELS WITH H/E VARIABLES
    barrelSigLead = (((inputValuesSigLead[np.abs(varValuesSigLead[:,10]) < 1.4442]).T)[0:14]).T #without Pt/M
#    barrelSigLead = (((inputValuesSigLead[np.abs(varValuesSigLead[:,10]) < 1.4442]).T)[0:15]).T #With Pt/M
    barrelSigLeadVars = varValuesSigLead[np.abs(varValuesSigLead[:,10]) < 1.4442]

    barrelSigSub = (((inputValuesSigSub[np.abs(varValuesSigSub[:,10]) < 1.4442]).T)[0:14]).T #without Pt/M
#    barrelSigSub = (((inputValuesSigSub[np.abs(varValuesSigSub[:,10]) < 1.4442]).T)[0:15]).T #With Pt/M
    barrelSigSubVars = varValuesSigSub[np.abs(varValuesSigSub[:,10]) < 1.4442]

    barrelPrompt = (((inputValuesPrompt[np.abs(varValuesPrompt[:,10]) < 1.4442]).T)[0:14]).T #without Pt/M
#    barrelPrompt = (((inputValuesPrompt[np.abs(varValuesPrompt[:,10]) < 1.4442]).T)[0:15]).T #With Pt/M
    barrelPromptVars = varValuesPrompt[np.abs(varValuesPrompt[:,10]) < 1.4442]

    barrelFake = (((inputValuesFake[np.abs(varValuesFake[:,10]) < 1.4442]).T)[0:14]).T #without Pt/M
#    barrelFake = (((inputValuesFake[np.abs(varValuesFake[:,10]) < 1.4442]).T)[0:15]).T #With Pt/M
    barrelFakeVars = varValuesFake[np.abs(varValuesFake[:,10]) < 1.4442]

    endcapSigLead = inputValuesSigLead[np.abs(varValuesSigLead[:,10]) > 1.556]
    endcapSigLeadVars = varValuesSigLead[np.abs(varValuesSigLead[:,10]) > 1.556]

    endcapSigSub = inputValuesSigSub[np.abs(varValuesSigSub[:,10]) > 1.556]
    endcapSigSubVars = varValuesSigSub[np.abs(varValuesSigSub[:,10]) > 1.556]

    endcapPrompt = inputValuesPrompt[np.abs(varValuesPrompt[:,10]) > 1.556]
    endcapPromptVars = varValuesPrompt[np.abs(varValuesPrompt[:,10]) > 1.556]

    endcapFake = inputValuesFake[np.abs(varValuesFake[:,10]) > 1.556]
    endcapFakeVars = varValuesFake[np.abs(varValuesFake[:,10]) > 1.556]

#    FOR MODELS WITHOUT H/E VARIABLES
#    barrelSigLead = (((inputValuesSigLead[np.abs(varValuesSigLead[:,10]) < 1.4442]).T)[0:12]).T
#    barrelSigLeadVars = varValuesSigLead[np.abs(varValuesSigLead[:,10]) < 1.4442]
#
#    barrelSigSub = (((inputValuesSigSub[np.abs(varValuesSigSub[:,10]) < 1.4442]).T)[0:12]).T
#    barrelSigSubVars = varValuesSigSub[np.abs(varValuesSigSub[:,10]) < 1.4442]
#
#    barrelPrompt = (((inputValuesPrompt[np.abs(varValuesPrompt[:,10]) < 1.4442]).T)[0:12]).T
#    barrelPromptVars = varValuesPrompt[np.abs(varValuesPrompt[:,10]) < 1.4442]
#
#    barrelFake = (((inputValuesFake[np.abs(varValuesFake[:,10]) < 1.4442]).T)[0:12]).T
#    barrelFakeVars = varValuesFake[np.abs(varValuesFake[:,10]) < 1.4442]
#
#    endcapSigLead = (((inputValuesSigLead[np.abs(varValuesSigLead[:,10]) > 1.556]).T)[0:14]).T
#    endcapSigLeadVars = varValuesSigLead[np.abs(varValuesSigLead[:,10]) > 1.556]
#
#    endcapSigSub = (((inputValuesSigSub[np.abs(varValuesSigSub[:,10]) > 1.556]).T)[0:14]).T
#    endcapSigSubVars = varValuesSigSub[np.abs(varValuesSigSub[:,10]) > 1.556]
#
#    endcapPrompt = (((inputValuesPrompt[np.abs(varValuesPrompt[:,10]) > 1.556]).T)[0:14]).T
#    endcapPromptVars = varValuesPrompt[np.abs(varValuesPrompt[:,10]) > 1.556]
#
#    endcapFake = (((inputValuesFake[np.abs(varValuesFake[:,10]) > 1.556]).T)[0:14]).T
#    endcapFakeVars = varValuesFake[np.abs(varValuesFake[:,10]) > 1.556]

    startTime = time.time()
   
    modelFileB = open(modelNamesB[i],'rb')
    modelB = pickle.load(modelFileB)
    
    yPredSigLeadB = np.stack((modelB.predict_proba(barrelSigLead)[:,1],barrelSigLeadVars[:,23]),axis=1)
    yPredSigSubB = np.stack((modelB.predict_proba(barrelSigSub)[:,1],barrelSigSubVars[:,23]),axis=1)
    yPredPromptB = np.stack((modelB.predict_proba(barrelPrompt)[:,1],barrelPromptVars[:,23]),axis=1)
    yPredFakeB = np.stack((modelB.predict_proba(barrelFake)[:,1],barrelFakeVars[:,23]),axis=1)
    
    
    modelFileE = open(modelNamesE[i],'rb')
    modelE = pickle.load(modelFileE)

    yPredSigLeadE = np.stack((modelE.predict_proba(endcapSigLead)[:,1],endcapSigLeadVars[:,23]),axis=1)
    yPredSigSubE = np.stack((modelE.predict_proba(endcapSigSub)[:,1],endcapSigSubVars[:,23]),axis=1)
    yPredPromptE = np.stack((modelE.predict_proba(endcapPrompt)[:,1],endcapPromptVars[:,23]),axis=1)
    yPredFakeE = np.stack((modelE.predict_proba(endcapFake)[:,1],endcapFakeVars[:,23]),axis=1)
    
    yPredSigLeadTotal = np.row_stack((yPredSigLeadB,yPredSigLeadE))
    yPredSigLeadTotal[np.argsort(yPredSigLeadTotal[:,1])]
    yPredSigSubTotal = np.row_stack((yPredSigSubB,yPredSigSubE))
    yPredSigSubTotal[np.argsort(yPredSigSubTotal[:,1])]
    yPredPromptTotal = np.row_stack((yPredPromptB,yPredPromptE))
    yPredPromptTotal[np.argsort(yPredPromptTotal[:,1])]
    yPredFakeTotal = np.row_stack((yPredFakeB,yPredFakeE))
    yPredFakeTotal[np.argsort(yPredFakeTotal[:,1])]

#    auc_val_test = metrics.roc_auc_score(target_values,y_pred,sample_weight=train_weights)
#    print auc_val_test
    
    endTime = time.time()
    timeSpent = endTime - startTime
    timeSpentMin = timeSpent/60.0
    print 'Run Time(min) =  {0:.4f}'.format(timeSpentMin)

    xgbScoreSigLead = array('f',[0.0])
    weightSigLead = array('f',[0.0])
    varValsSigLead = array('f',23*[0.0])
    
    xgbScoreSigSub = array('f',[0.0])
    weightSigSub = array('f',[0.0])
    varValsSigSub = array('f',23*[0.0])

    xgbScorePrompt = array('f',[0.0])
    weightPrompt = array('f',[0.0])
    varValsPrompt = array('f',23*[0.0])
    
    xgbScoreFake = array('f',[0.0])
    weightFake = array('f',[0.0])
    varValsFake = array('f',23*[0.0])

    tSigLead = TTree("sigLead", "Signal Lead Photons")
    tSigLead.Branch( 'xgbScore', xgbScoreSigLead, 'xgbScore/F' )
    tSigLead.Branch( 'weight', weightSigLead, 'weight/F' )
    tSigLead.Branch( 'varVals', varValsSigLead, 'varVals[23]/F' )

    tSigSub = TTree("sigSub", "Signal SubLead Photons")
    tSigSub.Branch( 'xgbScore', xgbScoreSigSub, 'xgbScore/F' )
    tSigSub.Branch( 'weight', weightSigSub, 'weight/F' )
    tSigSub.Branch( 'varVals', varValsSigSub, 'varVals[23]/F' )

    tPrompt = TTree("promptPhotons", "GJet Prompt")
    tPrompt.Branch( 'xgbScore', xgbScorePrompt, 'xgbScore/F' )
    tPrompt.Branch( 'weight', weightPrompt, 'weight/F' )
    tPrompt.Branch( 'varVals', varValsPrompt, 'varVals[23]/F' )
    
    tFake = TTree("fakePhotons", "GJet Fake")
    tFake.Branch( 'xgbScore', xgbScoreFake, 'xgbScore/F' )
    tFake.Branch( 'weight', weightFake, 'weight/F' )
    tFake.Branch( 'varVals', varValsFake, 'varVals[23]/F' )

    xgbScoreSigLeadArray = []
    weightSigLeadArray = []
    xgbScoreSigSubArray = []
    weightSigSubArray = []

    xgbScorePromptArray = []
    weightPromptArray = []
   
    xgbScoreFakeArray = []
    weightFakeArray = []

#    sigSubLen = len(yPredSigLead)
#    for n in range (0,sigSubLen):
#        if targetValuesSigSub[n] == 1:
#            xgbScoreSigArray.append(y_pred[n])
#            weightSigArray.append(orig_weights[n])
#        if target_values[n] == 0:
#            xgbScoreBkgArray.append(y_pred[n])
#            weightBkgArray.append(orig_weights[n])
#
    sigLeadLen = len(yPredSigLeadTotal)
    for n in range(0,sigLeadLen):
        xgbScoreSigLead[0] = yPredSigLeadTotal[n][0]
        weightSigLead[0] = origWeightsSigLead[0][n]
        for i in range(0,23):
            varValsSigLead[i] = varValuesSigLead[n][i]
        tSigLead.Fill()
    
    sigSubLen = len(yPredSigSubTotal)
    for n in range(0,sigSubLen):
        xgbScoreSigSub[0] = yPredSigSubTotal[n][0]
        weightSigSub[0] = origWeightsSigSub[0][n]
        for i in range(0,23):
            varValsSigSub[i] = varValuesSigSub[n][i]
        tSigSub.Fill()
    
    promptLen = len(yPredPromptTotal)
    for n in range(0,promptLen):
        xgbScorePrompt[0] = yPredPromptTotal[n][0]
        weightPrompt[0] = origWeightsPrompt[0][n]
        for i in range(0,23):
            varValsPrompt[i] = varValuesPrompt[n][i]
        tPrompt.Fill()
    
    fakeLen = len(yPredFakeTotal)

    for n in range(0,fakeLen):
        xgbScoreFake[0] = yPredFakeTotal[n][0]
        weightFake[0] = origWeightsFake[n]
        for i in range(0,23):
            varValsFake[i] = varValuesFake[n][i]
        tFake.Fill()
#
#    for n in range(0,bkgLen):
#        xgbScoreBkg[0] = xgbScoreBkgArray[n]
#        weightBkg[0] = weightBkgArray[n]
#        for i in range(0,22):
#            varValsBkg[i] = varValsBkgIn[n][i]
#        tFake.Fill()
#
    outFileRoot.Write()
    outFileRoot.Close()
#













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
