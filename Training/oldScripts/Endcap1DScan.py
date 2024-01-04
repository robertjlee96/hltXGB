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
import EndcapPtUtils
import EndcapTMVAUtils
#import oldXGBoost2tmva
import time,pickle
from tqdm import tqdm
import sys
from tempfile import TemporaryFile
import time
from array import array
#import joblib
#import gc
gROOT.ProcessLine(
                "struct myStruct {\
                Double_t     yVal;\
                Double_t     yPred;\
                };" );

from ROOT import myStruct

# for sklearn, see
np.random.seed(1337)

nEst = 1500
fin = uproot.open("../NTuples/outNTuples/GGH13andData_RelaxedV24_TrainingCuts_0502_Train.root")

prompt = fin['sigTree']
fake = fin['bkgTree']
## for barrel
geometry_selection = lambda tree: np.abs(tree.array('eta')) > 1.556
#dptCut = lambda tree: abs((tree.array('genPt')-tree.array('pt'))/tree.array('genPt')) < 0.1
ptCut = lambda tree: tree.array('et') > 14.25
massCut = lambda tree: tree.array('mass') > 90.0


input_values, target_values, orig_weights, train_weights, pt, scEta, input_vars = EndcapPtUtils.load_file(fin,geometry_selection,ptCut,massCut)

print "input_vars", input_vars
# ### split into training and test set
X_train, X_test, w_train, dummy, y_train, y_test, dummy, w_test, pt_train, pt_test, scEta_train, scEta_test = train_test_split(input_values,train_weights,target_values,orig_weights,pt,scEta,test_size=0.25)

#MD = 20
#nModels = 5
#mdVals = [13,15,17,19,20]
#lrVals = [0.1,0.09,0.07,0.05,0.045]
#mdStrings = ["13","15","17","19","20"]
#lrStrings = ["1","09","07","05","045"]
nModels = 3
mdVals = [9,13,17]
lrVals = [0.1,0.07,0.5]
mdStrings = ["9","13","17"]
lrStrings = ["10","07","05"]
#lrVals = [0.030,0.06]
#lrVals = [0.06,0.07]
#lrVals = [0.045]
#LR = 0.045

#for MD in mdVals:
for i in range(0,nModels):
    nNow = 0
    MD = mdVals[i]
    LR = lrVals[i]
    fNameGen = "endcapOut/Model_NewGGH_MD" + mdStrings[i] + "LR" + lrStrings[i] + "_M90_Endcap_0504"
    fName = fNameGen + ".npz"
    outFile = open(fName,"w")
    fNameRoot = fNameGen + ".root"
    outFileRoot = TFile(fNameRoot,'recreate')

    xgBoostParams = []
    aucTest = []
    aucTrain = []
    lossTest = []
    lossTrain = []
    mvaValsTest = []
    mvaValsTrain = []
    bestVals = []
    fprVals = []
    tprVals = []
    timeVals = []

    startTime = time.time()

    model = XGBClassifier(max_depth = MD, learning_rate = LR, subsample = 0.9, max_delta_step = 0.0, min_child_weight = 0.0, reg_alpha = 0.0, reg_lambda = 4.0, colsample_bytree = 0.65, verbosity = 1, n_estimators=nEst, nthread = 16, feature_names = input_vars)

    boostParamsTmp = []
    boostParamsTmp.append(model.get_xgb_params()["max_depth"])
    boostParamsTmp.append(model.get_xgb_params()["learning_rate"])
    boostParamsTmp.append(model.get_xgb_params()["gamma"])
    boostParamsTmp.append(model.get_xgb_params()["min_child_weight"])
    boostParamsTmp.append(model.get_xgb_params()["max_delta_step"])
    boostParamsTmp.append(model.get_xgb_params()["colsample_bytree"])
    boostParamsTmp.append(model.get_xgb_params()["subsample"])
    boostParamsTmp.append(model.get_xgb_params()["reg_alpha"])
    boostParamsTmp.append(model.get_xgb_params()["reg_lambda"])
    xgBoostParams.append(boostParamsTmp)
            
    eval_set = [(X_train, y_train), (X_test, y_test)]

    model.fit(X_train, y_train, sample_weight = w_train, eval_metric=["auc","logloss"], eval_set=eval_set, early_stopping_rounds=50, verbose=True)

#    model_name = "Model_MD18_LR" + LR + "_M95PTM25_OnlyPFPairs_UL2017_0113"
    model_fname = fNameGen + ".model"
    model_xmlName = fNameGen + ".xml"

    input_variables = []
    for name in input_vars:#model.get_booster().feature_names:
        input_variables.append((name, 'F'))
    #    print "name in config file = {}".format(name)

    #MoreVarsTmvaUtils.convert_model(model.get_booster().get_dump(), input_variables = input_variables, output_xml=model_xmlName)

    #gc.collect()
    pickle.dump(model, open(model_fname, "wb"))
    #joblib.dump(model, model_fname)
     
    ###FIRST, EVALUATE AUC
    y_train_pred = model.predict_proba(X_train)[:,1]
    y_test_pred = model.predict_proba(X_test)[:,1]

    mvaValsTest.append(y_test_pred)
    mvaValsTrain.append(y_train_pred)
        
    fpr, tpr, thresholds = roc_curve(y_test, y_test_pred,sample_weight=w_test)
    fprVals.append(fpr)
    tprVals.append(tpr)
    auc_val_test = metrics.auc(fpr, tpr)

    #NOW, write fpr and tpr to a root file
    xgbScoreTrainSig = array('f',[0.0])
    xgbScoreTestSig = array('f',[0.0])
#    etaSig = array('f',[0.0])
#    nEgsSig = array('f',[0.0])
#    etSig = array('f',[0.0])
#    massSig = array('f',[0.0])
#    passFailStdSig = array('i',[0])

    tSigTrain = TTree("sigTreeTrain", "sigTreeTrain")
    tSigTrain.Branch( 'xgbScore', xgbScoreTrainSig, 'xgbScore/F' )

    tSigTest = TTree("sigTreeTest", "sigTreeTest")
    tSigTest.Branch( 'xgbScore', xgbScoreTestSig, 'xgbScore/F' )
#    tSig.Branch( 'eta', etaSig, 'eta/F' )
#    tSig.Branch( 'nEgs', nEgsSig, 'nEgs/F' )
#    tSig.Branch( 'et', etSig, 'et/F' )
#    tSig.Branch( 'mass', massSig, 'mass/F' )
#    tSig.Branch( 'passFailStd', passFailStdSig, 'passFailStd/I' )
    
    xgbScoreTrainBkg = array('f',[0.0])
    xgbScoreTestBkg = array('f',[0.0])
#    etaBkg = array('f',[0.0])
#    nEgsBkg = array('f',[0.0])
#    etBkg = array('f',[0.0])
#    massBkg = array('f',[0.0])
#    passFailStdBkg = array('i',[0])

    tBkgTrain = TTree("bkgTreeTest", "bkgTreeTest")
    tBkgTrain.Branch( 'xgbScore', xgbScoreTrainBkg, 'xgbScore/F' )
    tBkgTest = TTree("bkgTreeTest", "bkgTreeTest")
    tBkgTest.Branch( 'xgbScore', xgbScoreTestBkg, 'xgbScore/F' )
#    tBkg.Branch( 'eta', etaBkg, 'eta/F' )
#    tBkg.Branch( 'nEgs', nEgsBkg, 'nEgs/F' )
#    tBkg.Branch( 'et', etBkg, 'et/F' )
#    tBkg.Branch( 'mass', massBkg, 'mass/F' )
#    tBkg.Branch( 'passFailStd', passFailStdBkg, 'passFailStd/I' )

    for n in range(0,len(y_train)):
        if y_train[n] == 1:
            xgbScoreTrainSig[0] = y_train_pred[n]
            tSigTrain.Fill()
        if y_train[n] == 0:
            xgbScoreTrainBkg[0] = y_train_pred[n]
            tBkgTrain.Fill()

    for n in range(0,len(y_test)):
        if y_test[n] == 1:
            xgbScoreTestSig[0] = y_test_pred[n]
            tSigTest.Fill()
        if y_test[n] == 0:
            xgbScoreTestBkg[0] = y_test_pred[n]
            tBkgTest.Fill()

    outFileRoot.Write()
    outFileRoot.Close()

    ###NEXT, EVALUATE LOG LOSS
    results = model.evals_result() ##NO TEST REWEIGHTING
                
    aucTest.append(results['validation_1']['auc'])
    aucTrain.append(results['validation_0']['auc'])
    lossTest.append(results['validation_1']['logloss'])
    lossTrain.append(results['validation_0']['logloss'])
            
    bestRound = model.best_iteration

    bestValsTmp = [bestRound,results['validation_0']['auc'][bestRound],results['validation_1']['logloss'][bestRound],results['validation_0']['logloss'][bestRound]]
                
    bestVals.append(bestValsTmp)
                
    print 'Done Training {0:.0f}'.format(nNow)
        
    endTime = time.time()
    timeSpent = endTime - startTime
    timeSpentHR = timeSpent/3600.0
    print 'Run Time(hr) =  {0:.4f}'.format(timeSpentHR)
    timeVals.append(timeSpent)

    EndcapTMVAUtils.variable_importance(model,input_variables,fNameGen)

        
    del results
    del model

    np.savez(outFile, xgBoostParams = xgBoostParams, aucTest = aucTest, aucTrain = aucTrain, lossTest = lossTest, lossTrain = lossTrain, mvaValsTest = mvaValsTest, mvaValsTrain = mvaValsTrain, bestVals = bestVals, fprVals = fprVals, tprVals = tprVals, timeVals = timeVals)


