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

inFileName ="../../CMSSW_12_4_8/src/allPhotonsNTuples/GGH13andDataB_DiphotonValidation_0724.root"
modelNamesB = []
modelNamesE = []
fNames = []

mNamesGen =  []

mNamesGen.append("0728/M9LR15_GGH13andData_M60_0728")

for mName in mNamesGen:
    modelNamesB.append("../Training/barrelOut/" + mName + "_Barrel.model")
    modelNamesE.append("../Training/endcapOut/" + mName + "_Endcap.model")
    fNames.append("validationNTuples/" + mName + "Validated")

## for barrel
geometry_selection = lambda tree: np.abs(tree.array('eta')) < 2.5
ptCut = lambda tree: tree.array('et') > 0.0
#massCut = lambda tree: tree.array('mass') > 60.0
#massCut = lambda tree: tree.array('mass') > 0.0

#inputValuesSig, inputValuesBkg, varValuesSig, varValuesBkg, inputVarNames, varNames = validationUtils.load_file(fin,geometry_selection,ptCut)
inputValuesSig, inputValuesBkg, varValuesSig, varValuesBkg, inputVarNames, varNames, nEventsSig, nEventsBkg = validationUtils.load_file(inFileName,geometry_selection,ptCut)

#inputValuesSigLead, inputValuesSigSub, inputValuesBkgLead, inputValuesBkgSub, targetValuesSigLead, targetValuesSigSub, targetValuesBkgLead,targetValuesBkgSub, origWeightsSigLead, origWeightsSigSub,origWeightsBkgLead,origWeightsBkgSub, varValuesSigLead, varValuesSigSub, varValuesBkgLead, varValuesBkgSub, inputVarNames, varNames = diphotonUtils.load_file(fin,geometry_selection,ptCut)

for i in range(0,len(fNames)):
    fNameRoot = fNames[i] + ".root"
    outFileRoot = TFile(fNameRoot,'recreate')

    inputBarrelSig = inputValuesSig.T[np.abs(inputValuesSig[7,:]) < 1.5]
    varsBarrelSig = varValuesSig.T[np.abs(varValuesSig[7,:]) < 1.5]
    
    inputBarrelBkg = inputValuesBkg.T[np.abs(inputValuesBkg[7,:]) < 1.5]
    varsBarrelBkg = varValuesBkg.T[np.abs(varValuesBkg[7,:]) < 1.5]

    modelFileB = open(modelNamesB[i],'rb')
    modelB = pickle.load(modelFileB)
    
    inputEndcapSig = inputValuesSig.T[np.abs(inputValuesSig[7,:]) > 1.5]
    varsEndcapSig = varValuesSig.T[np.abs(varValuesSig[7,:]) > 1.5]
    
    inputEndcapBkg = inputValuesBkg.T[np.abs(inputValuesBkg[7,:]) > 1.5]
    varsEndcapBkg = varValuesBkg.T[np.abs(varValuesBkg[7,:]) > 1.5]

    modelFileE = open(modelNamesE[i],'rb')
    modelE = pickle.load(modelFileE)

    startTime = time.time()

    #21 is the event index determined before splitting. Allows Lead and Sub to be recombined after evaluation
    #So XGB results are an array with the score result and the event index to be recombined into diphoton event
    yPredSigB = np.stack((modelB.predict_proba(inputBarrelSig)[:,1],varsBarrelSig[:,21]),axis=1)
    yPredBkgB = np.stack((modelB.predict_proba(inputBarrelBkg)[:,1],varsBarrelBkg[:,21]),axis=1)

    yPredSigE = np.stack((modelE.predict_proba(inputEndcapSig)[:,1],varsEndcapSig[:,21]),axis=1)
    yPredBkgE = np.stack((modelE.predict_proba(inputEndcapBkg)[:,1],varsEndcapBkg[:,21]),axis=1)

    endTime = time.time()
    timeSpent = endTime - startTime
    print 'Run Time(sec) =  {0:.4f}'.format(timeSpent)


    yPredSigTotal = np.row_stack((yPredSigB,yPredSigE))
    yPredSigTotal = yPredSigTotal[np.argsort(yPredSigTotal[:,1])]


    yPredBkgTotal = np.row_stack((yPredBkgB,yPredBkgE))
    yPredBkgTotal = yPredBkgTotal[np.argsort(yPredBkgTotal[:,1])]
    
    varsSigTotal = np.row_stack((varsBarrelSig,varsEndcapSig))
    varsSigTotal = varsSigTotal[np.argsort(varsSigTotal[:,21])]

    varsBkgTotal = np.row_stack((varsBarrelBkg,varsEndcapBkg))
    varsBkgTotal = varsBkgTotal[np.argsort(varsBkgTotal[:,21])]
    
    maxPhotons = 50

#FIRST DO SIGNAL
    nEventSig = array('i',[0])
    nEgsSig = array('i',[0])
    nEgsPassingSig = array('i',[0])
    xgbScoreSig = array('f',[0.0]*maxPhotons)
    rawEnergySig = array('f',[0.0]*maxPhotons)
    r9HLTSig = array('f',[0.0]*maxPhotons)
    sigmaIEtaIEtaSig = array('f',[0.0]*maxPhotons)
    etaWidthSig = array('f',[0.0]*maxPhotons)
    phiWidthSig = array('f',[0.0]*maxPhotons)
    s4Sig = array('f',[0.0]*maxPhotons)
    trkIsoPhoSig = array('f',[0.0]*maxPhotons)
    etaSig = array('f',[0.0]*maxPhotons)
    hOvrESig = array('f',[0.0]*maxPhotons)
    ecalPFIsoSig = array('f',[0.0]*maxPhotons)
    etSig = array('f',[0.0]*maxPhotons)
    energySig = array('f',[0.0]*maxPhotons)
    phiSig = array('f',[0.0]*maxPhotons)
    massSig = array('f',[0.0])

    triggerBitsSig = array('i',[0])
    passFailStdSig = array('i',[0])
    passFailL1SingleSig = array('i',[0])
    passFailL1DoubleSig = array('i',[0])

    tSig = TTree("sigTree", "sigTree")
    tSig.Branch( 'nEvent', nEventSig, 'nEvent/I' )
    tSig.Branch( 'nEgs', nEgsSig, 'nEgs/I' )
    tSig.Branch( 'nEgsPassing', nEgsPassingSig, 'nEgsPassing/I' )
    tSig.Branch( 'xgbScore', xgbScoreSig, 'xgbScore[nEgs]/F' )
    tSig.Branch( 'rawEnergy', rawEnergySig, 'rawEnergy[nEgs]/F' )
    tSig.Branch( 'r9HLT', r9HLTSig, 'r9HLT[nEgs]/F' )
    tSig.Branch( 'sigmaIEtaIEta', sigmaIEtaIEtaSig, 'sigmaIEtaIEta[nEgs]/F' )
    tSig.Branch( 'etaWidth', etaWidthSig, 'etaWidth[nEgs]/F' )
    tSig.Branch( 'phiWidth', phiWidthSig, 'phiWidth[nEgs]/F' )
    tSig.Branch( 's4', s4Sig, 's4[nEgs]/F' )
    tSig.Branch( 'trkIsoPho', trkIsoPhoSig, 'trkIsoPho[nEgs]/F' )
    tSig.Branch( 'eta', etaSig, 'eta[nEgs]/F' )
    tSig.Branch( 'hOvrE', hOvrESig, 'hOvrE[nEgs]/F' )
    tSig.Branch( 'ecalPFIso', ecalPFIsoSig, 'ecalPFIso[nEgs]/F' )
    tSig.Branch( 'et', etSig, 'et[nEgs]/F' )
    tSig.Branch( 'energy', energySig, 'energy[nEgs]/F' )
    tSig.Branch( 'phi', phiSig, 'phi[nEgs]/F')
    tSig.Branch( 'mass', massSig, 'mass/F' )
    tSig.Branch( 'triggerBits', triggerBitsSig, 'triggerBits/I' )
    tSig.Branch( 'passFailStd', passFailStdSig, 'passFailStd/I' )
    tSig.Branch( 'passFailL1Single', passFailL1SingleSig, 'passFailL1Single/I' )
    tSig.Branch( 'passFailL1Double', passFailL1DoubleSig, 'passFailL1Double/I' )

    #Loop over all events
    photonIndexSig = 0

    #print varsSigTotal
    
    nTotalSig = 0
    nWithZeroSig = 0
    nWithOneSig = 0
    nWithTwoOrMoreSig = 0
    
    #First loop goes over all events
    for n in range(0,nEventsSig):
        nEventSig[0] = int(varsSigTotal[photonIndexSig][19])
        nEgsSig[0] = int(varsSigTotal[photonIndexSig][14])
        nEgsPassingSig[0] = int(varsSigTotal[photonIndexSig][20])

        nTotalSig += 1
        if nEgsSig[0] == 0:
            nWithZeroSig += 1
        elif nEgsSig[0] == 1:
            nWithOneSig += 1
        elif nEgsSig[0] > 1:
            nWithTwoOrMoreSig += 1
        
        massSig[0] = varsSigTotal[photonIndexSig][13]
        passFailStdSig[0] = int(varsSigTotal[photonIndexSig][15])
        passFailL1SingleSig[0] = int(varsSigTotal[photonIndexSig][16])
        passFailL1DoubleSig[0] = int(varsSigTotal[photonIndexSig][17])
        triggerBitsSig[0] = int(varsSigTotal[photonIndexSig][18])
                
        #Second event loops over all photons in said event
        for m in range(0,nEgsSig[0]):
            xgbScoreSig[m] = yPredSigTotal[photonIndexSig][0]
            rawEnergySig[m] = varsSigTotal[photonIndexSig][0]
            r9HLTSig[m] = varsSigTotal[photonIndexSig][1]
            sigmaIEtaIEtaSig[m] = varsSigTotal[photonIndexSig][2]
            etaWidthSig[m] = varsSigTotal[photonIndexSig][3]
            phiWidthSig[m] = varsSigTotal[photonIndexSig][4]
            s4Sig[m] = varsSigTotal[photonIndexSig][5]
            trkIsoPhoSig[m] = varsSigTotal[photonIndexSig][6]
            etaSig[m] = varsSigTotal[photonIndexSig][7]
            hOvrESig[m] = varsSigTotal[photonIndexSig][8]
            ecalPFIsoSig[m] = varsSigTotal[photonIndexSig][9]
            etSig[m] = varsSigTotal[photonIndexSig][10]
            energySig[m] = varsSigTotal[photonIndexSig][11]
            phiSig[m] = varsSigTotal[photonIndexSig][12]
            photonIndexSig += 1

        tSig.Fill()

#    print nTotalSig
#    print nWithZeroSig
#    print nWithOneSig
#    print nWithTwoOrMoreSig

##NOW DO BACKGROUND
    nEventBkg = array('i',[0])
    nEgsBkg = array('i',[0])
    nEgsPassingBkg = array('i',[0])
    xgbScoreBkg = array('f',[0.0]*maxPhotons)
    rawEnergyBkg = array('f',[0.0]*maxPhotons)
    r9HLTBkg = array('f',[0.0]*maxPhotons)
    sigmaIEtaIEtaBkg = array('f',[0.0]*maxPhotons)
    etaWidthBkg = array('f',[0.0]*maxPhotons)
    phiWidthBkg = array('f',[0.0]*maxPhotons)
    s4Bkg = array('f',[0.0]*maxPhotons)
    trkIsoPhoBkg = array('f',[0.0]*maxPhotons)
    etaBkg = array('f',[0.0]*maxPhotons)
    hOvrEBkg = array('f',[0.0]*maxPhotons)
    ecalPFIsoBkg = array('f',[0.0]*maxPhotons)
    etBkg = array('f',[0.0]*maxPhotons)
    energyBkg = array('f',[0.0]*maxPhotons)
    phiBkg = array('f',[0.0]*maxPhotons)
    massBkg = array('f',[0.0])

    triggerBitsBkg = array('i',[0])
    passFailStdBkg = array('i',[0])
    passFailL1SingleBkg = array('i',[0])
    passFailL1DoubleBkg = array('i',[0])

    tBkg = TTree("bkgTree", "bkgTree")
    tBkg.Branch( 'nEvent', nEventBkg, 'nEvent/I' )
    tBkg.Branch( 'nEgs', nEgsBkg, 'nEgs/I' )
    tBkg.Branch( 'nEgsPassing', nEgsPassingBkg, 'nEgsPassing/I' )
    tBkg.Branch( 'xgbScore', xgbScoreBkg, 'xgbScore[nEgs]/F' )
    tBkg.Branch( 'rawEnergy', rawEnergyBkg, 'rawEnergy[nEgs]/F' )
    tBkg.Branch( 'r9HLT', r9HLTBkg, 'r9HLT[nEgs]/F' )
    tBkg.Branch( 'sigmaIEtaIEta', sigmaIEtaIEtaBkg, 'sigmaIEtaIEta[nEgs]/F' )
    tBkg.Branch( 'etaWidth', etaWidthBkg, 'etaWidth[nEgs]/F' )
    tBkg.Branch( 'phiWidth', phiWidthBkg, 'phiWidth[nEgs]/F' )
    tBkg.Branch( 's4', s4Bkg, 's4[nEgs]/F' )
    tBkg.Branch( 'trkIsoPho', trkIsoPhoBkg, 'trkIsoPho[nEgs]/F' )
    tBkg.Branch( 'eta', etaBkg, 'eta[nEgs]/F' )
    tBkg.Branch( 'hOvrE', hOvrEBkg, 'hOvrE[nEgs]/F' )
    tBkg.Branch( 'ecalPFIso', ecalPFIsoBkg, 'ecalPFIso[nEgs]/F' )
    tBkg.Branch( 'et', etBkg, 'et[nEgs]/F' )
    tBkg.Branch( 'energy', energyBkg, 'energy[nEgs]/F' )
    tBkg.Branch( 'phi', phiBkg, 'phi[nEgs]/F')
    tBkg.Branch( 'mass', massBkg, 'mass/F' )
    tBkg.Branch( 'triggerBits', triggerBitsBkg, 'triggerBits/I' )
    tBkg.Branch( 'passFailStd', passFailStdBkg, 'passFailStd/I' )
    tBkg.Branch( 'passFailL1Single', passFailL1SingleBkg, 'passFailL1Single/I' )
    tBkg.Branch( 'passFailL1Double', passFailL1DoubleBkg, 'passFailL1Double/I' )


    nTotalBkg = 0
    nWithZeroBkg = 0
    nWithOneBkg = 0
    nWithTwoOrMoreBkg = 0

    #Loop over all events
    photonIndexBkg = 0
    #First loop goes over all events
    for n in range(0,nEventsBkg):
        nEventBkg[0] = int(varsBkgTotal[photonIndexBkg][19])
        nEgsBkg[0] = int(varsBkgTotal[photonIndexBkg][14])
        nEgsPassingBkg[0] = int(varsBkgTotal[photonIndexBkg][20])
        nTotalBkg += 1
        if nEgsBkg[0] == 0:
            nWithZeroBkg += 1
        elif nEgsBkg[0] == 1:
            nWithOneBkg += 1
        elif nEgsBkg[0] > 1:
            nWithTwoOrMoreBkg += 1
        massBkg[0] = varsBkgTotal[photonIndexBkg][13]
        passFailStdBkg[0] = int(varsBkgTotal[photonIndexBkg][15])
        passFailL1SingleBkg[0] = int(varsBkgTotal[photonIndexBkg][16])
        passFailL1DoubleBkg[0] = int(varsBkgTotal[photonIndexBkg][17])
        triggerBitsBkg[0] = int(varsBkgTotal[photonIndexBkg][18])
                
        #Second event loops over all photons in said event
        for m in range(0,nEgsBkg[0]):
            xgbScoreBkg[m] = yPredBkgTotal[photonIndexBkg][0]
            rawEnergyBkg[m] = varsBkgTotal[photonIndexBkg][0]
            r9HLTBkg[m] = varsBkgTotal[photonIndexBkg][1]
            sigmaIEtaIEtaBkg[m] = varsBkgTotal[photonIndexBkg][2]
            etaWidthBkg[m] = varsBkgTotal[photonIndexBkg][3]
            phiWidthBkg[m] = varsBkgTotal[photonIndexBkg][4]
            s4Bkg[m] = varsBkgTotal[photonIndexBkg][5]
            trkIsoPhoBkg[m] = varsBkgTotal[photonIndexBkg][6]
            etaBkg[m] = varsBkgTotal[photonIndexBkg][7]
            hOvrEBkg[m] = varsBkgTotal[photonIndexBkg][8]
            ecalPFIsoBkg[m] = varsBkgTotal[photonIndexBkg][9]
            etBkg[m] = varsBkgTotal[photonIndexBkg][10]
            energyBkg[m] = varsBkgTotal[photonIndexBkg][11]
            phiBkg[m] = varsBkgTotal[photonIndexBkg][12]
            photonIndexBkg += 1
                        
        tBkg.Fill()
    
#    print nTotalBkg
#    print nWithZeroBkg
#    print nWithOneBkg
#    print nWithTwoOrMoreBkg

    outFileRoot.Write()
    outFileRoot.Close()
