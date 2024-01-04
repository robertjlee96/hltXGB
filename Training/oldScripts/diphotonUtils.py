#!/usr/bin/env python

import numpy as np
import uproot
from ROOT import TFile, TTree, TChain
import array

inputVars = [
    'rawEnergy',
    'r9HLT',
    'sigmaIEtaIEta',
    'etaWidth',
    'phiWidth',
    's4',
    'trkIsoPho',
    'eta',
    'hOvrE',
    'ecalPFIso'
]

extraVars = ['et','energy']
singleVars = ['mass','nEgs','passFailStd','passFailL1Single','passFailL1Double','triggerBits']

#----------------------------------------------------------------------

#def load_file(input_file, geoSelection = None, ptCuts = None, dptCuts = None):
def load_file(inputFile, geoSelection = None, ptCuts = None, massCuts = None):


    """input_file should be a uproot object corresponding to a ROOT file

    :return: input_values, target_values, orig_weights, train_weights, pt, scEta, input_var_names
    """
    inputValuesSig = []
    inputValuesSigLead = []
    inputValuesSigSub = []
    
    targetValuesSig = []

    inputValuesBkg = []
    inputValuesBkgLead = []
    inputValuesBkgSub = []
    
    targetValuesBkg = []


    # names of variables used as BDT input
    inputVarNames = []

    # original weights without pt/eta reweighting
    # we can use these also for evaluation
    origWeightsSig = []
    origWeightsBkg = []

    varValuesSig = []
    varValuesSigLead = []
    varValuesSigSub = []
    
    varValuesBkg = []
    varValuesBkgLead = []
    varValuesBkgSub = []
    
    varNames = []

    isFirstVar = True

    for varname in inputVars + extraVars + singleVars + ['index']:

        thisValuesSig = []
        thisValuesBkg = []

        isFirstProc = True

        for treeName, label in [
            ('sigTree', 1),
            ('bkgTree', 0)
        ]:

            tree = inputFile[treeName]
            
            #mask = []
            
            mask = massCuts(tree)
            #print mask
            
#            for i in range(0,len(tree.array('eta'))):
#                if tree.array('mass')[i][0] > 60:
#                    mask.append(True)
#                else:
#                    mask.append(False)

            if not mask is None:
                indices = mask
            else:
                indices = np.ones(len(tree.array(varname)), dtype = 'bool')
                
           
            if varname != 'index' and varname not in singleVars:
                if label == 0:
                    #varValuesBkg.append([(item[0],item[1]) for item in tree.array(varname)[indices]])
                    varValuesBkgLead.append([item[0] for item in tree.array(varname)[indices]])
                    varValuesBkgSub.append([item[1] for item in tree.array(varname)[indices]])
                if label == 1:
                    #varValuesSig.append([(item[0],item[1]) for item in tree.array(varname)[indices]])
                    varValuesSigLead.append([item[0] for item in tree.array(varname)[indices]])
                    varValuesSigSub.append([item[1] for item in tree.array(varname)[indices]])


                if isFirstProc:
                    varNames.append(varname)
#
                if varname in inputVars:
                    # BDT input variable
                    if label == 0:
                        #inputValuesBkg.append([item[0],item[1] for item in tree.array(varname)[indices]])
                        inputValuesBkgLead.append([item[0] for item in tree.array(varname)[indices]])
                        inputValuesBkgSub.append([item[1] for item in tree.array(varname)[indices]])
                    if label == 1:
                        #inputValuesSig.append([item[0],item[1] for item in tree.array(varname)[indices]])
                        inputValuesSigLead.append([item[0] for item in tree.array(varname)[indices]])
                        inputValuesSigSub.append([item[1] for item in tree.array(varname)[indices]])

                    if isFirstProc:
                        inputVarNames.append(varname)
                        
            elif varname in singleVars:
                if label == 0:
                    #varValuesBkg.append([(item[0],item[1]) for item in tree.array(varname)[indices]])
                    varValuesBkgLead.append([item for item in tree.array(varname)[indices]])
                    varValuesBkgSub.append([item for item in tree.array(varname)[indices]])
                if label == 1:
                    #varValuesSig.append([(item[0],item[1]) for item in tree.array(varname)[indices]])
                    varValuesSigLead.append([item for item in tree.array(varname)[indices]])
                    varValuesSigSub.append([item for item in tree.array(varname)[indices]])
                        
            elif varname == 'index':
                if label == 0:
                    indicesBkg = range(0,len(tree.array('eta')[indices]))
                    #indicesArray = [val for pair in zip(indicesBkg, indicesBkg) for val in pair]
                    varValuesBkgLead.append(indicesBkg)
                    varValuesBkgSub.append(indicesBkg)
                if label == 1:
                    indicesSig = range(0,len(tree.array('eta')[indices]))
                    #indicesArray = [val for pair in zip(indicesSig, indicesSig) for val in pair]
                    varValuesSigLead.append(indicesSig)
                    varValuesSigSub.append(indicesSig)
                

            # append target values and weights
            if isFirstVar:
#                #thisWeights =  tree.array('weight')[indices]
                thisWeights =  np.ones(len(mask))[indices]
                if label == 0:
                    targetValuesBkg.append(np.ones(len(inputValuesBkgLead[-1])) * 1)
                    origWeightsBkg.append(thisWeights)
                    
                if label == 1:
                    targetValuesSig.append(np.ones(len(inputValuesSigLead[-1])) * 1)
                    origWeightsSig.append(thisWeights)
    
            isFirstProc = False

        if isFirstVar:
            if label == 0:
                targetValuesBkg = np.hstack(targetValuesBkg)
                origWeightsBkg = np.hstack(origWeightsBkg)
           
            if label == 1:
                targetValuesSig = np.hstack(targetValuesSig)
                origWeightsSig = np.hstack(origWeightsSig)

        isFirstVar = False



    inputValuesSigLead = np.vstack(inputValuesSigLead)
    varValuesSigLead = np.vstack(varValuesSigLead)

    inputValuesBkgLead = np.vstack(inputValuesBkgLead)
    varValuesBkgLead = np.vstack(varValuesBkgLead)


    inputValuesSig = np.concatenate((inputValuesSigLead,inputValuesSigSub), axis=1)
    inputValuesBkg = np.concatenate((inputValuesBkgLead,inputValuesBkgSub), axis=1)
    
    
    varValuesSig = np.concatenate((varValuesSigLead,varValuesSigSub), axis=1)
    varValuesBkg = np.concatenate((varValuesBkgLead,varValuesBkgSub), axis=1)
    

#    inputValuesSigSub = np.array(np.vstack(inputValuesSigSub[0]).T)
#    varValuesSigSub = np.array(np.vstack(varValuesSigSub[0]).T)
#
#    inputValuesBkgLead = np.array(np.vstack(inputValuesBkgLead[0]).T)
#    varValuesBkgLead = np.array(np.vstack(varValuesBkgLead[0]).T)
#    inputValuesBkgSub = np.array(np.vstack(inputValuesBkgSub[0]).T)
#    varValuesBkgSub = np.array(np.vstack(varValuesBkgSub[0]).T)

    varNames = np.hstack(varNames)

    return inputValuesSig, inputValuesBkg, targetValuesSig, targetValuesBkg, origWeightsSig, origWeightsBkg, varValuesSig, varValuesBkg, inputVarNames, varNames
