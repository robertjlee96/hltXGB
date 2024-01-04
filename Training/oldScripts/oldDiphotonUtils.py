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
    inputValuesSigLead = []
    targetValuesSigLead = []
    inputValuesSigSub = []
    targetValuesSigSub = []

    inputValuesBkgLead = []
    targetValuesBkgLead = []
    inputValuesBkgSub = []
    targetValuesBkgSub = []

    # names of variables used as BDT input
    inputVarNames = []

    # original weights without pt/eta reweighting
    # we can use these also for evaluation
    origWeightsSigLead = []
    origWeightsSigSub = []
    origWeightsBkgLead = []
    origWeightsBkgSub = []

    varValuesSigLead = []
    varValuesSigSub = []
    varValuesBkgLead = []
    varValuesBkgSub = []
    varNames = []

    isFirstVar = True

    for varname in inputVars + extraVars + singleVars + ['index']:

        thisValuesSigLead = []
        thisValuesSigSub = []

        thisValuesBkgLead = []
        thisValuesBkgSub = []

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
                
           
            if varname not in singleVars and varname != 'index':
                if label == 0:
                    varValuesBkgLead.append([item[0] for item in tree.array(varname)[indices]])
                    varValuesBkgSub.append([item[1] for item in tree.array(varname)[indices]])
                if label == 1:
                    varValuesSigLead.append([item[0] for item in tree.array(varname)[indices]])
                    varValuesSigSub.append([item[1] for item in tree.array(varname)[indices]])

                if isFirstProc:
                    varNames.append(varname)
#
                if varname in inputVars:
                    # BDT input variable
                    if label == 0:
                        inputValuesBkgLead.append([item[0] for item in tree.array(varname)[indices]])
                        inputValuesBkgSub.append([item[1] for item in tree.array(varname)[indices]])
                    if label == 1:
                        inputValuesSigLead.append([item[0] for item in tree.array(varname)[indices]])
                        inputValuesSigSub.append([item[1] for item in tree.array(varname)[indices]])

                    if isFirstProc:
                        inputVarNames.append(varname)

            elif varname in singleVars:
                if label == 0:
                    varValuesBkgLead.append(tree.array(varname)[indices])
                    varValuesBkgSub.append(tree.array(varname)[indices])
                if label == 1:
                    varValuesSigLead.append(tree.array(varname)[indices])
                    varValuesSigSub.append(tree.array(varname)[indices])
                        
            elif varname == 'index':
                if label == 0:
                    varValuesBkgLead.append(range(0,len(tree.array('eta')[indices])))
                    varValuesBkgSub.append(range(0,len(tree.array('eta')[indices])))
                if label == 1:
                    varValuesSigLead.append(range(0,len(tree.array('eta')[indices])))
                    varValuesSigSub.append(range(0,len(tree.array('eta')[indices])))
                

            # append target values and weights
            if isFirstVar:
#                #thisWeights =  tree.array('weight')[indices]
                thisWeights =  np.ones(len(mask))[indices]
                if label == 0:
                    targetValuesBkgLead.append(np.ones(len(inputValuesBkgLead[-1])) * 1)
                    origWeightsBkgLead.append(thisWeights)
                    targetValuesBkgSub.append(np.ones(len(inputValuesBkgSub[-1])) * 1)
                    origWeightsBkgSub.append(thisWeights)
                    
                if label == 1:
                    targetValuesSigLead.append(np.ones(len(inputValuesSigLead[-1])) * 1)
                    origWeightsSigLead.append(thisWeights)
                    targetValuesSigSub.append(np.ones(len(inputValuesSigSub[-1])) * 1)
                    origWeightsSigSub.append(thisWeights)
            isFirstProc = False

        if isFirstVar:
            if label == 0:
                targetValuesBkgLead = np.hstack(targetValuesBkgLead)
                origWeightsBkgLead = np.hstack(origWeightsBkgLead)
                targetValuesBkgSub = np.hstack(targetValuesBkgSub)
                origWeightsBkgSub = np.hstack(origWeightsBkgSub)
            if label == 1:
                targetValuesSigLead = np.hstack(targetValuesSigLead)
                origWeightsSigLead = np.hstack(origWeightsSigLead)
                targetValuesSigSub = np.hstack(targetValuesSigSub)
                origWeightsSigSub = np.hstack(origWeightsSigSub)

        isFirstVar = False

    inputValuesSigLead = np.vstack(inputValuesSigLead)
    varValuesSigLead = np.vstack(varValuesSigLead)
    inputValuesSigSub = np.vstack(inputValuesSigSub)
    varValuesSigSub = np.vstack(varValuesSigSub)
    
    inputValuesBkgLead = np.vstack(inputValuesBkgLead)
    varValuesBkgLead = np.vstack(varValuesBkgLead)
    inputValuesBkgSub = np.vstack(inputValuesBkgSub)
    varValuesBkgSub = np.vstack(varValuesBkgSub)


#    inputValuesSigSub = np.array(np.vstack(inputValuesSigSub[0]).T)
#    varValuesSigSub = np.array(np.vstack(varValuesSigSub[0]).T)
#
#    inputValuesBkgLead = np.array(np.vstack(inputValuesBkgLead[0]).T)
#    varValuesBkgLead = np.array(np.vstack(varValuesBkgLead[0]).T)
#    inputValuesBkgSub = np.array(np.vstack(inputValuesBkgSub[0]).T)
#    varValuesBkgSub = np.array(np.vstack(varValuesBkgSub[0]).T)

    varNames = np.hstack(varNames)

    return inputValuesSigLead, inputValuesSigSub, inputValuesBkgLead, inputValuesBkgSub, targetValuesSigLead, targetValuesSigSub, targetValuesBkgLead,targetValuesBkgSub, origWeightsSigLead, origWeightsSigSub,origWeightsBkgLead,origWeightsBkgSub, varValuesSigLead, varValuesSigSub, varValuesBkgLead, varValuesBkgSub, inputVarNames, varNames
