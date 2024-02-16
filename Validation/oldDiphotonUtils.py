#!/usr/bin/env python

import numpy as np
import uproot
from ROOT import TFile, TTree, TChain
import array

barrelVars = [
    'SCRawE',
       'r9',
        'sigmaIetaIeta',
    'etaWidth',
    'phiWidth',
    'covIEtaIPhi',
    's4',
    'phoIso03',
    'chgIsoWrtChosenVtx',
    'chgIsoWrtWorstVtx',
    'scEta',
    'rho',
    'hadTowOverEm',
    'hadronicOverEm'
]
endcapVars = list(barrelVars) + [
    'esEffSigmaRR',
    'esEnergyOverRawE',

]

extraVars = ['pt','hggMass','trkSumPtHollowConeDR03','passSPresel','passDPresel','scEtaSecond','index']

#----------------------------------------------------------------------

#def load_file(input_file, geoSelection = None, ptCuts = None, dptCuts = None):
def load_file(inputFile, geoSelection = None, ptCuts = None):
    """input_file should be a uproot object corresponding to a ROOT file

    :return: input_values, target_values, orig_weights, train_weights, pt, scEta, input_var_names
    """
    inputValuesSigLead = []
    targetValuesSigLead = []
    inputValuesSigSub = []
    targetValuesSigSub = []
    
    inputValuesPrompt = []
    targetValuesPrompt = []
    
    inputValuesFake = []
    targetValuesFake = []

    # names of variables used as BDT input
    inputVarNames = []

    # original weights without pt/eta reweighting
    # we can use these also for evaluation
    origWeightsSigLead = []
    origWeightsSigSub = []
    origWeightsPrompt = []
    origWeightsFake = []

    varValuesSigLead = []
    varValuesSigSub = []
    varValuesPrompt = []
    varValuesFake = []
    varNames = []

    isFirstVar = True

#    for varname in barrelVars + ['esEffSigmaRR','esEnergyOverRawE'] + extraVars: #With Pt/M In model
    for varname in barrelVars + ['ptOvrM','esEffSigmaRR','esEnergyOverRawE'] + extraVars: #Without Pt/M in Model

        thisValuesSigLead = []
        thisValuesSigSub = []
        
        thisValuesPrompt = []
        thisValuesFake = []
        
        isFirstProc = True

        for treeName, label in [
            ('sigPhotonsLead', 0),
            ('sigPhotonsSub', 1),
            ('promptPhotons', 2),
            ('fakePhotons', 3)
        ]:

            tree = inputFile[treeName]
            
            geoMask = geoSelection(tree)
            ptMask = ptCuts(tree)

            mask = []
#            for i in range(0,len(geoMask)):
#                if geoMask[i] == True and ptMask[i] == True:
#                    mask.append(True)
#                else:
#                    mask.append(False)
            for i in range(0,len(geoMask)):
                    mask.append(True)

            if not mask is None:
                indices = mask
              
            else:
                indices = np.ones(len(tree.array(varname)), dtype = 'bool')
            if varname != 'index':
                if label == 0:
                    varValuesSigLead.append(tree.array(varname)[indices])
                if label == 1:
                    varValuesSigSub.append(tree.array(varname)[indices])
                if label == 2:
                    varValuesPrompt.append(tree.array(varname)[indices])
                if label == 3:
                    varValuesFake.append(tree.array(varname)[indices])
                
                if isFirstProc:
                    varNames.append(varname)

                if varname in endcapVars:
                    # BDT input variable
                    if label == 0:
                        inputValuesSigLead.append(tree.array(varname)[indices])
                    if label == 1:
                        inputValuesSigSub.append(tree.array(varname)[indices])
                    if label == 2:
                        inputValuesPrompt.append(tree.array(varname)[indices])
                    if label == 3:
                        inputValuesFake.append(tree.array(varname)[indices])
    #                if label == 0:
    #                    thisValuesSigLead.append(tree.array(varname)[indices])
    #                if label == 1:
    #                    thisValuesSigSub.append(tree.array(varname)[indices])
    #                if label == 2:
    #                    thisValuesPrompt.append(tree.array(varname)[indices])
    #                if label == 3:
    #                    thisValuesFake.append(tree.array(varname)[indices])
                    
                    if isFirstProc:
                        inputVarNames.append(varname)
            else:
                    if label == 0:
                        varValuesSigLead.append(range(0,len(tree.array('event')[indices])))
                    if label == 1:
                        varValuesSigSub.append(range(0,len(tree.array('event')[indices])))
                    if label == 2:
                        varValuesPrompt.append(range(0,len(tree.array('event')[indices])))
                    if label == 3:
                        varValuesFake.append(range(0,len(tree.array('event')[indices])))

            # append target values and weights
            if isFirstVar:
                thisWeights =  tree.array('weight')[indices]
                if label == 0:
                    targetValuesSigLead.append(np.ones(len(inputValuesSigLead[-1])) * 1)
                    origWeightsSigLead.append(thisWeights)
                if label == 1:
                    targetValuesSigSub.append(np.ones(len(inputValuesSigSub[-1])) * 1)
                    origWeightsSigSub.append(thisWeights)
                if label == 2:
                    targetValuesPrompt.append(np.ones(len(inputValuesPrompt[-1])) * 1)
                    origWeightsPrompt.append(thisWeights)
                if label == 3:
                    targetValuesFake.append(np.zeros(len(inputValuesFake[-1])))
                    origWeightsFake.append(thisWeights)

            isFirstProc = False

#            if thisValuesSigLead and label == 0:
#                print "len(thisValuesSigLead)",len(thisValuesSigLead)
#                inputValuesSigLead.append(np.hstack(thisValuesSigLead))
#            if label == 1 and thisValuesSigSub:
#                inputValuesSigSub.append(np.hstack(thisValuesSigSub))
#            if label == 2 and thisValuesPrompt:
#                inputValuesPrompt.append(np.hstack(thisValuesPrompt))
#            if label == 3 and thisValuesFake:
#                inputValuesFake.append(np.hstack(thisValuesFake))

        if isFirstVar:
            if label == 0:
                targetValuesSigLead = np.hstack(targetValuesSigLead)
                origWeightsSigLead = np.hstack(origWeightsSigLead)
            if label == 1:
                targetValuesSigSub = np.hstack(targetValuesSigSub)
                origWeightsSigSub = np.hstack(origWeightsSigSub)
            if label == 2:
                targetValuesPrompt = np.hstack(targetValuesPrompt)
                origWeightsPrompt = np.hstack(origWeightsPrompt)
            if label == 3:
                targetValuesFake = np.hstack(targetValuesFake)
                origWeightsFake = np.hstack(origWeightsFake)

        isFirstVar = False

    
    inputValuesSigLead = np.array(np.vstack(inputValuesSigLead).T)
    varValuesSigLead = np.array(np.vstack(varValuesSigLead).T)
    
    inputValuesSigSub = np.array(np.vstack(inputValuesSigSub).T)
    varValuesSigSub = np.array(np.vstack(varValuesSigSub).T)
    
    inputValuesPrompt = np.array(np.vstack(inputValuesPrompt).T)
    varValuesPrompt = np.array(np.vstack(varValuesPrompt).T)
    
    inputValuesFake = np.array(np.vstack(inputValuesFake).T)
    varValuesFake = np.array(np.vstack(varValuesFake).T)
   
    varNames = np.hstack(varNames)


    return inputValuesSigLead, inputValuesSigSub, inputValuesPrompt, inputValuesFake, targetValuesSigLead, targetValuesSigSub, targetValuesPrompt,targetValuesFake, origWeightsSigLead, origWeightsSigSub,origWeightsPrompt,origWeightsFake, varValuesSigLead, varValuesSigSub, varValuesPrompt, varValuesFake, inputVarNames, varNames
