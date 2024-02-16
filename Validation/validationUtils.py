#!/usr/bin/env python

import numpy as np
import uproot
from ROOT import TFile, TTree, TChain
import array

#barrel_vars = [
#    'SCRawE',
#       'r9',
#        'sigmaIetaIeta',
#    'etaWidth',
#    'phiWidth',
#    'covIEtaIPhi',
#    's4',
#    'phoIso03',
#    'scEta',
#    'rho',
#    'hadronicOverEm'
#]
#endcap_vars = list(barrel_vars) + [
#    'esEffSigmaRR',
#    'esEnergyOverRawE',
#
#]
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


extraVars = ['nEgs','et','energy','mass','passFailStd','passFailL1Single','passFailL1Double']


#----------------------------------------------------------------------

#def load_file(input_file, geoSelection = None, ptCuts = None, dptCuts = None):
def load_file(input_file, geoSelection = None, ptCuts = None, massCut = None):
    """input_file should be a uproot object corresponding to a ROOT file

    :return: input_values, target_values, orig_weights, train_weights, pt, scEta, input_var_names
    """
    inputValuesSig = []
    inputValuesBkg = []

    targetValuesSig = []
    targetValuesBkg = []

    # names of variables used as BDT input
    inputVarNames = []

    # original weights without pt/eta reweighting
    # we can use these also for evaluation
    origWeightsSig = []
    origWeightsBkg = []

    # weights which are pt/eta reweighted for prompt photons
    # and equal to orig_weights for fake photons

    varValuesSig = []
    varValuesBkg = []
    varNames = []

    is_first_var = True

    for varname in inputVars + extraVars:

        thisValuesSig = []
        thisValuesBkg = []

        is_first_proc = True

        for tree_name, label in [
            ('sigTree', 1),
            ('bkgTree', 0)
        ]:

            tree = input_file[tree_name]

            
            geoMask = geoSelection(tree)
            ptMask = ptCuts(tree)
            massMask = massCut(tree)

            mask = []
            if label == 0:
                for i in range(0,len(geoMask)):
                    if geoMask[i] == True:# and massMask[i] == True:
                        mask.append(True)
                    else:
                        mask.append(False)

                if not mask is None:
                    indices = mask

            if label == 1:
                for i in range(0,len(geoMask)):
                    if geoMask[i] == True:# and massMask[i] == True:
                        mask.append(True)
                    else:
                        mask.append(False)

                if not mask is None:
                    indices = mask

#            if not geoSelection is None:
#                indices = geoSelection(tree)
              
            #else:
                #indices = np.ones(len(tree.array(varname)), dtype = 'bool')
            if varname != 'dummy':
                if label == 0:
                    varValuesBkg.append(tree.array(varname)[indices])
                if label == 1:
                    varValuesSig.append(tree.array(varname)[indices])

                if is_first_proc:
                    varNames.append(varname)
                    
                if varname in inputVars:
                    if label == 0:
                        inputValuesBkg.append(tree.array(varname)[indices])
                    if label == 1:
                        inputValuesSig.append(tree.array(varname)[indices])
                    
                    if is_first_proc:
                        inputVarNames.append(varname)
                
            else:
                this_values.append(tree.array(varname)[indices])
                if is_first_proc:
                    input_var_names.append(varname)

            # append target values and weights
            if is_first_var:
                this_weights =  np.ones(len(mask))[indices]
                if label == 0:
                    targetValuesBkg.append(np.ones(len(inputValuesBkg[-1])) * label)
                    origWeightsBkg.append(this_weights)
                if label == 1:
                    targetValuesSig.append(np.ones(len(inputValuesSig[-1])) * label)
                    origWeightsSig.append(this_weights)
                    
            is_first_proc = False

        # end of loop over processes

        if is_first_var:
            if label == 0:
                targetValuesBkg = np.hstack(targetValuesBkg)
                origWeightsBkg = np.hstack(origWeightsBkg)
            if label == 1:
                targetValuesSig = np.hstack(targetValuesSig)
                origWeightsSig = np.hstack(origWeightsSig)

        is_first_var = False


    inputValuesSig = np.vstack(inputValuesSig).T
    varValuesSig = np.vstack(varValuesSig).T
    
    inputValuesBkg = np.vstack(inputValuesBkg).T
    varValuesBkg = np.vstack(varValuesBkg).T
    
    varNames = np.hstack(varNames)

    return inputValuesSig, inputValuesBkg, targetValuesSig, targetValuesBkg, origWeightsSig, origWeightsBkg, varValuesSig, varValuesBkg, inputVarNames, varNames
