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
barrel_vars = [
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
endcap_vars = list(barrel_vars)


#----------------------------------------------------------------------

#def load_file(input_file, geoSelection = None, ptCuts = None, dptCuts = None):
def load_file(input_file, geoSelection = None, ptCuts = None, massCut = None):
    """input_file should be a uproot object corresponding to a ROOT file

    :return: input_values, target_values, orig_weights, train_weights, pt, scEta, input_var_names
    """
    input_values = []

    target_values = []

    # names of variables used as BDT input
    input_var_names = []

    # original weights without pt/eta reweighting
    # we can use these also for evaluation
    orig_weights = []

    # weights which are pt/eta reweighted for prompt photons
    # and equal to orig_weights for fake photons
    train_weights = []

    pt_values = []
    scEta_values = []

    is_first_var = True

    for varname in barrel_vars + ['et']:

        this_values = []

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
                    if geoMask[i] == True:
                        mask.append(True)
                    else:
                        mask.append(False)

                if not mask is None:
                    indices = mask

            if label == 1:
                for i in range(0,len(geoMask)):
                    if geoMask[i] == True and massMask[i] == True:
                        mask.append(True)
                    else:
                        mask.append(False)

                if not mask is None:
                    indices = mask

#            if not geoSelection is None:
#                indices = geoSelection(tree)
              
            #else:
                #indices = np.ones(len(tree.array(varname)), dtype = 'bool')
                
            if varname == 'et':
                pt_values.append(tree.array(varname)[indices])
            elif varname == 'eta':
                scEta_values.append(tree.array(varname)[indices])
                this_values.append(tree.array(varname)[indices])
                if is_first_proc:
                    input_var_names.append(varname)
            else:
                # BDT input variable
                this_values.append(tree.array(varname)[indices])

                if is_first_proc:
                    input_var_names.append(varname)

            # append target values and weights
            if is_first_var:
                target_values.append(np.ones(len(this_values[-1])) * label)
                this_weights =  np.ones(len(mask))[indices]

#                this_weights =  tree.array('weight')[indices]

                orig_weights.append(this_weights)
#
#                if label == 1:
#                    # eta/pt reweighting is only for signal
#                    this_weights = this_weights * tree.array('PtvsEtaWeight')[indices]
#
#                train_weights.append(this_weights)
                train_weights.append(this_weights)

            is_first_proc = False

        # end of loop over processes

        if this_values:
            input_values.append(np.hstack(this_values))

        if is_first_var:
            target_values = np.hstack(target_values)
            orig_weights = np.hstack(orig_weights)
            train_weights = np.hstack(train_weights)

        is_first_var = False


    input_values = np.vstack(input_values).T
    pt_values = np.hstack(pt_values)
    scEta_values = np.hstack(scEta_values)

    return input_values, target_values, orig_weights, train_weights, pt_values, scEta_values, input_var_names
