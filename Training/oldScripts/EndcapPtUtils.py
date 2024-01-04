#!/usr/bin/env python

import numpy as np
import uproot
from ROOT import TFile, TTree, TChain
import array

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

    for varname in endcap_vars + ['et']:

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

#def saveScores(inputFile, outputFileName, xgbScores, geoSelectionLead = None, geoSelectionSub = None, geoSelectionBkg = None):
def saveScores(inputFile, inputFileName, outputFileName, xgbScores, geoSelection = None, ptCuts = None):

    inFilepyRoot = TFile.Open(inputFileName)
    outputFile  = TFile(outputFileName, 'recreate')
   
   
    #treesIn = [inFilepyRoot.Get('ggh_125_13TeV_GluGluHToGG_M125_$SYST'),inFilepyRoot.Get('fakePhotons')]
    #treesOut = [treesIn[0].CloneTree(0),treesIn[1].CloneTree(0)]
   
    treeIn1 = inFilepyRoot.Get('promptPhotons')
    treeOut1 = treeIn1.CloneTree(0)
    treeOut1.SetName("promptPhotons")
    
    treeName1 = 'promptPhotons'
    tree1 = inputFile[treeName1]
    #print "number of events in {} = {}".format(treeName1,treeIn1.GetEntries())
   
    geoMaskPrompt = geoSelection(tree1)
    ptMaskPrompt = ptCuts(tree1)
    
    maskPrompt = []
    for i in range(0,len(geoMaskPrompt)):
        if geoMaskPrompt[i] == True and ptMaskPrompt[i] == True:
            maskPrompt.append(True)
        else:
            maskPrompt.append(False)
    
    xgbValPrompt = array.array("d", [0.])
    xgbScorePrompt = treeOut1.Branch("xgbScore",xgbValPrompt,"xgbScore/D")
    
    lenGeoMaskTrue = 0
    for maskVal in maskPrompt:
        if maskVal == True:
            lenGeoMaskTrue += 1
    j = 0
    for i in range (0,treeIn1.GetEntries()):
        #if geoMaskPrompt[i] == True:
        if maskPrompt[i] == True:
            #treeRoot.GetEntry(i)
            treeIn1.GetEntry(i)
            xgbValPrompt[0] = float(xgbScores[0][j])
            treeOut1.Fill()
            #outTree.Fill()
            j += 1

    #treeOut1.Write()
    
    treeIn2 = inFilepyRoot.Get('fakePhotons')
    treeOut2 = treeIn2.CloneTree(0)
    
    treeName2 = 'fakePhotons'
    tree2 = inputFile[treeName2]
   # print "number of events in {} = {}".format(treeName2,treeIn2.GetEntries())

    geoMaskFake = geoSelection(tree2)
    ptMaskFake = ptCuts(tree2)
    
    maskFake = []
    for i in range(0,len(geoMaskFake)):
        if geoMaskFake[i] == True and ptMaskFake[i] == True:
            maskFake.append(True)
        else:
            maskFake.append(False)
    

    xgbValFake = array.array("d", [0.])
    xgbScoreFake = treeOut2.Branch("xgbScore",xgbValFake,"xgbScore/D")
    #xgbBkg = outTree.Branch("xgbScoreBkg",xgbBkgVal,"xgbScoreBkg/D")
    j = 0
    #for i in range (0,treeRoot.GetEntries()):
    for i in range (0,treeIn2.GetEntries()):
       # if geoMaskFake[i] == True:
        if maskFake[i] == True:
            xgbValFake[0] = float(xgbScores[1][j])
            #treeRoot.GetEntry(i)
            treeIn2.GetEntry(i)
            treeOut2.Fill()
            #outTree.Fill()
            j += 1
    #treeOut2.Write()
   
    #outTree.Write()
    
    #treesOut[sample].Write()
    #outTree.Delete()
    outputFile.Write()
    #soutputFile.Close()

    return
