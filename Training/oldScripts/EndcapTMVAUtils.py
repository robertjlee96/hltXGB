import re
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.cElementTree as ET
regex_float_pattern = r'[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'

def build_tree(xgtree, base_xml_element, var_indices,var_list):
    parent_element_dict = {'0':base_xml_element}
    pos_dict = {'0':'s'}
    
    varListGenNames = ['f0','f1','f2','f3','f4','f5','f6','f7','f8','f9','f10','f11','f12','f13','f14','f15']
    
#    varListBothNames = {}
#
#    for i in range(0,len(varListGenNames)):
#        var = varListGenNames[i]
#        varListBothNames[var] = var_list
        
    for line in xgtree.split('\n'):
        if not line: continue
        if ':leaf=' in line:
            #leaf node
            result = re.match(r'(\t*)(\d+):leaf=({0})$'.format(regex_float_pattern), line)
            if not result:
                print(line)
            depth = result.group(1).count('\t')
            inode = result.group(2)
            res = result.group(3)
            node_elementTree = ET.SubElement(parent_element_dict[inode], "Node", pos=str(pos_dict[inode]),
                                             depth=str(depth), NCoef="0", IVar="-1", Cut="0.0e+00", cType="1", res=str(res), rms="0.0e+00", purity="0.0e+00", nType="-99")
        else:
            #\t\t3:[var_topcand_mass<138.19] yes=7,no=8,missing=7
            result = re.match(r'(\t*)([0-9]+):\[(?P<var>.+)<(?P<cut>{0})\]\syes=(?P<yes>\d+),no=(?P<no>\d+)'.format(regex_float_pattern),line)
            if not result:
                print(line)
            depth = result.group(1).count('\t')
            inode = result.group(2)
            var = result.group('var')
            cut = result.group('cut')
            lnode = result.group('yes')
            rnode = result.group('no')
            pos_dict[lnode] = 'l'
            pos_dict[rnode] = 'r'
#            indexTmp = var_indicies.index(var)
#            newName = var_indices[
#            print "var_indices = {}".format(var_indices)
#            print "var_list = {}".format(var_list)
#            print "var = {}".format(var)
            varIndex = varListGenNames.index(var)
#            print "varIndex = {}".format(varIndex)
            varCorr = var_list[varIndex][0]
#            print "varCorr = {}".format(varCorr)
            #print "indexTmp = {}".format(var_indices[var])
        
            node_elementTree = ET.SubElement(parent_element_dict[inode], "Node",
                pos=str(pos_dict[inode]), depth=str(depth),
                NCoef="0", IVar=str(var_indices[varCorr]),
                Cut=str(cut), cType="1", res="0.0e+00", rms="0.0e+00", purity="0.0e+00", nType="0")
            parent_element_dict[lnode] = node_elementTree
            parent_element_dict[rnode] = node_elementTree
            
def convert_model(model, input_variables, output_xml):
    NTrees = len(model)
    var_list = input_variables
#    print "var_list Within Convert Model = {}".format(var_list)
    var_indices = {}
    
    # <MethodSetup>
    MethodSetup = ET.Element("MethodSetup", Method="BDT::BDT")

    # <Variables>
    Variables = ET.SubElement(MethodSetup, "Variables", NVar=str(len(var_list)))
    for ind, val in enumerate(var_list):
#        print "ind, val = {}, {}".format(ind,val)
        name = val[0]
        var_type = val[1]
        var_indices[name] = ind
        Variable = ET.SubElement(Variables, "Variable", VarIndex=str(ind), Type=val[1],
            Expression=name, Label=name, Title=name, Unit="", Internal=name,
            Min="0.0e+00", Max="0.0e+00")
#    print "var_indices = {}".format(var_indices)

    # <GeneralInfo>
    GeneralInfo = ET.SubElement(MethodSetup, "GeneralInfo")
    Info_Creator = ET.SubElement(GeneralInfo, "Info", name="Creator", value="xgboost2TMVA")
    Info_AnalysisType = ET.SubElement(GeneralInfo, "Info", name="AnalysisType", value="Classification")

    # <Options>
    Options = ET.SubElement(MethodSetup, "Options")
    Option_NodePurityLimit = ET.SubElement(Options, "Option", name="NodePurityLimit", modified="No").text = "5.00e-01"
    Option_BoostType = ET.SubElement(Options, "Option", name="BoostType", modified="Yes").text = "Grad"
    
    # <Weights>
    Weights = ET.SubElement(MethodSetup, "Weights", NTrees=str(NTrees), AnalysisType="1")
    
    for itree in range(NTrees):
        BinaryTree = ET.SubElement(Weights, "BinaryTree", type="DecisionTree", boostWeight="1.0e+00", itree=str(itree))
        build_tree(model[itree], BinaryTree, var_indices,var_list)
        
    tree = ET.ElementTree(MethodSetup)
    tree.write(output_xml)
    # format it with 'xmlli

def variable_importance(classifier, features, fileNameIn):
#    print "mylist", features[:]
    indices      = np.argsort(classifier.feature_importances_)[::-1]
    importances   = classifier.feature_importances_

#    for f in range(len(features)):
#        print("%d. feature %s %d (%f)" % (f + 1, features[indices[f]], indices[f], importances[indices[f]]))

    fileName = fileNameIn + "_FeatureImporance.pdf"
    plt.figure()
    plt.title("Feature importances")
    plt.bar(range(len(features)), importances[indices],
       color="r",  align="center")
    plt.xticks(range(len(features)), indices)
    plt.xlim([-1, len(features)])
    plt.savefig(fileName)
