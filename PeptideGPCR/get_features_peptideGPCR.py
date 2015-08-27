#!/usr/bin/python
# Extract features for PeptideGPCR.
# Songpeng Zu /zusongpeng@gmail.com/

#-- Import module
import sys
import os
sys.path.append(os.path.abspath("../FeatureExtraction"))

import ComFeatureExtract as cfe
import pandas as pd
from pandas import DataFrame, Series

# Import rdkit package for chemical feature extraction.
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
# Import 2D Pharmacophore Fingerprints
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
# Import Fingerprints.
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs, Torsions
# Import Descriptor Calculation
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

#sys.path.append("/usr/share/doc/rdkit-doc/Book/")

#-- Extract phychem and fingerprints for Peptide-GPCR.
def get_phychem_from_CPI_file(CPI_file,phychem_file,outfile):
    """
    Get the compounds' phy-chem characters from phychem_file.
    """
    CPIs = pd.read_table(CPI_file, sep="\t")
    known_phychem_chembl = pd.read_table(phychem_file,sep="\t",header=None)
    nms_phychem = [x[0] for x in cfe.Descriptors._descList]
    nms_phychem.remove("MolWt")
    nms_phychem.extend(['smiles','compound'])
    known_phychem_chembl.columns = nms_phychem
    phychem_in_CPIs = DataFrame.merge(CPIs,known_phychem_chembl,how='inner')
    phychem_in_CPIs.to_csv(outfile,mode="w",header=True, sep="\t",index=False)

def get_molsmi_from_CPI_file(CPI_file,ChEMBL2Smil_file):
    ChEMBL2Smil = pd.read_table(ChEMBL2Smil_file,header=None)
    ChEMBL2Smil.columns = ['compound','smiles']
    cpi_data = pd.read_table(CPI_file,sep="\t")
    cpi_smile = DatFrame.merge(cpi_data,ChEMBL2Smil,how='inner')
    molsmitmp = [Chem.MolFromSmiles(x) for x in cpi_smile['smiles']]
    i = 0
    molsmi = []
    for x in molsmitmp:
        if x is not None:
            x.SetProp("_Name", cpi_smile['compound'][i])
            molsmi.append(x)
        i += 1
    return(molsmi)

def get_fingerprint_from_DataFrame(chem_smile,fpfunc):
    molsmitmp = [Chem.MolFromSmiles(x) for x in chem_smile['smiles']]
    i = 0
    molsmi = []
    for x in molsmitmp:
        if x is not None:
            x.SetProp("_Name",chem_smile['compound'].loc[i])
            molsmi.append(x)
        i += 1
    fps = [fpfunc(x) for x in molsmi]
    fpsmat = np.matrix(fps)
    del fps
    df = DataFrame(fpsmat,index = [x.GetProp("_Name") for x in molsmi])
    del fpsmat
    df.insert(0,'chembl',df.index)
    df.insert(1,'smiles',[Chem.MolToSmiles(x) for x in molsmi])
    #df[1,'SMILES'] = [Chem.MolToSmiles(x) for x in molsmi]
    #df[0,'CHEMBL'] = df.index
    return(df)

def get_fingerprint_from_CPI_file(CPI_file,ChEMBL2Smil_file,fpfunc,outfilepath):
    """
    Get different kinds of fingerprints for Peptide-GPCR from RDKit.
    """
    ChEMBL2Smil = pd.read_table(ChEMBL2Smil_file,header=None)
    ChEMBL2Smil.columns = ['compound','smiles']
    cpi_data = pd.read_table(CPI_file,sep="\t")
    cpi_smile = DataFrame.merge(cpi_data,ChEMBL2Smil,how='inner')
    del cpi_data
    df = get_fingerprint_from_DataFrame(cpi_smile,fpfunc)
    df.to_csv(outfilepath,header=True,sep="\t",index=False)
#    chunk = 2000 # default chunksize
#    chunker = pd.read_table(CPI_file,sep="\t",chunksize=chunk)
#    for piece in chunker:
#        #piece.columns = ['protein','compound','affinity']
#        chembl_with_smiles = DataFrame.merge(piece,ChEMBL2Smil,how='inner')
#        df = get_fingerprint_from_DataFrame(chembl_with_smiles,fpfunc)
#        df.to_csv(outfilepath,mode="a",header=True,sep="\t",index=False)

def get_Morgan_from_CPI_file(CPI_file,ChEBML2Smil_file,outfilepath,radius=2):
    """
    Special function for Morgan finerprints. Set useFeatures=True, nBits=1024.
    User should set the radius, whose default value is 2 (like 4 in ECFP).
    """
    ChEMBL2Smil = pd.read_table(ChEBML2Smil_file,header=None)
    ChEMBL2Smil.columns = ['compound','smiles']
    cpi_data = pd.read_table(CPI_file, sep="\t")
    cpi_smile = DataFrame.merge(cpi_data,ChEMBL2Smil,how='inner')
    del cpi_data
    molsmitmp = [Chem.MolFromSmiles(x) for x in cpi_smile['smiles']]
    i = 0
    molsmi = []
    for x in molsmitmp:
        if x is not None:
            x.SetProp("_Name",cpi_smile['compound'][i])
            molsmi.append(x)
        i += 1
    fps = [AllChem.GetMorganFingerprintAsBitVect(x,radius,useFeatures=True,nBits=1024)
           for x in molsmi]
    fpsmat = np.matrix(fps)
    del fps
    df = DataFrame(fpsmat,index=[x.GetProp("_Name") for x in molsmi])
    del fpsmat
    df.insert(0,'chembl',df.index)
    df.insert(0,'smiles',[Chem.MolToSmiles(x) for x in molsmi])
    df.to_csv(outfilepath,header=True,sep="\t",index=False)

def get_topology_from_CPI_file(CPI_file, ChEMBL2Smile_file,outfilepath):
    """
    Get the topology fingerprints. Default output by RDKit is ExplicitBitVect. (A special class in RDKit)
    Here we transform it into general array, then matrix.
    """
    molsmi = get_molsmi_from_CPI_file(CPI_file, ChEMBL2Smile_file)
    fps = [FingerprintMols.FingerprintMol(x) for x in molsmi]
    list_array = [np.fromstring(DataStructs.cDataStructs.ExplicitBitVect.ToBitString(x),'u1') - ord('0')
                  for x in fps]
    fpsmat = np.matrix(np.array(list_array))




def main():
    # Input file path
    common_path = "/home/zusongpeng/lab/TransferCPIs/QSARMulT/"

    CPI_file = "CPIs_Peptide GPCR.txt" # Data of Peptide-GPCR.
    CPI_fullpath = common_path + "OriginalData/" + CPI_file

    phychem_file = "phychemRDKit_smiles2chembl.txt"
    phychem_fullpath = common_path + "FeatureExtraction/" + phychem_file
    chembl2smiles_file = "chembl2smilesStripSalts_tab"
    ChEMBL2Smil_fullpath = common_path+"FeatureExtraction/"+chembl2smiles_file

    #-- Get the phy-chem
    out_file = "phychem_compounds_PeptideGPCR"
    out_file_fullpath = common_path + "PeptideGPCR/" + out_file
    get_phychem_from_CPI_file(CPI_fullpath,phychem_fullpath, out_file_fullpathi,chunk=1000)

    #-- Get the fingerprints.
    out_file_fp_tag = "_fingerprint_PeptideGPCR"
    # MACCSKey
    macc_pre = "MACCSKeys"
    out_macc = common_path + "PeptideGPCR/" + macc_pre + out_file_fp_tag
    get_fingerprint_from_CPI_file(CPI_fullpath,ChEMBL2Smil_fullpath,
                                  MACCSkeys.GenMACCSKeys, out_macc)
    # TOPOLOGY
    topo_pre = "TOPOLOGY"
    out_topo = common_path + "PeptideGPCR/" + topo_pre + out_file_fp_tag
    get_fingerprint_from_CPI_file(CPI_fullpath, ChEMBL2Smil_fullpath,
                                  FingerprintMols.FingerprintMol,out_topo)
    # Atompair is huge, and we ignore it here.
#    atompair_pre = "ATOMPair"
#    out_atom = common_path + "PeptideGPCR/" + atompair_pre + out_file_fp_tag
#    get_fingerprint_from_CPI_file(CPI_fullpath, ChEMBL2Smil_fullpath,
#                                  Pairs.GetAtomPairFingerprintAsBitVect,
#                                  out_atom)
    # Torsions is huge, and we ignore it here.
#    torsions_pre = "TORSIONS"
#    out_torsions = common_path + "PeptideGPCR/" + torsions_pre + out_file_fp_tag
#    get_fingerprint_from_CPI_file(CPI_fullpath, ChEMBL2Smil_fullpath,
#                                  Torsions.GetTopologicalTorsionFingerprintAsIntVect,
#                                  out_torsions)
    # Morgan
    morgan_pre = "Morgan"
    out_morgan = common_path + "PeptideGPCR/" + morgan_pre + out_file_fp_tag
    radius = 2 # Set the radius for Morgan.
    get_Morgan_from_CPI_file(CPI_fullpath, ChEMBL2Smil_fullpath,
                             out_morgan,radius)

if __name__ == "__main__":
    main()
