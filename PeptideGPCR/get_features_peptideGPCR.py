#!/usr/bin/python
# Extract features for PeptideGPCR.
# Songpeng Zu /zusongpeng@gmail.com/

#-- Import module
import sys
import os
sys.path.append(os.path.abspath("../FeatureExtraction"))

# Need to reload these ?
import ComFeatureExtract as cfe
import pandas as pd
from pandas import DataFrame, Series

#-- Import Peptide-GPCR related Compounds
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

def main():
    common_path = "/home/zusongpeng/lab/TransferCPIs/QSARMulT/"
    CPI_file = "CPIs_Peptide GPCR.txt"
    CPI_fullpath = common_path + "OriginalData/" + CPI_file
    phychem_file = "phychemRDKit_smiles2chembl.txt"
    phychem_fullpath = common_path + "FeatureExtraction/" + phychem_file
    out_file = "phychem_compounds_PeptideGPCR"
    out_file_fullpath = common_path + "PeptideGPCR/" + out_file
    get_phychem_from_CPI_file(CPI_fullpath,phychem_fullpath, out_file_fullpath)

if __name__ == "__main__":
    main()
