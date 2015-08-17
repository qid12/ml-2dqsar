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
    known_phychem_chembl = pd.read_table(phychem_file,sep="\t")



def main():
    common_path = "~/lab/TransferCPIs/QSARMulT/"
    CPI_file = "CPIs_Peptide GPCR.txt"
    CPI_fullpath = common_path + "OriginalData/" + CPI_file
    phychem_file = "phychemRDKit_smiles2chembl.txt"
    phychem_fullpath = common_path + "FeatureExtraction/" + phychem_file
    out_file = "phychem_compounds_PeptideGPCR"
    out_file_fullpath = common_path + "PeptideGPCR" + out_file