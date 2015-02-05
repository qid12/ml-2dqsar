#!/usr/bin/python
# Extract Compounds' Features by RDKit.
# Songpeng Zu /zusongpeng@gmail.com/

#-- Import packages

# Import numpy and pandas for data structures in python
import numpy as np
from pandas import DataFrame, Series
import pandas as pd

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
# Import urllib2 to get the the compounds from CHEMBL directly.
import urllib2
## Import the regular expression in python.
#import re

#-- Define a global variable to restrict the lines to read at one time.
MAXLINES = 1000

#-- Definitions of functions.
def readFILE(filen):
    with open(filen) as f:
        compounds = f.readlines()
    f.close()
    return(compounds)

def getSMILESfromCHEMBL(chemblid):
    """Generate the SMILES format of a compound recorded as CHEMBL ID.
    The record is like "CHEMBLE116"

    Keyword arguments:
    chemblid: the compound's CHEMBL ID.

    Return: a string represented the SMILES format.
    """
    webstr = "https://www.ebi.ac.uk/chembldb/download_helper/getmol/"
    comid = chemblid.replace("CHEMBL")
    cominfor = urllib2.urlopen(webstr+comid)
    d = cominfor.read()
    m = Chem.MolFromMolBlock(cominfor.read())
    return(Chem.MolToSmiles(m))

def phychemDataFrame(chempandas,smicol):
    """ Generate the physicochemical properties of the compounds.
    The compounds are stored in the DataFrame Structure defined by pandas.

    Keyword arguments:
    chempandas: the compounds stored in DataFrame, which contain the name and SMILES as columns.
    smicol: the column number of SMILES in the DataFrame.

    Return: a DataFrame of the compounds merging chempadas and the phychem by columns. 
    If None is detected given a SMILES-like string, it would be not deleted. 
    """
    assert chempandas.shape[0] <= MAXLINES
    molsmitmp = [Chem.MolFromSmiles(x) for x in chempandas.iloc[:,smicol]]
    molsmi = [x for x in molsmitmp if x is not None]
    nms = [x[0] for x in Descriptors._descList]
    nms.remove("MolWt")
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
    descrs = [calc.CalcDescriptors(x) for x in molsmi]
    descrsmat = np.matrix(descrs)
    df = DataFrame(descrsmat, index = [Chem.MolToSmiles(x) for x in molsmi], columns=nms)
    df['SMILES'] = df.index
    return(df)

def phychemFileStream(inputfile,smicol,outputfile,chunknum=1000):
    """ Generate the physicochemical properties of the compounds from files.
    The compounds are recorded as the SMILES format.

    Keyword argument:
    inputfile: the file of compounds in which each row a compound.
    smicol: the column number of SMILES in the DataFrame.
    outputfile: the file of results in which each row a compound and each column a physiclchemical feature.
    chunknum: number of lines are read once a time due to the file might be quite large.
              default is 1000, should be smaller than MAXLINES.
    
    Return: outputfile containing these properties.
    """
    assert chunknum <= MAXLINES
    chunker = pd.read_table("chemble2cid2smiles.txt",sep="\t",chunksize=chunknum)
    nms = [x[0] for x in Descriptors._descList]
    nms.remove("MolWt")
#    title = "\t".join(nms.append("SMILES"))+"\n"
#    f = open(outputfile,"w")
#    f.write(title)
#    f.close()
    for piece in chunker:
        df = phychemDataFrame(piece,smicol)
#       dfnew = df.append(Series(list(df.columns.values),index=list(df.columns.values)),ignore_index=True)
        df.to_csv(outputfile,mode="a",header=False,sep="\t",index=False)
    return(0)

def main(inputfile, smicol,outputfile):
    phychemFileStream(inputfile,smicol,outputfile)


if __name__ == "__main__":
    inputfile = "chemble2cid2smiles.txt"
    outputfile = "phychemRDKit.txt"
    smicol = 2
    main(inputfile,smicol,outputfile)
    
